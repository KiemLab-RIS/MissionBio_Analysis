//
//  ReadFastq.cpp
//  mbFilter
//
//  Created by mark enstrom on 1/15/25.
//

#include "ReadFastq.hpp"
//--------------------------------------------------------------------------------
//
// Function to get the complement of a DNA base
//
//--------------------------------------------------------------------------------

char dnacomp(char ch) {
    if (ch == 'A') return 'T';
    if (ch == 'T') return 'A';
    if (ch == 'C') return 'G';
    if (ch == 'G') return 'C';
    return ch;
}
//--------------------------------------------------------------------------------
//
// Function to get the reverse complement of a DNA sequence
//
//--------------------------------------------------------------------------------

std::string revcomp(const std::string& s) {
    std::string result;
    result.reserve(s.size());

    for (int i = (int)s.size() - 1; i >= 0; --i) {
        result.push_back(dnacomp(s[i]));
    }

    return result;
}


//--------------------------------------------------------------------------------
//
// checkQuality: how many bases are below Q
//
//--------------------------------------------------------------------------------
int checkQuality(const std::string qualSeq) {
    //
    // count how many bases are below low quality
    //
    int totalQ = 0;
    for (char c :qualSeq) {
        int ordinalValue = static_cast<int>(c);
        int ordinalBase = static_cast<int>('!');
        int q = ordinalValue - ordinalBase;
        totalQ += q;
    }
    double aveQ = (double)totalQ / (double)(qualSeq.length());
    int iQ = (int)floor(aveQ);
    return iQ;
}

//-------------------------------------------------------------
//
// setup
//
//-------------------------------------------------------------
ReadFastq::ReadFastq() {

}
//-------------------------------------------------------------
//
// read
//
//----------------------------------------------------------
void ReadFastq::read(Params params,CELLDATA *pData,std::string fileR1,std::string fileR2) {
    NWAlign nwAlign;
    std::cout << "Read from " << fileR1 << "\n";
    //
    // open
    //
    std::ifstream fastqFileR1(fileR1);
    if (!fastqFileR1.is_open()) {
        std::cerr << "Error opening file " << fileR1 << std::endl;
        return;
    }
    //
    // Open the FASTQ file
    //
    std::ifstream fastqFileR2(fileR2);
    if (!fastqFileR2.is_open()) {
        std::cerr << "Error opening file " << fileR2  << std::endl;
        return;
    }
    //
    // read both files together
    //
    std::string lineR1;
    std::string lineR2;
    std::string identifier;
    std::string sequenceR1;
    std::string sequenceR2;
    int failQ = 0;
    int failOffset = 0;
    int failID = 0;
    int failStartSeq = 0;
    int failHBG1 = 0;
    int failHBG2 = 0;
    int failCD33 = 0;
    //    // Define the map to store sequence data
    std::map<int,int> buildStat;
    std::map<std::string,int> badStartSeq;
    int count = 0;
    try {
        while (1) {
            //
            // read1
            //
            std::string barcode = "";
            std::getline(fastqFileR1,lineR1);
            if (lineR1 == "") {
                std::cout << "End of file reached\n";
                break;
            }
            count += 1;
            std::getline(fastqFileR1,sequenceR1);  // Line 2: Sequence
            std::getline(fastqFileR1,lineR1);      // Line 3: Plus line (+) (ignored)
            std::getline(fastqFileR1,lineR1);      // Line 4: Quality
            int aveQ1 = checkQuality(lineR1);
            buildStat[aveQ1] += 1;
            //
            // read2
            //
            std::getline(fastqFileR2,lineR2);      // Line 3: Plus line (+)
            std::getline(fastqFileR2,sequenceR2);  // Line 2: Sequence
            std::getline(fastqFileR2,lineR2);      // Line 3: Plus line (+)
            std::getline(fastqFileR2,lineR2);      // Line 4: Quality
            int aveQ2 = checkQuality(lineR2);
            //
            // filter bad reads
            //
            if ((aveQ1 < 22) | (aveQ2 < 22)) {
                failQ ++;
                continue;
            }
            //

            //
            // missionBio constant seq
            //
            std::string test2 = sequenceR1.substr(9+offset+9,15);
            if (test2 == seq1Const2) {
                std::string seqData1 = sequenceR1.substr(9+offset+9+15);
                std::string seqData2 = revcomp(sequenceR2);
                std::pair<std::string, int> result = params.identifySequence(seqData1);
                std::string ampID = result.first;
                int ampIndex = result.second;

		if (ampIndex == -1) {
		  // nostart primer found
		  failID ++;
		} else if (ampIndex == -2) {
		  // primer found but bad start sequence
	          failStartSeq++;
		  badStartSeq[ampID]++;
		} else {
                    //
                    // save a list of barcodes
                    // save ampID for non-edit sequences
                    // save full sequence for edit IDs
                    //
                    if ((ampID == geCD33) | (ampID == geHBG1) | (ampID == geHBG2)) {
		        //
			// try to build full sequence from R1/R2
			//
                        std::string full = params.buildSeq(seqData1, seqData2, ampID);
                        if (full != "") {
                            //std::cout << "good " << ampID << "\n";
                            pData->cellEdit[barcode][ampID].push_back(full);
                            pData->barcodes[barcode] += 1;
                            pData->cells[barcode][ampID] += 1;
                        } else {
			    //
			    // a fail here will make the ampID count and the edit seq count different
			    //
                            //std::cout << "fail build sequence " << ampID << "\n";
                            if (ampID == geHBG1) {
                                failHBG1++;
                            } else if (ampID == geHBG2) {
                                failHBG2++;
                            } else if (ampID == geCD33) {
			        failCD33++;
			    }

                            failID++;
                        }
                    } else {
                      pData->barcodes[barcode] += 1;
                      pData->cells[barcode][ampID] += 1;
		    }
                    //
                }
            }

            //
            // track
            //
            if ((count % 1000000) == 0) {
                std::cout << count << " cells = " << pData->cells.size() <<  "\n";
            }
//            if (count > 1000000) {
//              std::cout << "Break on line count\n";
//              break;
//            }
        }
    } catch (int x) {
        std::cout << "Exception Caught \n";
    }
    std::cout << "Done reading " << count << "\n";
    std::cout << pData->barcodes.size() << " barcodes\n";
    std::cout << pData->cells.size() << " cells\n";
    std::cout << pData->cellEdit.size() << " edit cells\n";
    std::cout << "reads        = " << count << "\n";
    std::cout << "failQ        = " << failQ << "\n";
    std::cout << "failOffset   = " << failOffset << "\n";
    std::cout << "failID       = " << failID << "\n";
    std::cout << "failStartSeq = " << failStartSeq << "\n";
    std::cout << "failHBG1     = " << failHBG1 << "\n";
    std::cout << "failHBG2     = " << failHBG2 << "\n";
    std::cout << "failCD33     = " << failCD33 << "\n";
    for (auto pair:badStartSeq){
    	std::cout << pair.first << "  " << pair.second << "\n";
    }
    //
    // stats on reads
    //
    //    for (auto pair:buildStat) {
    //        std::cout << pair.first << "  " << pair.second << "\n";
    //    }
    return;
}
