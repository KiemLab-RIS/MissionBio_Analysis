//
//  Params.cpp
//  mbFilter
//
//  Created by mark enstrom on 1/15/25.
//
#include <string>
#include <iostream>
#include <fstream>
#include "Params.hpp"
//-------------------------------------------------------------
//
// split by tabs
//
//-------------------------------------------------------------
std::vector<std::string>
splitByTabs(const std::string &str) {
    std::vector<std::string> result;
    size_t start = 0;
    size_t end = str.find('\t');

    while (end != std::string::npos) {
        result.push_back(str.substr(start, end - start));
        start = end + 1;
        end = str.find('\t', start);
    }

    // Add the last substring after the final tab
    result.push_back(str.substr(start));

    return result;
}
//-------------------------------------------------------------
//
// look up index from string
//
//-------------------------------------------------------------
int Params::getID(const std::string &ampID) {
    for (int i = 0; i < (int)paramData.size(); i++) {
        auto p = paramData[i];
        if (p->ampID == ampID) {
            return(i);
        }
    }
    return(-1);
}

//-------------------------------------------------------------
//
// read primer file
// <primer format>
// AmpID    chr    fwd_seq    rev_seq    amplicon_seq    Targets    Genes
//
// full id = ampID + targets
//
//-------------------------------------------------------------
bool Params::readPrimerFile(std::string filename) {
   
    // Open the FASTQ file
    std::ifstream primerFile(filename);
    if (!primerFile.is_open()) {
        std::cerr << "Error opening file " << filename << std::endl;
        return false;
    }

    std::string line;
    // skip header
    std::getline(primerFile, line);
    if (line == "") {
        return true;
    }
    while (1) {
        //
        // read1
        //
        std::getline(primerFile, line);
        if (line == "") {
            break;
        }
        //
        // break line by tabs
        //
        std::vector<std::string> items = splitByTabs(line);
        PENTRY p = new(ENTRY);
        p->ampID    = items[0];
        p->chr      = items[1];
        p->fwdSeq   = items[2];
        p->revSeq   = items[3];
        p->amplicon = items[4];
        p->target   = items[5];
        p->gene     = items[6];
        
        paramData.push_back(p);
    }
    
    
    std::cout << paramData.size() << " items in param file\n";
    
    return true;
}

//-------------------------------------------------------------
//
// must be same len, allow 'N' for second param
//
//-------------------------------------------------------------
int stringDiff(const std::string &s1, const std::string &s2,int maxError) {
    if (s1.size() != s2.size()) {
        return -1;
    }
    int error = 0;
    for (int index = 0; index < (int)s1.size(); index++) {
        if (s1[index] != s2[index]) {
            if (s2[index] != 'N') {
                error += 1;
                if (error > maxError) {
                    return (error);
                }
            }
        }
    }
    return(error);
}
//-------------------------------------------------------------
//
// find a match of start of seq to start primer in params
//
//-------------------------------------------------------------

std::pair<std::string,int> Params::identifySequence(std::string seq) {
    for (size_t paramIndex = 0; paramIndex < paramData.size(); paramIndex++) {
        PENTRY p = paramData[paramIndex];
        size_t sFwd = p->fwdSeq.length();
        std::string ts = seq.substr(0,sFwd);
        int m = stringDiff(p->fwdSeq, ts, 2);
	// allow 1 error in primer
        if (m <= 1) {
            //
            // check amplicon score
            //
            std::string s1;
            std::string s2;
            
            if (p->amplicon.length() > seq.length()) {
                s1 = p->amplicon.substr(0,seq.length());
                s2 = seq;
            } else {
                s1 = p->amplicon;
                s2 = seq.substr(0,p->amplicon.length());
            }
            int a = stringDiff(s1, s2, 20);
	    // allow 19 errors in fwd seq comp with amplicon
            if (a > 19) {
                return std::make_pair("",-1);
//                int start = 0;
//                int l   = 100;
//                while (start+l < (int)s1.length()) {
//                    std::cout << "\n";
//                    std::cout << s1.substr(start,l) << "\n";
//                    std::cout << s2.substr(start,l) << "\n";
//                    std::cout << "\n";
//                    start += l;
//                }
            }
            return std::make_pair(p->ampID,paramIndex);
        }
    }
    return std::make_pair("",-1);
}

//-------------------------------------------------------------
//
//  Function to build the combined sequence
//
//-------------------------------------------------------------
// Function to compare two sequences and return the number of errors
std::pair<int, std::string> compSeq(const std::string& s1, const std::string& s2) {
    int errors = 0;
    std::string comparison;
    for (size_t i = 0; i < s1.size() && i < s2.size(); ++i) {
        if ((s1[i] != 'N') & (s2[i] != 'N') & (s1[i] != s2[i])) {
            errors++;
            comparison += 'N';  // Mark mismatch with 'N'
        } else {
            comparison += s1[i];
        }
    }
    return {errors, comparison};
}
//-------------------------------------------------------------
//
//  Function to build the combined sequence
//
//-------------------------------------------------------------
std::string
Params::buildSeq(const std::string& fwd_seq, const std::string& seq2,const std::string& ampID) {
    // Get the amplicon from the primers table
    int tableID = getID(ampID);
    if (tableID == -1) {return "";}
    std::string amp = paramData[tableID]->amplicon;
    int overlap = (int)fwd_seq.size() + (int)seq2.size() - (int)amp.size();
    std::string bSeq;

    if (overlap >= 0) {
        // There is overlap
        bSeq = fwd_seq + seq2.substr(overlap);

        // Check if overlap bases are equal
        int x1 = (int)fwd_seq.size() - overlap;
//        int x2 = overlap;
//        int x3 = (int)seq2.size() - overlap;

        std::string t1 = fwd_seq.substr(x1);
        std::string t2 = seq2.substr(0, overlap);

        auto [e, ts] = compSeq(t1, t2);
        if (e >= 4) {

            //std::cout << "-------------------------------------\n";
            //std::cout << " fail of " << ampID << "\n";
	    //std::cout << "e = " << e << "\n";
	    //std::cout << "t1 " << t1 << "\n";
	    //std::cout << "t2 " << t2 << "\n";
	    //std::cout << "ts " << ts << "\n";
	    //std::cout << "x1 " << x1 << "\n";
	    //std::cout << "fw " << fwd_seq << "\n";
	    //std::cout << "s2 " << seq2 << "\n";
	    //std::cout << "am " << amp << "\n";
	    //std::cout << "ov " << overlap << "\n";
            return {""};
        } else {
	    //std::cout << "----- good " << ampID << "-----\n";
	}
    } else {
        // No overlap, must fill with 'N'
        int fill = -overlap;
        bSeq = fwd_seq + std::string(fill, 'N') + seq2;
    }

    // Check for errors in the final sequence
    auto [errors, disp] = compSeq(amp, bSeq);
    
    if (errors < 30) {
        return bSeq;
    } else {
        //std::cout << "-------------------------------------\n";
        //std::cout << " fail of " << ampID << "\n";
        //std::cout << amp.substr(0,100) << "\n";
        //std::cout << bSeq.substr(0,100) << "\n";
        //std::cout << disp.substr(0,100) << "\n";
        //std::cout << "\n";
        //std::cout << amp.substr(100,100) << "\n";
        //std::cout << bSeq.substr(100,100) << "\n";
        //std::cout << disp.substr(100,100) << "\n";
        return "";
    }
}




bool META::readMetaFile(std::string filename) {
    std::string line;
    std::ifstream metaFile(filename);
    if (!metaFile.is_open()) {
        std::cerr << "Error opening file " << filename << std::endl;
        return false;
    }
    while (1) {
        //
        // read1
        //
        std::getline(metaFile, line);
        if (line == "") {
            break;
        }
        //
        // break line by tabs
        //
        std::vector<std::string> read1 = splitByTabs(line);
        std::getline(metaFile, line);
        if (line == "") {
            std::cout << "bad metafile format\n";
            return false;
        }
        std::vector<std::string> read2 = splitByTabs(line);
        if ( (read1[0] == "D") & (read1[1] == "1") &
             (read2[0] == "D") & (read2[1] == "2")) {
            vec_DNA_R1.push_back(read1[2]);
            vec_DNA_R2.push_back(read2[2]);
        } else if ((read1[0] == "P") & (read1[1] == "1") &
                   (read2[0] == "P") & (read2[1] == "2")) {
            vec_PROTEIN_R1.push_back(read1[2]);
            vec_PROTEIN_R2.push_back(read2[2]);

        }
    }
    return true;
}
