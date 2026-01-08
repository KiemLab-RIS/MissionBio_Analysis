//
//  Output.cpp
//  mbFilter
//
//  Created by mark enstrom on 1/15/25.
//

#include "Output.hpp"
#include <algorithm>

//-------------------------------------------------------------
//
// vector compare functions
//
//-------------------------------------------------------------
bool cmp(std::pair<std::string, int>& a,
         std::pair<std::string, int>& b)
{
    return a.second > b.second;
}
bool cmp1(std::pair<std::string,std::map<std::string,int>>& a,
          std::pair<std::string,std::map<std::string,int>>& b)
{
    return a.second.size() > b.second.size();
}


void Output::write(CELLDATA *pData, const std::string &filename) {
    std::ofstream file(filename);
    std::cout << "cells size = " << pData->cells.size() << std::endl;
    //
    // sort barcodes by captures
    //
    std::vector<std::pair<std::string,int>> i1;
    for (auto entry:pData->barcodes) {
        if (entry.second > 1000) {
            i1.push_back(entry);
        }
    }
    std::cout << "i1 size = " << i1.size() << std::endl;
    sort(i1.begin(),i1.end(),cmp);
    
    for (auto pair:i1) {
        std::string cellBarcode = pair.first;
        int count = pair.second;
        file << "#" << cellBarcode << "\t" << count << "\n";
        //
        // cellID data
        //
        std::map<std::string,int> capID = pData->cells[cellBarcode];
        for (auto cellPair:capID) {
            std::string ampID = cellPair.first;
            int count = cellPair.second;
            file << cellBarcode << "\t" << ampID << "\t" << count << "\n";
        }
        //
        // celledit data
        //
//        std::map<std::string,std::vector<std::string>> editMap = pData->cellEdit[cellBarcode];
//        for (auto editPair:editMap) {
//            std::string ampID = editPair.first;
//            for (auto s:editPair.second) {
//                file << cellBarcode << "\t" << ampID << "\t" << s << "\n";
//            }
//        }
    }
    file.close();
}
//
// CD33 has real and pseudo rate
//
std::pair<int,int> cd33EditCount(std::vector<std::string> seqVec) {
    if (seqVec.size() == 0) {return(std::pair<int,int>(0,0));}
    int edited = 0;
    int total = 0;
    for (auto seq:seqVec) {
        // check not pseudogene
        if (seq[29] == 'T') {
          total += 1;
          //std::cout << "CD33   " << seq[184] << " " << seq[186] << "\n";
          if ((seq[111] != 'A') | (seq[109] != 'A') ) {
              edited += 1;
          }
	}
    }
    return(std::pair<int,int>(edited, total));
    
}
std::pair<int,int> cd33PseudoEditCount(std::vector<std::string> seqVec) {
    if (seqVec.size() == 0) {return(std::pair<int,int>(0,0));}
    int edited = 0;
    int total = 0;
    for (auto seq:seqVec) {
        // check IS pseudogene
        if (seq[29] != 'T') {
          total += 1;
          //std::cout << "CD33   " << seq[184] << " " << seq[186] << "\n";
          if ((seq[111] != 'A') | (seq[109] != 'A') ) {
              edited += 1;
          }
	}
    }
    return(std::pair<int,int>(edited, total));
}
//
//
//
int cd117EditCount(std::vector<std::string> seqVec) {
    if (seqVec.size() == 0) {return 0;}
    int edited = 0;
    int total = 0;
    for (auto seq:seqVec) {
        total += 1;
        //std::cout << "CD117   " << seq[183] << " " << seq[184] << " " << seq[186] << seq[188] << " " << "\n";
        if ((seq[183] == 'C') | (seq[184] == 'C') | (seq[186] == 'C') | (seq[188] == 'C')) {
            edited += 1;
        }
    }
    return(edited);
}

double hbg1EditCount(std::vector<std::string> seqVec) {
    if (seqVec.size() == 0) {return 0;}
    int edited = 0;
    int total = 0;
    for (auto seq:seqVec) {
        total += 1;
        //std::cout << "HBG1   " << seq[265] << " " << seq[268] << " " << seq[269] << seq[271] << " " << "\n";
        if ((seq[265] == 'G') | (seq[268] == 'G') | (seq[269] == 'G')  | (seq[271] == 'G')) {
            edited += 1;
        }
    }
    return(edited);
}

double hbg2EditCount(std::vector<std::string> seqVec) {
    if (seqVec.size() == 0) {return 0;}
    int edited = 0;
    int total = 0;
    for (auto seq:seqVec) {
        total += 1;
        //std::cout << "HBG2   " << seq[266] << " " << seq[269] << " " << seq[270] << seq[272] << " " << "\n";
        if ((seq[266] == 'G') | (seq[269] == 'G') | (seq[270] == 'G') | (seq[272] == 'G')) {
            edited += 1;
        }
    }
    return(edited);
}
//
// abBarcodeVec and abNameVec are constants defined in ReadProtein
//
//
// class CELLDATA {
// public:
//     std::map<std::string, int> barcodes;
//     std::map<std::string,std::map<std::string,int>> cells;   //  just ampID count
//     std::map<std::string,std::map<std::string,int>> scaddenMap;
//     std::map<std::string,std::map<std::string,std::vector<std::string>>> cellEdit;  // full sequence data
// };
//
// the cellEdit count may not match the ampID count if the ReadFastq  could not build the sequence 
//
void Output::bigTable(CELLDATA *pData,Params &pr, const std::string &filename) {
    //
    // go through ampID entries for each barcpde
    //
    std::ofstream file(filename);
    std::cout << "cells size = " << pData->cells.size() << std::endl;
    //
    // header...ampID then edit
    //
    file << "CellBarcode\tnCap\t";
    for (int i = 0; i < (int)pr.items(); i++) {
        std::string ampID = pr.amp(i);
        file << ampID << "\t";
    }
    //
    // header for edit/not edit counts
    //
    file << "CD117_e\tCD117_t\tCD33_e\tCD33_t\tCD33P_e\tCD33P_t\tHBG1_e\tHBG1_t\tHBG2_e\tHBG2_t";
    //
    // proteins
    //
    for (auto abID:abNameVec) {
        file << "\t" << abID;
    }
    file << "\n";
    //
    // sort barcodes by captures
    //
    std::vector<std::pair<std::string,int>> i1;
    for (auto entry:pData->barcodes) {
        if (entry.second > 10) {
            i1.push_back(entry);
        }
    }
    std::cout << "cells above limit = " << i1.size() << std::endl;
    sort(i1.begin(),i1.end(),cmp);
    //
    // for each cell barcode find the amp data
    // then cell edit data
    //
    for (auto pair:i1) {
        std::string cellBarcode = pair.first;
        int count = pair.second;
        file << cellBarcode << "\t" << count << "\t";
        //
        // cellID data
        //
        // for each ampID
        //
        std::map<std::string,int> capID = pData->cells[cellBarcode];
        for (int i = 0; i < (int)pr.items(); i++) {
            std::string ampID = pr.amp(i);
            auto count = capID[ampID];
            file << count << "\t";
        }
        //
        // cell edit data for cd117,hbg1,hbg2
        //
        auto cellEditMap = pData->cellEdit[cellBarcode];
        auto seqVec = cellEditMap["CO-819_AMP29"];
        int editCount = cd117EditCount(seqVec);
        file << editCount << "\t" << seqVec.size() << "\t";
        //std::cout << "cd117 seq size = " << seqVec.size() << " r = " << r << "\n";
	//
	// CD33 good
	//
	seqVec = cellEditMap["CO-819_AMP47"];
	//std::cout << "AMP47 CD33 sequece count = " << seqVec.size() << "\n";
	auto cd33e = cd33EditCount(seqVec);
        file << cd33e.first << "\t" << cd33e.second << "\t";
	//std::cout << "CD33 edit count = " << cd33e.first << " out of " << cd33e.second << "\n";
	//
	// CD33 Pseudo
	//
	auto cd33pe = cd33PseudoEditCount(seqVec);
	//std::cout << "CD33 Pseudo edit count = " << cd33pe.first << " out of " << cd33pe.second << "\n";
        file << cd33pe.first << "\t" << cd33pe.second << "\t";
        //
        //
        //
        seqVec = cellEditMap["HbG_AMP1"];
        editCount = hbg1EditCount(seqVec);
        file << editCount << "\t" << seqVec.size() << "\t";
        //
        //
        //
        seqVec = cellEditMap["HbG_AMP2"];
        editCount = hbg2EditCount(seqVec);
        file << editCount << "\t" << seqVec.size(); // note not \t here for AB loop
        //
        // proteins
        //
        
        std::string cellType = "NA";
        
        for (auto abID:abNameVec) {
            auto count = capID[abID];
            if (abID == "pCellType") {
                file << "\t" << cellType;
            } else {
                file << "\t" << count;
            }
        }
        //
        // scadden barcodes
        //
        
        auto scaddenCellMap = pData->scaddenMap[cellBarcode];
        
        if (cellBarcode== "TTCACTGGA_CTTCTATGA") {
            std::cout << "found TTCACTGGA_CTTCTATGA l = " << scaddenCellMap.size() << "\n";
        }
        for (auto pair:scaddenCellMap) {
            file << "\t" << pair.first << "\t" << pair.second;
        }
        
        file << "\n";

    }
    file.close();
}
