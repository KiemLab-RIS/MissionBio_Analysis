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
        if (entry.second > 1000) { // 1000 for full size
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
        std::map<std::string,std::vector<std::string>> editMap = pData->cellEdit[cellBarcode];
        for (auto editPair:editMap) {
            std::string ampID = editPair.first;
            for (auto s:editPair.second) {
                file << cellBarcode << "\t" << ampID << "\t" << s << "\n";
            }
        }
    }
    file.close();
}
//
// edit counte
//
std::pair<int,int> calculateEditCount(std::vector<std::string> seqVec,std::string gene,EditMap &edits) {
    if (seqVec.size() == 0) {return(std::pair<int,int>(0,0));}
    int edited = 0;
    int total = 0;
    for (auto seq:seqVec) {
    	total += 1;
	auto entryVector = edits[gene];
	for (auto entry:entryVector) {
	     int location = entry.location;
	     char base = entry.base;
	     if (seq[location] == base){
	       edited += 1;
	       break;  // if even one base is edited, sequence is edited
	     }
	}
    }
    return(std::pair<int,int>(edited, total));
}
//
// edit detail
//
void editDetail(
		std::vector<std::string> seqVec,
		std::string gene,
		EditMap &edits,
		std::ofstream &ofile,
		std::string bc) {
    //
    // recode details of edit sites
    //
    if (seqVec.size() == 0) {return;}
    std::map<std::string,int> seqMap;
    for (auto seq:seqVec) {
        std::string so = gene + ",";
	auto entryVector = edits[gene];
	for (auto entry:entryVector) {
	     int location = entry.location;
             so = so + seq[location];
	}
	seqMap[so] += 1;
    }
    for (auto pair:seqMap) {
      ofile << bc << "," << pair.first << "," << pair.second << "\n";
    }
    return;
}
//
// the cellEdit count may not match the ampID count if the ReadFastq  could not build the sequence
//
void Output::bigTable(CELLDATA *pData,Params &pr, EditMap &edits,const std::string &filename) {
    //
    // go through ampID entries for each barcpde
    //
    std::ofstream file(filename);
    size_t lastDot = filename.find_last_of('.');
    std::string seqFileName = filename.substr(0, lastDot) + "_seq.txt";
    std::ofstream seqFile(seqFileName);
    seqFile << "#cell sequence file\n";
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
	seqFile << "#\n";
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
        
	//
	// CD33 good
	//
	seqVec = cellEditMap["CO-819_AMP47"];
	//
	// separate CD33 from pseudo
	//
	std::vector<std::string> cd33Vec;
	std::vector<std::string> cd33PVec;
	for (auto seq:seqVec) {
          if (seq[29] == 'T') {
	    cd33Vec.push_back(seq);
	  } else {
	    cd33PVec.push_back(seq);
	  }
	}
	//std::cout << "AMP47 CD33 sequece count = " << seqVec.size() << "\n";
        auto cd33e = calculateEditCount(cd33Vec,"CD33",edits);
	//auto cd33e = cd33EditCount(seqVec,edits);
        file << cd33e.first << "\t" << cd33e.second << "\t";
        editDetail(cd33Vec,"CD33",edits,seqFile,cellBarcode);
	//std::cout << "CD33 edit count = " << cd33e.first << " out of " << cd33e.second << "\n";
	//
	// CD33 Pseudo
	//
	//auto cd33pe = cd33PseudoEditCount(seqVec);
        auto cd33pe = calculateEditCount(cd33PVec,"CD33P",edits);
	//std::cout << "CD33 Pseudo edit count = " << cd33pe.first << " out of " << cd33pe.second << "\n";
        file << cd33pe.first << "\t" << cd33pe.second << "\t";
	editDetail(cd33PVec,"CD33P",edits,seqFile,cellBarcode);

        //
        //
        //
        seqVec = cellEditMap["HbG_AMP1"];
        auto hbg1e = calculateEditCount(seqVec,"HBG1",edits);
        file << hbg1e.first << "\t" << hbg1e.second << "\t";
	editDetail(seqVec,"HBG1",edits,seqFile,cellBarcode);
        //
        //
        //
        seqVec = cellEditMap["HbG_AMP2"];
        auto hbg2e = calculateEditCount(seqVec,"HBG2",edits);
        file << hbg2e.first << "\t" << hbg2e.second; // note not \t here for AB loop
	editDetail(seqVec,"HBG2",edits,seqFile,cellBarcode);

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

        file << "\n";
    }
    file.close();
}
