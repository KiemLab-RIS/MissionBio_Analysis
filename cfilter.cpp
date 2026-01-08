#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include<algorithm>
#include <stdexcept>
#include "NWAlign.hpp"


////-------------------------------------------------------------
////
//// must be same len, allow 'N' for second param
////
////-------------------------------------------------------------
//int stringDiff(const std::string &s1, const std::string &s2,int maxError) {
//  int error = 0;
//  for (int index = 0; index < (int)s1.size(); index++) {
//    if (s1[index] != s2[index]) {
//      if (s2[index] != 'N') {
//        error += 1;
//	if (error > maxError) {
//	  return (error);
//	}
//      }
//    }
//  }
//  return(error);
//}
////-------------------------------------------------------------
////
//// look for one of the sequence offset constants
////
////-------------------------------------------------------------
//int findOffset(const std::string &seq) {
//  int loc = seq.find(offc17);
//  if (loc == 9) return 17;
//
//  loc = seq.find(offc16);
//  if (loc == 9) return 16;
//
//  loc = seq.find(offc15);
//  if (loc == 9) return 15;
//
//  loc = seq.find(offc14);
//  if (loc == 9) return 14;
//
//  return(-1);
//}
////-------------------------------------------------------------
////
//// vector compare functions
////
////-------------------------------------------------------------
//bool cmp(std::pair<std::string, int>& a,
//	 std::pair<std::string, int>& b)
//{
//	return a.second > b.second;
//}
//bool cmp1(std::pair<std::string,std::map<std::string,int>>& a,
//		  std::pair<std::string,std::map<std::string,int>>& b)
//{
//	return a.second.size() > b.second.size();
//}
////-------------------------------------------------------------
////
//// main
////
////-------------------------------------------------------------
//int main(int argc, char* argv[]) {
//    if (argc != 5) {
//      std::cout << "usage: filter file_R1.fastq  file_R2.fastq outfile bcfile\n";
//      return(-1);
//    }
//    //
//    // read primer data
//    //
//    std::vector<std::vector<std::string>> primerData = readPrimerFile();
//    if (primerData.size() == 0) return(0);
//    //
//    // fastq files
//    //
//    std::string fileR1 = argv[1];
//    std::string fileR2 = argv[2];
//    std::string outFile = argv[3];
//    std::string errorFile = outFile.substr(0,outFile.size()-4) + "_error.csv";
//    std::string bcFile = argv[4];
//    // Define the map to store sequence data
//    std::map<std::string, int> barcodes;
//    std::map<std::string,std::map<std::string,int>> cells;
//    //
//    // Open the FASTQ file
//    //
//    std::ifstream fastqFileR1(fileR1);
//    if (!fastqFileR1.is_open()) {
//        std::cerr << "Error opening file " << fileR1 << std::endl;
//        return 1;
//    }
//    //
//    // Open the FASTQ file
//    //
//    std::ifstream fastqFileR2(fileR2);
//    if (!fastqFileR2.is_open()) {
//        std::cerr << "Error opening file " << fileR2  << std::endl;
//        return 1;
//    }
//    //
//    // read both files together
//    //
//    std::string lineR1;
//    std::string lineR2;
//    std::string identifier;
//    std::string sequenceR1;
//    std::string sequenceR2;
//    int count = 0;
//    try {
//      while (1) {
//   	//
//   	// read1
//   	//
//        std::string barcode = ""; 
//  	std::getline(fastqFileR1, lineR1);
//  	if (lineR1 == "") {
//	  std::cout << "End of file reached\n";
//	  break;
//	}
//        std::getline(fastqFileR1, sequenceR1);  // Line 2: Sequence
//        std::getline(fastqFileR1, lineR1);      // Line 3: Plus line (+)
//        std::getline(fastqFileR1, lineR1);      // Line 4: Quality (ignored)
//        int offset = findOffset(sequenceR1);
//	if ((int)sequenceR1.size() < 9 + offset + 9 + 15) {
//	  std::cout << "erorr in R1 "  << sequenceR1.size() << "\n";
//	  break;
//	}
// 	//
// 	// read2
// 	//
//        std::getline(fastqFileR2, lineR2);      // Line 3: Plus line (+)
//        std::getline(fastqFileR2, sequenceR2);  // Line 2: Sequence
//        std::getline(fastqFileR2, lineR2);      // Line 3: Plus line (+)
//        std::getline(fastqFileR2, lineR2);      // Line 4: Quality (ignored)
//	//
//	// if good MissionBio read then save
//	//
// 	if (offset != -1) {
//    	    barcode = sequenceR1.substr(0,9) + "_" + sequenceR1.substr(9 + offset, 9);
//	    std::string test2 = sequenceR1.substr(9+offset+9,15);
//	    if (test2 == seq1Const2) {
//    	    	barcodes[barcode] += 1;
//		std::string seqData = sequenceR1.substr(9+offset+9+15) + '\t' + sequenceR2;
//		cells[barcode][seqData] += 1;
//	    }
//	}
//	//
//	// track
//	//
//	count += 1;
//	if ((count % 100000) == 0) {
//	  std::cout << count << " cells = " << cells.size() <<  "\n";
//	}
//	if (count > 8000000000) {
//	  std::cout << "Break on line count\n";
//	  break;
//	}
//      }
//    } catch (int x) {
//        std::cout << "Exception Caught \n";
//    }
//    std::cout << "Done reading " << count << "\n";
//    fastqFileR1.close();
//    fastqFileR2.close();
//
//    // Display the sequences stored in the map (for verification)
//    std::cout << "barcode size = " << barcodes.size() << "\n";
//    std::cout << "cells size = " << cells.size() << "\n";
//
//
//    std::ofstream file(outFile);
//    std::ofstream fileB(bcFile);
//    std::ofstream fileE(errorFile);
//    // 
//    // sourt cells by number of dna seq captured
//    //
//    std::cout << "Sort\n";
//    std::vector<std::pair<std::string,std::map<std::string,int>>> t1;
//    std::vector<std::pair<std::string,std::map<std::string,int>>> t2;
//    for (auto entry:cells) {
//      if (entry.second.size() > 20) {
//        t1.push_back(entry);
//      }
//    }
//    sort(t1.begin(),t1.end(),cmp1);
//    t2 = t1;
//    //
//    //
//    // 
//    int indexCount = 0;
//    int cellCount = t2.size();
//    
//    for (auto entry:t2) {
//      std::string barcode = entry.first;
//      int bcCount = entry.second.size();
//      if (bcCount > 100) {
//        fileB << barcode << "," << bcCount << "\n"; 
//        //file << "# cell\n";
//        for (auto pair:entry.second) {
//	   std::string seqData = pair.first;
//	   file << barcode << "\t" << seqData << "\n";
//        }
//      }
//      if ((indexCount % 1000) == 0) {
//        std::cout << indexCount << " of " << cellCount << "\n";
//      }
//      indexCount += 1;
//    }
//
//    file.close();
//    fileB.close();
//    fileE.close();
//
//    return 0;
//}
