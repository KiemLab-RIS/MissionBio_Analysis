//
//  ReadProtein.cpp
//  mbFilter
//
//  Created by mark enstrom on 1/16/25.
//

#include "ReadProtein.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include<algorithm>
#include <stdexcept>



//
// must be same len
//
int stringDiff(const std::string &s1, const std::string &s2) {
  int error = 0;
  for (int index = 0; index < (int)s1.size(); index++) {
    if (s1[index] != s2[index]) error += 1;
  }
  return(error);
}
//
//
//
std::string findAntibody(std::string &seq) {
    // exact match
    for (int index = 0; index < (int)abBarcodeVec.size(); index++) {
        if (seq == abBarcodeVec[index]) {
            return(abNameVec[index]);
        }
    }
    // error correct 1
    int best = (int)seq.size();
    int bestIndex = 0;
    for (int index = 0; index < (int)abBarcodeVec.size(); index++) {
        int mismatch = stringDiff(seq,abBarcodeVec[index]);
        if (mismatch < best) {
            best = mismatch;
            bestIndex = index;
        }
    }
    if (best <= 2) {
        return(abNameVec[bestIndex]);
    }
    //std::cout << "-----ba match error----" << best << "\n";
    //std::cout << seq << "\n";
    //std::cout << abBarcodeVec[bestIndex] << "\n";
    //std::cout << abNameVec[bestIndex] << "\n";
    
    return("UNMATCH_" + seq);
}
//--------------------------------------------------------------------
//
// Function to find all matches allowing one(two) mismatch
//
//
//--------------------------------------------------------------------

size_t findPatternWithMismatch(const std::string& text, const std::string& pattern) {
    size_t patternLength = pattern.length();
    size_t textLength = text.length();
    // Sliding window approach
    for (size_t i = 0; i <= textLength - patternLength; ++i) {
        int mismatchCount = 0;
        // Compare the current window with the pattern
        for (size_t j = 0; j < patternLength; ++j) {
            if (text[i + j] != pattern[j]) {
                ++mismatchCount;
                if (mismatchCount > 2) break; // Stop early if more than one mismatch
            }
        }

        // If at most one mismatch, store the index
        if (mismatchCount <= 2) {
            return (size_t)i;
        }
    }
    return std::string::npos;
}


//
// main
//
void ReadProtein::read(CELLDATA *pData,std::string fileR1,std::string fileR2) {
    std::cout << "Read Protein from " << fileR1 << "\n";
    //
    // Open the FASTQ file
    //
    std::ifstream fastqFileR1(fileR1);
    if (!fastqFileR1.is_open()) {
        std::cerr << "Error opening file " << fileR1 << std::endl;
        return;
    }
    // Open the FASTQ file
    std::ifstream fastqFileR2(fileR2);
    if (!fastqFileR2.is_open()) {
        std::cerr << "Error opening file " << fileR2 << std::endl;
        return;
    }
    std::string lineR1;
    std::string lineR2;
    std::string identifier;
    std::string sequenceR1;
    std::string sequenceR2;
    int count = 0;
    int badOffset = 0;
    int badLoCap = 0;
    int badAB = 0;
    try {
        while (1) {
            //
            // read1
            //
            std::string barcode = "";
            std::getline(fastqFileR1, lineR1);
            if (lineR1 == "") {
                break;
            }
            std::getline(fastqFileR1, sequenceR1);  // Line 2: Sequence
            std::getline(fastqFileR1, lineR1);      // Line 3: Plus line (+)
            std::getline(fastqFileR1, lineR1);      // Line 4: Quality (ignored)            
            //
            // read2
            //
            std::getline(fastqFileR2, lineR2);      // Line 3: Plus line (+)
            std::getline(fastqFileR2, sequenceR2);  // Line 2: Sequence
            std::getline(fastqFileR2, lineR2);      // Line 3: Plus line (+)
            std::getline(fastqFileR2, lineR2);      // Line 4: Quality (ignored)
            
            count += 1;
            if ((count % 1000000) == 0) {
                std::cout << count << "\n";
            }
            
            int offset = findOffset(sequenceR1);
            
            if (offset == -1) {
                badOffset +=1;
                continue;
            }
            
            if ((int)sequenceR1.size() < 9 + offset + 9) {
                std::cout << "erorr in R1 "  << sequenceR1.size() << "\n";
                break;
            }

            barcode = sequenceR1.substr(0,9) + "_" + sequenceR1.substr(9 + offset, 9);
            pData->barcodes[barcode] += 1;
            //
            // quick location of antibody cap sequence
            //
            size_t loCap = sequenceR2.find(antibodyCapSeq);
            
            if (loCap == std::string::npos) {
                //
                // allow one mismatch error
                //
                loCap = findPatternWithMismatch(sequenceR2,antibodyCapSeq);
                if (loCap == std::string::npos) {
                    //std::cout << "R2 " << sequenceR2 << " " << sequenceR2.length() << "\n";
                    //std::string sp(38, ' ');
                    //std::cout << "CapSeq " << sp << antibodyCapSeq << "\n";
                    badLoCap++;
                    continue;
                }
            }
            
            if ((loCap >= 15) && (barcode != "")) {
                if (sequenceR2.size() < loCap) {
                    std::cout << "error in seqr2 " << loCap << " " << sequenceR2.size() << "\n";
                    break;
                }
                std::string ab = sequenceR2.substr(loCap-15,15);
                std::string abMap = findAntibody(ab);
                pData->cells[barcode][abMap] += 1;
                if (abMap[0] == 'U') {
                    badAB ++;
                }
            }
        }
    } catch (int x) {
        std::cout << "Exception Caught \n";
    }
    std::cout << "Done reading " << count << "\n";
    std::cout << "bad Offset = " << badOffset << "\n";
    std::cout << "bad loCap  = " << badLoCap << "\n";
    std::cout << "bad AB     = " << badAB << "\n";
    std::cout << "barcodes   = " << pData->barcodes.size() << "\n";
    std::cout << "cells      = " << pData->cells.size() << "\n";
    
    fastqFileR1.close();
    fastqFileR2.close();
}
//    // Display the sequences stored in the map (for verification)
//    std::cout << "map size = " << barcodes.size() << "\n";
//    // for (const auto &entry : barcodes) {
//    //      if (entry.second > 10) {
//    //         std::cout << "Identifier: " << entry.first << "\nSequence: " << entry.second << "\n" << std::endl;
//    //    }
//    //}
//    //
//    //
//    //  output
//    //
//    //
//    //
//    // summarize AB
//    //
//    
//    std::vector<std::pair<std::string,std::vector<std::string>>> cellList;
//    
//    for (const auto &entry : cells) {
//        cellList.push_back(entry);
//    }
//    
//    std::sort(cellList.begin(),cellList.end(),cmp1);
//    
//    
//    std::ofstream fileB(outBC);
//    
//    for (const auto &entry : cellList) {
//        if (entry.second.size() > 1) {
//            fileB << entry.first << "," << entry.second.size() << "\n";
//        }
//    }
//    fileB.close();
//    
//    
//    std::ofstream file(outFile);
//    
//    for (const auto &entry : cellList) {
//        //
//        // summarize matching AB sequence in a map
//        // then convert map to vector
//        // sort vector for AB counts
//        //
//        std::map<std::string,int> abSum;
//        int totalAB = 0;
//        for (const auto ab:entry.second) {
//            abSum[ab] += 1;
//            totalAB += 1;
//        }
//        //
//        // at least 10 antibody capture
//        //
//        if (totalAB < 10) {
//            continue;
//        }
//        
//        std::vector<std::pair<std::string,int>> fred;
//        
//        for (auto pair:abSum) {
//            fred.push_back(pair);
//        }
//        
//        std::sort(fred.begin(),fred.end(),cmp);
//        file << "#Cell\n";
//        for (const auto abEntry:fred) {
//            if (abEntry.second > 1) {
//                file << entry.first << "\t" << abEntry.first << "\t" << abEntry.second << "\n";
//            }
//        }
//    }
//    
//    file.close();
//    
//    return 0;
//}
