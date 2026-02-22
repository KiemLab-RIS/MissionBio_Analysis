//
//  Whitelist.cpp
//  mbFilter
//
//  Created by mark enstrom on 1/18/25.
//

#include "Whitelist.hpp"


//---------------------------------------------------------------------
//
// file format
//
// celltag,count
// GGTTGTTC,0
// ATACCTTT,0
// TCGAAAGC,0
// TAGTCTAT,0
//
//
//---------------------------------------------------------------------
Whitelist::Whitelist(const std::string& whitelistFile) {
    std::string line;
    std::ifstream istream(whitelistFile);
    if (!istream.is_open()) {
        std::cerr << "Error opening file " << whitelistFile << std::endl;
    }
    //
    // skip header
    //
    std::getline(istream, line);
    //
    // read
    //
    
    while (std::getline(istream, line)) {
        //
        // read1
        //
        barcodes.push_back(line.substr(0,8));
    }
    std::cout << barcodes.size() << " whitelist barcodes read from file\n";
}
//---------------------------------------------------------------------
//
// try to error-correct scadden barcodes
//
//
//---------------------------------------------------------------------
bool Whitelist::findExact(std::string bc) {
    for (auto s:barcodes) {
        if (bc == s) {
            return true;
        }
    }
    return false;
}

//
// must be same len
//
int Whitelist::stringDiff(const std::string &s1, const std::string &s2) {
  int error = 0;
  for (int index = 0; index < (int)s1.size(); index++) {
    if (s1[index] != s2[index]) error += 1;
  }
  return(error);
}

std::string Whitelist::findClose(std::string bc) {
    for (auto s:barcodes) {
        if (stringDiff(bc,s) == 1) {
            //std::cout << "find match " << s << "\n";
            return s;
        }
    }
    return "";
}
//---------------------------------------------------------------------
//
// try to error-correct scadden barcodes
// build a new map to contain good barcodes and corrected
// new barcodes
//
//------------------------------------------------------------------


void Whitelist::errorCorrect(std::map<std::string, int>& scadMap) {
    if (scadMap.size() <= 1) return;
    // get a list of keys
    std::vector<std::string> keys;
    for (auto pair:scadMap) {
        keys.push_back(pair.first);
    }
    // loop through keys
    for (auto barcode:keys) {
        int count = scadMap[barcode];
        // exact match in whitelist
        if (!findExact(barcode)) {
            std::cout << barcode << " mismatch \n";
            std::string mBarcode = findClose(barcode);
            if (mBarcode != "") {
                std::cout << "repair " << mBarcode << " from " << barcode << "\n";
                // add corrected
                scadMap[mBarcode] += count;
                // erase old
                scadMap.erase(barcode);
            }
        }
    }
    return;
}
//-------------------------------------------------------------
//
// vector compare functions
//
//-------------------------------------------------------------
bool cmpV(std::pair<std::string, int>& a,
         std::pair<std::string, int>& b)
{
    return a.second > b.second;
}
//-------------------------------------------------------------
//
// string comp
//
//-------------------------------------------------------------se
bool differsByOneBase(const std::string& seq1, const std::string& seq2) {
    if (seq1.length() != seq2.length()) return false;

    int diffCount = 0;
    for (size_t i = 0; i < seq1.length(); ++i) {
        if (seq1[i] != seq2[i]) {
            ++diffCount;
            if (diffCount > 1) return false;
        }
    }
    return diffCount == 1;
}

void getSortedKeys(std::vector<std::string>& keys,std::map<std::string, int>& scadMap) {
    
    std::vector<std::pair<std::string,int>> i1;
    
    for (auto entry:scadMap) {
        i1.push_back(entry);
    }
    sort(i1.begin(),i1.end(),cmpV);
    for (const auto& pair : i1) {
        keys.push_back(pair.first);
    }
}

//-------------------------------------------------------------
//
// vector compare functions
//
//-------------------------------------------------------------

void combineSequences(std::map<std::string, int>& dnaMap) {
    bool merged = true;

    // Keep merging until no more sequences can be combined
    while (merged) {
        merged = false;
        std::vector<std::string> keys; // Store the current keys to iterate
        getSortedKeys(keys,dnaMap);
        for (size_t i = 0; i < keys.size(); ++i) {
            for (size_t j = i + 1; j < keys.size(); ++j) {
                if (differsByOneBase(keys[i], keys[j])) {
                    // Merge counts
                    std::cout << "merge " << keys[i] << " " << keys[j] << "\n";
                    dnaMap[keys[i]] += dnaMap[keys[j]];
                    // Remove the second sequence from the map
                    dnaMap.erase(keys[j]);
                    merged = true;
                    break; // Break out of inner loop
                }
            }
            if (merged) break; // Break out of outer loop
        }
    }
}
//-------------------------------------------------------------------
//
// combine barcodes with up to one base pair difference then
// also correct 1 base to whitelist
//
//------------------------------------------------------------------
void Whitelist::ScaddenCorrectToWhite(CELLDATA *pData) {
    std::cout << " check " << pData->scaddenMap.size() << " barcodes\n";
    for (auto& pair:pData->scaddenMap) {
        auto &innerMap = pair.second;
        
        std::string cellBarcode = pair.first;
        std::cout << "-------------correct Cell " << cellBarcode << "---------\n";
        
        //std::map<std::string,int> scadMap = pair.second;
        
        for (auto pr:innerMap) {
            std::cout << "   b.c.   " << pr.first << " " << pr.second << "\n";
        }
        
        combineSequences(innerMap);
        std::cout << "map contains " << innerMap.size() << " entries\n";

        errorCorrect(innerMap);
        
        for (auto pr:pData->scaddenMap[cellBarcode]) {
            std::cout << "   a.c.   " << pr.first << " " << pr.second << "\n";
        }
    }
}

