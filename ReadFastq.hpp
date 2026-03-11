//
//  ReadFastq.hpp
//  mbFilter
//
//  Created by mark enstrom on 1/15/25.
//

#ifndef ReadFastq_hpp
#define ReadFastq_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include "Params.hpp"
#include <string>
#include "math.h"
#include "NWAlign.hpp"

//
// 
//
int findOffset(const std::string &seq);

const std::string geCD33  = "CO-819_AMP47";
const std::string geHBG1  = "HbG_AMP1";
const std::string geHBG2  = "HbG_AMP2";


//
// data storage
//
//
//
//
class CELLDATA {
public:
    std::map<std::string, int> barcodes;
    std::map<std::string,std::map<std::string,int>> cells;
    std::map<std::string,std::map<std::string,std::vector<std::string>>> cellEdit;

    void setScaddenBarcodeMap(std::string cellBarcode,std::map<std::string,int> *psMap) {
        for (auto pair:scaddenMap[cellBarcode]) {
            std::cout << "BEFORE " << pair.first << " " << pair.second << "\n";
        }
        scaddenMap[cellBarcode].clear();
        for (auto pair:(*psMap)) {
            scaddenMap[cellBarcode][pair.first] = pair.second;
        }
        for (auto pair:scaddenMap[cellBarcode]) {
            std::cout << "REPAIRED " << pair.first << " " << pair.second << "\n";
        }
    }
};
//
//
//
class ReadFastq {
private:

    std::string fileR1;
    std::string fileR2;
    std::map<std::string, int> barcodes;
    std::map<std::string,std::map<std::string,int>> cells;
    std::map<std::string,std::map<std::string,int>> cellEdit;
public:
    ReadFastq();
    void read(Params params,CELLDATA *pData,std::string r1,std::string r2);

};

#endif /* ReadFastq_hpp */
