//
//  ReadProtein.hpp
//  mbFilter
//
//  Created by mark enstrom on 1/16/25.
//

#ifndef ReadProtein_hpp
#define ReadProtein_hpp

#include <stdio.h>
#include <stdio.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include "ReadFastq.hpp"
#include "Params.hpp"

//
// fastq constants
//
const std::string antibodyCapSeq = "GAGCCGATCTAGTATCTCAGT";
//
// antibody barcodes
//
const std::string CD8    = "GCGCAACTTGATGAT";
const std::string CD11b  = "GACAAGTGATCTGCA";
const std::string CD14   = "TCTCAGACCTCCGTA";
const std::string CD16   = "AAGTTCACTCTTTGC";
const std::string CD20   = "TTCTGGGTCCCTAGA";
const std::string CD90   = "GCATTGTACGATTCA";
const std::string CD117  = "AGACTAATAGCTGAC";
const std::string CD34   = "GATTCCACTCTCCGA";
const std::string CD4    = "AACGCTTGAGATGCC";
const std::string Custom = "CGACGTAACGTAACT";

const std::vector<std::string> abBarcodeVec = {CD8,CD11b,CD14,CD16,CD20,CD90,CD117,CD34,CD4,Custom};
const std::vector<std::string> abNameVec = {"CD8","CD11b","CD14","CD16","CD20","CD90","CD117","CD34","CD4","Custom","pCellType"};

class ReadProtein {
public:
    void read(CELLDATA *pData,std::string fileR1,std::string fileR2);
};




#endif /* ReadProtein_hpp */

