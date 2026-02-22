//
//  main.cpp
//  mbFilter
//
//  Created by mark enstrom on 1/15/25.
//

#include <iostream>
#include "NWAlign.hpp"
#include "Params.hpp"
#include "ReadFastq.hpp"
#include "ReadProtein.hpp"
#include "Output.hpp"
#include "Whitelist.hpp"

int main(int argc, const char * argv[]) {
    // insert code here...
    if (argc != 6) {
        std::cout << "usage: filter outfile param Metafile whitelist\n";
        return(-1);
    }
    std::string paramName = argv[1];
    std::string outFile = argv[2];
    std::string metaFile = argv[3];
    std::string whiteListFile = argv[4];
    std::string editFile = argv[5];
    
    Whitelist whitelist(whiteListFile);
    EditMap edits = read_gene_csv(editFile); 
    
    META meta;
    if (!meta.readMetaFile(metaFile)) {
        std::cout << "error reading metafile\n";
        return(-1);
    }
    //
    // set up params
    //
    Params params(paramName);
    //
    // read fastq files
    //
    CELLDATA *pData = new CELLDATA();
    
    ReadFastq rfq;
    
    ReadProtein rpr;
    //
    // read protein files
    //
    for (size_t i = 0; i < meta.sizePROTEIN(); i++) {
        auto r1 = meta.protein_r1(i);
        auto r2 = meta.protein_r2(i);
        rpr.read(pData, r1, r2);
    }
    //
    // read dna files
    //
    for (int i = 0; i < (int)meta.sizeDNA(); i++) {
        auto r1 = meta.dna_r1(i);
        auto r2 = meta.dna_r2(i);
        rfq.read(params,pData,r1,r2);
    }
    //
    // whitelist correct
    //
    whitelist.ScaddenCorrectToWhite(pData);
    //
    // classify external proteins
    //
    
    //
    // output
    //
    Output out;
    out.write(pData, outFile + "_cell.tsv");
    out.bigTable(pData,params,edits,outFile + "_table.tsv");
    return 0;
}
