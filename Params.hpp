//
//  Params.hpp
//  mbFilter
//
//  Created by mark enstrom on 1/15/25.
//

#ifndef Params_hpp
#define Params_hpp

#include <stdio.h>
#include <stdio.h>
#include <iostream>
#include <vector>

typedef struct _ENTRY {
    std::string ampID;
    std::string chr;
    std::string fwdSeq;
    std::string revSeq;
    std::string amplicon;
    std::string target;
    std::string gene;
} ENTRY,*PENTRY;

std::vector<std::string> splitByTabs(const std::string &str);

class Params {
private:
    std::vector<PENTRY> paramData;

    int getID(const std::string &ampID);
public:
    bool readPrimerFile(std::string filename);
    Params(std::string filename) {
        readPrimerFile(filename);
    }
    ~Params() {}
    size_t items() {return paramData.size();}
    // not validated
    std::string amp(int index) {return paramData[index]->ampID;}
    std::string ref(int index) {return paramData[index]->amplicon;}
    //
    //
    //
    std::pair<std::string,int> identifySequence(std::string seq);
    
    std::string buildSeq(
      const std::string& fwd_seq,
      const std::string& seq2,
      const std::string& ampID);
};


class META {
public:
    bool readMetaFile(std::string filename);
    
    std::vector<std::string> vec_DNA_R1;
    std::vector<std::string> vec_DNA_R2;
    
    std::vector<std::string> vec_PROTEIN_R1;
    std::vector<std::string> vec_PROTEIN_R2;
    
    
    size_t sizeDNA() {return vec_DNA_R1.size();}
    size_t sizePROTEIN() {return vec_PROTEIN_R1.size();}
    
    std::string dna_r1(size_t index) {
        if (index < vec_DNA_R1.size()) {
            return(vec_DNA_R1[index]);
        }
        return("");
    }
    std::string dna_r2(size_t index) {
        if (index < vec_DNA_R2.size()) {
            return(vec_DNA_R2[index]);
        }
        return("");
    }
    std::string protein_r1(size_t index) {
        if (index < vec_PROTEIN_R1.size()) {
            return(vec_PROTEIN_R1[index]);
        }
        return("");
    }
    std::string protein_r2(size_t index) {
        if (index < vec_PROTEIN_R2.size()) {
            return(vec_PROTEIN_R2[index]);
        }
        return("");
    }
    
    
};

#endif /* Params_hpp */
