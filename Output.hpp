//
//  Output.hpp
//  mbFilter
//
//  Created by mark enstrom on 1/15/25.
//

#ifndef Output_hpp
#define Output_hpp

#include <stdio.h>
#include <stdio.h>
#include <stdio.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include "Params.hpp"
#include "ReadFastq.hpp"
#include "ReadProtein.hpp"



class Output {
public:
    void write(CELLDATA *pData, const std::string &filename);
    void bigTable(CELLDATA *pData,Params &pr,EditMap &edits,const std::string &filename);
    std::string predict1(std::map<std::string,int>& capID);
    std::string predictCellType(std::map<std::string,int>& capID);
};
#endif /* Output_hpp */
