//
//  Whitelist.hpp
//  mbFilter
//
//  Created by mark enstrom on 1/18/25.
//

#ifndef Whitelist_hpp
#define Whitelist_hpp


#include <stdio.h>
#include <stdio.h>
#include <stdio.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include "ReadFastq.hpp"


class Whitelist {
private:
    std::vector<std::string> barcodes;
    void errorCorrect(std::map<std::string, int>& scadMap);
    bool findExact(std::string bc);
    std::string findClose(std::string bc);
    int stringDiff(const std::string &s1, const std::string &s2);
public:
    void ScaddenCorrectToWhite(CELLDATA *pData);
    Whitelist(const std::string& whitelistFile);

};







#endif /* Whitelist_hpp */
