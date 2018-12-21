/* 
Copyright 2016 Simon Weis and Philipp Schoenhoefer

This file is part of Pomelo.

Pomelo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Pomelo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Pomelo.  If not, see <http://www.gnu.org/licenses/>.

The development of Pomelo took place at the Friedrich-Alexander University of Erlangen and was funded by the German Research Foundation (DFG) Forschergruppe FOR1548 "Geometry and Physics of Spatial Random Systems" (GPSRS). 
*/
#ifndef CSPLITSTRING_H_INCLUDED
#define CSPLITSTRING_H_INCLUDED

#include <string>
#include <vector>
#include <iostream>

class cSplitString : public std::string
{
    std::vector<std::string> flds;
public:
    cSplitString(char* s) : std::string(s) { };
    std::vector<std::string>& split(char delim, int rep=0)
    {
        if (!flds.empty()) flds.clear();  // empty vector if necessary
        std::string work = data();
        std::string buf = "";
        unsigned int i = 0;
        while (i < work.length())
        {
            if (work[i] != delim)
                buf += work[i];
            else if (rep == 1)
            {
                flds.push_back(buf);
                buf = "";
            }
            else if (buf.length() > 0)
            {
                flds.push_back(buf);
                buf = "";
            }
            i++;
        }
        if (!buf.empty())
            flds.push_back(buf);
        return flds;
    }
};


#endif // CSPLITSTRING_H_INCLUDED
