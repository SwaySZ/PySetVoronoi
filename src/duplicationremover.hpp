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
#ifndef DUPLICATIONREMOVER_H_1234567
#define DUPLICATIONREMOVER_H_1234567

#include <vector>
#include "pointpattern.hpp"

class duplicationremover
{
public:
    duplicationremover(unsigned int _x, unsigned int _y, unsigned int _z) : x(_x), y(_y), z(_z)
    {
        list.resize(x*y*z);
    };

    // read points from a pointpattern
    void addPoints ( pointpattern& porig, bool useCellIDs = false)
    {
        std::cout << "\tadding N= " <<porig.points.size() <<  " points for duplication check" << std::endl;
        for (  auto it = porig.points.begin();
                it != porig.points.end();
                ++it)
        {
            if (useCellIDs)
            {
                //std::cout << "adding point with label : " << it->l << " in cellID " << it->cellID << "at (" << it->x << " " << it->y << " " << it->z << std::endl;
                addpoint(it->x, it->y, it->z, it->l, it->cellID);
            }
            else
            {
                //std::cout << "adding point with label : " << it->l << " at " << it->x << " " << it->y << " " << it->z << std::endl;
                addpoint(it->x, it->y, it->z, it->l);
            }
        }
    }

    // TODO periodic Boundaries
    void setboundaries (double xmi, double xma, double ymi, double yma, double zmi, double zma)
    {
        xmin = xmi;
        ymin = ymi;
        zmin = zmi;
        xmax = xma;
        ymax = yma;
        zmax = zma;
    }

    void addpoint(double dx, double dy, double dz, unsigned long long l, long cellID = -1)
    {
        if (dx > xmax || dx < xmin)
        {
            std::cerr << "Warning: point with label: " << l << " out of X boundaries, adding in outmost subcell" << std::endl;
        }
        if (dy > ymax || dy < ymin)
        {
            std::cerr << "Warning: point with label: " << l << " out of Y boundaries, adding in outmost subcell" << std::endl;
        }
        if (dz > zmax || dz < zmin)
        {
            std::cerr << "Warning: point with label: " << l << " out of Z boundaries, adding in outmost subcell" << std::endl;
        }

        //std::cout << "adding particle " << dx << " " << dy << " "  << dz << " " << zmax  << std::endl;
        double sx = (xmax-xmin)/x;
        double sy = (ymax-ymin)/y;
        double sz = (zmax-zmin)/z;

        int cx = static_cast<int>((dx - xmin)/sx);
        int cy = static_cast<int>((dy - ymin)/sy);
        int cz = static_cast<int>((dz - zmin)/sz);


        // if strange index has been calculated, shift it to the valid range
        while (cx < 0 ) cx++;
        while (cy < 0 ) cy++;
        while (cz < 0 ) cz++;

        while (static_cast<unsigned int>(cx) >= x ) cx--;
        while (static_cast<unsigned int>(cy) >= y ) cy--;
        while (static_cast<unsigned int>(cz) >= z ) cz--;

        unsigned int index = getindex(cx,cy,cz);
        if (cellID == -1)
        {
            list[index].addpoint(l, dx,dy,dz);
        }
        else
        {
            list[index].addpointForCell(dx,dy,dz,l,cellID,cellID);
        }

    }

    // store all points to passed parameter p
    void getallPoints ( pointpattern& p)
    {
        p.points.clear();
        unsigned int i = 0;
        // loop over all subcells
        for (
            auto it = list.begin();
            it != list.end();
            ++it)
        {
            // ... in each subcell, loop over all points
            for (
                auto itp = (*it).points.begin();
                itp != (*it).points.end();
                ++itp)
            {
                // add them to p
                p.addpoint(itp->l, itp->x, itp->y, itp->z);
                i++;
            }
        }
        //std::cout << "found " << i << std::endl;

    };


    // remove any duplicated points
    void removeduplicates(double epsilon)
    {
        std::cout << "\tremoving duplicates within one box" << std::endl;
        for (
            auto it = list.begin();
            it != list.end();
            ++it)
        {
            (*it).removeduplicates(epsilon);
        }

        std::cout << "\tremoving duplicates between neighboring boxes" << std::endl;
        for (unsigned int cx = 0; cx != x; ++cx)
            for (unsigned int cy = 0; cy != y; ++cy)
                for (unsigned int cz = 0; cz != z; ++cz)
                {
                    for (int i = -1; i != 1; ++i)
                        for (int j = -1; j != 1; ++j)
                            for (int k = -1; k != 1; ++k)
                            {
                                if (i == 0 && j == 0 && k == 0) continue;
                                int ox = cx + i;
                                int oy = cy + j;
                                int oz = cz + k;
                                if (ox < 0 ||static_cast<unsigned int >(ox) >= x ) continue;
                                if (oy < 0 ||static_cast<unsigned int >(oy) >= y ) continue;
                                if (oz < 0 ||static_cast<unsigned int >(oz) >= z ) continue;

                                int index1 = getindex(cx,cy,cz);
                                int index2 = getindex(ox,oy,oz);
                                list[index1].removeduplicates(epsilon, list[index2]);
                            }

                }

        for (
            auto it = list.begin();
            it != list.end();
            ++it)
        {

            for (auto it2 = (*it).indexShift.begin();
                    it2 != (*it).indexShift.end();
                    ++it2)
            {
                indexShift[(*it2).first] = (*it2).second;
            }
        }
    };

    unsigned int getindex (unsigned int cx, unsigned int cy, unsigned int cz)
    {
        return cx + x * cy + x*y*cz;
    }

    std::map<unsigned int, long> indexShift;


private:
    std::vector<pointpattern> list;

    unsigned int x;
    unsigned int y;
    unsigned int z;

    double xmin, xmax, ymin, ymax, zmin, zmax;
};

#endif
