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
#include <iostream>
#include "pointpattern.hpp"


void pointpattern::addpointForCell(double x, double y, double z, int l, long cf, long cC)
{
    point p (x,y,z,l, cf,cC);
    points.push_back(p);
}

void pointpattern::addpoint(int l, double x, double y, double z)
{
    //std::cout << "adding point " << l << " " << x << " " << y << " " << z << std::endl;
    point p (x,y,z,l);
    points.push_back(p);
}

void pointpattern::print ()
{
    std::cout << "number of points " << points.size() << std::endl;

}

void pointpattern::removeduplicates (double epsilon)
{
    if(points.empty())
    {
        return;
    }
    std::vector<point> newpoints;

    for(unsigned int i = 0; i != points.size(); ++i)
    {
        if(i%1000==0)
        {
            //std::cout << i << " / " <<  points.size() << "\n";
        }
        bool addthis =true;
        point p1 = points[i];
        point p2(0,0,0,-1,-1);
        for(unsigned int j = i+1; j != points.size(); ++j)
        {
            if (j >= points.size()) break;
            p2 = points[j];
            //if(p2.cellID == -1 || p2.l == -1) continue;
            if (p2.l == 0) std::cout << "particle with label 0 detected" << std::endl;
            if (checkdistancecloserthan(p1, p2, epsilon))
            {
                if (p1.cellID == p2.cellID && p1.l != -1)
                {
                    //std::cout << "removing point " << p1.l <<" for cell " << p1.cellID << std::endl;
                    //std::cout << " point "<<   p1.l << " and point " << p2.l << " too close together" << "\n";
                    addthis = false;
                    break;
                }
            }
        }
        if (addthis)
        {
            unsigned int l = p1.l;
            indexShift[l] = -1;
            //std::cout << "adding" << std::endl;
            newpoints.push_back(p1);
        }
        else
        {
            unsigned int l1 = p1.l;
            unsigned int l2 = p2.l;
            indexShift[l1] = l2;
        }
    }
    points = newpoints;
}

void pointpattern::removeduplicates (double epsilon, pointpattern& p)
{
    if(points.empty()) return;
    std::vector<point> newpoints;
    for(unsigned int i = 0; i != points.size(); ++i)
    {
        if(i%1000==0)
        {
            //std::cout << i << " / " <<  points.size() << "\n";
        }
        bool addthis =true;
        point p1 = points[i];
        point p2(0,0,0,-1, -1);
        for(unsigned int j = 0; j != p.points.size(); ++j)
        {
            p2 = p.points[j];
            //if(p2.cellID == -1 || p2.l == -1) continue;
            if (p2.l == 0) std::cout << "particle with label 0 detected" << std::endl;
            if (checkdistancecloserthan(p1, p2, epsilon))
            {
                if (p1.cellID == p2.cellID && p1.l != -1)
                {
                    //std::cout << "removing point " << p1.l <<" for cell " << p1.cellID << std::endl;
                    //std::cout << " point "<<   p1.l << " and point " << p2.l << " too close together" << "\n";
                    addthis = false;
                    break;
                }
                //addthis = false;
                //break;
            }
        }
        if (addthis)
        {
            //std::cout << "adding" << std::endl;
            newpoints.push_back(p1);
            unsigned int l = p1.l;
            indexShift[l] = -1;
        }
        else
        {
            //std::cout << "removing point" << std::endl;
            unsigned int l1 = p1.l;
            unsigned int l2 = p2.l;
            indexShift[l1] = l2;
        }
    }
    points = newpoints;
}
