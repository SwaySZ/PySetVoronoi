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
#ifndef POSTPROCESSING_H_123456_GUARD
#define POSTPROCESSING_H_123456_GUARD

#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <map>
#include <algorithm>

inline void DoPostProcessing(std::string folder)
{
    std::string filepath = folder + "custom.dat";
    std::ifstream infile(filepath);
    
    std::string line;

    // faceshistogram
    std::map<unsigned long, long> histogramFaces;
    std::vector<double> vecVolume;

    while(std::getline(infile, line))
    {
        std::istringstream iss(line);
        unsigned long l = 0;
        unsigned long f = 0;
        double v = 0.0;
        
        if (!(iss >>l >> f >> v))
        {
            break;  //error
        }

        vecVolume.push_back(v);

        if (!histogramFaces.count(f))
        {
            histogramFaces.insert(std::pair<unsigned int,long >(f,0));
        }
        histogramFaces.find(f)->second +=1.0;
    };
    
    // calculate volume histogram


    // the following Code to calculate volume distributions does not seem to work right no, so it stays commented until it will work
    /*double minVol = *std::min_element(vecVolume.begin(), vecVolume.end());
    double maxVol = *std::max_element(vecVolume.begin(), vecVolume.end());
    double count = vecVolume.size();


    double binCount = count/5;
    if (binCount >= 100)
        binCount = 100;

    double dv = (maxVol - minVol)/(binCount + 1);

    std::vector<long> histogramVolume (binCount+1, 0);
    for (double v : vecVolume)
    {
        unsigned int idx = static_cast<unsigned int>((v -minVol)/dv);
        histogramVolume[idx] += 1;
    }
*/
    // save output 
    
    {
        std::ofstream outFaces (folder+"faces.stat");
        outFaces << "#NumberOfFaces Occurrence" << std::endl;

        for (auto it = histogramFaces.begin(); it != histogramFaces.end(); ++it)
        {
            outFaces << it->first << " " << it->second << std::endl;
        }
        outFaces.close();
    }
/*
    {
        std::ofstream outVolumes(folder+"volumes.stat");
        outVolumes << "#volume  occurrence" << std::endl;
        for (unsigned int i = 0; i != histogramVolume.size(); ++i)
        {
            double v = minVol + i * dv;
            double c = histogramVolume[i];

            outVolumes << v << " " << c << std::endl;
        }
    }*/
}

#endif
