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

#include<fstream>
#include<iostream>

#include<sstream>
#include<string>
#include "fileloader.hpp"

#include "csplitstring.hpp"
#include <stdio.h>
#include <cstdlib>
bool fileloader::read(std::string filename, std::vector<particleparameterset>& set)
{
    std::ifstream infile;
    infile.open(filename.data(),std::ifstream::in);
    if (infile.fail())
    {
        std::cout << "Cannot load file " << filename << std::endl;
        return false;
    }
#pragma GCC diagnostic ignored "-Wwrite-strings"
    cSplitString line("");
    unsigned int linesloaded = 0;
    //std::getline(infile, line);
    while (std::getline(infile, line))
    {
        if(line.find("#")!=std::string::npos) continue; // ignore comment lines

        std::vector<std::string> thisparams = line.split('\t');//changed by Sway

        particleparameterset p;
        for(int i = 0;i<thisparams.size();++i){

        	double d = atof(thisparams.at(i).c_str());
        	p.push_back(d);
        }

        linesloaded++;
        set.push_back(p);

    }
    std::cout << "Lines loaded: " << linesloaded << std::endl << std::endl;

    infile.close();
    return true;
}
void fileloader::readWall(std::string filename, double& xmin,double& xmax,double& ymin,double& ymax,double& zmin,double& zmax)
{
    std::ifstream infile;
    infile.open(filename.data(), std::ifstream::in);
    if (infile.fail())
    {
        std::cout << "Cannot load file " << filename << std::endl;
        return;
    }
#pragma GCC diagnostic ignored "-Wwrite-strings"
    cSplitString line("");
    unsigned int linesloaded = 0;
    //std::getline(infile, line);
    std::cout<<"Warning! The wallfile should include six lines in order with xmin,xmax,ymin,ymax,zmin,zmax. Each line has three components."<<std::endl;
    std::vector< std::vector<double> > wallpos;    
    while (std::getline(infile, line))
    {
        if(line.find("#")!=std::string::npos) continue; // ignore comment lines

        std::vector<std::string> xyzstring = line.split('\t');//
        //std::cout<<xyzstring[0]<<std::endl;
        std::vector<double> xyz;
        //xyz point
        for(int i =0;i<xyzstring.size();i++){
        	double d = atof(xyzstring.at(i).c_str());
        	xyz.push_back(d);
        }

        //check the legality of xyz data
        if(3 > xyz.size())
        {
            std::cout << "Warning!! The data at Line " << linesloaded <<" in the dataset has wrong format. I regected it for robust running."<< std::endl << std::endl;
        }else
        {
            wallpos.push_back(xyz);
        }
        linesloaded++;

    }
    //check if the wallpos has all data
    if (wallpos.size()!=6){
        std::cout<<"Error: the wallpos has not six-line data."<<wallpos.size()<<std::endl;
    }else{
        //get container size in order
        xmin = wallpos[0][0];
        xmax = wallpos[1][0];
        ymin = wallpos[2][1];
        ymax = wallpos[3][1];
        zmin = wallpos[4][2];
        zmax = wallpos[5][2];
    }
    std::cout << "Wall position was imported successfully!" << std::endl << std::endl;

    infile.close();
}
void fileloader::readRawParticle(std::string filename, pointpattern& pp, unsigned int particleid)
{
    std::ifstream infile;
    infile.open(filename.data(),std::ifstream::in);
    if (infile.fail())
    {
        std::cout << "Cannot load file " << filename << std::endl;
        return;
    }
#pragma GCC diagnostic ignored "-Wwrite-strings"
    cSplitString line("");
    unsigned int linesloaded = 0;
    std::getline(infile, line);
    while (std::getline(infile, line))
    {
        if(line.find("#")!=std::string::npos) continue; // ignore comment lines

        std::vector<std::string> xyzstring = line.split('\t');//
        std::vector<double> xyz;
        //xyz point
        for (int i =0;i<xyzstring.size();i++){
        	double d = atof(xyzstring.at(i).c_str());
        	xyz.push_back(d);
        }

        //check the legality of xyz data
        if(3 > xyz.size())
        {
            std::cout << "Warning!! The data at Line " << linesloaded <<" in the dataset has wrong format. I regected it for robust running."<< std::endl << std::endl;
        }else
        {
            pp.addpoint(particleid,xyz[0],xyz[1],xyz[2]);
        }
        linesloaded++;
       

    }
    std::cout << "Lines loaded: " << linesloaded << std::endl << std::endl;

    infile.close();
}
