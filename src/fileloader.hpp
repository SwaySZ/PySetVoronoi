/* 
Copyright 2016 Simon Weis and Philipp Schoenhoefer
Sway Zhao 2016
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
#ifndef FILELOADER_GUARD
#define FILELOADER_GUARD

#include <vector>
#include <string>
#include "pointpattern.hpp"
#include "particleparameterset.hpp"


class fileloader
{
public:
    bool read(std::string filename, std::vector<particleparameterset>& set);
    void readWall(std::string filename, double& xmin,double& xmax,double& ymin,double& ymax,double& zmin,double& zmax);
    void readRawParticle(std::string filename, pointpattern& pp, unsigned int id);//xyz file read
};


#endif
