#ifndef __CELL_MACHINE__
#define __CELL_MACHINE__
/*
 * =====================================================================================
 *
 *       Filename:  CellMachine.hpp
 *
 *    Description:  CellMachine is designed for handling the Voronoi cell of a single
 *                  particle surrounded by certain particles.
 *
 *        Version:  1.0
 *        Created:  12/22/2018
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  zhswee (zhswee@gmail.com)
 *   Organization: South China University of Technology
 *
 * =====================================================================================
 */
 #include <iostream>
 #include <fstream>
 #include <map>
 #include "config.hpp"
 #include <limits>
 #include <cmath>
 #include <sys/stat.h>
 #include <voro++.hh>
 #include "fileloader.hpp"
 #include "pointpattern.hpp"
 #include "duplicationremover.hpp"
 #include "polywriter.hpp"
 #include "postprocessing.hpp"
 #include "csplitstring.hpp"
 #include "parAttr.hpp"

class CellMachine{//this class is designed to handle the computation of a single particle's cell.
private:
  //configuration
	std::string in_folder;//folder for input data
  std::string out_folder;//folder for output data
  std::string wallFile;//wall position file; just for box wall
	int nx, ny, nz;
  double xmin,xmax,ymin,ymax,zmin,zmax;//box size
  double wall_xmin,wall_xmax,wall_ymin,wall_ymax,wall_zmin,wall_zmax;//the global wall size
  bool xpbc, ypbc,zpbc;
	bool savereduced, removeduplicate, withboundary, savevtk, cellVTK, savepov, savepoly;
	voro::pre_container *pcon;//container
  unsigned long long pid;//point id
  //std::map < unsigned long long, int> labelidmap;//no need to use map, just vector is ok
  std::vector<int> labelidmap;
  double scale;//scale the data to avoid numerical errors when processing voro++. Voro++'s nplnae throws bugs here.
  double boxScale;//scale box to consider more points for tessellation
  unsigned int blockMem;//memory for each block to be allocated by dufault
  //pointpattern *pp;
  //particleAttr pAttr;
  double rr;//resolution used to clip facets of a cell
  double delta;//shrink the box to get gaps between points and box-walls
  //cell properties
  int cid;//cell id, also particle id
  double cellVolume;
  double cellSurfaceArea;
  double cellNormalTensor[6];//six components of normalTensor(1-6)
  double cellNormalAreaTensor[6];//six components of normalAreaTensor(1-6)
public:
	CellMachine(std::string input_folder,std::string output_folder);
  CellMachine();
	~CellMachine(){
    if(pcon){
      //pcon->clear();
      //delete pcon;
    }
  };
	void pushPoints(particleAttr& pAttr);//push points to the box for voronoi tessellation
	void writeGlobal();//write data to files with info from all particles
	void writeLocal(polywriter *pw);//write data to files with info from just this particle
	void processing();
  void getFaceVerticesOfFace( std::vector<int>& f, unsigned int k, std::vector<unsigned int>& vertexList);
  void readWall(std::string filename);//read the boundar of the global wall
  void readParticle(std::string filename, bool flag, int particleId);//xyz file read
  bool isInBox(double x, double y, double z);//test if a point is inside the box
  void initial();//initialize a cell machine
  void reset();//reset the cell machine for another cell computation if needed.
  //
  double get_scale(){return scale;}
  void set_scale(double sc){scale = sc;}
  double get_boxScale(){return boxScale;}
  void set_boxScale(double sc){boxScale = sc;}
  void set_cellVTK(bool cv){cellVTK=cv;}
  bool get_cellVTK(){return cellVTK;}
  void set_wallFile(std::string pf){wallFile = pf;}
  std::string get_wallFile(){return wallFile;}
};
#endif
