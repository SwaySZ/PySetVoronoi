/*
 * =====================================================================================
 *
 *       Filename:  CellFactory.hpp
 *
 *    Description:  a wrapper of Class CellMachine exposed in Python.
 *
 *        Version:  1.0
 *        Created:  12/22/2018
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Shiwei Zhao (zhswee@gmail.com)
 *   Organization: The Hong Kong University of Science and Technology
 *                 South China University of Technology
 *
 * =====================================================================================
 */
#ifndef __CELL_FACTORY__
#define __CELL_FACTORY__
#include "Superquadrics.hpp"
#include "CellMachine.hpp"
#include "pointpattern.hpp"
#include "duplicationremover.hpp"
#include "particleparameterset.hpp"
#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>

#define CF_OPENMP//using openMP, to do
#define CF_DEBUG_no

class CellFactory{
private:
  std::string in_folder;
  std::string out_folder;
  std::string posFile;//position file of particles
  std::vector<particleAttr> parAttrlist;//store particle attributes
  double parShrink;//shrink a particle by a small gap (parShrink) during generating its point cloud.
  //double w_slices;//slicing
  //double h_slices;
  bool cellVTK;//output cells' vtk files
  //configuration
  double searchRadius;//neighbour lists would be built within a sphere with a radius of seachRadius
public:
  CellFactory();
  ~CellFactory(){};
  void genPointClouds(double w_slices,double h_slices);
  void neighborSearch(void);
  void processing(void);
  void autoWorkFlow(void);
  bool pointCloud_Superquadric(unsigned int id, std::string outfile, double scaledist,std::vector<double> &set, int w_slices, int h_slices, double& area, double& volume, particleAttr& pattr);
  //
  void set_infolder(std::string _infolder){in_folder = _infolder;}
  std::string get_infolder(){return in_folder;}
  void set_outfolder(std::string _outfolder){out_folder = _outfolder;}
  std::string get_outfolder(){return out_folder;}
  void set_searchRadius(double sr){searchRadius = sr;}
  double get_searchRadius(){return searchRadius;}
  void set_posFile(std::string pf){posFile = pf;}
  std::string get_posFile(){return posFile;}
  void set_cellVTK(bool cv){cellVTK=cv;}
  bool get_cellVTK(){return cellVTK;}
  void set_parShrink(double ps){parShrink = ps;}
  double get_parShrink(){return parShrink;}

};
#endif
