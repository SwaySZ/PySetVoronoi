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
#include "config.hpp"
#include "Superquadrics.hpp"
#include "CellMachine.hpp"
#include "pointpattern.hpp"
#include "duplicationremover.hpp"
#include "particleparameterset.hpp"
#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

class CellFactory{
private:
  std::string in_folder;
  std::string out_folder;
  std::string posFile;//position file of particles
  std::string wallFile;//wall position file
  std::vector<particleAttr> parAttrlist;//store particle attributes
  double parShrink;//shrink a particle by a small gap (parShrink) during generating its point cloud.
  double boxScale;
  double scale;
  //double w_slices;//slicing
  //double h_slices;
  bool cellVTK;//output cells' vtk files
  bool cellPOV;
  std::vector<int> visualized_ids;//id list for particle cells needed visualizing via vtk or pov.
  //configuration
  double searchRadius;//neighbour lists would be built within a sphere with a radius of seachRadius
  unsigned int threadNum;
public:
  CellFactory();
  ~CellFactory(){};
  void genPointClouds(double w_slices,double h_slices);
  void neighborSearch(void);
  void processing(void);
  void processingOne(unsigned int num);
  void autoWorkFlow(void);
  bool checkCreateFolder(std::string target);
  bool isInVisualIds(int id);
  bool pointCloud_Superquadric(unsigned int id, std::string outfile, double scaledist,std::vector<double> &set, int w_slices, int h_slices, double& area, double& volume, particleAttr& pattr);
  //
  void set_infolder(std::string _infolder){in_folder = _infolder;std::cout<<"setting infolder"<<std::endl;checkCreateFolder(in_folder);parAttrlist.clear();}
  std::string get_infolder(){return in_folder;}
  void set_outfolder(std::string _outfolder){out_folder = _outfolder;checkCreateFolder(out_folder);checkCreateFolder(out_folder+"/tmp/");parAttrlist.clear();}
  std::string get_outfolder(){return out_folder;}
  void set_searchRadius(double sr){searchRadius = sr;}
  double get_searchRadius(){return searchRadius;}
  void set_posFile(std::string pf){posFile = pf;parAttrlist.clear();}
  std::string get_posFile(){return posFile;}
  void set_wallFile(std::string pf){wallFile = pf;parAttrlist.clear();}
  std::string get_wallFile(){return wallFile;}
  void set_cellVTK(bool cv){cellVTK=cv;}
  bool get_cellVTK(){return cellVTK;}
  void set_cellPOV(bool cv){cellPOV=cv;}
  bool get_cellPOV(){return cellPOV;}
  void set_visual_ids(std::vector<int> cv){visualized_ids=cv;}
  std::vector<int> get_visual_ids(){return visualized_ids;}
  void set_parShrink(double ps){parShrink = ps;}
  double get_parShrink(){return parShrink;}
  double get_scale(){return scale;}
  void set_scale(double sc){scale = sc;}
  double get_boxScale(){return boxScale;}
  void set_boxScale(double sc){boxScale = sc;}
  unsigned int get_threadNum(){return threadNum;}
  void set_threadNum(int tn){threadNum = std::min(abs(tn),omp_get_max_threads());}
};
#endif
