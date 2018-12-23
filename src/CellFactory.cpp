/*
 * =====================================================================================
 *
 *       Filename:  CellFactory.cpp
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
#include<boost/python.hpp>
#include "Superquadrics.hpp"
#include "CellMachine.hpp"
#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>

#define CF_OPENMP//using openMP, to do
#define CF_DEBUG_no

namespace py = boost::python;
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
CellFactory::CellFactory(){
  searchRadius = 4.0;
  parShrink = 0.01e-3;
  //w_slices = 20;
  //h_slices = 20;
  cellVTK = false;
}
void CellFactory::genPointClouds(double w_slices=20,double h_slices = 20){
  std::vector<particleparameterset> setlist;
//judge whether the file can be open or not
  fileloader loader;
  bool flag = loader.read(posFile,setlist);
  if (!flag){std::cerr<<"Error: the file of particle positons has not been loaded correctly!"<<std::endl;}//if the current posfile can not be found, then skipt it.
  unsigned int ids = 0;
  //preprocessing data
  std::ofstream fp1;
  std::string path1 = out_folder + "/parproperties.txt";
  fp1.open(path1.c_str(),std::ios::out);
  fp1 << "#particleID particleVolume particleSurfaceArea" << std::endl;
  for(auto it = setlist.begin(); it != setlist.end(); ++it )
  {
    // read one particle from the position file
    particleparameterset set = (*it);
    double area=0,volume=0;
    std::string outfile = in_folder + "/"+std::to_string(ids)+".dat";//store point clouds of particles
    particleAttr pa;
    pointCloud_Superquadric(ids, outfile, parShrink,set.parameter,w_slices,h_slices,area,volume,pa);
    parAttrlist.push_back(pa);
    fp1 << ids << "\t"<<volume<<"\t" << area << std::endl;
    ids ++;
  }
  std::cout << "point cloud created!" << std::endl;
}

void CellFactory::neighborSearch(void){
  for(int i = 0 ;i< parAttrlist.size();i++){//mayby using pointer is faster
//		particleAttr p1;
    double dist = 0;
    for(int j = 0 ; j < parAttrlist.size();j++){
      //FIXME:using AABB to speed up?
      if(j == i) continue;//
      particleAttr p2  = parAttrlist.at(j);
      dist = sqrt(pow((parAttrlist.at(i).centerx-p2.centerx),2)+pow((parAttrlist.at(i).centery-p2.centery),2)+pow((parAttrlist.at(i).centerz-p2.centerz),2));
      if(parAttrlist.at(i).radius * searchRadius > dist ){
        //p2 close to p1
        parAttrlist.at(i).surroundedID.push_back(p2.ID); // find all particle close to p1 not including current particle itself;
      }
    }
  }
}

void CellFactory::processing(void){
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  //here we can use openMP (or MPI) for parallel computation
  for(int i = 0 ;i< parAttrlist.size();i++){
  //#endif
    CellMachine CM = CellMachine(in_folder,out_folder);
    CM.set_cellVTK(cellVTK);
    //loading point clouds from local files
    CM.pushPoints(parAttrlist[i]);
    //comupute cells
    CM.processing();
  }
  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count();
  std::cout<<"Time ellapsed is:"<<duration<<" seconds"<<std::endl;
}

void CellFactory::autoWorkFlow(void){
  #ifdef CF_DEBUG
  std::cout << "DEBUG: accessfile() called by process " << ::getpid() << " (parent: " << ::getppid() << ")" << std::endl;
	std::this_thread::sleep_for(std::chrono::seconds(10));//sleep for gdb debug
  #endif
  //generate point clouds
  genPointClouds();
  //get neighbor list
  neighborSearch();
  //process
  processing();
}


BOOST_PYTHON_MODULE(cellfactory){

  py::scope().attr("__doc__")="CellFactory is a wrapper for handling Set Voronoi Tessellation.";

	py::docstring_options docopt;
	docopt.enable_all();
	docopt.disable_cpp_signatures();
  /*
  py::class_<parAttr>("parAttr","struct for particle attributes.")
    .add_property("scale",&CellMachine::get_scale,&CellMachine::set_scale,"Scale up data to avoid numerical issues when conducting tessellation.")
    ;
    */
  py::class_<CellFactory>("CellFactory","CellFactory.")
		.add_property("infolder",&CellFactory::get_infolder,&CellFactory::set_infolder,"directory where the data would be input.")
    .add_property("outfolder",&CellFactory::get_outfolder,&CellFactory::set_outfolder,"directory where the data would be output.")
    .add_property("searchRadius",&CellFactory::get_searchRadius,&CellFactory::set_searchRadius,"neighbour lists would be built within a sphere with a radius of seachRadius")
    .add_property("posFile",&CellFactory::get_posFile,&CellFactory::set_posFile,"path of the file of particle positions.")
    .add_property("cellVTK",&CellFactory::get_cellVTK,&CellFactory::set_cellVTK,"whether to write vtk files for all cells.")
    .add_property("parShrink",&CellFactory::get_parShrink,&CellFactory::set_parShrink,"shrink a particle by a small gap (parShrink) during generating its point cloud.")
    .def("processing",&CellFactory::processing,"Processing tessellation.")
    .def("genPointClouds",&CellFactory::genPointClouds,(py::arg("w_slices")=20,py::arg("h_slices")=20),"initialize a cell machine")
    .def("neighborSearch",&CellFactory::neighborSearch,"reset the cell machine for another cell computation if needed.")
    .def("autoWorkFlow",&CellFactory::autoWorkFlow,"")
		;

	py::class_<CellMachine>("CellMachine","CellMachine builds a single Voronoi cell of a given particle surrounded by others.")
		.add_property("scale",&CellMachine::get_scale,&CellMachine::set_scale,"Scale up data to avoid numerical issues when conducting tessellation.")
.add_property("cellVTK",&CellMachine::get_cellVTK,&CellMachine::set_cellVTK,"whether to write vtk files for all cells.")
    .def("processing",&CellMachine::processing,"Processing tessellation.")
    .def("initial",&CellMachine::initial,"initialize a cell machine")
    .def("reset",&CellMachine::reset,"reset the cell machine for another cell computation if needed.")
		;
}
