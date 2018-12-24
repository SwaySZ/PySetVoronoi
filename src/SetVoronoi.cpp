/*
 * =====================================================================================
 *
 *       Filename:  SetVoronoi.cpp
 *
 *    Description:  a wrapper of CellFactory and CellMachine exposed in Python.
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
#include "CellMachine.hpp"
#include "CellFactory.hpp"
#include<boost/python.hpp>

namespace py = boost::python;

BOOST_PYTHON_MODULE(setvoronoi){

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
    .def("processingOne",&CellFactory::processingOne,"process only one particle.")
		;

	py::class_<CellMachine>("CellMachine","CellMachine builds a single Voronoi cell of a given particle surrounded by others.")
		.add_property("scale",&CellMachine::get_scale,&CellMachine::set_scale,"Scale up data to avoid numerical issues when conducting tessellation.")
    .add_property("cellVTK",&CellMachine::get_cellVTK,&CellMachine::set_cellVTK,"whether to write vtk files for all cells.")
    .def("processing",&CellMachine::processing,"Processing tessellation.")
    .def("initial",&CellMachine::initial,"initialize a cell machine")
    .def("reset",&CellMachine::reset,"reset the cell machine for another cell computation if needed.")
		;
}
