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
#include "../config.hpp"
#include "../CellMachine.hpp"
#include "../CellFactory.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h> //Some automatic conversions
#include <pybind11/stl_bind.h>
#include <pybind11/iostream.h> //Capturing standard output from ostream

namespace py = pybind11;

  PYBIND11_MODULE(setvoronoi, m) {

  
  py::class_<CellFactory>(m,"CellFactory")
    .def(py::init<>())
		.def_property("infolder",&CellFactory::get_infolder,&CellFactory::set_infolder,"directory where the data would be input.")
    .def_property("outfolder",&CellFactory::get_outfolder,&CellFactory::set_outfolder,"directory where the data would be output.")
    .def_property("searchRadius",&CellFactory::get_searchRadius,&CellFactory::set_searchRadius,"neighbour lists would be built within a sphere with a radius of seachRadius")
    .def_property("posFile",&CellFactory::get_posFile,&CellFactory::set_posFile,"path of the file of particle positions.")
    .def_property("wallFile",&CellFactory::get_wallFile,&CellFactory::set_wallFile,"path of the file of wall positions.")
    .def_property("cellVTK",&CellFactory::get_cellVTK,&CellFactory::set_cellVTK,"whether to write vtk files for all cells.")
    .def_property("cellPOV",&CellFactory::get_cellPOV,&CellFactory::set_cellPOV,"whether to write pov files for all cells.")
    .def_property("visualized_ids",&CellFactory::get_visual_ids,&CellFactory::set_visual_ids,"for the cells with specified ids to be visualized.")
    .def_property("scale",&CellFactory::get_scale,&CellFactory::set_scale,"")
    .def_property("boxScale",&CellFactory::get_boxScale,&CellFactory::set_boxScale,"")
    .def_property("verbose",&CellFactory::get_verbose,&CellFactory::set_verbose,"")
    .def_property("parShrink",&CellFactory::get_parShrink,&CellFactory::set_parShrink,"shrink a particle by a small gap (parShrink) during generating its point cloud.")
    .def_property("threadNum",&CellFactory::get_threadNum,&CellFactory::set_threadNum,"the thread number for openMP parallization.")
    .def("processing",&CellFactory::processing,"Processing tessellation.")
    .def("genPointClouds",&CellFactory::genPointClouds,py::arg("w_slices")=20,py::arg("h_slices")=20,"initialize a cell machine")
    .def("neighborSearch",&CellFactory::neighborSearch,"reset the cell machine for another cell computation if needed.")
    .def("autoWorkFlow",&CellFactory::autoWorkFlow,"")
    .def("processingOne",&CellFactory::processingOne,"process only one particle.")
		;

	py::class_<CellMachine>(m,"CellMachine","CellMachine builds a single Voronoi cell of a given particle surrounded by others.")
    .def(py::init<>())
		.def_property("scale",&CellMachine::get_scale,&CellMachine::set_scale,"Scale up data to avoid numerical issues when conducting tessellation.")
    .def_property("cellVTK",&CellMachine::get_cellVTK,&CellMachine::set_cellVTK,"whether to write vtk files for all cells.")
    .def("processing",&CellMachine::processing,"Processing tessellation.")
    .def("initial",&CellMachine::initial,"initialize a cell machine")
    .def("reset",&CellMachine::reset,"reset the cell machine for another cell computation if needed.")
		;
}
