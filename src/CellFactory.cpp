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
#include "CellMachine.hpp"

namespace py = boost::python;

BOOST_PYTHON_MODULE(cellfactory){
	//SUDODEM_SET_DOCSTRING_OPTS;
  py::scope().attr("__doc__")="CellFactory is a wrapper for handling Set Voronoi Tessellation.";

	py::docstring_options docopt;
	docopt.enable_all();
	docopt.disable_cpp_signatures();

	//py::def("View",createView,"Create a new 3d view.");


	py::class_<CellMachine>("CellMachine","CellMachine builds a single Voronoi cell of a given particle surrounded by others.")
		.add_property("scale",&CellMachine::get_scale,&CellMachine::set_scale,"Scale up data to avoid numerical issues when conducting tessellation.")
		.def("processing",&CellMachine::processing,"Processing tessellation.")
		;
}
