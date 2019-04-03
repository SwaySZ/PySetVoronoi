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
#include "config.hpp"
#include "CellMachine.hpp"
#include "CellFactory.hpp"
#include<boost/python.hpp>
#include <boost/foreach.hpp>//BOOST_FOREACH

namespace py = boost::python;

/*** c++-list to python-list */
template<typename containedType>
struct custom_vector_to_list{
	static PyObject* convert(const std::vector<containedType>& v){
		boost::python::list ret;
    	BOOST_FOREACH(const containedType& e, v) ret.append(e);
		return boost::python::incref(ret.ptr());
	}
};
template<typename containedType>
struct custom_vector_from_seq{
	custom_vector_from_seq(){ boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id<std::vector<containedType> >()); }
	static void* convertible(PyObject* obj_ptr){
		// the second condition is important, for some reason otherwise there were attempted conversions of Body to list which failed afterwards.
		if(!PySequence_Check(obj_ptr) || !PyObject_HasAttrString(obj_ptr,"__len__")) return 0;
		return obj_ptr;
	}
	static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data){
		 void* storage=((boost::python::converter::rvalue_from_python_storage<std::vector<containedType> >*)(data))->storage.bytes;
		 new (storage) std::vector<containedType>();
		 std::vector<containedType>* v=(std::vector<containedType>*)(storage);
		 int l=PySequence_Size(obj_ptr); if(l<0) abort(); /*std::cerr<<"l="<<l<<"; "<<typeid(containedType).name()<<std::endl;*/ v->reserve(l); for(int i=0; i<l; i++) { v->push_back(boost::python::extract<containedType>(PySequence_GetItem(obj_ptr,i))); }
		 data->convertible=storage;
	}
};


BOOST_PYTHON_MODULE(setvoronoi){

  py::scope().attr("__doc__")="CellFactory is a wrapper for handling Set Voronoi Tessellation.";

	py::docstring_options docopt;
	docopt.enable_all();
	docopt.disable_cpp_signatures();


  #define VECTOR_SEQ_CONV(Type) custom_vector_from_seq<Type>();  boost::python::to_python_converter<std::vector<Type>, custom_vector_to_list<Type> >();
		VECTOR_SEQ_CONV(int);
		//VECTOR_SEQ_CONV(bool);
  #undef VECTOR_SEQ_CONV


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
    .add_property("wallFile",&CellFactory::get_wallFile,&CellFactory::set_wallFile,"path of the file of wall positions.")
    .add_property("cellVTK",&CellFactory::get_cellVTK,&CellFactory::set_cellVTK,"whether to write vtk files for all cells.")
    .add_property("cellPOV",&CellFactory::get_cellPOV,&CellFactory::set_cellPOV,"whether to write pov files for all cells.")
    .add_property("visualized_ids",&CellFactory::get_visual_ids,&CellFactory::set_visual_ids,"for the cells with specified ids to be visualized.")
    .add_property("scale",&CellFactory::get_scale,&CellFactory::set_scale,"")
    .add_property("boxScale",&CellFactory::get_boxScale,&CellFactory::set_boxScale,"")
    .add_property("parShrink",&CellFactory::get_parShrink,&CellFactory::set_parShrink,"shrink a particle by a small gap (parShrink) during generating its point cloud.")
    .add_property("threadNum",&CellFactory::get_threadNum,&CellFactory::set_threadNum,"the thread number for openMP parallization.")
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
