project(setvoronoi CXX)
cmake_minimum_required(VERSION 3.10)

##########   PYBIND MODULE ##########

# Use just one of these:
# GCC/clang:
#set(PYBIND11_CPP_STANDARD -std=c++11)
#set(PYBIND11_CPP_STANDARD -std=c++14)
#set(PYBIND11_CPP_STANDARD -std=c++1z) # Experimental C++17 support
# MSVC:
#set(PYBIND11_CPP_STANDARD /std:c++14)
#set(PYBIND11_CPP_STANDARD /std:c++latest) # Enables some MSVC C++17 features

# Temporarily, ignore the python library
add_subdirectory(${CMAKE_SOURCE_DIR}/../thirdparty/pybind11 ${CMAKE_BINARY_DIR}/pybind11build)  # or find_package(pybind11)


## here use the CMAKE_INSTALL_PREFIX
## default cmake install prefix
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/../install" CACHE PATH "This is default path" FORCE)
endif()

#user defined cmake install prefix
GET_FILENAME_COMPONENT(CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX} ABSOLUTE)
MESSAGE(STATUS "setvoronoi will be installed to ${CMAKE_INSTALL_PREFIX}")

if(APPLE)
  cmake_policy(SET CMP0042 NEW)
  set(CMAKE_C_COMPILER clang)
  set(CMAKE_CXX_COMPILER clang++)
  set(CMAKE_CXX_STANDARD 17)
  if(CMAKE_C_COMPILER_ID MATCHES "Clang\$")
      set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
      set(OpenMP_C_LIB_NAMES "omp")
      set(OpenMP_omp_LIBRARY /usr/local/opt/libomp/lib/libomp.dylib)
  endif()

  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang\$")
      set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
      set(OpenMP_CXX_LIB_NAMES "omp")
      set(OpenMP_omp_LIBRARY /usr/local/opt/libomp/lib/libomp.dylib)
  endif()

endif()

IF (LINUX)
  SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS}  -save-temps -ftrack-macro-expansion=0 -fopenmp")
ENDIF()

set(DEBUG OFF)
IF (DEBUG)
	SET(CMAKE_VERBOSE_MAKEFILE 1)
    SET(CMAKE_BUILD_TYPE Debug)
ELSE (DEBUG)
	SET(CMAKE_BUILD_TYPE Release)
    SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fPIC -O2")
ENDIF (DEBUG)
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++17 -Wno-deprecated-declarations")




find_package(OpenMP REQUIRED)
MESSAGE(STATUS "OpenMP FOUND: ${OpenMP_FOUND} version: ${OpenMP_VERSION}")

#set include and lib for openMP on MacOS
IF(APPLE)
  INCLUDE_DIRECTORIES("/usr/local/opt/libomp/include") #for omp.h
  #LINK_DIRECTORIES("/usr/local/opt/libomp/lib") #for libomp.dylib
ENDIF()


## install path
SET(INSTALL_LIBPATH "${CMAKE_INSTALL_PREFIX}/lib")

SET(DEFAULT ON CACHE INTERNAL "Default value for enabled by default options")

## ADD HEADER FILES
SET(EIGEN_INCLUDE "${CMAKE_SOURCE_DIR}/../thirdparty/eigen3.3.5/")
SET(VORO_INCLUDE "${CMAKE_SOURCE_DIR}/../thirdparty/voro++/src/")
INCLUDE_DIRECTORIES(${EIGEN_INCLUDE})
INCLUDE_DIRECTORIES(${VORO_INCLUDE})

## ADD SRC FILES
FILE(GLOB SRC_LIB "*.cpp")

FILE(GLOB VORO_SRC_LIB  "${CMAKE_SOURCE_DIR}/../thirdparty/voro++/src/voro++.cc")
ADD_LIBRARY(voro++ SHARED ${VORO_SRC_LIB})
SET_TARGET_PROPERTIES(voro++ PROPERTIES PREFIX "" INSTALL_RPATH "$ORIGIN;")
TARGET_LINK_LIBRARIES(voro++)
INSTALL(TARGETS voro++ DESTINATION ${INSTALL_LIBPATH})

## ADD LIBRARY
add_library(setvoronoi MODULE ${CMAKE_CURRENT_SOURCE_DIR}/pymodule/SetVoronoi.cpp)
ADD_LIBRARY(setvoronoilib SHARED ${SRC_LIB})
SET_TARGET_PROPERTIES(setvoronoilib  PROPERTIES INSTALL_RPATH "$ORIGIN" )

IF(LINUX)
  #SET_TARGET_PROPERTIES(setvoronoilib  PROPERTIES LINK_FLAGS "-Wl,--as-needed" INSTALL_RPATH "$ORIGIN" )
ENDIF()
IF(APPLE)
  SET_TARGET_PROPERTIES(setvoronoilib  PROPERTIES LINK_FLAGS "-Wl" INSTALL_RPATH "$ORIGIN" )
ENDIF()


## LINK
TARGET_LINK_LIBRARIES(setvoronoilib  PRIVATE OpenMP::OpenMP_CXX voro++)

TARGET_LINK_LIBRARIES(setvoronoi PRIVATE pybind11::module OpenMP::OpenMP_CXX setvoronoilib) 

IF(MINGW)
  SET_TARGET_PROPERTIES(setvoronoi PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}" SUFFIX "${PYTHON_MODULE_EXTENSION}" INSTALL_RPATH "${MPM_INSTALL_LIBPATH};${DEM_INSTALL_LIBPATH};$ORIGIN;")
ELSE(MINGW)
  SET_TARGET_PROPERTIES(setvoronoi PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}" INSTALL_RPATH "${INSTALL_LIBPATH};$ORIGIN;")
ENDIF()
## INSTALL
INSTALL(TARGETS setvoronoilib DESTINATION ${INSTALL_LIBPATH})
INSTALL(TARGETS setvoronoi DESTINATION ${INSTALL_LIBPATH})


