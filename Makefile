#CXX = clang++ -Wall -Wextra -O3  -I./lib/Lua#-std=c++11
CXX = g++ -Wextra -g -O3 -std=c++11 -fopenmp -fPIC#-I./lib/Lua # -I/scratch/softwares/gcc4.9/include
#CXXVORO = g++ -std=c++11 -g -O3
# Relative include and library paths for compilation of the examples
E_INC=-I./thirdparty/eigen3.3.5
#voro++
VORO_INC=-I./thirdparty/voro++/src/
#boost
BOOST_INC=-I/home/swayzhao/software/DEM/3rdlib/boost167/include
BOOST_LIB=-L/home/swayzhao/software/DEM/3rdlib/boost167/lib
#python 2.7
PYTHON_INC=-I/usr/include/python2.7
MAIN = src/SetVoronoi


MKDIR_P = mkdir -p
OUT_DIR = ./obj ./lib
.PHONY: directories

all: directories obj/setvoronoi.o LINK

directories: ${OUT_DIR}

${OUT_DIR}:
	${MKDIR_P} ${OUT_DIR}


obj/mathbase.o: src/mathbase.cpp
	$(CXX) -c $(E_INC) -o obj/mathbase.o src/mathbase.cpp
obj/voro.o: thirdparty/voro++/src/voro++.*
	$(CXX) -c  -o  obj/voro.o thirdparty/voro++/src/voro++.cc

obj/pointpattern.o: src/pointpattern.*
	$(CXX) -c -o obj/pointpattern.o src/pointpattern.cpp
obj/cellmachine.o: src/CellMachine.*
		$(CXX) $(E_INC) $(VORO_INC) -c  -o obj/cellmachine.o src/CellMachine.cpp
obj/superquadric.o: src/Superquadrics.cpp
	$(CXX)  -c  $(E_INC) -o obj/superquadric.o src/Superquadrics.cpp
obj/cellfactory.o: src/CellFactory.cpp
	$(CXX)  -c  $(E_INC) $(VORO_INC) $(BOOST_INC) -o obj/cellfactory.o src/CellFactory.cpp
obj/setvoronoi.o: src/SetVoronoi.cpp
	$(CXX)  -c  $(E_INC) $(PYTHON_INC) $(BOOST_INC) $(VORO_INC) -o obj/setvoronoi.o src/SetVoronoi.cpp
LINK: obj/setvoronoi.o obj/cellfactory.o obj/voro.o obj/pointpattern.o obj/cellmachine.o obj/superquadric.o obj/mathbase.o
	$(CXX) obj/setvoronoi.o obj/cellfactory.o obj/voro.o obj/pointpattern.o obj/cellmachine.o obj/superquadric.o obj/mathbase.o -o lib/setvoronoi.so -shared `pkg-config python --cflags` $(BOOST_LIB) -lboost_python27 -lboost_filesystem -lm -ldl  -lpthread

%o: %.cpp
	$(CXX) -o $*.o -c $<

clean:
	rm obj/*
	rm lib/*
