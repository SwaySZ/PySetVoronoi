#CXX = clang++ -Wall -Wextra -O3  -I./lib/Lua#-std=c++11
CXX = g++ -Wall -Wextra -g -O3 -std=c++11  -fPIC#-I./lib/Lua # -I/scratch/softwares/gcc4.9/include
#CXXVORO = g++ -std=c++11 -g -O3
# Relative include and library paths for compilation of the examples
E_INC=-I./lib/eigen3.3.5
#voro++
VORO_INC=-I./lib/voro++/src/
#boost
BOOST_INC=-I/home/swayzhao/software/DEM/3rdlib/boost167/include
BOOST_LIB=-L/home/swayzhao/software/DEM/3rdlib/boost167/lib
#python 2.7
PYTHON_INC=-I/usr/include/python2.7
MAIN = src/CellFactory


all: obj/cellfactory.o LINK

obj/voro.o: lib/voro++/src/voro++.*
	$(CXX) -c  -o  obj/voro.o lib/voro++/src/voro++.cc

obj/fileloader.o: src/fileloader.*
	$(CXX) -c -o obj/fileloader.o src/fileloader.cpp

obj/pointpattern.o: src/pointpattern.*
	$(CXX) -c -o obj/pointpattern.o src/pointpattern.cpp
obj/cellmachine.o: src/CellMachine.*
		$(CXX) $(VORO_INC) -c  -o obj/cellmachine.o src/CellMachine.cpp

obj/cellfactory.o: src/CellFactory.cpp
	$(CXX)  -c  $(E_INC) $(PYTHON_INC) $(BOOST_INC) $(VORO_INC) -o obj/cellfactory.o src/CellFactory.cpp

LINK: obj/cellfactory.o obj/voro.o obj/fileloader.o obj/pointpattern.o obj/cellmachine.o
	$(CXX) obj/cellfactory.o obj/voro.o obj/fileloader.o obj/pointpattern.o obj/cellmachine.o -o bin/cellfactory.so -shared `pkg-config python --cflags` $(BOOST_LIB) -lboost_python27 -lm -ldl



%o: %.cpp
	$(CXX) -o $*.o -c $<

clean:
	rm obj/*
	rm bin/pomelo
