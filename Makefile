#CXX = clang++ -Wall -Wextra -O3  -I./lib/Lua#-std=c++11
CXX = g++ -Wall -Wextra -g -O3 -std=c++11 -I./lib/Lua # -I/scratch/softwares/gcc4.9/include
CXXVORO = g++ -std=c++11 -g -O3
# Relative include and library paths for compilation of the examples
E_INC=-I./lib/eigen-3.2.3

MAIN = src/main

all: obj/main.o LINK 

obj/voro.o: lib/voro++/src/voro++.*
	$(CXXVORO) -c -o obj/voro.o lib/voro++/src/voro++.cc

obj/fileloader.o: src/fileloader.*
	$(CXX) -c -o obj/fileloader.o src/fileloader.cpp

obj/pointpattern.o: src/pointpattern.*
	$(CXX) -c -o obj/pointpattern.o src/pointpattern.cpp
#obj/ray.o: src/ray.*
#	$(CXX) -c -o obj/ray.o src/ray.cpp

#obj/raytracer.o: obj/ray.o src/raytracer.*
#	$(CXX) -c -o obj/raytracer.o src/raytracer.cpp

obj/main.o: src/main.cpp  src/duplicationremover.hpp src/polywriter.hpp src/postprocessing.hpp src/Superquadrics.hpp src/process.hpp
	$(CXX)  -c $(E_INC) -o obj/main.o src/main.cpp 

LINK: obj/main.o obj/voro.o obj/fileloader.o obj/pointpattern.o
	$(CXX) obj/main.o obj/voro.o obj/fileloader.o obj/pointpattern.o -o bin/pomelo -L./lib/Lua/ -llua -lm -ldl



%o: %.cpp 
	$(CXX) -o $*.o -c $<

clean:
	rm obj/*
	rm bin/pomelo
