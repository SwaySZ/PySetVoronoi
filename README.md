![](./doc/cover.png)
# PySetVoronoi

A Python interface for Set Voronoi Tessellation of poly-superellipsoids and general point clouds from irregular particles, where the kenerl is written in C++.

## License

PySetVoronoi is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PySetVoronoi is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the code.  If not, see <http://www.gnu.org/licenses/>.

> **Note: The source codes of the thirdparty libraries (Eigen, Voro++, Pybind11) have been presented here, and their licenses are distributed with the source codes.**

## Features

- Voronoi tesselatio of particles with arbitrary shapes
- Basic geometric info of individual cells
- Common Minkowski tensors for cell shape
- Support vtk and pov files for individual cells
- Support parallel computation with OpenMP
- Easy setting in Python

## Supported platforms

- [x] Windows
- [x] Linux
- [x] Mac

## Sample usage

```py
import sys,os,os.path
sys.path.insert(1, '../../install/lib')
#os.environ['OMP_NUM_THREADS']=str(2)
import setvoronoi as sv
#print cf.__doc__
mycf = sv.CellFactory() #establish a factory (class) to handle the computation
mycf.infolder = "./input" #this folder contains point-cloud data for each particle (will be generated)
mycf.outfolder = "./output" #this folder contains data after computation
mycf.posFile = "./Particles.dat" #particle data with info such as position, orientation, etc
mycf.wallFile = "./Walls.dat" #wall data with positions of each walll for a cubic container
mycf.cellVTK = True #yield vtk files for cells
mycf.cellPOV = True #yield pov files for cells
mycf.scale = 1000 #the parameter used to scale up the data during computation (due to a bug in vtk)
mycf.boxScale = 2.0 #the parameter used to scale up the AABB box of a given particle
mycf.parShrink = 0.1e-3 #shrink particles inward to avoid contact particles (with intersection in DEM)
mycf.threadNum = 2 #threads in OpenMP
mycf.visualized_ids = [0,1,2] #id list of particle/cell that will be visualized by vtk/pov, empty for all.
#you can execute it step by step
mycf.genPointClouds(w_slices=30,h_slices=20)#point-cloud generation. here you can put your raw data
mycf.neighborSearch()
#mycf.processing()#processing all particles
pid = 0 #as a demonstrate, we calculate only a single particle with id = 0.
mycf.processingOne(pid)#processing only one particle with id of pid
#or you can conduct an automatic work flow
#mycf.autoWorkFlow() #this line will execute all processes starting from point-cloud generation.
```

## Compile & Install

Your compiler should support the C++17 standard.

You can create a fresh folder like 'build' at the root folder of the source
 code, and then compile the code according to the following brief instructions.

cd PySetVoronoi/

mkdir build

cd build

cmake ../src

make

make install

> **Note: By default, 'make install' will install the compiled libraries into the folder named "install" next to 'build'.**

Then, you can run the examples in your terminal, e.g.,
```sh
python3 testsp.py
```

## Cite this work

If you plan to publish your research work using this code, please consider citing one of the following publications:

- Zhang, C., Zhao, S., Zhao, J., Zhou, X., 2021. **Three-dimensional Voronoi analysis of realistic grain packing: An XCT assisted set Voronoi tessellation framework**. _Powder Technology_ 379, 251–264. https://doi.org/10.1016/j.powtec.2020.10.054

- Zhao, S., Zhao, J., Guo, N., 2020. **Universality of internal structure characteristics in granular media under shear**. _Phys. Rev. E_ 101, 012906. https://doi.org/10.1103/PhysRevE.101.012906

- Zhao, S., Evans, T.M., Zhou, X., 2018. **Three-dimensional Voronoi analysis of monodisperse ellipsoids during triaxial shear**. _Powder Technology_ 323, 323–336. https://doi.org/10.1016/j.powtec.2017.10.023
