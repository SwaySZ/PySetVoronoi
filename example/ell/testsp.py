import sys,os,os.path
sys.path.insert(1, '../../install/lib')

import setvoronoi as sv
print sv.__doc__
mycf = sv.CellFactory()
mycf.infolder = "./input"
mycf.outfolder = "./output"
mycf.posFile = "./Particles.dat"
mycf.wallFile = "./Walls.dat"
mycf.threadNum = 2#how many threads used
mycf.cellVTK = True
mycf.scale = 1.0
mycf.boxScale =2.0
mycf.parShrink = 0.1
pid = 0
#you can execute it step by step
mycf.genPointClouds(w_slices=50,h_slices=40)#here you can put your raw data
mycf.neighborSearch()
#mycf.processing()
mycf.processingOne(pid)
#or you can conduct an automatic work flow
#mycf.autoWorkFlow()
