import sys,os,os.path
sys.path.insert(1, '/home/swayzhao/software/DEM/EPomelo2/bin')

import setvoronoi as sv
#print cf.__doc__
mycf = sv.CellFactory()
mycf.infolder = "./input"
mycf.outfolder = "./output"
mycf.posFile = "./Particles.dat"
mycf.wallFile = "./Walls.dat"
mycf.cellVTK = True
mycf.threadNum = 2
mycf.scale = 1000
mycf.boxScale = 3.0
mycf.parShrink = 0.1e-3
mycf.searchRadius = 4.0
#you can execute it step by step
mycf.genPointClouds(w_slices=50,h_slices=40)#here you can put your raw data
mycf.neighborSearch()
mycf.processing()#processing all particles
pid = 0
#mycf.processingOne(pid)#processing only one particle with id of pid
#or you can conduct an automatic work flow
#mycf.autoWorkFlow()

