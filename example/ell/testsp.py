import sys,os,os.path
sys.path.insert(1, '/home/swayzhao/software/DEM/EPomelo2/bin')

import cellfactory as cf
print cf.__doc__
mycf = cf.CellFactory()
mycf.infolder = "./input"
mycf.outfolder = "./output"
mycf.posFile = "./input/Particles.dat"
mycf.cellVTK = True
mycf.scale = 1.0
mycf.parShrink = 0.01
#you can execute it step by step
mycf.genPointClouds(w_slices=30,h_slices=30)#here you can put your raw data
mycf.neighborSearch()
mycf.processing()
#or you can conduct an automatic work flow
#mycf.autoWorkFlow()
