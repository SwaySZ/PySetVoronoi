import sys,os,os.path
sys.path.insert(1, '../../install/lib')
#os.environ['OMP_NUM_THREADS']=str(2)
import setvoronoi as sv
#print cf.__doc__
mycf = sv.CellFactory()
mycf.infolder = "./input"
mycf.outfolder = "./output"
mycf.posFile = "./Particles.dat"
mycf.wallFile = "./Walls.dat"
mycf.cellVTK = True
mycf.cellPOV = True
mycf.scale = 1000
mycf.boxScale = 2.0
mycf.parShrink = 0.1e-3
mycf.threadNum = 2
mycf.visualized_ids = [0,1,2]
#you can execute it step by step
mycf.genPointClouds(w_slices=30,h_slices=20)#here you can put your raw data
mycf.neighborSearch()
#mycf.processing()#processing all particles
pid = 0
mycf.processingOne(pid)#processing only one particle with id of pid
#or you can conduct an automatic work flow
#mycf.autoWorkFlow()
