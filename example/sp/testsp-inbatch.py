import sys,os,os.path
sys.path.insert(1, '../../install/lib')
#os.environ['OMP_NUM_THREADS']=str(2)
import setvoronoi as sv
#print cf.__doc__
mycf = sv.CellFactory()
#mycf.infolder = "./input"
#mycf.outfolder = "./output"
#mycf.posFile = "./Particles.dat"
#mycf.wallFile = "./Walls.dat"
mycf.cellVTK = True
mycf.cellPOV = False#True
mycf.scale = 1000
mycf.boxScale = 2.0
mycf.parShrink = 0.1e-3
mycf.threadNum = 2
mycf.visualized_ids = [1,2]
#you can execute it step by step
parents = ["./case1","./case2"]
file_lists = ['D_D_init', 'D_D_26', 'D_D_45', 'D_U_26', 'D_U_45', 'M_D_init', 'M_D_26', 'M_D_45', 'M_U_26', 'M_U_45','L_D_init', 'L_D_26', 'L_D_45']
for par in parents:
   mycf.infolder = par + "/input"
   mycf.outfolder = par + "/output"
   mycf.posFile = par + "/"+prefix+"_Particles.dat"
   mycf.wallFile = par + "/"+prefix+"_Walls.dat"
   mycf.genPointClouds(w_slices=30,h_slices=20)#here you can put your raw data
   mycf.neighborSearch()
   mycf.processing()#processing all particles
#pid = 0
#mycf.processingOne(pid)#processing only one particle with id of pid
#or you can conduct an automatic work flow
#mycf.autoWorkFlow()
