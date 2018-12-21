--process multi-jobs in batch
isBatch = false
DatasetDir = "./sphere/"

--folder for output files
outputFolder = "./sphere/output"

--process only one job
positionfile = "./sphere/sphere0.dat"
wallfile = "./sphere/Walls0.dat"
RawData = false


cellVTK = true
cellIDsfile = "./sphere/ids.dat"


w_slices = 20
h_slices = 20

nx = 120
ny = 120
nz = 120

scale = 0.05

epsilon = 1e-6
rr=0.01
boundary = "none"
withboundary = true
postprocessing = false
savepoly = true
savereduced = false
savevtk=false
savesurface = true
savepov=true
outputfile="data"
--output on console
outputSteps = 100
