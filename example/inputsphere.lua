--process multi-jobs in batch
isBatch = false
DatasetDir = "./input/"

--folder for output files
outputFolder = "./output"
inputFolder = "./input"

--process only one job
positionfile = "./input/Particles.dat"
wallfile = "./input/walls.dat"
RawData = false


cellVTK = true
cellIDsfile = "./ids.dat"


w_slices = 30
h_slices = 30

nx = 120
ny = 120
nz = 120

scale = 0.01e-3

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
