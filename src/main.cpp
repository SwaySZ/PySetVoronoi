#include <iostream>
#include <fstream>
#include <map>
//#include "process.hpp"
#include <limits>
#include <cmath>
#include <sys/stat.h>
#include "include.hpp"
#include "fileloader.hpp"
#include "pointpattern.hpp"
#include "duplicationremover.hpp"
#include "polywriter.hpp"
#include "postprocessing.hpp"
#include "Superquadrics.hpp"
#include <chrono>
#include <thread>

#include <string>
#include <vector>
//#define POLYDATA 1
extern "C"{
#include <lualib.h>
//#include <lua2.2/lauxlib.h>
#include <lua.h>
}
std::string version = "0.5.0";

using namespace sel;
using namespace voro;
struct particleAttr{
		int ID;
		double centerx;
		double centery;
		double centerz;
		double xmin,xmax,ymin,ymax,zmin,zmax;
		double xrange,yrange,zrange;
		double radius;
		std::vector<int> surroundedID;
};
std::string char2string(const char* charstr){
        std::string s(charstr);
        return s;
}


const char* con_fname(int i,std::string path,std::string prefix, std::string suffix){
        //std::string str = "contactinfo-";
        //std::cout << str << std::endl;
        std::stringstream out;
        std::string s;
        out << i;
        s = out.str();
        prefix = path+prefix + s+suffix;
	std::cout << "file prefix: " << prefix << std::endl;
        return prefix.c_str();
}
void calculate_Attributes(std::vector<pointpattern>& pplist,std::vector<particleAttr>& parAttrlist,double radiation_radius ,double epsilon){
	for(unsigned int i = 0;i<pplist.size();i++){
//		std::cout<<"number of particle list:"<< pplist.size()<<std::endl;

		pointpattern pp = pplist.at(i);
		if(pp.points.empty()){
			continue;
		}

		double xsum = 0;
		double ysum = 0;
		double zsum = 0 ;
		double xmax = 0 , ymax = 0 , zmax = 0 ;
		double xmin = 1e10 , ymin = 1e10,zmin = 1e10;
		particleAttr pta;

		for(unsigned int j = 0 ;j<pp.points.size();j++){
			point p = pp.points.at(j);
			pta.ID = p.l;
			xsum += p.x;
			ysum += p.y;
			zsum += p.z;
			if(p.x > xmax )xmax = p.x;
			if(p.y > ymax )ymax = p.y;
			if(p.z > zmax) zmax = p.z;
			if(p.x < xmin) xmin = p.x;
			if(p.y < ymin) ymin = p.y;
			if(p.z < zmin) zmin = p.z;
		}
		pta.centerx = xsum / pp.points.size();
		pta.centery = ysum / pp.points.size();
		pta.centerz = zsum / pp.points.size();
		pta.xmax = xmax;
		pta.xmin = xmin;
		pta.ymax = ymax;
		pta.ymin = ymin;
		pta.zmax = zmax;
		pta.zmin = zmin;
		pta.xrange = xmax-xmin;
		pta.yrange = ymax-ymin;
		pta.zrange = zmax-zmin;
		//max range
		double maxrange = pta.xrange;
		if(pta.yrange > maxrange) maxrange = pta.yrange;
		if(pta.zrange > maxrange) maxrange = pta.zrange;
		pta.radius = maxrange/2;
		parAttrlist.push_back(pta);
	}
	//search for particles within the specific range and get the ID of the surrounding particles

	//using Exhaustive method
	int particlenum = parAttrlist.size();
	for(int i = 0 ;i< particlenum;i++){
//		particleAttr p1;
		double dist = 0;
		for(int j = 0 ; j < particlenum;j++){
//			if(j == i) continue;//
			particleAttr p2  = parAttrlist.at(j);
			dist = sqrt(pow((parAttrlist.at(i).centerx-p2.centerx),2)+pow((parAttrlist.at(i).centery-p2.centery),2)+pow((parAttrlist.at(i).centerz-p2.centerz),2));

			if(parAttrlist.at(i).radius * radiation_radius > dist ){
				//p2 close to p1
				parAttrlist.at(i).surroundedID.push_back(p2.ID); // find all particle close to p1 not including current particle itself;
			}
//			std::cout<<"dist:"<<dist<<"\tp1.radius:"<<p1.radius * radiation_radius - dist<<std::endl;
		}
//		std::cout<<p1.radius<<"\t"<<p1.surroundedID.size()<<std::endl;
//		for(int k = 0 ;k<p1.surroundedID.size();k++){
//			std::cout<<p1.surroundedID.at(k)<<"\t";
//		}
//		std::cout<<std::endl;
	}
//	std::cout<<parAttrlist.at(0).surroundedID.size()<<std::endl;
}

void  getFaceVerticesOfFace( std::vector<int>& f, unsigned int k, std::vector<unsigned int>& vertexList)
{
    vertexList.clear();

    // iterate through face-vertices vector bracketed
    // (number of vertices for face 1) vertex 1 ID face 1, vertex 2 ID face 1 , ... (number of vertices for face 2, ...
    unsigned long long index = 0;
    // we are at face k, so we have to iterate through the face vertices vector up to k
    for (unsigned long long cc = 0; cc <= k; cc++)
    {

        unsigned long long b = f[index];    // how many vertices does the current face (index) have?
        // iterate index "number of vertices of this face" forward
        for (unsigned long long bb = 1; bb <= b; bb++)
        {
            index++;
            // if we have found the correct face, get all the vertices for this face and save them
            if (cc == k)
            {
                //std::cout << "\t" << f[index] << std::endl;
                int vertexindex = f[index];
                vertexList.push_back(vertexindex);
            }
        }
        index++;
    }
}

std::map < unsigned long long, double > cellVolumes;//
std::map < unsigned long long, double > cellSurfaceArea;//surface area of each cell
std::map < unsigned long long, double > particleVolumes;//
std::map < unsigned long long, double > particleSurfaceArea;//
std::map < unsigned long long, std::vector<double> > cellNormalTensor;//six components of normal tensor of a cell surface, i.e., n11,n12,n13,n22,n23,n33.
std::map < unsigned long long, std::vector<double> > cellNormalAreaTensor;//six components of normal by area tensor of a cell surface, i.e., n11*a,n12*a,n13*a,n22*a,n23*a,n33*a.
fileloader loader;
void local_voronoi(State& state,std::vector<pointpattern>& pplist,std::string wallfile,double radiation_radius ,double epsilon,std::string sequenceNum)
{
//	std::vector<particleAttr> parAttrlist;

//we need to know xmin, xmax, ymin, ymax, zmin, zmax, nx, ny, nz, xpbc, ypbc, zpbc, 32
 int nx, ny, nz;
    nx = state["nx"];
    ny = state["ny"];
    nz = state["nz"];
    //import wall position
    double xmin1,xmax1,ymin1,ymax1,zmin1,zmax1;
    if (wallfile.empty()){
        std::cout<<"Wall file is not specified. Container size is specified manually."<<std::endl;
        xmin1 = state["xmin"];
        ymin1 = state["ymin"];
        zmin1 = state["zmin"];
        xmax1 = state["xmax"];
        ymax1 = state["ymax"];
        zmax1 = state["zmax"];
    }else{
        std::cout<<"Import wall info..."<<std::endl;
        loader.readWall(wallfile, xmin1,xmax1,ymin1,ymax1,zmin1,zmax1);
    }
     double xmin = xmin1;
     double xmax = xmax1;
     double ymin = ymin1;
     double ymax = ymax1;
     double zmin = zmin1;
     double zmax = zmax1;

     std::cout<<xmin<<"\t"<<xmax<<"\t"<<ymin<<"\t"<<ymax<<"\t"<<zmin<<"\t"<<zmax<<std::endl;
	//xpbc , ypbc , zpbc what are their meanings?
    bool xpbc = false;
    bool ypbc = false;
    bool zpbc = false;

    std::string boundary = state["boundary"];
	//周期性边界？？
    if(boundary == "periodic")
    {
        std::cout << "boundary condition mode 'periodic' selected." ;
        xpbc = state["xpbc"];
        ypbc = state["ypbc"];
        zpbc = state["zpbc"];
        std::cout << "x: " << xpbc << "\ny: " << ypbc << "\nz: " << zpbc << std::endl;
    }
    else if (boundary == "none")
    {
        std::cout << "boundary condition mode 'none' selected." << std::endl;
    }
    else
    {
        std::cerr << "bondary condition mode " << boundary << " not known" << std::endl;
    }
    std::cout << std::endl;


 container con(xmin, xmax, ymin, ymax, zmin, zmax, nx, ny, nz, xpbc, ypbc, zpbc, 32);

 std::vector<particleAttr> parAttrlist;
 calculate_Attributes(pplist, parAttrlist,radiation_radius ,epsilon);
 std::cout<<"number of particle list:"<< pplist.size()<<std::endl;
 std::cout<<"number of particle Attributes  list:"<< parAttrlist.size()<<std::endl;
 //now we get Attributes of particles
 //we just add several particles to the container and at last return only one cell;
 //std::vector<polywriter> pwlist;
 // print out pointpattern file
 std::string folder = state["outputFolder"];
 for(unsigned int id = 0 ; id <parAttrlist.size();id++){
	 int currentparticleID = id;
	 std::cout << "creating label id map:           currentparticleID: "<<currentparticleID ;
	 std::map < unsigned long long, unsigned long long > labelidmap;
	 unsigned long long pid = 0;//point id
	 particleAttr parAttr = parAttrlist.at(id);
	 std::vector<int> surroundedID = parAttr.surroundedID;
	 std::cout<<"    surroundedID size: "<<surroundedID.size()<< "   ";



	 for(unsigned int j = 0 ;j<surroundedID.size();j++){
		 int particleID = surroundedID.at(j);
		 pointpattern pp = pplist.at(particleID);
		       for(    auto it = pp.points.begin();
		               it != pp.points.end();
		               ++it)
		       {
		           con.put(pid, it->x, it->y, it->z);
		           unsigned int l = it->l;
		           labelidmap[pid] = l;
		           ++pid;
		       }
	 }
	 std::cout<<std::endl;
	 polywriter pw;
	 unsigned long long numberofpoints = labelidmap.size();
	 // merge voronoi cells to set voronoi diagram
	 std::cout << "merge voronoi cells :    ";
	 // loop over all voronoi cells
	 c_loop_all cla(con);
	 // cell currently worked on
	  unsigned long long status = 0;
	  // counter for process output
	  double outputSteps = state["outputSteps"];
	  double tenpercentSteps = 1/outputSteps*static_cast<double>(numberofpoints);
	  std::cout<<"number of points : " << numberofpoints ;
	  double target = tenpercentSteps;
	  if(cla.start())
	  {
		  std::cout << "started\n" << std::flush;
		         do
		         {
		             voronoicell_neighbor c;
		             status++;
		             if(con.compute_cell(c,cla))
		             {
		                 if ( status >= target)
		                 {
		                     target += tenpercentSteps;
		                     std::cout << static_cast<int>(static_cast<double>(status)/static_cast<double>(numberofpoints)*outputSteps) << " \%%\t" << std::flush;
		                 }
		                 //std::cout << "computed"  << std::endl;
		                 double xc = 0;
		                 double yc = 0;
		                 double zc = 0;
		                 // Get the position of the current particle under consideration
		                 cla.pos(xc,yc,zc);

		                 unsigned int id = cla.pid();
		                 unsigned long long l = labelidmap[id];
		                 //cell volume
		                 cellVolumes[l] += c.volume();
		                 std::vector<int> f; // list of face vertices (bracketed, as ID)
		                 c.face_vertices(f);

		                 std::vector<double> vertices;   // all vertices for this cell
		                 c.vertices(xc,yc,zc, vertices);

		                 std::vector<int> w; // neighbors of faces
		                 c.neighbors(w);
		                 // for this cell, loop over all faces and get the corresponding neighbors
		                 for (unsigned long long k = 0; k != w.size(); ++k)
		                 {
		                     // compare if id for this cell and the face-neighbor is the same
		                     int n = w[k];   // ID of neighbor cell
		                     if (labelidmap[n] == l)
		                     {
		                         // discard this neighbour/face, since they have the same id
		                         // std::cout << "discarding face " << l << " " << labelidmap[n] << std::endl;
		                     }
		                     else
		                     {
		                         std::vector<unsigned int> facevertexlist;
		                         getFaceVerticesOfFace(f, k, facevertexlist);
		                         std::vector<double> positionlist;
		                         for (
		                             auto it = facevertexlist.begin();
		                             it != facevertexlist.end();
		                             ++it)
		                         {
		                             unsigned int vertexindex = (*it);
		                             double x = vertices[vertexindex*3];
		                             double y = vertices[vertexindex*3+1];
		                             double z = vertices[vertexindex*3+2];
		                             positionlist.push_back(x);
		                             positionlist.push_back(y);
		                             positionlist.push_back(z);
		                         }
		                          //we only add faces belonging to the current particle
		 							pw.addface(positionlist, l);

		                         //#endif
		 			//calculate surface area
		                         double ux, uy, uz, vx,vy,vz, wx,wy,wz,area;
		                         area = 0.0;
		                         for (unsigned int i=1;i < positionlist.size()/3-1;i++)
		                         {
		                             ux = positionlist[3*i] - positionlist[0];
		                             uy = positionlist[3*i+1] - positionlist[1];
		                             uz = positionlist[3*i+2] - positionlist[2];
		                             vx = positionlist[3*i+3] - positionlist[0];
		                             vy = positionlist[3*i+4] - positionlist[1];
		                             vz = positionlist[3*i+5] - positionlist[2];
		                             wx = uy*vz - uz*vy;
		                             wy = uz*vx - ux*vz;
		                             wz = ux*vy - uy*vx;
		                             double area_tmp1 = wx*wx+wy*wy+wz*wz;//squarenorm of the vector (wx,wy,wz)
		                             double area_tmp2 = sqrt(area_tmp1);//
		                             //area += sqrt(wx*wx+wy*wy+wz*wz);
		                             //if(i<2){//just get a normal vector of the first triangle
		                                 //normalize the vector

																		 //std::cout<<"l="<<l<<std::endl;
		                                 cellNormalTensor[l][0] += wx*wx/area_tmp1;//n11
		                                 cellNormalTensor[l][1] += wx*wy/area_tmp1;//n12
		                                 cellNormalTensor[l][2] += wx*wz/area_tmp1;//n13
		                                 cellNormalTensor[l][3] += wy*wy/area_tmp1;//n22
		                                 cellNormalTensor[l][4] += wy*wz/area_tmp1;//n23
		                                 cellNormalTensor[l][5] += wz*wz/area_tmp1;//n33
		                                 //
		                                 cellNormalAreaTensor[l][0] += 0.5*wx*wx/area_tmp2;//n11
		                                 cellNormalAreaTensor[l][1] += 0.5*wx*wy/area_tmp2;//n12
		                                 cellNormalAreaTensor[l][2] += 0.5*wx*wz/area_tmp2;//n13
		                                 cellNormalAreaTensor[l][3] += 0.5*wy*wy/area_tmp2;//n22
		                                 cellNormalAreaTensor[l][4] += 0.5*wy*wz/area_tmp2;//n23
		                                 cellNormalAreaTensor[l][5] += 0.5*wz*wz/area_tmp2;//n33
		                             //}

		                             area += area_tmp2;
		                         }
		                         cellSurfaceArea[l] += 0.5*area;
		                         //caculate the unit normal vector of the facet

		                     }
		                 }

		             }
		         }
		         while (cla.inc());
	  }
	  //pwlist.push_back(pw);//we should dump all data out as soon as possible to save memory
	  std::cout<<std::endl;
	  std::cout<<"writing data.."<<std::endl;//std::cout<< "add polywriter successfully! size of pwlist:	" << pwlist.size()<<std::endl;
		//output
		//output polydata if applicable
	//	if(state["isBatch"] == false)//output polydata when running one job
	//	{
			// save point pattern output
			if(state["savereduced"] == true)
			{
					pw.savePointPatternForGnuplot(folder + "/reduced.xyz");
			}

			std::cout << std::endl;
			//remove duplicates and label back indices
			if(state["removeduplicate"] == true){
				pw.removeduplicates(epsilon, xmin, xmax, ymin, ymax, zmin, zmax);
			}
			std::cout << std::endl;
			bool withboundary = true;
			withboundary = state["withboundary"];
			//write poly vtk file
			if(state["savevtk"] == true)
			{
				double rr = state["rr"];
					pw.savePolyVTK(folder + "/cell.vtk",xmin,xmax,ymin,ymax,zmin,zmax,rr,withboundary);
			}
			//pw.saveOnePolyVTK("onecell.vtk",14);//for test
			//write vtk files of individual cells
			if(state["cellVTK"] == true)
			{
				std::string cellIDsfile = state["cellIDsfile"];
				if(cellIDsfile.empty())
					{
							std::cerr << "The file of cellIDs is not valid." << std::endl;
							//return -1;
					}
				pw.outputCellVTK(cellIDsfile,folder+"/onecell");
			}

			//write pov file
			if(state["savepov"] == true)
			{
				double rr = state["rr"];
					pw.savePolyPOV(folder + "/cell.pov",xmin,xmax,ymin,ymax,zmin,zmax,rr);
			}
			// Write poly file for karambola
			if(state["savepoly"] == true)
			{
	//		    std::cout << "writing poly file" << std::endl;
					std::ofstream file;
					std::string polypath =folder + "/cell.poly";
					file.open(polypath.c_str(),std::ios::out);
					file << pw;
					file.close();
			}

 std::cout << std::endl;
 con.clear();
 }
//}
//void output(State& state,std::string posfile,std::string sequenceNum,std::string folder,std::vector<pointpattern> pplist,std::vector<polywriter>pwlist
//		){
    if (state["savesurface"] == true)
    {
        std::cout << "save point pattern file" << std::endl;
        std::ofstream file;
        std::string surfacepath = folder + "/surface_triangulation.xyz";
        file.open(surfacepath.c_str(),std::ios::out);
        for(unsigned int i = 0 ;i<pplist.size();i++){
        	file << pplist.at(i);
        }
        file.close();
    }

    //output local porosity
    std::string volareafile = state["outputfile"];
    std::ofstream fp2;
    std::string path2 = folder + "/"+volareafile + sequenceNum+".txt";
    fp2.open(path2.c_str(),std::ios::out);
    //FILE *fp2=fopen(folder + "/"+volareafile+".txt","w");
    fp2 << "#cellVolume particleVolume cellSurfaceArea particleSurfaceArea normalTensor(1-6) normalAreaTensor(1-6)" << std::endl;
    //fprintf(fp2,"#cellVolume particleVolume cellSurfaceArea particleSurfaceArea\n");
    for(std::map < unsigned long long, double >::iterator  iter = cellVolumes.begin(); iter != cellVolumes.end(); iter++) {
        // iter->first = key
        // iter->second = value
        if(state["RawData"]){
            //to do
            //just print cell volumes
            //fprintf(fp2,"%g \n",cellVolumes[iter->first]);
            fp2 <<  std::setprecision(12) << cellVolumes[iter->first] << std::endl;

        }else{//print porosity
            //porosity = 1.0 - particleVolumes[iter->first]/cellVolumes[iter->first];
            //fprintf(fp2,"%g %g %g %g\n",cellVolumes[iter->first],particleVolumes[iter->first],cellSurfaceArea[iter->first],particleSurfaceArea[iter->first]);
            std::vector<double> cn = cellNormalTensor[iter->first];
            std::vector<double> cna = cellNormalAreaTensor[iter->first];
            fp2 <<  std::setprecision(12) << cellVolumes[iter->first] << " " << particleVolumes[iter->first] << " " << cellSurfaceArea[iter->first]<< " " <<particleSurfaceArea[iter->first]
            << " " <<cn[0]<< " " <<cn[1]<< " " <<cn[2]<< " " <<cn[3]<< " " <<cn[4]<< " " <<cn[5]
            << " " <<cna[0]<< " " <<cna[1]<< " " <<cna[2]<< " " <<cna[3]<< " " <<cna[4]<< " " <<cna[5]
             << std::endl;

        }
    }
    //fclose(fp2);
    fp2.close();


//	}
    //#endif


}

void process_voronoi(State& state,std::string posfile,std::string wallfile,std::string sequenceNum){

    std::string folder = state["outputFolder"];
    if(folder.empty())
    {
        std::cerr << "outfilepath is not valid" << std::endl;
        return;
    }
    //std::cout<<"wallfile: "<<wallfile.empty()<<std::endl;
    //std::string readfile = state["readfile"];

    std::cout << "Parsing Position File... \nWorking on " << posfile << std::endl;

    // pp contains the triangulation of the particle surfaces


    std::vector<pointpattern> pplist; //all particles

    if(state["RawData"]){//directly input surface points of particles
        std::string posDatasetDir = state["posDatasetDir"];
        std::string posFilePrefix = state["posFilePrefix"];
        std::string posFileSuffix = state["posFileSuffix"];
        const unsigned int particleNum = state["particleNum"];
        for(unsigned int i=0;i<particleNum;i++){//reading particles one by one
            //
        	pointpattern pp;
            loader.readRawParticle(con_fname(i, posDatasetDir, posFilePrefix, posFileSuffix), pp, i);
	    std::cout<< posDatasetDir << posFilePrefix << i << posFileSuffix <<std::endl;

        }


    }else{
        // read particle parameters and positions
        std::vector<particleparameterset> setlist;
		//judge whether the file can be open or not
        bool flag = loader.read(posfile,setlist);
        if (!flag){return;}//if the current posfile can not be found, then skipt it.
        std::cout << "Creating Surface Triangulation... " << std::flush;

        const int w_slices = state["w_slices"];
        const int h_slices = state["h_slices"];
        const double scale =  state["scale"];

        // scope for the readstate to ensure it won't lack out to anything else
        {
            // create a readstate that translates the particle parameters to surface shapes
            //State readstate {true};
            //readstate["pointpattern"].SetClass<pointpattern> ("addpoint", &pointpattern::addpoint );
            //readstate.Load(readfile);

            unsigned int ids = 0;
            //by chris

            for(auto it = setlist.begin(); it != setlist.end(); ++it )
            {
				pointpattern pp;
                // read one particle from the position file
                particleparameterset set = (*it);
                // put all parameters for this particle to the lua readstate
                //for(unsigned int i = 0; i != set.parameter.size(); ++i)
                //{   std::cout << set.get(i) << std::endl;
                    //readstate["s"][i] = set.get(i);
               // }

                    cellVolumes[ids] = 0.0;
                    cellSurfaceArea[ids] = 0.0;
                    double area=0,volume=0;

                    extractSuperquadric(pp, ids, scale,set.parameter,w_slices,h_slices,area,volume);

                    std::cout << "points created: " << pp.points.size() << std::endl << std::endl;
                    //store all particle in a list pplist
                    pplist.push_back(pp);
                    particleVolumes[ids] = volume;
                    particleSurfaceArea[ids] = area;
                    //surface tensors
                    std::vector<double> normalComp={0.,0.,0.,0.,0.,0.};//six components
                    std::vector<double> normalAreaComp={0.,0.,0.,0.,0.,0.};//six components
                    cellNormalTensor[ids] = normalComp;
                    cellNormalAreaTensor[ids] = normalAreaComp;
										ids ++;//particle id from 0. FIXME:We need read the particle id from the local file. Here the ids are sequencially numbered.
            }
        }
    }//end parameterized particle surface reading

    std::cout << "finished!" << std::endl;
    const double epsilon = state["epsilon"];
    std::cout << std::endl;

    std::cout << "\nDone! :) "<< std::endl;

//call local_voronoi
    double radiation_radius = 4;
local_voronoi(state,pplist,wallfile,radiation_radius,epsilon,sequenceNum);

}
int main (int argc, char* argv[])
{
	std::cout << "DEBUG: accessfile() called by process " << ::getpid() << " (parent: " << ::getppid() << ")" << std::endl;
	std::this_thread::sleep_for(std::chrono::seconds(10));//sleep for gdb debug
    // command line argument parsing
    if(argc != 2 )
    {
        std::cerr << "Commandline parameters not correct .... aborting "  << std::endl;
        std::cerr << std::endl <<  "Use pomelo this way:\n\t ./pomelo [path-to-lua-file] "  << std::endl;
        return -1;
    }

    std::cout << "\nP O M E L O\n\nVersion " << version << "\nCopyright (C) 2016\nSimon Weis and Philipp Schoenhoefer\n\n";
    std::cout << "The current version is enhanced by Sway Zhao.\n";
    std::cout << "pomelo home:\t\t http://www.theorie1.physik.uni-erlangen.de/" << std::endl << std::endl << std::endl;


    const std::string filename = argv[1];
    //std::string folder = argv[2];


    // lua state for the global parameter file
    //Not used in the enhanced version
    State state {true};
    if(!state.Load(filename))
    {
        std::cerr << "error loading lua parameter file: " << filename << std::endl;
        return -1;
    }
    std::string folder = state["outputFolder"];
    std::cout << "command line arguments parsed:\nLUA Parameter File: " << filename << "\noutfolder: " << folder << std::endl << std::endl;

    // command line sanity check
    /*
    if(folder.empty())
    {
        std::cerr << "outfilepath is not valid" << std::endl;
        return -1;
    }
    */
    // sanitize folder path and create folder
    //char lastCharOfFolder = *folder.rbegin();
    //if (lastCharOfFolder != '/')
    //    folder += '/';
    //mkdir(folder.c_str(),0755);
    bool isBatch = state["isBatch"];
    if (isBatch){
        std::string DatasetDir = state["DatasetDir"];
        std::string posFilePrefix = state["posFilePrefix"];
        std::string posFileSuffix = state["posFileSuffix"];
        std::string wallFilePrefix = state["wallFilePrefix"];
        std::string wallFileSuffix = state["wallFileSuffix"];
        const int jobNumInterval = state["jobNumInterval"];
        const int jobNum = state["jobNum"];
        if(DatasetDir.empty())
        {
            std::cerr << "outfilepath is not valid" << std::endl;
            return -1;
        }
         std::string posfile, wallfile;
        for(int i=0;i<jobNum;i++){//process jobs in bach
            //
            std::cout<<"Job "<<i+1<<" of "<<jobNum<< " begins..."<<std::endl;
            posfile = DatasetDir + posFilePrefix + std::to_string(i*jobNumInterval) + posFileSuffix;
            wallfile = DatasetDir + wallFilePrefix + std::to_string(i*jobNumInterval) + wallFileSuffix;
	//process_voronoi(state,posfile,wallfile,std::to_string(i*jobNpostprocessingumInterval));
            process_voronoi(state,posfile,wallfile,std::to_string(i*jobNumInterval));
            std::cout<<"Job "<<i+1<<" of "<<jobNum<< " finished!"<<std::endl;
        }
    }else{//process only one job
        // parse global parameters from lua file
        std::string posfile = state["positionfile"];
        std::string wallfile = state["wallfile"];
		//the most important step is the "process_voronoi"

        process_voronoi(state,posfile,wallfile,char2string(""));
    }
    return 0;
}
