/*
 * =====================================================================================
 *
 *       Filename:  CellFactory.cpp
 *
 *    Description:  a wrapper of Class CellMachine exposed in Python.
 *
 *        Version:  1.0
 *        Created:  12/22/2018
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Shiwei Zhao (swzhao (at) scut.edu.cn)
 *   Organization: The Hong Kong University of Science and Technology
 *                 South China University of Technology
 *
 * =====================================================================================
 */

#include "Superquadrics.hpp"
#include "CellFactory.hpp"
#include "CellMachine.hpp"
#include <pthread.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>
#include <experimental/filesystem> //use experimental for gcc below 8.0

#include <iterator>
#include <regex>

namespace fs = std::experimental::filesystem;

CellFactory::CellFactory(){
  searchRadius = 4.0;
  parShrink = 0.01e-3;
  scale = 1000.0;
  boxScale = 2.0;
  verbose = 7;//output all info
  //w_slices = 20;
  //h_slices = 20;
  cellVTK = false;
  cellPOV = false;
  visualized_ids.clear();
  threadNum = 2;
  //check and create folders
  //checkCreateFolder(in_folder);
  //checkCreateFolder(out_folder);
  //checkCreateFolder(out_folder+"/tmp/");

}
bool CellFactory::checkCreateFolder(std::string target){
  if(target.length() == 0){std::cout<<"Directory is not created!"<<std::endl;return false;}
  if(! fs::exists(target)){// the target folder doest not exist, and create it.
    bool ret = fs::create_directory(target);
    if(ret){std::cout<<"Folder "<<target<<" is created."<<std::endl;}
    else{std::cout<<"Failed to create Folder "<<target<<std::endl;}
    return ret;
  }
  else{
    std::cout<<"Folder "<<target<<" exists!"<<std::endl;
    return true;
  }
}

bool CellFactory::isInVisualIds(int id){
  bool flag = false;
  if (visualized_ids.size()==0) return true;
  for (std::vector<int>::iterator it = visualized_ids.begin() ; it != visualized_ids.end(); ++it){
    if (id == (*it)) flag = true;
  }
  return flag;
}
void CellFactory::genPointClouds(double w_slices=20,double h_slices = 20){
  std::vector<std::vector<double> > setlist;
  bool flag = read(setlist);
  if (!flag){std::cerr<<"Error: the file of particle positons has not been loaded correctly!"<<std::endl;}//if the current posfile can not be found, then skipt it.
  unsigned int ids = 0;
  //preprocessing data
  std::ofstream fp1;
  std::string path1 = out_folder + "/parproperties.txt";
  fp1.open(path1.c_str(),std::ios::out);
  fp1 << "#particleID particleVolume particleSurfaceArea" << std::endl;
  for(auto it = setlist.begin(); it != setlist.end(); ++it )
  {
    // read one particle from the position file
    std::vector<double> set = (*it);
    double area=0,volume=0;
    std::string outfile = in_folder + "/"+std::to_string(ids)+".dat";//store point clouds of particles
    particleAttr pa;
    bool ret;
    ret = pointCloud_Superquadric(ids, outfile, parShrink,set,w_slices,h_slices,area,volume,pa);
    if(ret){
      parAttrlist.push_back(pa);
      fp1 << ids << "\t"<<volume<<"\t" << area << std::endl;
      ids ++;
    }else{
      //warning
    }
  }
  std::cout << "point cloud created!" << std::endl;
}

bool CellFactory::read( std::vector<std::vector<double> >& set)
{
    std::ifstream infile;
    infile.open(posFile.data(),std::ifstream::in);
    if (infile.fail())
    {
        std::cout << "Cannot load file " << posFile << std::endl;
        return false;
    }
#pragma GCC diagnostic ignored "-Wwrite-strings"
    std::string line("");
    unsigned int linesloaded = 0;
    //std::getline(infile, line);
    while (std::getline(infile, line))
    {
        if(line.find("#")!=std::string::npos) continue; // ignore comment lines

        std::vector<double> p;

        const std::regex ws_re("\\t+"); // table
        auto line_begin = std::sregex_token_iterator(line.begin(), line.end(), ws_re, -1);
        auto line_end = std::sregex_token_iterator();
        
        for (auto it = line_begin; it != line_end; ++it ){
            double d = atof((*it).str().c_str());
            p.push_back(d);
        }
        linesloaded++;
        set.push_back(p);

    }
    std::cout << "Lines loaded: " << linesloaded << std::endl << std::endl;

    infile.close();
    return true;
}


void CellFactory::readRawParticle(std::string filename, pointpattern& pp, unsigned int particleid)
{
    std::ifstream infile;
    infile.open(filename.data(),std::ifstream::in);
    if (infile.fail())
    {
        std::cout << "Cannot load file " << filename << std::endl;
        return;
    }
#pragma GCC diagnostic ignored "-Wwrite-strings"
    std::string line("");
    unsigned int linesloaded = 0;
    std::getline(infile, line);
    while (std::getline(infile, line))
    {
        if(line.find("#")!=std::string::npos) continue; // ignore comment lines

        std::vector<double> xyz;
        //xyz point
        const std::regex ws_re("\\t+"); // table
        auto line_begin = std::sregex_token_iterator(line.begin(), line.end(), ws_re, -1);
        auto line_end = std::sregex_token_iterator();
        
        for (auto it = line_begin; it != line_end; ++it ){
            double d = atof((*it).str().c_str());
            xyz.push_back(d);
        }
        //check the legality of xyz data
        if(3 > xyz.size())
        {
            std::cout << "Warning!! The data at Line " << linesloaded <<" in the dataset has wrong format. I regected it for robust running."<< std::endl << std::endl;
        }else
        {
            pp.addpoint(particleid,xyz[0],xyz[1],xyz[2]);
        }
        linesloaded++;


    }
    std::cout << "Lines loaded: " << linesloaded << std::endl << std::endl;

    infile.close();
}


void CellFactory::neighborSearch(void){
  //#ifdef CF_OPENMP
  //  #pragma omp parallel for schedule(dynamic,1) num_threads(threadNum)
  //#endif
  for(int i = 0 ;i< parAttrlist.size();i++){//mayby using pointer is faster
//		particleAttr p1;
    double dist = 0;
    for(int j = 0 ; j < parAttrlist.size();j++){
      //FIXME:using AABB to speed up?
      if(j == i) continue;//
      particleAttr p2  = parAttrlist.at(j);
      dist = sqrt(pow((parAttrlist.at(i).centerx-p2.centerx),2)+pow((parAttrlist.at(i).centery-p2.centery),2)+pow((parAttrlist.at(i).centerz-p2.centerz),2));
      if(parAttrlist.at(i).radius * searchRadius > dist ){
        //p2 close to p1
        parAttrlist.at(i).surroundedID.push_back(p2.ID); // find all particle close to p1 not including current particle itself;
      }
    }
  }
}

void CellFactory::processing(void){
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  //here we can use openMP (or MPI) for parallel computation
  #ifdef CF_OPENMP
  //omp_set_dynamic(0);     // Explicitly disable dynamic teams
  //omp_set_num_threads(2); // Use 4 threads for all consecutive parallel regions
  //#pragma omp parallel
  //std::cout<<"max number "<<omp_get_max_threads()<<" thread(s)"<<std::endl;
  //omp_set_num_threads(2);
  //setenv("OMP_NUM_THREADS", "1", 2);
    #pragma omp parallel for schedule(dynamic,1) num_threads(threadNum)
  #endif
    for(int i = 0 ;i< parAttrlist.size();i++){
    //#endif
      CellMachine CM = CellMachine(in_folder,out_folder);
      if(isInVisualIds(parAttrlist[i].ID)){//output all cells for empty visual list, or check it the specified cell should be output.
        CM.set_cellVTK(cellVTK);
        CM.set_cellPOV(cellPOV);
      }
      CM.set_scale(scale);
      CM.set_verbose(verbose);
      CM.set_boxScale(boxScale);
      CM.set_wallFile(wallFile);
      //loading point clouds from local files
      CM.pushPoints(parAttrlist[i]);
      //comupute cells
      CM.processing();
    }
  //}
  //clear temporary files
  #ifdef CF_OPENMP
  std::ofstream outfile;
  std::string outpath = out_folder+"/cellProperties.dat";
  outfile.open(outpath.c_str());
  outfile << "#id cellVolume cellSurfaceArea normalizedNormalAreaTensor(n11,n12,n13,n22,n23,n33)" << std::endl;
  for(int i = 0 ;i< parAttrlist.size();i++){
    std::ifstream infile;
    std::string path = out_folder + "/tmp/"+std::to_string(i)+"cellProperties.tmp";
    infile.open(path.c_str(),std::ios::in);
    if (infile.fail()){std::cout << "Cannot load file "<< std::endl;continue;}
    #pragma GCC diagnostic ignored "-Wwrite-strings"
    std::string line("");
    std::getline(infile, line);
    outfile << line <<std::endl;;
    /*while (std::getline(infile, line))
    {
        if(line.find("#")!=std::string::npos) continue; // ignore comment lines
        outfile << line;
    }*/
    //std::cout << "Lines loaded: " << linesloaded << std::endl << std::endl;
    infile.close();
  }
  outfile.close();
  std::uintmax_t n = fs::remove_all(out_folder+"/tmp/");
  std::cout << "Deleted " << n << " files or directories\n";
  #endif
  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count();
  std::cout<<"Time ellapsed is:"<<duration<<" seconds"<<std::endl;
}

void CellFactory::processingOne(unsigned int pid){
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  //here we can use openMP (or MPI) for parallel computation
  if(pid>=parAttrlist.size()){std::cerr << "the particle id is not valid!" << std::endl;return;}
  //#endif
    CellMachine CM = CellMachine(in_folder,out_folder);
    if(isInVisualIds(parAttrlist[pid].ID)){//output all cells for empty visual list, or check it the specified cell should be output.
      CM.set_cellVTK(cellVTK);
      CM.set_cellPOV(cellPOV);
    }
    CM.set_scale(scale);
    CM.set_verbose(verbose);
    CM.set_boxScale(boxScale);
    CM.set_wallFile(wallFile);
    //loading point clouds from local files
    CM.pushPoints(parAttrlist[pid]);
    //comupute cells
    CM.processing();
  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count();
  std::cout<<"Time ellapsed is:"<<duration<<" seconds"<<std::endl;
}

void CellFactory::autoWorkFlow(void){
  #ifdef CF_DEBUG
  //std::cout << "DEBUG: accessfile() called by process " << ::getpid() << " (parent: " << ::getppid() << ")" << std::endl;
	//std::this_thread::sleep_for(std::chrono::seconds(10));//sleep for gdb debug
  #endif
  //generate point clouds
  genPointClouds();
  //get neighbor list
  neighborSearch();
  //process
  processing();
}

bool CellFactory::pointCloud_Superquadric(unsigned int id, std::string outfile, double scaledist,std::vector<double> &set, int w_slices, int h_slices, double& area, double& volume, particleAttr& pattr){
				double rx1,ry1,rz1,rx2,ry2,rz2;//,rmin;
        double eps1,eps2;
        Vector3r Position;
        Quaternionr Ori;
				bool polysuper = false;
        bool sphere = false;
				if(set.size()==15){
					polysuper = true;
				}else if(set.size()==12){
					polysuper = false;
        }else if(set.size()==4){//sphere
          polysuper = false;
          sphere = true;
				}else{
					std::cerr<<"Error: the column number of the position file is not correct!"<<std::endl;
					return false;
				}
        if(sphere){//processing a sphere
          rx1 = set[0];
          Position = Vector3r(set[1], set[2], set[3]);
          pattr.ID = id;
  				pattr.centerx = Position[0];
  				pattr.centery = Position[1];
  				pattr.centerz = Position[2];
  				double xmin(1e10),xmax(-1e10),ymin(1e10),ymax(-1e10),zmin(1e10),zmax(-1e10);
  				pattr.radius = rx1;
          int i,j,w=w_slices,h=h_slices;
          double a=0.0,b=0.0,phi0,phi1;
          double hStep=M_PI/(h-1);
          double wStep=2*M_PI/w;
          Quaternionr Orisp;
          //random orientation
          double heading = double(rand())/RAND_MAX*2.0*M_PI;//rotate around z
          double attitude = double(rand())/RAND_MAX*2.0*M_PI; //y
          double bank = double(rand())/RAND_MAX*2.0*M_PI; //y
          double c1 = cos(heading/2);
          double s1 = sin(heading/2);
          double c2 = cos(attitude/2);
          double s2 = sin(attitude/2);
          double c3 = cos(bank/2);
          double s3 = sin(bank/2);
          double c1c2 = c1*c2;
          double s1s2 = s1*s2;

          Orisp.w() =c1c2*c3 - s1s2*s3;
          Orisp.x() =c1c2*s3 + s1s2*c3;
          Orisp.y() =s1*c2*c3 + c1*s2*s3;
          Orisp.z() =c1*s2*c3 - s1*c2*s3;

          Matrix3r rot = Orisp.toRotationMatrix();
          Vector3r p,n;
  				std::ofstream fp;
          double shrinkR = rx1 - scaledist;
          if(shrinkR<=0){shrinkR = 0.5*rx1;}//scaledist is too large
  				fp.open(outfile.c_str(),std::ios::out|std::ios::binary);//write data into a binary file instead of a text file for a better performance in reading.
  				//pointpattern pp;
  		    //caution:the two polar points should be degenerated.
          //std::ofstream fp3;
          //std::string outfile1 = "./test0.xyz";
          //fp3.open(outfile1.c_str(),std::ios::out|std::ios::app);
          for(a=hStep,i=0;i<h-2;i++,a+=hStep)
          {
            for(b=0.0,j=0;j<w;j++,b+=wStep)
            {     phi0 = b;
                  phi1 = a-M_PI_2l;

            //get surface point
          	double xyz[3];
            Vector3r p = shrinkR*Vector3r(cos(phi0)*cos(phi1),sin(phi0)*cos(phi1),sin(phi1));
            p = rot*p;
          	xyz[0] = p[0] + Position[0];
          	xyz[1] = p[1] + Position[1];
          	xyz[2] = p[2] + Position[2];

          	//p = Position + Vector3r(x,y,z);

              //pp.addpoint(0,p(0),p(1),p(2));
            //if(id==0) fp3 <<std::scientific<< p(0)<<"\t" << p(1)<<"\t" << p(2) <<std::endl;
          	//fp <<std::scientific<< p(0)<<"\t" << p(1)<<"\t" << p(2) <<std::endl;
            fp.write( (char*)&xyz, sizeof xyz);

           }
          }
  		    //two polar points
          for(a=0.0,i=0;i<2;i++,a+=M_PI)
          {
            phi0 = 0;
          	phi1 = a-M_PI_2l;

            //get surface point
            double xyz[3];
            Vector3r p = rot*Vector3r(0,0,shrinkR*sin(phi1));
            xyz[0] = p[0] + Position[0];
            xyz[1] = p[1] + Position[1];
            xyz[2] = p[2] + Position[2];
            //p = Position + Vector3r(x,y,z);
            //pp.addpoint(0,p(0),p(1),p(2));
            //if(id==0) fp3 <<std::scientific<< p(0)<<"\t" << p(1)<<"\t" << p(2) <<std::endl;
          	//fp <<std::scientific<< p(0)<<"\t" << p(1)<<"\t" << p(2) <<std::endl;
            fp.write( (char*)&xyz, sizeof xyz);
            //if (id==0) std::cout<<std::scientific<< xyz[0]<<"\t" << xyz[1]<<"\t" << xyz[2] <<std::endl;
          }
          //we could use explicit function to get the AABB
          pattr.xmin = Position[0] - rx1;
          pattr.xmax = Position[0] + rx1;
          pattr.ymin = Position[1] - rx1;
          pattr.ymax = Position[1] + rx1;
          pattr.zmin = Position[2] - rx1;
          pattr.zmax = Position[2] + rx1;
          // it may be not necessary to store xrange, yrange and zrange.
          pattr.xrange = rx1;
          pattr.yrange = rx1;
          pattr.zrange = rx1;
          //fp3.close();//test
  				fp.close();
  		      //pattr: updated over

        	volume = 4.0/3.0*M_PI*pow(rx1,3);     //volume = 2a1a2a3ep1ep2B(ep1/2+1,ep1)B(ep2/2,ep2/2) see the reference
            //surface area
          area = 4.0*M_PI*pow(rx1,2);//2500 discretized points on the surface
        	return true;
        }//sphere processing ends


				if(polysuper){
					rx1 = set[0];
	        rx2 = set[1];
	        ry1 = set[2];
					ry2 = set[3];
	        rz1 = set[4];
	        rz2 = set[5];
	        eps1 = set[6];
	        eps2 = set[7];//
	        Position = Vector3r(set[8], set[9], set[10]);
					Ori.w() = set[11];
					Ori.x() = set[12];
					Ori.y() = set[13];
					Ori.z() = set[14];
				}else{
					rx1 = rx2 = set[0];
	        ry1 = ry2 = set[1];
	        rz1 = rz2 = set[2];
	        eps1 = set[3];
	        eps2 = set[4];//
	        Position = Vector3r(set[5], set[6], set[7]);
					Ori.w() = set[8];
					Ori.x() = set[9];
					Ori.y() = set[10];
					Ori.z() = set[11];
				}
				PolySuperellipsoid PS = PolySuperellipsoid(rx1,rx2,ry1,ry2,rz1,rz2,eps1,eps2);

        Matrix3r A = Ori.toRotationMatrix();
				//update particle attribution
				pattr.ID = id;
				pattr.centerx = Position[0];
				pattr.centery = Position[1];
				pattr.centerz = Position[2];
				double xmin(1e10),xmax(-1e10),ymin(1e10),ymax(-1e10),zmin(1e10),zmax(-1e10);

        //find rmin
        //rmin = (rx<ry)?rx:ry;
        //rmin = (rmin<rz)?rmin:rz;
				//find rmax
				pattr.radius = PS.getr_max();
        int i,j,w=w_slices,h=h_slices;
        double a=0.0,b=0.0,phi0,phi1;
        double hStep=M_PI/(h-1);
        double wStep=2*M_PI/w;
				//double scale = 0.95;

        Vector3r p,n;
        //std::cout<<"scale1="<<scale<<std::endl;
        //double scaledist=0.0;
        //scaledist = scale*rmin;
				std::ofstream fp;
				fp.open(outfile.c_str(),std::ios::out|std::ios::binary);
				//pointpattern pp;
		//caution:the two polar points should be degenerated.
        for(a=hStep,i=0;i<h-2;i++,a+=hStep)
        {
          for(b=0.0,j=0;j<w;j++,b+=wStep)
          {     phi0 = b;
                phi1 = a-M_PI_2l;

	        //get surface point
					p = PS.getSurfaceMC(Vector2r(phi0,phi1));
						//get out-ward normal at the surface point
					n = PS.getNormal(Vector2r(phi0,phi1));
            n.normalize();
            //std::cout<<"pnorm="<<p.norm()<<std::endl;
            //std::cout<<"rmin="<<rmin<<std::endl;
            //scale = 1.0 - scaledist/p.norm();
            //std::cout<<"scale="<<scale<<std::endl;
            //p = p*scale;
            p = p - n*scaledist;
            p = A*p + Position;
            double xyz[3];
            xyz[0] = p(0);xyz[1] = p(1);xyz[2] = p(2);
            fp.write( (char*)&xyz, sizeof xyz);
						//fp <<std::scientific<< p(0)<<"\t" << p(1)<<"\t" << p(2) <<std::endl;
						if(p(0)<xmin) xmin = p(0);
						if(p(1)<ymin) ymin = p(1);
						if(p(2)<zmin) zmin = p(2);
						if(p(0)>xmax) xmax = p(0);
						if(p(1)>ymax) ymax = p(1);
						if(p(2)>zmax) zmax = p(2);
         }
      }
		//two polar points
		  for(a=0.0,i=0;i<2;i++,a+=M_PI)
      {
		    phi0 = 0;
      	phi1 = a-M_PI_2l;

        //get surface point
				p = PS.getSurfaceMC(Vector2r(phi0,phi1));
          //get out-ward normal at the surface point
				n = PS.getNormal(Vector2r(phi0,phi1));
        n.normalize();
        //std::cout<<"pnorm="<<p.norm()<<std::endl;
        //std::cout<<"rmin="<<rmin<<std::endl;
        //scale = 1.0 - scaledist/p.norm();
        //std::cout<<"scale="<<scale<<std::endl;
        //p = p*scale;
        p = p - n*scaledist;
        p = A*p + Position;
        double xyz[3];
        xyz[0] = p(0);xyz[1] = p(1);xyz[2] = p(2);
        fp.write( (char*)&xyz, sizeof xyz);
        //pp.addpoint(0,p(0),p(1),p(2));
				//fp <<std::scientific<< p(0)<<"\t" << p(1)<<"\t" << p(2) <<std::endl;
				if(p(0)<xmin) xmin = p(0);
				if(p(1)<ymin) ymin = p(1);
				if(p(2)<zmin) zmin = p(2);
				if(p(0)>xmax) xmax = p(0);
				if(p(1)>ymax) ymax = p(1);
				if(p(2)>zmax) zmax = p(2);
      }
				 //we could use explicit function to get the AABB
				 pattr.xmin = xmin;
				 pattr.xmax = xmax;
				 pattr.ymin = ymin;
				 pattr.ymax = ymax;
				 pattr.zmin = zmin;
				 pattr.zmax = zmax;
				 // it may be not necessary to store xrange, yrange and zrange.
				 pattr.xrange = xmax - xmin;
				 pattr.yrange = ymax - ymin;
				 pattr.zrange = zmax - zmin;
         //pattr.nearBoundary = false;//by default
				 /*std::cout << "polywriter: remove duplicates" << std::endl;
				 duplicationremover d(16,16,16);
				 d.setboundaries(xmin, xmax, ymin, ymax, zmin, zmax);
				 std::cout << "\tadding"<<pp.points.size()<<" points" << std::endl;
				 d.addPoints(pp, false);
				 std::cout << "\tremoving duplicates" << std::endl;
				 double epsilon = 1e-6;
				 d.removeduplicates(epsilon);
				 d.getallPoints(pp);
				 std::cout << "\tget back "<<pp.points.size()<<" points" << std::endl;
				 //write a file
				 for(unsigned int i = 0;i<pp.points.size();i++){
					 fp << pp.points.at(i).x<<"\t" << pp.points.at(i).y<<"\t" << pp.points.at(i).z<<std::endl;
        }*/
				fp.close();
		//pattr: updated over

	volume = PS.getVolume();     //volume = 2a1a2a3ep1ep2B(ep1/2+1,ep1)B(ep2/2,ep2/2) see the reference
    //surface area
  area = PS.getSurfaceArea(50,50);//2500 discretized points on the surface
	return true;
}
