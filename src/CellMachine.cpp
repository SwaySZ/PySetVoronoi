/*
 * =====================================================================================
 *
 *       Filename:  CellMachine.cpp
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  12/22/2018
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Shiwei Zhao (zhswee@gmail.com)
 *   Organization: The Hong Kong University of Science and Technology
 *                 South China University of Technology
 *
 * =====================================================================================
 */
#include "CellMachine.hpp"
#include "duplicationremover.hpp"

void CellMachine::initial(){
  cid = -1;
  pid = 0;
  rr = 0.01;
  cellVolume = 0.0;
  cellSurfaceArea = 0.0;
  for(int i=0;i<6;i++) cellNormalTensor[i]=0.0;
  for(int i=0;i<6;i++) cellNormalAreaTensor[i]=0.0;
  labelidmap.clear();
}
CellMachine::CellMachine(){
	xpbc = ypbc = zpbc = false;//using non-periodic boundary by default
  savereduced =false;
  removeduplicate = false;
  withboundary = false;
  savevtk = false;
  cellVTK = false;
  savepov = false;
  savepoly = false;
  delta = 0.1e-3;//
  scale = 1000.0;//scale up 3 order of magnitude. Thus, an input unit of meter yields an output unit of millimeter. However, for the sake of consistence, we scale down the results for writing out.
  pcon = NULL;
  //pp = NULL;
  initial();
  blockMem = 32;
  nx = ny = nz =10;//we try to guess a better value for each of them
	//to do,
}
CellMachine::CellMachine(std::string input_folder,std::string output_folder){
	xpbc = ypbc = zpbc = false;//using non-periodic boundary by default
  in_folder = input_folder;
  out_folder = output_folder;
  wallFile = input_folder+"Walls.dat";
  savereduced =false;
  removeduplicate = false;
  withboundary = false;
  savevtk = false;
  cellVTK = false;
  savepov = false;
  savepoly = false;
  delta = 0.1e-3;//
  boxScale = 2.0;
  scale = 1000.0;//scale up 3 order of magnitude. Thus, an input unit of meter yields an output unit of millimeter. However, for the sake of consistence, we scale down the results for writing out.
  pcon = NULL;
  //pp = NULL;
  initial();
  blockMem = 32;
  nx = ny = nz =10;//we try to guess a better value for each of them
	//to do,
}
void CellMachine::reset(){
  initial();
}
void  CellMachine::getFaceVerticesOfFace( std::vector<int>& f, unsigned int k, std::vector<unsigned int>& vertexList)
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
//test a point is inside the box of this cell machine
bool CellMachine::isInBox(double x, double y, double z){
  bool ret = (x<xmin)||(x>xmax)||(y<ymin)||(y>ymax)||(z<zmin)||(z>zmax);
  return !ret;
}

void CellMachine::readWall(std::string filename)
{
    std::ifstream infile;
    infile.open(filename.data(), std::ifstream::in);
    if (infile.fail())
    {
        std::cout << "Cannot load file " << filename << std::endl;
        return;
    }
#pragma GCC diagnostic ignored "-Wwrite-strings"
    cSplitString line("");
    unsigned int linesloaded = 0;
    //std::getline(infile, line);
    //std::cout<<"Warning! The wallfile should include six lines in order with xmin,xmax,ymin,ymax,zmin,zmax. Each line has three components."<<std::endl;
    std::vector< std::vector<double> > wallpos;
    while (std::getline(infile, line))
    {
        if(line.find("#")!=std::string::npos) continue; // ignore comment lines

        std::vector<std::string> xyzstring = line.split('\t');//
        //std::cout<<xyzstring[0]<<std::endl;
        std::vector<double> xyz;
        //xyz point
        for(unsigned int i =0;i<xyzstring.size();i++){
        	double d = atof(xyzstring.at(i).c_str());
        	xyz.push_back(d);
        }

        //check the legality of xyz data
        if(3 > xyz.size())
        {
            std::cout << "Warning!! The data at Line " << linesloaded <<" in the dataset has wrong format. I regected it for robust running."<< std::endl << std::endl;
        }else
        {
            wallpos.push_back(xyz);
        }
        linesloaded++;

    }
    //check if the wallpos has all data
    if (wallpos.size()!=6){
        std::cout<<"Error: the wallpos has not six-line data."<<wallpos.size()<<std::endl;
    }else{
        //get container size in order
        wall_xmin = wallpos[0][0];
        wall_xmax = wallpos[1][0];
        wall_ymin = wallpos[2][1];
        wall_ymax = wallpos[3][1];
        wall_zmin = wallpos[4][2];
        wall_zmax = wallpos[5][2];
    }
    //std::cout << "Wall position was imported successfully!" << std::endl << std::endl;

    infile.close();
}
void CellMachine::readParticle(std::string filename, bool flag, int particleId)
{   //flag: if true, check inside or outside the box
    std::ifstream infile;
    infile.open(filename.data(),std::ifstream::in);
    //testing
    #ifdef DEBUG_CM
    std::ofstream fp2;
    std::string outfile = "./test.xyz";
    std::cout<<"particleID="<<particleId<<filename<<std::endl;
    if(!flag){
      fp2.open(outfile.c_str(),std::ios::out);
    }else{
      fp2.open(outfile.c_str(),std::ios::out|std::ios::app);
    }

    #endif
    if (infile.fail())
    {
        std::cout << "Cannot load file " << filename << std::endl;
        return;
    }
#pragma GCC diagnostic ignored "-Wwrite-strings"
    cSplitString line("");
    unsigned int linesloaded = 0;
    std::getline(infile, line);
    while (std::getline(infile, line))
    {
        if(line.find("#")!=std::string::npos) continue; // ignore comment lines

        std::vector<std::string> xyzstring = line.split('\t');//
        std::vector<double> xyz;
        //xyz point
        for (unsigned int i =0;i<xyzstring.size();i++){
        	double d = atof(xyzstring.at(i).c_str());
        	xyz.push_back(d);
        }

        //check the legality of xyz data
        if(3 > xyz.size())
        {
            std::cout << "Warning!! The data at Line " << linesloaded <<" in the dataset has wrong format. I regected it for robust running."<< std::endl << std::endl;
        }else
        {
            //pp.addpoint(particleid,);
            //std::cout<<"hre"<<con<<std::endl;
            if(flag){
              if(isInBox(xyz[0],xyz[1],xyz[2])){
                pcon->put(pid, xyz[0]*scale,xyz[1]*scale,xyz[2]*scale);
                //pp->addpoint(particleId,xyz[0],xyz[1],xyz[2]);
                labelidmap.push_back(particleId);
                #ifdef DEBUG_CM
                fp2 << xyz[0]<<"\t" << xyz[1]<<"\t" << xyz[2] <<std::endl;
                #endif
                ++pid;
              }
            }else{
              pcon->put(pid, xyz[0]*scale,xyz[1]*scale,xyz[2]*scale);
              //pp->addpoint(particleId,xyz[0],xyz[1],xyz[2]);
              labelidmap.push_back(particleId);
              #ifdef DEBUG_CM
              fp2 << xyz[0]<<"\t" << xyz[1]<<"\t" << xyz[2] <<std::endl;
              #endif
              ++pid;
            }

        }
        linesloaded++;


    }
    //std::cout << "Lines loaded: " << linesloaded << std::endl << std::endl;

    infile.close();
    //fp2.close();
}

void CellMachine::pushPoints(particleAttr& pAttr){
  cid = pAttr.ID;
  //read wall boundary
  readWall(wallFile);
  //cellVTK
  /*bool ret = ((pAttr.centerx-pAttr.radius)<wall_xmin)||((pAttr.centerx+pAttr.radius)>wall_xmax)||
             ((pAttr.centery-pAttr.radius)<wall_ymin)||((pAttr.centery+pAttr.radius)>wall_ymax)||
             ((pAttr.centerz-pAttr.radius)<wall_zmin)||((pAttr.centerz+pAttr.radius)>wall_zmax);
  if(ret){cellVTK = true;}else{cellVTK=false;}
  */
  //define the box size by considering the global wall
  std::cout<<"scale="<<scale<<" boxScale="<<boxScale<<std::endl;
  xmin = std::max(pAttr.centerx - pAttr.xrange*boxScale, wall_xmin);
  xmax = std::min(pAttr.centerx + pAttr.xrange*boxScale, wall_xmax);
  ymin = std::max(pAttr.centery - pAttr.yrange*boxScale, wall_ymin);
  ymax = std::min(pAttr.centery + pAttr.yrange*boxScale, wall_ymax);
  zmin = std::max(pAttr.centerz - pAttr.zrange*boxScale, wall_zmin);
  zmax = std::min(pAttr.centerz + pAttr.zrange*boxScale, wall_zmax);
  //create a voro container
  //std::cout<<"creating a container..."<<std::endl;
  //std::cout<<pAttr.centerx<<" "<<pAttr.centery<<" "<<pAttr.centerz<<std::endl;
  //std::cout<<"ymin"<<pAttr.centery - pAttr.yrange*boxScale<<std::endl;
  //std::cout<<pAttr.yrange<<pAttr.xrange<<pAttr.zrange<<std::endl;
  //std::cout<<"wallboxsize="<<wall_xmin<<" "<<wall_xmax<<" "<<wall_ymin<<" "<<wall_ymax<<" "<<wall_zmin<<" "<<wall_zmax<<std::endl;
  //std::cout<<"boxsize="<<xmin<<" "<<xmax<<" "<<ymin<<" "<<ymax<<" "<<zmin<<" "<<zmax<<std::endl;
  pcon=new voro::pre_container(xmin*scale, xmax*scale, ymin*scale, ymax*scale, zmin*scale, zmax*scale, xpbc, ypbc, zpbc);
  //con = new voro::container(xmin, xmax, ymin, ymax, zmin, zmax, nx, ny, nz, xpbc, ypbc, zpbc, 32);
  //pp = new pointpattern();
  //std::cout<<"con="<<con<<std::endl;
  //con = con_tmp;
  std::vector<int> surroundedID = pAttr.surroundedID;
	std::cout << "currentparticleID: "<<cid<<"    surroundedID size: "<<surroundedID.size()<<std::endl;
  readParticle(in_folder + "/"+std::to_string(cid)+".dat", false, cid);//read the particle itself
  //read the other surrounding particles
  for(std::vector<int>::iterator it = surroundedID.begin();it!=surroundedID.end();it++){
    //std::cout<<"spid="<<*it<<std::endl;
		readParticle(in_folder + "/"+std::to_string((*it))+".dat", true, (*it));
	}

  /*std::cout << "polywriter: remove duplicates" << std::endl;
  duplicationremover d(16,16,16);
  d.setboundaries(xmin, xmax, ymin, ymax, zmin, zmax);
  std::cout << "\tadding"<<pp->points.size()<<" points" << std::endl;
  d.addPoints(*pp, false);
  std::cout << "\tremoving duplicates" << std::endl;
  double epsilon = 1e-5;
  d.removeduplicates(epsilon);
  d.getallPoints(*pp);
  std::cout << "\tget back "<<pp->points.size()<<" points" << std::endl;

  for(unsigned int i = 0;i<pp->points.size();i++){
    pcon->put(pid, pp->points.at(i).x*1000,pp->points.at(i).y*1000,pp->points.at(i).z*1000);
    ++pid;
    labelidmap.push_back(pp->points.at(i).l);
 }
 delete pp;
 pp = NULL;
 */
  pcon->guess_optimal(nx,ny,nz);
  std::cout<<"nx="<<nx<<" ny="<<ny<<" nz="<<nz<<std::endl;
}
void CellMachine::writeGlobal(){
			std::ofstream fp;
      #ifdef CF_OPENMP
      //std::string path = out_folder + "/"+std::to_string(cid)+"cellProperties.dat";
      std::string path = out_folder +"/tmp/"+std::to_string(cid)+"cellProperties.tmp";
      fp.open(path.c_str(),std::ios::out);
      //fp << "#id cellVolume cellSurfaceArea normalTensor(1-6) normalAreaTensor(1-6)" << std::endl;
      #else
			std::string path = out_folder + "/cellProperties.dat";
      if(cid==0){
        fp.open(path.c_str(),std::ios::out);
        fp << "#id cellVolume cellSurfaceArea normalAreaTensor(n11,n12,n13,n22,n23,n33)" << std::endl;
      }else{
        fp.open(path.c_str(),std::ios::out| std::ios::app);
      }
      #endif
      /*
			fp <<  cid<<" "<<std::scientific << cellVolume/pow(scale,3) << " " << cellSurfaceArea/pow(scale,2)
							<< " " <<cellNormalTensor[0]<< " " <<cellNormalTensor[1]<< " " <<cellNormalTensor[2]
              << " " <<cellNormalTensor[3]<< " " <<cellNormalTensor[4]<< " " <<cellNormalTensor[5]
							<< " " <<cellNormalAreaTensor[0]<< " " <<cellNormalAreaTensor[1]<< " " <<cellNormalAreaTensor[2]
              << " " <<cellNormalAreaTensor[3]<< " " <<cellNormalAreaTensor[4]<< " " <<cellNormalAreaTensor[5]
							<< std::endl;
      */
      fp <<  cid<<" "<<std::scientific << cellVolume/pow(scale,3) << " " << cellSurfaceArea/pow(scale,2)
							<< " " <<cellNormalAreaTensor[0]<< " " <<cellNormalAreaTensor[1]<< " " <<cellNormalAreaTensor[2]
              << " " <<cellNormalAreaTensor[3]<< " " <<cellNormalAreaTensor[4]<< " " <<cellNormalAreaTensor[5]
							<< std::endl;

			fp.close();
}
void CellMachine::writeLocal(polywriter *pw){
	std::cout<<"writing cell data ..."<<std::endl;
  	// save point pattern output
		/*if(savereduced)
		{
				pw->savePointPatternForGnuplot(out_folder + "/reduced.xyz");
		}*/

		std::cout << std::endl;
		//remove duplicates and label back indices
		/*if(removeduplicate){
			pw->removeduplicates(epsilon, xmin, xmax, ymin, ymax, zmin, zmax);
		}*/
		std::cout << std::endl;
		//write poly vtk file
		/*if(savevtk)
		{
				pw->savePolyVTK(out_folder + "/cell.vtk",xmin,xmax,ymin,ymax,zmin,zmax,rr,withboundary);
		}*/
		//pw.saveOnePolyVTK("onecell.vtk",14);//for test
		//write vtk files of individual cells
		if(cellVTK)
		{
        pw->saveOnePolyVTKnew(out_folder+"/"+std::to_string(cid)+".vtk", cid);
		}

		//write pov file
		if(savepov)
		{
				pw->savePolyPOV(out_folder + "/cell.pov",xmin,xmax,ymin,ymax,zmin,zmax,rr);
		}
		// Write poly file for karambola
		if(savepoly)
		{
//		    std::cout << "writing poly file" << std::endl;
				std::ofstream file;
				std::string polypath =out_folder + "/cell.poly";
				file.open(polypath.c_str(),std::ios::out);
				file << (*pw);
				file.close();
		}
}
//comupute a single Voronoi cell based on data of the corresponding particle and the surrounding ones.
void CellMachine::processing(){
  if(cid<0){std::cout<<"Using parAttr to initialize the Cell Machine"<<std::endl;}
    if(nx*ny*nz>12000){//the number is set by experience
      blockMem = 64;
    }
    voro::container con(xmin*scale, xmax*scale, ymin*scale, ymax*scale, zmin*scale, zmax*scale, nx, ny, nz, xpbc, ypbc, zpbc, blockMem);
    pcon->setup(con);
    delete pcon;
    pcon = NULL;
		polywriter *pw = new polywriter();
		unsigned long long numberofpoints = labelidmap.size();
		// merge voronoi cells to set voronoi diagram
		// loop over all voronoi cells
		voro::c_loop_all cla(con);
		// cell currently worked on
		 unsigned long long status = 0;
		 // counter for process output
		 double outputSteps = 100;
		 double tenpercentSteps = 10/outputSteps*static_cast<double>(numberofpoints);
		 std::cout << "merge voronoi cells of "<< numberofpoints<<" points..."<<std::endl;
		 double target = tenpercentSteps;
		 if(cla.start())
		 {
			//std::cout << "started\n" << std::flush;
			do
			{
					voro::voronoicell_neighbor c;
					status++;
					if(con.compute_cell(c,cla))
					{
							if ( status >= target)
							{
									target += tenpercentSteps;
									std::cout << static_cast<int>(static_cast<double>(status)/static_cast<double>(numberofpoints)*outputSteps) << " \%\t" << std::flush;
							}
							//std::cout << "computed"  << std::endl;
							double xc = 0;
							double yc = 0;
							double zc = 0;
							// Get the position of the current particle under consideration
							cla.pos(xc,yc,zc);

							unsigned int id = cla.pid();
							unsigned int l = labelidmap[id];
              //std::cout<<"cellID="<<cid<<"pid="<<id<<"label="<<l<<std::endl;
              if (l!=cid) continue;//we just need only one cell
							//cell volume
							cellVolume += c.volume();
							std::vector<int> f; // list of face vertices (bracketed, as ID)
							c.face_vertices(f);

							std::vector<double> vertices;   // all vertices for this cell
							c.vertices(xc,yc,zc, vertices);

							std::vector<int> w; // neighbors of faces
							c.neighbors(w);
							// for this cell, loop over all faces and get the corresponding neighbors
							for (unsigned int k = 0; k != w.size(); ++k)
							{
									// compare if id for this cell and the face-neighbor is the same
									int n = w[k];   // ID of neighbor cell
                  //std::cout<<"n="<<n<<"size="<<labelidmap.size()<<std::endl;
                  //not sure why n could be negative
                  if(n<0) continue;
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
                      if(cellVTK){//FIXME:should be another flag.
                        pw->addface(positionlist, l);
                      }
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
                          /*
													cellNormalTensor[0] += wx*wx/area_tmp1;//n11
													cellNormalTensor[1] += wx*wy/area_tmp1;//n12
													cellNormalTensor[2] += wx*wz/area_tmp1;//n13
													cellNormalTensor[3] += wy*wy/area_tmp1;//n22
													cellNormalTensor[4] += wy*wz/area_tmp1;//n23
													cellNormalTensor[5] += wz*wz/area_tmp1;//n33
                          */
													//just compute the product of normal and area
													cellNormalAreaTensor[0] += 0.5*wx*wx/area_tmp2;//n11
													cellNormalAreaTensor[1] += 0.5*wx*wy/area_tmp2;//n12
													cellNormalAreaTensor[2] += 0.5*wx*wz/area_tmp2;//n13
													cellNormalAreaTensor[3] += 0.5*wy*wy/area_tmp2;//n22
													cellNormalAreaTensor[4] += 0.5*wy*wz/area_tmp2;//n23
													cellNormalAreaTensor[5] += 0.5*wz*wz/area_tmp2;//n33
											//}

													area += area_tmp2;
											}
											cellSurfaceArea += 0.5*area;
											//caculate the unit normal vector of the facet

									}
							}

					}
			}
			while (cla.inc());
	}
	std::cout<<std::endl;
	writeLocal(pw);
  delete pw;
  writeGlobal();
	std::cout << std::endl;
	con.clear();
}
