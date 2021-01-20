
#ifndef __POLY_WRITER__
#define __POLY_WRITER__

#include <iomanip>
#include <vector>
#include <map>
#include <cmath>

#include <iterator>
#include <regex>

#include "duplicationremover.hpp"

class polywriter
{
public:
    void addface( std::vector<double> positionlist, unsigned int cellID)
    {
        unsigned int faceID = currentFaceLabel;
        currentFaceLabel++;
        faceCellMap[faceID] = cellID;


        std::vector<unsigned int> facevertexIDs;
        for(unsigned int i = 0; i != positionlist.size(); ++(++(++i)) )
        {
            double x = positionlist[i];
            double y = positionlist[i+1];
            double z = positionlist[i+2];
            unsigned int l = currentVertexLabel;

            p.addpointForCell(x, y, z, l, faceID, cellID);
            facevertexIDs.push_back(l);
            currentVertexLabel++;
        }

        faces[faceID] = facevertexIDs;
    };


    friend std::ostream& operator << (std::ostream &f, const polywriter& p)
    {
        f << "POINTS" << std::endl;
        f << std::fixed;
        for(auto it =  p.p.points.begin();
                it != p.p.points.end();
                ++it)
        {
           f << it->l << ":    " <<  std::setprecision(12) << it->x << " " << std::setprecision(12) << it-> y << " " << std::setprecision(12) << it->z<< std::endl;
        }

        f << "POLYS" <<  std::endl;
        for (
            auto it = p.faces.begin();
            it != p.faces.end();
            ++ it)
        {
            unsigned int faceID = it->first;
            unsigned int cellID = p.faceCellMap.at(faceID);
	        std::vector<unsigned int> testing;
            for (auto it2 = it->second.rbegin(); it2 != it->second.rend(); ++it2)
            {
		        bool doppelt = false;
		        for(unsigned int kk = 0; kk < testing.size(); kk++ ){
			        if(testing[kk] == (*it2) ){
				        doppelt = true;
				        break;
			        }
		        }
		        if(doppelt) continue;
                testing.push_back( (*it2) );
            }
	        if(2 < testing.size()){
	            f << it->first << ":    ";
	    	    for(unsigned int kk = 0; kk < testing.size(); kk++ ){
		            f << testing[kk] << " ";
	    	    }
            	f << "< c(0, 0, 0, " << cellID << ")" << std::endl;
	        }
	        else std::cout << testing.size() << " " << cellID << std::endl;
        }

        f << "END";
        return f;
    };


    void savePointPatternForGnuplot( std::string filename)
    {
        std::cout << "writing PointPattern file" << std::endl;
        std::ofstream file;
        file.open(filename.c_str(),std::ios::out);
        file >> p;
        file.close();
    }
    void savePolyVTK(std::string filename, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double rr,bool withboundary)
    {
//        std::cout << "writing Polydata into a VTK file" << std::endl;
        std::ofstream f;
        f.open(filename.c_str(),std::ios::out);
        //writing the file head
        f << "# vtk DataFile Version 2.0" << std::endl;
        f << "Sway data processing" << std::endl;
        f << "ASCII" << std::endl;
        f << "DATASET POLYDATA" << std::endl;
        f << "POINTS   "<<this->p.points.size() + 1 <<" float" << std::endl;//index is from 0 in VTK
        f << "0   0   0" << std::endl;//putting a point to take the place.


        f << std::fixed;
        for(auto it =  this->p.points.begin();
                it != this->p.points.end();
                ++it)
        {
           f <<  std::setprecision(12) << it->x << " " << std::setprecision(12) << it-> y << " " << std::setprecision(12) << it->z<< std::endl;
        }
        //polygon
        //std::cout<<xmin<<" "<<xmax<<" "<<ymin<<" "<<ymax<<" "<<zmin<<" "<<zmax<<" "<<rr<<std::endl;
        unsigned int face_num, enterty_num;
        std::vector<unsigned long> polydata;//store polygon data:face_num vertice1_id vertice2_id ...
        face_num = 0;
        enterty_num = 0;
        for (
            auto it = this->faces.begin();
            it != this->faces.end();
            ++ it)
        {
            unsigned int faceID = it->first;
            unsigned int cellID = this->faceCellMap.at(faceID);
	        std::vector<unsigned int> testing;
            bool boundary = true;//face on the boundry?
            double xx,yy,zz;

            for (auto it2 = it->second.rbegin(); it2 != it->second.rend(); ++it2)
            {
		        bool doppelt = false;
		        for(unsigned int kk = 0; kk < testing.size(); kk++ ){
			        if(testing[kk] == (*it2) ){
				        doppelt = true;
				        break;
			        }
		        }
		        if(doppelt) continue;
                testing.push_back( (*it2) );
            }
	        if(2 < testing.size()){
                //check the point is on the boundary or not
                /*
                count = 0;
                xx = (p.points[testing[0]-1].x+p.points[testing[1]-1].x+p.points[testing[2]-1].x)/3.0;
                yy = (p.points[testing[0]-1].y+p.points[testing[1]-1].y+p.points[testing[2]-1].y)/3.0;
                zz = (p.points[testing[0]-1].z+p.points[testing[1]-1].z+p.points[testing[2]-1].z)/3.0;
                if ((xx>xmax-rr)||(xx<xmin+rr)||(yy>ymax-rr)||(yy<ymin+rr)||(zz>zmax-rr)||(zz<zmin+rr)){
                    boundary = true;
                }*/
                if(!withboundary){
                    for(auto itt = testing.begin();itt != testing.end(); ++itt){
                        xx = p.points[(*itt)-1].x;
                        yy = p.points[(*itt)-1].y;
                        zz = p.points[(*itt)-1].z;

                        if (!((xx>xmax-rr)||(xx<xmin+rr)||(yy>ymax-rr)||(yy<ymin+rr)||(zz>zmax-rr)||(zz<zmin+rr))){
                                boundary = false;
                                break;

                        }
                    }
                }
                if (!withboundary && boundary){continue;}
	            face_num += 1;
	            enterty_num += testing.size() + 1;
                polydata.push_back(testing.size());
                for(unsigned int kk = 0; kk < testing.size(); kk++ ){
		            polydata.push_back(testing[kk]);
	    	    }

	        }
	        else { std::cout << testing.size() << " " << cellID << std::endl;}
        }

        f << "POLYGONS " << face_num <<" "<<enterty_num<< std::endl;
        for(unsigned long kk = 0; kk < polydata.size(); kk++ ){
            f << polydata[kk]<<"   ";
        }

    f.close();

    };

    void saveOnePolyVTK(std::string filename, unsigned int targetCellID)
    {
//        std::cout << "writing Polydata into a VTK file" << std::endl;
        std::ofstream f;
        f.open(filename.c_str(),std::ios::out);
        //writing the file head
        f << "# vtk DataFile Version 2.0" << std::endl;
        f << "Sway data processing" << std::endl;
        f << "ASCII" << std::endl;
        f << "DATASET POLYDATA" << std::endl;
        f << "POINTS   "<<this->p.points.size() + 1 <<" float" << std::endl;//index is from 0 in VTK
        f << "0   0   0" << std::endl;//putting a point to take the place.


        f << std::fixed;
        for(auto it =  this->p.points.begin();
                it != this->p.points.end();
                ++it)
        {
           f <<  std::setprecision(12) << it->x << " " << std::setprecision(12) << it-> y << " " << std::setprecision(12) << it->z<< std::endl;
        }
        //polygon
        //std::cout<<xmin<<" "<<xmax<<" "<<ymin<<" "<<ymax<<" "<<zmin<<" "<<zmax<<" "<<rr<<std::endl;
        unsigned int face_num, enterty_num;
        std::vector<unsigned long> polydata;//store polygon data:face_num vertice1_id vertice2_id ...
        face_num = 0;
        enterty_num = 0;
        for (
            auto it = this->faces.begin();
            it != this->faces.end();
            ++ it)
        {
            unsigned int faceID = it->first;
            unsigned int cellID = this->faceCellMap.at(faceID);
            if (cellID != targetCellID) continue;
	        std::vector<unsigned int> testing;



            for (auto it2 = it->second.rbegin(); it2 != it->second.rend(); ++it2)
            {
		        bool doppelt = false;
		        for(unsigned int kk = 0; kk < testing.size(); kk++ ){
			        if(testing[kk] == (*it2) ){
				        doppelt = true;
				        break;
			        }
		        }
		        if(doppelt) continue;
                testing.push_back( (*it2) );
            }
	        if(2 < testing.size()){


	            face_num += 1;
	            enterty_num += testing.size() + 1;
                polydata.push_back(testing.size());
                for(unsigned int kk = 0; kk < testing.size(); kk++ ){
		            polydata.push_back(testing[kk]);

	    	    }

	        }
	        else { std::cout << testing.size() << " " << cellID << std::endl;}
        }

        f << "POLYGONS " << face_num <<" "<<enterty_num<< std::endl;
        for(unsigned long kk = 0; kk < polydata.size(); kk++ ){
            f << polydata[kk]<<"   ";
        }
    f.close();

    };

    void saveOnePolyVTKnew(std::string filename, unsigned int targetCellID,double scale)
    {
//        std::cout << "writing Polydata into a VTK file" << std::endl;
        std::ofstream f;
        f.open(filename.c_str(),std::ios::out);
		//store points and polygons
		std::vector< point> vtkPoints;
		unsigned long numPoints = 0;
        unsigned int face_num, enterty_num;
        std::vector<unsigned long> polydata;//store polygon data:face_num vertice1_id vertice2_id ...
        face_num = 0;
        enterty_num = 0;
        for (
            auto it = this->faces.begin();
            it != this->faces.end();
            ++ it)
        {
            unsigned int faceID = it->first;

            unsigned int cellID = this->faceCellMap.at(faceID);
            if (cellID != targetCellID) continue;
	        std::vector<unsigned int> testing;



            for (auto it2 = it->second.rbegin(); it2 != it->second.rend(); ++it2)
            {
		        bool doppelt = false;
		        for(unsigned int kk = 0; kk < testing.size(); kk++ ){
			        if(testing[kk] == (*it2) ){
				        doppelt = true;
				        break;
			        }
		        }
		        if(doppelt) continue;
                testing.push_back( (*it2) );
            }
	        if(2 < testing.size()){


	            face_num += 1;
	            enterty_num += testing.size() + 1;
                polydata.push_back(testing.size());
                for(unsigned int kk = 0; kk < testing.size(); kk++ ){
		            //polydata.push_back(testing[kk]);
      					vtkPoints.push_back(this->p.points[testing[kk]-1]);
      					numPoints++;
                    polydata.push_back(numPoints-1);
	    	    }

	        }
	        else { std::cout << testing.size() << " " << cellID << std::endl;}
        }
        //writing the VTK file
        //writing the file head
        f << "# vtk DataFile Version 2.0" << std::endl;
        f << "Sway data processing" << std::endl;
        f << "ASCII" << std::endl;
        f << "DATASET POLYDATA" << std::endl;
        f << "POINTS   "<<vtkPoints.size() <<" float" << std::endl;//index is from 0 in VTK
		f << std::fixed;
        for(auto it =  vtkPoints.begin();
                it != vtkPoints.end();
                ++it)
        {
           f <<  std::setprecision(12) << it->x/scale << " "<< it-> y/scale << " " << it->z/scale<< std::endl;
        }
        //polygon
        f << "POLYGONS " << face_num <<" "<<enterty_num<< std::endl;
        for(unsigned long kk = 0; kk < polydata.size(); kk++ ){
            f << polydata[kk]<<"   ";
        }
    f.close();

    };
	void outputCellVTK(std::string filename,std::string outfileprefix,double scale)
	{
		std::ifstream infile;
		infile.open(filename.c_str(), std::ifstream::in);
		if (infile.fail())
		{
		    std::cout << "Cannot load file " << filename << std::endl;
		    return;
		}
	#pragma GCC diagnostic ignored "-Wwrite-strings"
		std::string line("");
		unsigned int linesloaded = 0;
		//std::getline(infile, line);
		std::cout<<"Warning!"<<std::endl;
		std::vector<int> targetCellIDs;
		while (std::getline(infile, line))
		{
		    if(line.find("#")!=std::string::npos) continue; // ignore comment lines
            //split
            
            const std::regex ws_re("\\s+"); // whitespace
            auto line_begin = std::sregex_token_iterator(line.begin(), line.end(), ws_re, -1);
            auto line_end = std::sregex_token_iterator();
            
            for (auto p = line_begin; p != line_end; ++p ){
                int d = std::stoi((*p));
                std::cout<<*p<<std::endl;
		        targetCellIDs.push_back(d);
            }
		    linesloaded++;

		}
		std::cout << "TargetCellIDs were imported successfully!" << std::endl << std::endl;
		infile.close();
		//write vtk files of all individual cells
		std::string suffix(".vtk");
		for(auto id = targetCellIDs.begin(); id != targetCellIDs.end(); ++id)
		{
			saveOnePolyVTKnew(outfileprefix+std::to_string(*id)+suffix, *id,scale);
		}
	}
  void saveOnePolyPOVnew(std::string filename, unsigned int targetCellID)
  {
//        std::cout << "writing Polydata into a VTK file" << std::endl;
      std::ofstream f;
      f.open(filename.c_str(),std::ios::out);
  //store points and polygons
  std::vector< point> vtkPoints;
  unsigned long numPoints = 0;
      unsigned int face_num, enterty_num;
      std::vector<unsigned long> polydata;//store polygon data:face_num vertice1_id vertice2_id ...
      face_num = 0;
      enterty_num = 0;
      for (
          auto it = this->faces.begin();
          it != this->faces.end();
          ++ it)
      {
          unsigned int faceID = it->first;

          unsigned int cellID = this->faceCellMap.at(faceID);
          if (cellID != targetCellID) continue;
        std::vector<unsigned int> testing;



          for (auto it2 = it->second.rbegin(); it2 != it->second.rend(); ++it2)
          {
          bool doppelt = false;
          for(unsigned int kk = 0; kk < testing.size(); kk++ ){
            if(testing[kk] == (*it2) ){
              doppelt = true;
              break;
            }
          }
          if(doppelt) continue;
              testing.push_back( (*it2) );
          }
        if(2 < testing.size()){


            face_num += 1;
            enterty_num += testing.size() + 1;
              polydata.push_back(testing.size());
              for(unsigned int kk = 0; kk < testing.size(); kk++ ){
              //polydata.push_back(testing[kk]);
              vtkPoints.push_back(this->p.points[testing[kk]-1]);
              numPoints++;
                  polydata.push_back(numPoints-1);
          }

          //draw a face as a polygon
          f<<"polygon {"<<testing.size()<<",";
          for(auto itt = testing.begin();itt != testing.end(); ++itt){
              f << " <"<<p.points[(*itt)-1].x<<","<<p.points[(*itt)-1].y<<","<<p.points[(*itt)-1].z<<">";
          }
          f<<"}\n";//f<<"\nscale 0.01}\n";//scaling using for avoiding "Possible Parse Error: Singular matrix in MInvers"
        }
        else { std::cout << testing.size() << " " << cellID << std::endl;}
      }

  f.close();

  };
    void savePolyPOV(std::string filename, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double rr)
    {
//        std::cout << "writing Pov file for postprocessing in Povray" << std::endl;
        std::ofstream f;
        f.open(filename.c_str(),std::ios::out);

        for (
            auto it = this->faces.begin();
            it != this->faces.end();
            ++ it)
        {
            unsigned int faceID = it->first;
            unsigned int cellID = this->faceCellMap.at(faceID);
	        std::vector<unsigned int> testing;
            bool boundary = true;//face on the boundry?
            double xx,yy,zz;

            for (auto it2 = it->second.rbegin(); it2 != it->second.rend(); ++it2)
            {
		        bool doppelt = false;
		        for(unsigned int kk = 0; kk < testing.size(); kk++ ){
			        if(testing[kk] == (*it2) ){
				        doppelt = true;
				        break;
			        }
		        }
		        if(doppelt) continue;
                testing.push_back( (*it2) );
            }
	        if(2 < testing.size()){
                //check the point is on the boundary or not

                for(auto itt = testing.begin();itt != testing.end(); ++itt){
                    xx = p.points[(*itt)-1].x;
                    yy = p.points[(*itt)-1].y;
                    zz = p.points[(*itt)-1].z;

                    if (!((xx>xmax-rr)||(xx<xmin+rr)||(yy>ymax-rr)||(yy<ymin+rr)||(zz>zmax-rr)||(zz<zmin+rr))){
                            boundary = false;
                            break;
                    }
                }
                if (boundary){continue;}
                //draw a face as a polygon
				f<<"polygon {"<<testing.size()<<",";
                for(auto itt = testing.begin();itt != testing.end(); ++itt){
                    xx = p.points[(*itt)-1].x;
                    yy = p.points[(*itt)-1].y;
                    zz = p.points[(*itt)-1].z;
                    f << " <"<<xx*100.0<<","<<yy*100.0<<","<<zz*100.0<<">";
                }
                f<<"\nscale 0.01}\n";//scaling using for avoiding "Possible Parse Error: Singular matrix in MInvers"

	        }
	        else { std::cout << testing.size() << " " << cellID << std::endl;}
        }

    f.close();

    };
/*	void save_cells_pov( std::string filename)
    {
        std::cout << "writing Pov file" << std::endl;
        std::ofstream file;
        file.open(filename);
        file >> p;
		//writing the file

        for (auto it = p.faces.rbegin();it != p.faces.rend();++ it)
        {
            unsigned int faceID = it->first;
            unsigned int cellID = p.faceCellMap.at(faceID);
			file << "// cell " << cellID << "\n";
	    	std::vector<unsigned int> testing;
            for (auto it2 = it->second.rbegin(); it2 != it->second.rend(); ++it2)
            {
			bool doppelt = false;
			for(unsigned int kk = 0; kk < testing.size(); kk++ ){
				if(testing[kk] == (*it2) ){
					doppelt = true;
					break;
				}
			}
			if(doppelt) continue;
		            testing.push_back( (*it2) );
		        }
			if(2 < testing.size()){
			    f << it->first << ":\t";
				for(unsigned int kk = 0; kk < testing.size(); kk++ ){
				//p.p.points[p.p.indexShift.at(testing[kk])];
				//file1 << "sphere{<"<<positionlist[i*3]<<","<<positionlist[i*3+1]<<","<<positionlist[i*3+2]<<">,r}\n";
				//index = (i+1)%3*3;
				//file1 << "cylinder{<"<<positionlist[i*3]<<","<<positionlist[i*3+1]<<","<<positionlist[i*3+2]<<">,<"<<positionlist[index]<<","<<positionlist[index+1]<<","<<positionlist[index+2]<<">,r}\n";

				f << testing[kk] << " ";
				}
		        	f << "< c(0, 0, 0, " << cellID << ")" << std::endl;
			}
			else std::cout << testing.size() << " " << cellID << std::endl;
		    }

		//

		//close the file
        file.close();
    }
        */
    void removeduplicates (double epsilon, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
    {
        std::cout << "polywriter: remove duplicates" << std::endl;
        duplicationremover d(16,16,16);
        d.setboundaries(xmin, xmax, ymin, ymax, zmin, zmax);
        std::cout << "\tadding points" << std::endl;
        d.addPoints(p, true);
        std::cout << "\tremoving duplicates" << std::endl;
        d.removeduplicates(epsilon);
        std::cout << "\tget back points" << std::endl;
        d.getallPoints(p);
        std::cout << "\tmatch back indices" << std::endl;
        rearrangeIndices(d.indexShift);

        //std::cout << "\torder indices" << std::endl;
        orderIndices();
    }

    void orderIndices ()
    {
        std::cout << "order indices" << std::endl;
        std::map<unsigned int, long> indexShift;
        unsigned int label = 1;
        for( auto it = p.points.begin(); it != p.points.end(); ++it)
        {
            indexShift[(*it).l] = label;
            (*it).l = label;
            label++;
        }
        rearrangeIndices(indexShift, false);
    }

    void rearrangeIndices(std::map<unsigned int, long>& indexShift, bool multipleTimes = true)
    {

        //std::cout << "### index map" << std::endl;
        unsigned int maxIndex = (*indexShift.rbegin()).first;
        for(unsigned int i = 1; i != maxIndex; ++i)
        {
            if(indexShift.find(i) == indexShift.end())
            {
                indexShift[i] = -1;
            }
            //std::cout << i << " " << indexShift[i] << std::endl;
        }

        std::cout << "\tPolywriter rearrange Indices"<< std::endl;
        unsigned int i = 0;
        double fivepercentSteps = 0.05*static_cast<double>(faces.size());
        double target = fivepercentSteps;
        for (auto it = faces.begin(); it != faces.end(); ++ it)
        {
            i++;
            if ( i >= target)
            {
                target += fivepercentSteps;
                std::cout << static_cast<int>(static_cast<double>(i)/static_cast<double>(faces.size())*100) << " \%\t"<< std::flush;
            }
            unsigned int j = 0;
            for (auto it2 = it->second.begin(); it2 != it->second.end(); ++it2)
            {
                j++;
                //std::cout << "\t\t\t" << j << " / " << it->second.size() << std::endl;
                unsigned int vertexIndex = (*it2);

                if( multipleTimes)
                {
                    while (indexShift[vertexIndex] != -1 )
                    {
                        //std::cout << "\t\t\t\t" << vertexIndex << std::endl;
                        vertexIndex = indexShift[vertexIndex];
                        if (indexShift.find(vertexIndex) == indexShift.end())
                        {
                            indexShift[vertexIndex] = -1;
                            break;
                        }
                    }
                    (*it2) = vertexIndex;
                }
                else
                {
                    (*it2) =indexShift[(*it2)];
                }
            }
        }
        std::cout << std::endl;
    }

    polywriter(){};
    ~polywriter(){};
    pointpattern p; // holds all the points
    std::map<unsigned int, unsigned int> faceCellMap;   // first is face id, second is cell id
    std::map<unsigned int, std::vector<unsigned int > > faces;
private:

    unsigned int currentVertexLabel = 1;
    unsigned int currentFaceLabel = 1;
};

#endif
