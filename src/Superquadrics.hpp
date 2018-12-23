#pragma once
/*
 * =====================================================================================
 *
 *       Filename:  superquadrics.h
 *
 *    Description:  superquadrics have been defined.
 *
 *        Version:  1.0
 *        Created:  07/22/2015 05:16:35 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  zhswee (zhswee@gmail.com)
 *   Organization: South China University of Technology
 *
 * =====================================================================================
 */
#include "pointpattern.hpp"
#include "duplicationremover.hpp"
#include "particleparameterset.hpp"
#include "parAttr.hpp"
#include <fstream>
#include <string>


#include <Eigen/Dense>

typedef Eigen::Quaternion<double> Quaternionr;
typedef Eigen::Matrix<double,3,3> Matrix3r;

using namespace Eigen;

#define M_PIl                3.141592653589793238462643383279502884L /* pi */
# define M_PI_2l        1.570796326794896619231321691639751442L /* pi/2 */
# define M_PI_4l        0.785398163397448309615660845819875721L /* pi/4 */


////////////////////////////////some auxiliary functions/////////////////////////////////////////////////////////////////
//double cot(double x){return tan(M_PI_2l-x);}
int Sign(double f){ if(f<0) return -1; if(f>0) return 1; return 0; }
/////////////////////////gamma func----referring to NR3.0
double gammln(const double xx) {
	int j;
	double x,tmp,y,ser;
	static const double cof[14]={57.1562356658629235,-59.5979603554754912,
	14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,
	.465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,
	-.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,
	.844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5};
	if (xx <= 0) throw("bad arg in gammln");
	y=x=xx;
	tmp = x+5.24218750000000000;
	tmp = (x+0.5)*log(tmp)-tmp;
	ser = 0.999999999999997092;
	for (j=0;j<14;j++) ser += cof[j]/++y;
	return tmp+log(2.5066282746310005*ser/x);
}
///////////////////////beta func
double beta(const double z, const double w) {
	return exp(gammln(z)+gammln(w)-gammln(z+w));
}

double triangleArea(Vector3d p1,Vector3d p2,Vector3d p3){
    Vector3d v12 = p2 -p1;
    Vector3d v13 = p3 - p1;
    Vector3d vcross = v12.cross(v13);
    return 0.5*vcross.norm();

}

double polygonArea(Vector3d p1,Vector3d p2,Vector3d p3,Vector3d p4){
   //calculate surface area
    double area = 0.0;
    area = triangleArea(p1,p2,p3) + triangleArea(p1,p3,p4);

    return area;

}

double getSurfaceArea(double rx,double ry,double rz, double eps1, double eps2,int w, int h){
////
    std::vector<Vector3d> vertices;//discretized vertices of a particle
    double area = 0.0;
    int i,j;
    double a=0.0,b=0.0,phi0,phi1;
    double hStep=M_PI/(h-1);
    double wStep=2*M_PI/w;

    Vector3d Surf;
    //surface discretization
    for(a=0.0,i=0;i<h;i++,a+=hStep){
      for(b=0.0,j=0;j<w;j++,b+=wStep)
      {
        phi0 = b;
        phi1 = a-M_PI_2l;

        //get surface point
        Surf(0) = Sign(cos(phi0))*rx*pow(fabs(cos(phi0)), eps1)*pow(fabs(cos(phi1)), eps2);
        Surf(1) = Sign(sin(phi0))*ry*pow(fabs(sin(phi0)), eps1)*pow(fabs(cos(phi1)), eps2);
        Surf(2) = Sign(sin(phi1))*rz*pow(fabs(sin(phi1)), eps2);
	    vertices.push_back(Surf);
      }
    }

    //polygons of surface slices

	for(i=0;i<h-1;i++)
	{
	    for(j=0;j<w-1;j++)
	    {

            area += polygonArea(vertices[i*w+j],vertices[i*w+j+1],vertices[(i+1)*w+j+1],vertices[(i+1)*w+j]);
	    }

        area += polygonArea(vertices[i*w+j],vertices[i*w],vertices[(i+1)*w],vertices[(i+1)*w+j]);

	}
    return area;
}

void extractSuperquadric(pointpattern& pp, unsigned int id, double scaledist,std::vector<double> &set, int w_slices, int h_slices, double& area, double& volume){
        double rx,ry,rz,rmin;
        double eps1,eps2;
        Vector3d Position;
        Quaternionr Ori;
        rx = set[0];
        ry = set[1];
        rz = set[2];
        eps1 = set[3];
        eps2 = set[4];//
        Position = Vector3d(set[5], set[6], set[7]);

        Ori.w() = set[8];
        Ori.x() = set[9];
        Ori.y() = set[10];
        Ori.z() = set[11];
        Matrix3r A = Ori.toRotationMatrix();

        //find rmin
        rmin = (rx<ry)?rx:ry;
        rmin = (rmin<rz)?rmin:rz;
        int i,j,w=w_slices,h=h_slices;
        double a=0.0,b=0.0,phi0,phi1;
        double hStep=M_PI/(h-1);
        double wStep=2*M_PI/w;
	//double scale = 0.95;

        Vector3d p,n;
        //std::cout<<"scale1="<<scale<<std::endl;
        //double scaledist=0.0;
        //scaledist = scale*rmin;
		//caution:the two polar points should be degenerated.
        for(a=hStep,i=0;i<h-2;i++,a+=hStep)
          {
          for(b=0.0,j=0;j<w;j++,b+=wStep)
          {     phi0 = b;
                phi1 = a-M_PI_2l;

	        //get surface point
	        p(0) = Sign(cos(phi0))*rx*pow(fabs(cos(phi0)), eps1)*pow(fabs(cos(phi1)), eps2);
	        p(1) = Sign(sin(phi0))*ry*pow(fabs(sin(phi0)), eps1)*pow(fabs(cos(phi1)), eps2);
	        p(2) = Sign(sin(phi1))*rz*pow(fabs(sin(phi1)), eps2);
            //get out-ward normal at the surface point

	        n(0) = Sign(cos(phi0)) /rx *pow(fabs(cos(phi0)), 2.-eps1)*pow(fabs(cos(phi1)), 2.-eps2);
	        n(1) = Sign(sin(phi0)) /ry *pow(fabs(sin(phi0)), 2.-eps1)*pow(fabs(cos(phi1)), 2.-eps2);
	        n(2) = Sign(sin(phi1)) /rz *pow(fabs(sin(phi1)), 2.-eps2);
            n.normalize();
            //std::cout<<"pnorm="<<p.norm()<<std::endl;
            //std::cout<<"rmin="<<rmin<<std::endl;
            //scale = 1.0 - scaledist/p.norm();
            //std::cout<<"scale="<<scale<<std::endl;
            //p = p*scale;
			//debuging
			{
			/*	if (id==8){
					std::cout<<p(0)<<" "<<p(1)<<" "<<p(2)<<std::endl;
				}*/
			}
            p = p - n*scaledist;
            p = A*p + Position;
            pp.addpoint(id,p(0),p(1),p(2));
         }
         }
		//two polar points
		  for(a=0.0,i=0;i<2;i++,a+=M_PI)
          {
		    phi0 = 0;
            phi1 = a-M_PI_2l;

	        //get surface point
	        p(0) = Sign(cos(phi0))*rx*pow(fabs(cos(phi0)), eps1)*pow(fabs(cos(phi1)), eps2);
	        p(1) = Sign(sin(phi0))*ry*pow(fabs(sin(phi0)), eps1)*pow(fabs(cos(phi1)), eps2);
	        p(2) = Sign(sin(phi1))*rz*pow(fabs(sin(phi1)), eps2);
            //get out-ward normal at the surface point

	        n(0) = Sign(cos(phi0)) /rx *pow(fabs(cos(phi0)), 2.-eps1)*pow(fabs(cos(phi1)), 2.-eps2);
	        n(1) = Sign(sin(phi0)) /ry *pow(fabs(sin(phi0)), 2.-eps1)*pow(fabs(cos(phi1)), 2.-eps2);
	        n(2) = Sign(sin(phi1)) /rz *pow(fabs(sin(phi1)), 2.-eps2);
            n.normalize();
            //std::cout<<"pnorm="<<p.norm()<<std::endl;
            //std::cout<<"rmin="<<rmin<<std::endl;
            //scale = 1.0 - scaledist/p.norm();
            //std::cout<<"scale="<<scale<<std::endl;
            //p = p*scale;
			//debuging
			{/*
				if (id==8){
					std::cout<<p(0)<<" "<<p(1)<<" "<<p(2)<<std::endl;
				}*/
			}
            p = p - n*scaledist;
            p = A*p + Position;
            pp.addpoint(id,p(0),p(1),p(2));

         }
    double alpha;
    alpha = rx*ry*rz*eps1*eps2;
	//beta1 = beta(1.5*eps2, 0.5*eps2)*beta(0.5*eps1,2.*eps1+1);
	//beta2 = beta(0.5*eps2, 0.5*eps2+1)*beta(1.5*eps1, eps1+1);
	volume = 2.*alpha*beta(eps1*0.5+1,eps1)*beta(eps2*0.5,eps2*0.5);     //volume = 2a1a2a3ep1ep2B(ep1/2+1,ep1)B(ep2/2,ep2/2) see the reference
    //surface area
    area = getSurfaceArea( rx, ry, rz,  eps1, eps2, 50,50);//2500 discretized points on the surface

}

void pointCloud_Superquadric(unsigned int id, std::string outfile, double scaledist,std::vector<double> &set, int w_slices, int h_slices, double& area, double& volume, particleAttr& pattr){
        double rx,ry,rz,rmin,rmax;
        double eps1,eps2;
        Vector3d Position;
        Quaternionr Ori;
        rx = set[0];
        ry = set[1];
        rz = set[2];
        eps1 = set[3];
        eps2 = set[4];//
        Position = Vector3d(set[5], set[6], set[7]);

        Ori.w() = set[8];
        Ori.x() = set[9];
        Ori.y() = set[10];
        Ori.z() = set[11];
        Matrix3r A = Ori.toRotationMatrix();
				//update particle attribution
				pattr.ID = id;
				pattr.centerx = set[5];
				pattr.centery = set[6];
				pattr.centerz = set[7];
				double xmin(1e10),xmax(-1e10),ymin(1e10),ymax(-1e10),zmin(1e10),zmax(-1e10);

        //find rmin
        rmin = (rx<ry)?rx:ry;
        rmin = (rmin<rz)?rmin:rz;
				//find rmax
				rmax = (rx>ry)?rx:ry;
				rmax = (rmax>rz)?rmax:rz;
				pattr.radius = rmax;
        int i,j,w=w_slices,h=h_slices;
        double a=0.0,b=0.0,phi0,phi1;
        double hStep=M_PI/(h-1);
        double wStep=2*M_PI/w;
	//double scale = 0.95;

        Vector3d p,n;
        //std::cout<<"scale1="<<scale<<std::endl;
        //double scaledist=0.0;
        //scaledist = scale*rmin;
				std::ofstream fp;
				fp.open(outfile.c_str(),std::ios::out);
				//pointpattern pp;
		//caution:the two polar points should be degenerated.
        for(a=hStep,i=0;i<h-2;i++,a+=hStep)
          {
          for(b=0.0,j=0;j<w;j++,b+=wStep)
          {     phi0 = b;
                phi1 = a-M_PI_2l;

	        //get surface point
	        p(0) = Sign(cos(phi0))*rx*pow(fabs(cos(phi0)), eps1)*pow(fabs(cos(phi1)), eps2);
	        p(1) = Sign(sin(phi0))*ry*pow(fabs(sin(phi0)), eps1)*pow(fabs(cos(phi1)), eps2);
	        p(2) = Sign(sin(phi1))*rz*pow(fabs(sin(phi1)), eps2);
            //get out-ward normal at the surface point

	        n(0) = Sign(cos(phi0)) /rx *pow(fabs(cos(phi0)), 2.-eps1)*pow(fabs(cos(phi1)), 2.-eps2);
	        n(1) = Sign(sin(phi0)) /ry *pow(fabs(sin(phi0)), 2.-eps1)*pow(fabs(cos(phi1)), 2.-eps2);
	        n(2) = Sign(sin(phi1)) /rz *pow(fabs(sin(phi1)), 2.-eps2);
            n.normalize();
            //std::cout<<"pnorm="<<p.norm()<<std::endl;
            //std::cout<<"rmin="<<rmin<<std::endl;
            //scale = 1.0 - scaledist/p.norm();
            //std::cout<<"scale="<<scale<<std::endl;
            //p = p*scale;
			//debuging
			{
			/*	if (id==8){
					std::cout<<p(0)<<" "<<p(1)<<" "<<p(2)<<std::endl;
				}*/
			}
            p = p - n*scaledist;
            p = A*p + Position;
            //pp.addpoint(0,p(0),p(1),p(2));
						fp <<std::scientific<< p(0)<<"\t" << p(1)<<"\t" << p(2) <<std::endl;
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
	        p(0) = Sign(cos(phi0))*rx*pow(fabs(cos(phi0)), eps1)*pow(fabs(cos(phi1)), eps2);
	        p(1) = Sign(sin(phi0))*ry*pow(fabs(sin(phi0)), eps1)*pow(fabs(cos(phi1)), eps2);
	        p(2) = Sign(sin(phi1))*rz*pow(fabs(sin(phi1)), eps2);
            //get out-ward normal at the surface point

	        n(0) = Sign(cos(phi0)) /rx *pow(fabs(cos(phi0)), 2.-eps1)*pow(fabs(cos(phi1)), 2.-eps2);
	        n(1) = Sign(sin(phi0)) /ry *pow(fabs(sin(phi0)), 2.-eps1)*pow(fabs(cos(phi1)), 2.-eps2);
	        n(2) = Sign(sin(phi1)) /rz *pow(fabs(sin(phi1)), 2.-eps2);
            n.normalize();
            //std::cout<<"pnorm="<<p.norm()<<std::endl;
            //std::cout<<"rmin="<<rmin<<std::endl;
            //scale = 1.0 - scaledist/p.norm();
            //std::cout<<"scale="<<scale<<std::endl;
            //p = p*scale;
            p = p - n*scaledist;
            p = A*p + Position;
            //pp.addpoint(0,p(0),p(1),p(2));
						fp <<std::scientific<< p(0)<<"\t" << p(1)<<"\t" << p(2) <<std::endl;
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
    double alpha;
    alpha = rx*ry*rz*eps1*eps2;
	//beta1 = beta(1.5*eps2, 0.5*eps2)*beta(0.5*eps1,2.*eps1+1);
	//beta2 = beta(0.5*eps2, 0.5*eps2+1)*beta(1.5*eps1, eps1+1);
	volume = 2.*alpha*beta(eps1*0.5+1,eps1)*beta(eps2*0.5,eps2*0.5);     //volume = 2a1a2a3ep1ep2B(ep1/2+1,ep1)B(ep2/2,ep2/2) see the reference
    //surface area
    area = getSurfaceArea( rx, ry, rz,  eps1, eps2, 50,50);//2500 discretized points on the surface

}
