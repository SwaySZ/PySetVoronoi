#ifndef __SUPERQUADRIC__
#define __SUPERQUADRIC__
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
#include "parAttr.hpp"
#include <fstream>
#include <string>


#include <Eigen/Dense>

typedef Eigen::Quaternion<double> Quaternionr;
typedef Eigen::Matrix<double,3,3> Matrix3r;
typedef Eigen::Matrix<double,2,1> Vector2r;
typedef Eigen::Matrix<double,6,1> Vector6r;

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


class PolySuperellipsoid{
	public:
		//constructor
		PolySuperellipsoid(Vector2r eps, Vector6r halflen);//{/*createIndex();*/ rxyz = Vector3r(x,y,z); eps = Vector2r(ep1,ep2); Initial();};
		virtual ~PolySuperellipsoid(){};
		Vector6r getrxyz(){return rxyz;}
		Vector6r getrxyz_ref(){return rxyz_ref;}
		//double geteps1(){return eps1;}
		//double geteps2(){return eps2;}
		Vector2r geteps(){return eps;}
		double getr_max(){return r_max;}
		double getr_max_ref(){return r_max_ref;}
		Vector3r getPosition(){return Position;}
		double getVolume(){Initial();return Volume;}
		Vector3r getMassCenter(){Initial();return massCenter;}
		Vector3r getMassCenter_ref(){return massCenter_ref;}
		Vector3r getInertia(){Initial();return Inertia;}
		Quaternionr getOrientation(){return Orientation;}

		void setPosition(Vector3r p){Position = p;}
		void setOrientation(Quaternionr Ori){
		Ori.normalize();Orientation = Ori;}
		void setRxyz(Vector6r rxyz_in){rxyz = rxyz_in;}
		void setMassCenter(Vector3r mc){massCenter = mc;}
		void setRmax(double rm){r_max = rm;}

		Vector3r getSurface(Vector2r phi) const;//at the local with respect to the geometric center
		Vector3r getSurfaceMC(Vector2r phi) const;//at the local with respect to the mass center
		Vector3r getNormal(Vector2r);
		double getSurfaceArea(int w, int h);
	protected:
		//double rx, ry, rz, eps1, eps2;//the main parameters of a superquadric
		//r_max:the maximum r in the three axes
	  double r_max;
		double ep1,ep2;
		Vector6r rxyz;
		Vector3r Position;			//the center position of a superquadric, i.e (x, y, z)
		double Volume;

		Vector3r massCenter;
		Quaternionr Orientation;	//orientation
};

PolySuperellipsoid::PolySuperellipsoid(double rx1_, double rx2_,double ry1_, double ry2_,double rz1_, double rz2_,double ep1,double ep2)
{

	eps1=ep1;
	eps2=ep2;
  rxyz = Vector6r(rx1_,rx2_,ry1_,ry2_,rz1_,rz2_);

        //
		//get the mass center
		std::vector<Vector3r> centers;
		//std::vector<double> vols;
		massCenter = Vector3r(0,0,0);
		Volume = 0.0;
		for(int k=0;k<2;k++){
			for(int j=0;j<2;j++){
				for(int i=0;i<2;i++){
					Vector3r center;
					double beta1,beta2,rx,ry,rz;
					//set rx,ry,rz
					rx = rxyz[i];
					ry = rxyz[j+2];
					rz = rxyz[k+4];
					beta1 = beta(0.5*eps1,eps1)*beta(0.5*eps2,1.5*eps2);
					beta2 = beta(0.5*eps1,0.5*eps1)*beta(0.5*eps2,eps2);
					center[0] = (1-2*i)*rx*3.0/4.0*beta1/beta2;
					center[1] = (1-2*j)*ry*3.0/4.0*beta1/beta2;
					center[2] = (1-2*k)*rz*3.0/4.0*beta(eps2,eps2)/beta(eps2,0.5*eps2);
					centers.push_back(center);
					double a,v;//some coefficients
					a = 0.125*rx*ry*rz*eps1*eps2;
					beta1 = beta(1.5*eps1, 0.5*eps1)*beta(0.5*eps2,2.*eps2+1);
					beta2 = beta(0.5*eps1, 0.5*eps1+1)*beta(1.5*eps2, eps2+1);
					v = 2.0*a*beta(eps1*0.5,eps1*0.5)*beta(eps2*0.5+1,eps2);     //volume = 2a1a2a3ep1ep2B(ep1/2+1,ep1)B(ep2/2,ep2/2) see the reference
					//vols.push_back(v);
					Volume += v;
					massCenter += center*v;
				}
			}
		}
		massCenter /= Volume;
	double rx,ry,rz,rx1,ry1,rz1,rx2,ry2,rz2;
	rx1 = rxyz[0]-massCenter[0];
	rx2 = rxyz[1]+massCenter[0];
	ry1 = rxyz[2]-massCenter[1];
	ry2 = rxyz[3]+massCenter[1];
	rz1 = rxyz[4]-massCenter[2];
	rz2 = rxyz[5]+massCenter[2];
	rx = max(rx1,rx2);
	ry = max(ry1,ry2);
	rz = max(rz1,rz2);
	r_max = (rx > ry) ? rx : ry;
	r_max = r_max > rz ? r_max : rz;
	if ((eps1 < 1.0) || (eps2 < 1.0)){
					r_max *= pow(2.0,0.5);
	}
	Orientation = Quaternionr::Identity();
	//rotation matrices
  rot_mat2local = (Orientation).conjugate().toRotationMatrix();//to particle's system
	rot_mat2global = (Orientation).toRotationMatrix();//to global system
	Position = Vector3r::Zero();

}

double PolySuperellipsoid::getSurfaceArea(int w, int h){
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
				Surf = getSurface(Vector2r(phi0,phi1));
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
//****************************************************************************************
Vector3r PolySuperellipsoid::getSurface(Vector2r phi) const//in local coordinates system
{
		assert(phi[0]>=0);
		double x,y,z;
		x = Mathr::Sign(cos(phi(0)))*pow(fabs(cos(phi(0))), eps[0])*pow(fabs(cos(phi(1))), eps[1]);
		y = Mathr::Sign(sin(phi(0)))*pow(fabs(sin(phi(0))), eps[0])*pow(fabs(cos(phi(1))), eps[1]);
		z = Mathr::Sign(sin(phi(1)))*pow(fabs(sin(phi(1))), eps[1]);
		x *= rxyz[(x>0?0:1)];
		y *= rxyz[(y>0?2:3)];
		z *= rxyz[(z>0?4:5)];
		return Vector3r(x,y,z);
}
Vector3r PolySuperellipsoid::getSurfaceMC(Vector2r phi) const//in local coordinates system with the origin at the mass center
{
		return getSurface(phi) - massCenter;
}

Vector3r PolySuperellipsoid::getNormal(Vector2r phi)
{
	Vector3r n;
	n(0) = Mathr::Sign(cos(phi(0))) *pow(fabs(cos(phi(0))), 2.-eps1)*pow(fabs(cos(phi(1))), 2.-eps2);
	n(1) = Mathr::Sign(sin(phi(0))) *pow(fabs(sin(phi(0))), 2.-eps1)*pow(fabs(cos(phi(1))), 2.-eps2);
	n(2) = Mathr::Sign(sin(phi(1))) *pow(fabs(sin(phi(1))), 2.-eps2);
	n(0) /= rxyz[(n(0)>0?0:1)];
	n(1) /= rxyz[(n(1)>0?2:3)];
	n(2) /= rxyz[(n(2)>0?4:5)];
	return n;
}
#endif
