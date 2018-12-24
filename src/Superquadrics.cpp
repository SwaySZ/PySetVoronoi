#include "Superquadrics.hpp"

PolySuperellipsoid::PolySuperellipsoid(double rx1_in, double rx2_in,double ry1_in, double ry2_in,double rz1_in, double rz2_in,double ep1,double ep2)
{

	eps1=ep1;
	eps2=ep2;
	//Vector2r a = Vector2r(ep1,ep2);

  rxyz <<rx1_in, rx2_in, ry1_in, ry2_in, rz1_in, rz2_in;

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
	rx = std::max(rx1,rx2);
	ry = std::max(ry1,ry2);
	rz = std::max(rz1,rz2);
	r_max = (rx > ry) ? rx : ry;
	r_max = r_max > rz ? r_max : rz;
	if ((eps1 < 1.0) || (eps2 < 1.0)){
					r_max *= pow(2.0,0.5);
	}
	Orientation = Quaternionr::Identity();
	Position = Vector3r::Zero();

}

double PolySuperellipsoid::getSurfaceArea(int w, int h){
////
    std::vector<Vector3r> vertices;//discretized vertices of a particle
    double area = 0.0;
    int i,j;
    double a=0.0,b=0.0,phi0,phi1;
    double hStep=M_PI/(h-1);
    double wStep=2*M_PI/w;

    Vector3r Surf;
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
		x = Sign(cos(phi(0)))*pow(fabs(cos(phi(0))), eps1)*pow(fabs(cos(phi(1))), eps2);
		y = Sign(sin(phi(0)))*pow(fabs(sin(phi(0))), eps1)*pow(fabs(cos(phi(1))), eps2);
		z = Sign(sin(phi(1)))*pow(fabs(sin(phi(1))), eps2);
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
	n(0) = Sign(cos(phi(0))) *pow(fabs(cos(phi(0))), 2.-eps1)*pow(fabs(cos(phi(1))), 2.-eps2);
	n(1) = Sign(sin(phi(0))) *pow(fabs(sin(phi(0))), 2.-eps1)*pow(fabs(cos(phi(1))), 2.-eps2);
	n(2) = Sign(sin(phi(1))) *pow(fabs(sin(phi(1))), 2.-eps2);
	n(0) /= rxyz[(n(0)>0?0:1)];
	n(1) /= rxyz[(n(1)>0?2:3)];
	n(2) /= rxyz[(n(2)>0?4:5)];
	return n;
}
