#include "mathbase.hpp"
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

double triangleArea(Vector3r p1,Vector3r p2,Vector3r p3){
    Vector3r v12 = p2 -p1;
    Vector3r v13 = p3 - p1;
    Vector3r vcross = v12.cross(v13);
    return 0.5*vcross.norm();

}

double polygonArea(Vector3r p1,Vector3r p2,Vector3r p3,Vector3r p4){
   //calculate surface area
    double area = 0.0;
    area = triangleArea(p1,p2,p3) + triangleArea(p1,p3,p4);

    return area;

}
