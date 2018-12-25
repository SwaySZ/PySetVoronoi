#ifndef __PAR_ATTR__
#define __PAR_ATTR__
#include <vector>
struct particleAttr{
		int ID;
		double centerx;
		double centery;
		double centerz;
		double xmin,xmax,ymin,ymax,zmin,zmax;
		double xrange,yrange,zrange;
		double radius;
		//bool nearBoundary;//is the particle at or near the boundary? This flag is used to output vtk files of boundary particles automatically.
		std::vector<int> surroundedID;
};
#endif
