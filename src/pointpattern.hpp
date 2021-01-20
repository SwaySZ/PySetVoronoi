/*
This file is copied from Pomelo by Simon Weis and Philipp Schoenhoefer 2016.
*/
#ifndef __POINT_PATTERN__
#define __POINT_PATTERN__

#include <iomanip>
#include <vector>
#include <map>
#include <cmath>
struct point
{
    point(double cx, double cy, double cz, int cl, long cf = -1, long cC = -1): x(cx), y (cy), z(cz), l(cl), faceID(cf), cellID(cC) {};
    double x, y, z;
    int l;
    long faceID, cellID;
};

inline bool checkdistancecloserthan (point& a, point& b, double e)
{
    double dx =fabs(a.x - b.x);
    double dy =fabs(a.y - b.y);
    double dz =fabs(a.z - b.z);

    return (dx < e && dy < e && dz < e);
};

class pointpattern
{
public:
    void addpoint(int l, double x, double y, double z);
    void addpointForCell(double x, double y, double z, int l, long cf, long cC);
    void print();
    void removeduplicates ( double epsilon);
    void removeduplicates ( double epsilon, pointpattern& p);
    std::vector<point> points;
    std::map<unsigned int, long> indexShift;    // first is particle ID, second is new Index or -1

    friend std::ostream& operator << (std::ostream &f, const pointpattern& p)
    {
        if(p.points.empty())
            return f;
        int oldl = p.points[0].l;


        for(unsigned int i = 0;i<p.points.size();i++){
        	if(oldl != p.points.at(i).l){
        		f << "\n\n";
        		oldl = p.points.at(i).l;
        	}
        	f << p.points.at(i).l << " " << std::setw(5) << p.points.at(i).x <<" " << std::setw(5) << p.points.at(i).y << " " << std::setw(5) << p.points.at(i).z << "\n";
        }

        return f;
    };

    friend std::ostream& operator >> (std::ostream &f, const pointpattern& p)
    {
        if(p.points.empty())
            return f;
	double xx = p.points[0].x;
	double yy = p.points[0].y;
	double zz = p.points[0].z;
        int oldf = p.points[0].faceID;
	unsigned int oldc = p.points[0].cellID;

        for(unsigned int i =0;i<p.points.size();i++){
        	if(oldf != p.points.at(i).faceID){
        		f << oldc << " " <<  std::setw(8)<< xx << " " << std::setw(8) << yy << " " << std::setw(8) << zz << "\n\n\n";
        	      oldf = p.points.at(i).faceID;
        	                oldc = p.points.at(i).cellID;
        	                xx = p.points.at(i).x;
        	                yy =p.points.at(i).y;
        	                zz =p.points.at(i).z;
        	}
        	f << p.points.at(i).cellID << " " <<  std::setw(8) << p.points.at(i).x << " " << std::setw(8) <<p.points.at(i).y << " " << std::setw(8) <<p.points.at(i).z << "\n";
        }
        f << oldc << " " <<  std::setw(8)<< xx << " " << std::setw(8) << yy << " " << std::setw(8) << zz;

        return f;
    };
};

#endif
