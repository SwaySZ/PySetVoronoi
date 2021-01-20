#pragma once
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
 *         Author:  zhswee (swzhao (at) scut.edu.cn)
 *   Organization: South China University of Technology
 *
 * =====================================================================================
 */
#include "parAttr.hpp"
#include <fstream>
#include <string>
#include "mathbase.hpp"


class PolySuperellipsoid{
	public:
		//constructor
		PolySuperellipsoid(double rx1_, double rx2_,double ry1_, double ry2_,double rz1_, double rz2_,double ep1,double ep2);
		virtual ~PolySuperellipsoid(){};
		Vector6r getrxyz(){return rxyz;}

		double getr_max(){return r_max;}
		Vector3r getPosition(){return Position;}

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
		double getVolume(){return Volume;}
		double getSurfaceArea(int w, int h);
	protected:
		//double rx, ry, rz, eps1, eps2;//the main parameters of a superquadric
		//r_max:the maximum r in the three axes
	  double r_max;
		double eps1,eps2;
		Vector6r rxyz;
		Vector3r Position;			//the center position of a superquadric, i.e (x, y, z)
		double Volume;

		Vector3r massCenter;
		Quaternionr Orientation;	//orientation
};
#endif
