#pragma once
#ifndef __MATH_BASE__
#define __MATH_BASE__

#include <Eigen/Dense>

typedef Eigen::Quaternion<double> Quaternionr;
typedef Eigen::Matrix<double,3,3> Matrix3r;
typedef Eigen::Matrix<double,2,1> Vector2r;
typedef Eigen::Matrix<double,3,1> Vector3r;
typedef Eigen::Matrix<double,6,1> Vector6r;

//using namespace Eigen;

#define M_PIl                3.141592653589793238462643383279502884L /* pi */
#define M_PI_2l        1.570796326794896619231321691639751442L /* pi/2 */
#define M_PI_4l        0.785398163397448309615660845819875721L /* pi/4 */

int Sign(double f);
/////////////////////////gamma func----referring to NR3.0
double gammln(const double xx);
///////////////////////beta func
double beta(const double z, const double w);

double triangleArea(Vector3r p1,Vector3r p2,Vector3r p3);

double polygonArea(Vector3r p1,Vector3r p2,Vector3r p3,Vector3r p4);
#endif
