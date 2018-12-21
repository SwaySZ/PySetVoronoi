/*
 * process.hpp
 *
 *  Created on: Dec 16, 2018
 *      Author: chris
 */

#ifndef SRC_PROCESS_HPP_
#define SRC_PROCESS_HPP_

#include "pointpattern.hpp"

#include <vector>
//define a new type
	struct particleAttr{
			int ID;
			double centerx;
			double centery;
			double centerz;
			double xmin,xmax,ymin,ymax,zmin,zmax;
			double xrange,yrange,zrange;
			double radius;
			std::vector<int> surroundedID;
	};

class process {

public:

	//calculate attributes of each particle
		void hello();
//	void calculate_Attributes(std::vector<pointpattern> & pplist,std::vector<particleAttr>& parAttrlist,double scale ,double epsilon);
		void calculate_Attributes(std::vector<pointpattern>& pplist,std::vector<particleAttr>& parAttrlist,double scale ,double epsilon);
};

#endif /* SRC_PROCESS_H_ */
