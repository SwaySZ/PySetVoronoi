/* 
Copyright 2016 Simon Weis and Philipp Schoenhoefer

This file is part of Pomelo.

Pomelo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Pomelo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Pomelo.  If not, see <http://www.gnu.org/licenses/>.

The development of Pomelo took place at the Friedrich-Alexander University of Erlangen and was funded by the German Research Foundation (DFG) Forschergruppe FOR1548 "Geometry and Physics of Spatial Random Systems" (GPSRS). 
*/
#ifdef __GNUC__
#pragma GCC system_header
#include "../lib/selene/include/selene.h"
#pragma GCC system_header
#include "../lib/voro++/src/voro++.hh"
#endif

// since selene produces a lot of ugly template warnings with gcc's -Wall and -Wextra, we will use this dirty hack and define selene as a system header for gcc. This will supress the warnings for this header only and print out warnings for all other issues in the code
