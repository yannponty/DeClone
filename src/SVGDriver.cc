/*  DeClone: A software for computing and analyzing ancestral adjacency scenarios.
 *  Copyright (C) 2015 Cedric Chauve, Yann Ponty, Ashok Rajaraman, Joao P.P. Zanetti
 *
 *  This file is part of DeClone.
 *  
 *  DeClone is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  DeClone is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with DeClone.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Contact: <yann.ponty@lix.polytechnique.fr>.
 *
 *
 *  DeClone uses the Quickhull algorithm implementation programmed by 
 *  Anatoly V. Tomilov. The code is available on <https://bitbucket.org/tomilov/quickhull/src/585267abb3a63794c04fc8325aa9ec9f726112ed/include/quickhull.hpp?at=master>.
 *
 *  Contact: <tomilovanatoliy@gmail.com>
 */

#include "SVGDriver.hh"
#include <fstream>

using namespace std;

ostream & operator<< (ostream & o, const Color & c)
{
	int rpc = (int) ((((double) c.r) ) * 100);
	int gpc = (int) ((((double) c.g) ) * 100);
	int bpc = (int) ((((double) c.b) ) * 100);
	return o << "rgb(" << rpc << "%, " << gpc << "%, " << bpc << "%)";
}

Point operator+ (const Point & p1, const Point & p2)
{
	return Point(p1.x+p2.x,p1.y+p2.y);
}

Point operator* (double k, const Point & p)
{
	return Point(k*p.x,k*p.y);
}

Point operator* (const Point & p,double k)
{
	return Point(k*p.x,k*p.y);
}
