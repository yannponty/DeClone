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


#include "DeClone-countcoopts.hh"

#include <cfloat>
#include <cmath>
#include <map>
#include <utility>
#include "RecTrees.hh"
#include "utils.hh"

#define RESULT_TYPE pair<double,long>
#define INF pair<double,long>(DBL_MAX, 0)
#define AdjGain pair<double,long>(adjacency_gain, 1)
#define AdjBreak pair<double,long>(adjacency_break, 1)
#define PLUS(a,b) pair<double,long>(a.first+b.first, a.second*b.second)
#define MIN(a,b,c,adj,g1,g2) pair<double,long>(min(a.first,b.first), (a.first==min(a.first,b.first)?a.second:0)+(b.first==min(a.first,b.first)?b.second:0))
#define ZERO pair<double,long>(0., 1)

#define RESCALING_FACTOR(a) pair<double,long>(0., (long)pow(scalingFactor,a))

#define allocateMatrix allocateMatrixCountCoopts
#define deleteMatrix deleteMatrixCountCoopts
#define computeMatrix computeMatrixCountCoopts

#include "DeCoDP.cc"

pair<double,long> countCooptimalAdjacencyTrees(RecTree *tree1, RecTree *tree2,  map<string, map<string,string> > & adjacencies)
{
    return MIN(computeMatrixCountCoopts(tree1,tree2,true, adjacencies),
	       computeMatrixCountCoopts(tree1,tree2,false, adjacencies),
	       "Root","","",false);

}

ostream& operator<<(ostream & o, const pair<double,long> & v)
{
     o <<"("<<v.first<<","<<v.second<<")";
     return o;
}
