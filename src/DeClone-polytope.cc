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

#include <cfloat>
#include <map>
#include <utility>
#include "RecTrees.hh"
#include "utils.hh"
#include "ConvexPolytope.hh"
//#include <fstream>
#include <iostream>
#include <cstdlib>

#define RESULT_TYPE Polytope
#define DIMENSION 2

#define allocateMatrix allocateMatrixPolytope
#define deleteMatrix deleteMatrixPolytope
#define computeMatrix computeMatrixPolytope

#define AdjGain Polytope(DIMENSION,adjgain)
#define AdjBreak Polytope(DIMENSION,adjbreak)
#define ZERO Polytope(DIMENSION,origin)
#define INF Polytope(DIMENSION)
#define RESCALING_FACTOR(a) ZERO 

#define PLUS(a,b) minkowski_addition(a,b) 
#define MIN(a,b,c,adj,g1,g2) convex_hull(a,b,adj)


const double adjgain[] = {1.,0.};
const double adjbreak[] = {0.,1.};
const double origin[] = {0.,0.};

//Polytope PARAMETER_POLYTOPE_OLD(const int dim,const double* point){
//  return Polytope(dim, point);
//}
//

Polytope minkowski_addition(Polytope p1, Polytope p2)
{
  return p1.minkovskiSum(p2);
}


Polytope convex_hull(Polytope p1, \
                     Polytope p2, \
                     const bool adj){
    return p1.convexSum(p2);
}


#include "DeCoDP.cc"


Polytope polycomputeValidAdjacencyTrees(RecTree *tree1, RecTree *tree2,  map<string, map<string,string> > & adjacencies)
{

    return MIN(MIN(INF,computeMatrixPolytope(tree1,tree2,true, adjacencies),"Root",true,"N/A","N/A"),
               computeMatrixPolytope(tree1,tree2,false, adjacencies),
  	           "Root",false,"N/A","N/A");

}
