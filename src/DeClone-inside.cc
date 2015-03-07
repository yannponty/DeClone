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


#include "DeClone-inside.hh"

#include <cfloat>
#include <map>
#include <utility>
#include "RecTrees.hh"
#include "utils.hh"


#define RESULT_TYPE double
#define INF 0.

#define AdjGain exp(-adjacency_gain/kT)
#define AdjBreak exp(-adjacency_break/kT)
#define PLUS(a,b) (a*b)
#define MIN(a,b,c,adj,g1,g2) a+b
#define ZERO 1.

#define RESCALING_FACTOR(a) pow(scalingFactor,a)

#define allocateMatrix allocateMatrixInside
#define deleteMatrix deleteMatrixInside 
#define computeMatrix computeMatrixInside

#include "DeCoDP.cc"

double computeInside(RecTree *tree1, RecTree *tree2,  map<string, map<string,string> > & adjacencies, double kTval)
{
  kT = kTval;
  return MIN(
      MIN(INF,computeMatrixInside(tree1,tree2,true, adjacencies),"Root",true,"",""),
	   	computeMatrixInside(tree1,tree2,false, adjacencies),"Root",false,"",""
	);
}
