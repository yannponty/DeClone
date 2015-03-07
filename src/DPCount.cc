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

// Dyn. Prog.

#include "DPCount.hh"
#include <iostream>

#define PLUS(a,b) (a*b)
#define MIN(a,b) (a+b)
#define ZERO 1.

#define TYPEDATA double
#define INFTY 0.
#define DUP_COST 1.
#define LOSS_COST 1.

#define allocateMatricesTASK allocateMatricesCount
#define fillMatricesTASK fillMatricesCount
#define deleteMatricesTASK deleteMatricesCount


#include "DPRaw.cc"

double countReconciliations(Tree * GeneTree, Tree * SpeciesTree)
{
	TYPEDATA** R = allocateMatricesTASK(GeneTree, SpeciesTree);
  fillMatricesTASK(GeneTree, SpeciesTree, R);
	deleteMatricesTASK(GeneTree, SpeciesTree, R);
	return R[GeneTree->getIndex()][SpeciesTree->getIndex()];
}
