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

#include "DPStochasticBacktrack.hh"
#include <iostream>
#include <cstdlib>

#include "math.h"

#define PLUS(a,b) (a*b)
#define MIN(a,b) (a+b)
#define ZERO 1

#define TYPEDATA double
#define INFTY 0
#define DUP_COST exp(-1./kT)
#define LOSS_COST exp(-1./kT)

#define allocateMatricesTASK allocateMatricesStochasticBacktrack
#define fillMatricesTASK fillMatricesStochasticBacktrack
#define deleteMatricesTASK deleteMatricesStochasticBacktrack


#include "DPRaw.cc"

EditTree * stochasticBacktrack(Tree * g, Tree * s,TYPEDATA** R);

vector<EditTree*> stochasticReconciliations(Tree * GeneTree, Tree * SpeciesTree, int numTrees)
{
	vector<EditTree*> result;
	TYPEDATA** R = allocateMatricesTASK(GeneTree, SpeciesTree);
  fillMatricesTASK(GeneTree, SpeciesTree, R);
  //srand((unsigned int) time(0));
  //srand();
  for(int i=0;i<numTrees;i++)
  {
		EditTree * scenario = stochasticBacktrack(GeneTree, SpeciesTree,R);
		result.push_back(scenario);
	}
	deleteMatricesTASK(GeneTree, SpeciesTree, R);
	return result;
}


EditTree * stochasticBacktrack(Tree * g, Tree * s,TYPEDATA** R)
{
	int i = g->getIndex();
	int j = s->getIndex();
	int ag = (g->getLeft()?  g->getLeft()->getIndex():-1);
	int bg = (g->getRight()? g->getRight()->getIndex():-1);
	int as = (s->getLeft()?  s->getLeft()->getIndex():-1);
	int bs = (s->getRight()? s->getRight()->getIndex():-1);
	Tree* agt = g->getLeft() ;
	Tree* bgt = g->getRight();
	Tree* ast = s->getLeft() ;
	Tree* bst = s->getRight();
	double r = (R[i][j]*rand())/(RAND_MAX+1.0); 
	if (g->isLeaf() && s->isLeaf())
	{				
			if (extentCompatible(g, s))
			{ 
				if (r<ZERO)
				{ return new EditTree(formatLabel("Match",s),MATCH_TYPE); }
			}
	}
	else
	{
		if ((ag!=-1) && (bg!=-1))
		{ 
			r -= PLUS(DUP_COST,PLUS(R[ag][j],R[bg][j]));
			if (r<0)
			{ 
				EditTree * rec1 = stochasticBacktrack(agt, s, R);	
				EditTree * rec2 = stochasticBacktrack(bgt, s, R);	
				return new EditTree(formatLabel("Dup",s),rec1,rec2,DUP_TYPE);
			}
		}
//		if ((ag!=-1) && (as!=-1) && (bg!=-1))
//		{ 
//			r -= PLUS(R[ag][as],R[bg][as]);
//			if (r<0) 
//			{ 
//				EditTree * rec1 = stochasticBacktrack(agt, ast, R);	
//				EditTree * rec2 = stochasticBacktrack(bgt, ast, R);	
//				return new EditTree(formatLabel("Spec_aa",s),rec1,rec2,SPEC_AA_TYPE);
//			}
//		}
		if ((ag!=-1) && (as!=-1) && (bg!=-1) && (bs!=-1)) 
		{ 
			r -= PLUS(R[ag][as],R[bg][bs]);
			if (r<0) 
			{ 
				EditTree * rec1 = stochasticBacktrack(agt, ast, R);	
				EditTree * rec2 = stochasticBacktrack(bgt, bst, R);	
				return new EditTree(formatLabel("Spec_ab",s),rec1,rec2,SPEC_AB_TYPE);
			}
		}
		if ((ag!=-1) && (as!=-1) && (bg!=-1) && (bs!=-1))
		{ 
			r -= PLUS(R[ag][bs],R[bg][as]);
			if (r<0) 
			{ 
				EditTree * rec1 = stochasticBacktrack(agt, bst, R);	
				EditTree * rec2 = stochasticBacktrack(bgt, ast, R);	
				return new EditTree(formatLabel("Spec_ba",s),rec1,rec2,SPEC_BA_TYPE);
			}
		}
//		if ((ag!=-1) && (bg!=-1) && (bs!=-1)) 
//		{ 
//			r -= PLUS(R[ag][bs],R[bg][bs]);
//			if (r<0) 
//			{ 
//				EditTree * rec1 = stochasticBacktrack(agt, bst, R);	
//				EditTree * rec2 = stochasticBacktrack(bgt, bst, R);	
//				return new EditTree(formatLabel("Spec_bb",s),rec1,rec2,SPEC_BB_TYPE);
//			}
//		}
		if ((as!=-1))
		{ 
			r -= PLUS(LOSS_COST,R[i][as]);
			if (r<0) 
			{ 
				EditTree * rec = stochasticBacktrack(g, ast, R);	
				return new EditTree(formatLabel("Loss_b",s),rec,NULL,LOSS_B_TYPE);
			}
		}
		if ((bs!=-1))
		{
			r -= PLUS(LOSS_COST,R[i][bs]);
			if (r<0) 
			{ 
				EditTree * rec = stochasticBacktrack(g, bst, R);	
				return new EditTree(formatLabel("Loss_a",s),NULL,rec,LOSS_A_TYPE);
			}
		}
	}
	cout << "Error: Could not backtrack for subtrees "<<g << " and "<<s<<endl;
	exit(2);
}


