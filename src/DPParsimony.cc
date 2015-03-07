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

#include "DPParsimony.hh"
#include <iostream>

#define PLUS(a,b) (a+b)
#define MIN(a,b) min(a,b)
#define ZERO 0
#define TYPEDATA double
#define INFTY DBL_MAX
#define DUP_COST 1.
#define LOSS_COST 1.

#define allocateMatricesTASK allocateMatricesMaxParsimony
#define fillMatricesTASK fillMatricesMaxParsimony
#define deleteMatricesTASK deleteMatricesMaxParsimony

#include "DPRaw.cc"

EditTree * backtrackMaxParsimony(Tree * g, Tree * s,TYPEDATA** R);


EditTree * computeMaxParsimony(Tree * GeneTree, Tree * SpeciesTree)
{
	TYPEDATA** R = allocateMatricesTASK(GeneTree, SpeciesTree);
  fillMatricesTASK(GeneTree, SpeciesTree, R);
  EditTree * t = NULL;
  if (R[GeneTree->getIndex()][SpeciesTree->getIndex()] != INFTY)
  {
		t = backtrackMaxParsimony(GeneTree, SpeciesTree, R);
	}
		deleteMatricesTASK(GeneTree, SpeciesTree, R);
		return t;
}


EditTree * backtrackMaxParsimony(Tree * g, Tree * s,TYPEDATA** R)
{
	int i = g->getIndex();
	int j = s->getIndex();
	TYPEDATA tmp = INFTY;
	int ag = (g->getLeft()?  g->getLeft()->getIndex():-1);
	int bg = (g->getRight()? g->getRight()->getIndex():-1);
	int as = (s->getLeft()?  s->getLeft()->getIndex():-1);
	int bs = (s->getRight()? s->getRight()->getIndex():-1);
	Tree* agt = g->getLeft() ;
	Tree* bgt = g->getRight();
	Tree* ast = s->getLeft() ;
	Tree* bst = s->getRight();
	if (g->isLeaf() && s->isLeaf())
	{				
			if (extentCompatible(g, s))
			{ 
				//cout << "Match("<<g->getLabel()<< "," << s->getLabel()<<")";
				return new EditTree(formatLabel("Match",s),MATCH_TYPE);
			}
	}
	else
	{
		if ((ag!=-1) && (bg!=-1) && (R[i][j] == PLUS(DUP_COST,PLUS(R[ag][j],R[bg][j]))))
		{ 
			//cout << "Dup_" <<s->getLabel()<< "(";
			EditTree * rec1 = backtrackMaxParsimony(agt, s, R);	
			EditTree * rec2 = backtrackMaxParsimony(bgt, s, R);	
			//cout << ")";
			return new EditTree(formatLabel("Dup",s),rec1,rec2,DUP_TYPE);
		}
//		if ((ag!=-1) && (as!=-1) && (bg!=-1) && (R[i][j] == PLUS(R[ag][as],R[bg][as])))
//		{ 
//			//cout << "Map_(a_g->a_s,b_g->a_s)" <<s->getLabel()<< "(";
//			EditTree * rec1 = backtrackMaxParsimony(agt, ast, R);	
//			EditTree * rec2 = backtrackMaxParsimony(bgt, ast, R);	
//			//cout << ")";
//			return new EditTree(formatLabel("Spec_aa",s),rec1,rec2,SPEC_AA_TYPE);
//		}
		if ((ag!=-1) && (as!=-1) && (bg!=-1) && (bs!=-1) && (R[i][j] == PLUS(R[ag][as],R[bg][bs]))) 
		{ 
			//cout << "Map_(a_g->a_s,b_g->b_s)" <<s->getLabel()<< "(";
			EditTree * rec1 = backtrackMaxParsimony(agt, ast, R);	
			EditTree * rec2 = backtrackMaxParsimony(bgt, bst, R);	
			//cout << ")";
			return new EditTree(formatLabel("Spec_ab",s),rec1,rec2,SPEC_AB_TYPE);
		}
		if ((ag!=-1) && (as!=-1) && (bg!=-1) && (bs!=-1) && (R[i][j] == PLUS(R[ag][bs],R[bg][as])))
		{ 
			//cout << "Map_(a_g->b_s,b_g->a_s)" <<s->getLabel()<< "(";
			EditTree * rec1 = backtrackMaxParsimony(agt, bst, R);	
			EditTree * rec2 = backtrackMaxParsimony(bgt, ast, R);	
			//cout << ")";
			return new EditTree(formatLabel("Spec_ba",s),rec1,rec2,SPEC_BA_TYPE);
		}
//		if ((ag!=-1) && (bg!=-1) && (bs!=-1) && (R[i][j] == PLUS(R[ag][bs],R[bg][bs]))) 
//		{ 
//			//cout << "Map_(a_g->b_s,b_g->b_s)" <<s->getLabel()<< "(";
//			EditTree * rec1 = backtrackMaxParsimony(agt, bst, R);	
//			EditTree * rec2 = backtrackMaxParsimony(bgt, bst, R);	
//			//cout << ")";
//			return new EditTree(formatLabel("Spec_bb",s),rec1,rec2,SPEC_BB_TYPE);
//		}
		if ((as!=-1) && (R[i][j] == PLUS(LOSS_COST,R[i][as])))
		{ 
			//cout << "Loss_b_s" <<s->getLabel()<< "(";
			EditTree * rec = backtrackMaxParsimony(g, ast, R);	
			//cout << ")";
			return new EditTree(formatLabel("Loss_b",s),rec,NULL,LOSS_B_TYPE);
		}
		if ((bs!=-1) && (R[i][j] == PLUS(LOSS_COST,R[i][bs])))
		{ 
			//cout << "Loss_a_s" <<s->getLabel()<< "(";
			EditTree * rec = backtrackMaxParsimony(g, bst, R);	
			//cout << ")";
			return new EditTree(formatLabel("Loss_a",s),NULL,rec,LOSS_A_TYPE);
		}
	}
	cerr << "Error: Could not backtrack for subtrees "<<g << " and "<<s<<endl;
	return NULL;
}


