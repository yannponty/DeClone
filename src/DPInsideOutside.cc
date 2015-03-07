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

#include <iostream>
#include "DPInsideOutside.hh"
#include "math.h"
#include "SVGDriver.hh"

#define PLUS(a,b) (a*b)
#define MIN(a,b) (a+b)
#define ZERO 1.

#define RT 1.0


#define TYPEDATA double
#define INFTY 0.
#define DUP_COST exp(-1./kT)
#define LOSS_COST exp(-1./kT)

#define allocateMatricesTASK allocateMatricesInsideOutside
#define fillMatricesTASK fillMatricesInsideOutside
#define deleteMatricesTASK deleteMatricesInsideOutside

#include "DPRaw.cc"

TYPEDATA** computeOutside(Tree * GeneTree, Tree * SpeciesTree,TYPEDATA** FwR);



double getPartitionFunction(Tree * GeneTree, Tree * SpeciesTree)
{
	TYPEDATA** FwR = allocateMatricesTASK(GeneTree, SpeciesTree);
	fillMatricesTASK(GeneTree, SpeciesTree, FwR);
	double res = 	FwR[GeneTree->getIndex()][SpeciesTree->getIndex()];
	deleteMatricesTASK(GeneTree, SpeciesTree, FwR);
	return res;
}

void getProbasDupl(Tree * GeneTree, Tree * SpeciesTree, double** probas)
{
	TYPEDATA** FwR = allocateMatricesTASK(GeneTree, SpeciesTree);
	fillMatricesTASK(GeneTree, SpeciesTree, FwR);
  //cerr << "Forward:"<< FwR[GeneTree->getIndex()][SpeciesTree->getIndex()]<<endl;
	TYPEDATA** BcR = computeOutside(GeneTree,SpeciesTree,FwR);
	vector<Tree*> DfoG = computeDepthFirstOrder(GeneTree);
	vector<Tree*> DfoS = computeDepthFirstOrder(SpeciesTree);
	for(int i=0;i<DfoG.size();i++) 
	{
		Tree * g = DfoG[i];
		for(int j=0;j<DfoS.size();j++) 
		{
			Tree * s = DfoS[j];
			int ag = (g->getLeft()?  g->getLeft()->getIndex():-1);
			int bg = (g->getRight()? g->getRight()->getIndex():-1);
			int as = (s->getLeft()?  s->getLeft()->getIndex():-1);
			int bs = (s->getRight()? s->getRight()->getIndex():-1);
			probas[i][j] = 0.;
			if (ag!=-1 && bg!=-1 )
			{ 
				probas[i][j] = (DUP_COST*BcR[i][j]*FwR[ag][j]*FwR[bg][j]) / FwR[GeneTree->getIndex()][SpeciesTree->getIndex()]; 
			}
		}
	}	
	deleteMatricesTASK(GeneTree, SpeciesTree, FwR);
	deleteMatricesTASK(GeneTree, SpeciesTree, BcR);
}

TYPEDATA** computeOutside(Tree * GeneTree, Tree * SpeciesTree,TYPEDATA** FwR)
{
	TYPEDATA** R = allocateMatricesTASK(GeneTree, SpeciesTree);
	vector<Tree*> DfoG = computeDepthFirstOrder(GeneTree);
	vector<Tree*> DfoS = computeDepthFirstOrder(SpeciesTree);
	for(int i=DfoG.size()-1;i>=0;i--) 
	{
		Tree * g = DfoG[i];
		for(int j=DfoS.size()-1;j>=0;j--) 
		{
			Tree * s = DfoS[j];
			if (g->isRoot() && s->isRoot())
			{				
					R[i][j] = ZERO;
			}
			else
			{
				TYPEDATA tmp = INFTY;
				bool sIsLeftChild = (s->getParent()?(s->getParent()->getLeft()==s):false);
				bool gIsLeftChild = (g->getParent()?(g->getParent()->getLeft()==g):false);
				int pg = (g->getParent()?  g->getParent()->getIndex():-1);
				int ps = (s->getParent()?  s->getParent()->getIndex():-1);
				int xg = (g->getParent()?
				  (g==g->getParent()->getLeft()?
					   g->getParent()->getRight()->getIndex()
					  :g->getParent()->getLeft()->getIndex())
					:-1);
				int xs = (s->getParent()?
				  (s==s->getParent()->getLeft()?
					   s->getParent()->getRight()->getIndex()
					  :s->getParent()->getLeft()->getIndex())
					:-1);
				if (pg!=-1)
				{ tmp = MIN(tmp,PLUS(DUP_COST,PLUS(R[pg][j],FwR[xg][j]))); }
//				if ((pg!=-1)&&(ps!=-1)&&(sIsLeftChild))
//				{
//					tmp = MIN(tmp,PLUS(R[pg][ps],FwR[xg][j])); 
//				}
				if ((pg!=-1)&&(ps!=-1)&&!(gIsLeftChild^sIsLeftChild))
				{
					tmp = MIN(tmp,PLUS(R[pg][ps],FwR[xg][xs])); 
				}
				if ((pg!=-1)&&(ps!=-1)&&(gIsLeftChild^sIsLeftChild))
				{
					tmp = MIN(tmp,PLUS(R[pg][ps],FwR[xg][xs])); 
				}
//				if ((pg!=-1)&&(ps!=-1)&&(!sIsLeftChild))
//				{
//					tmp = MIN(tmp,PLUS(R[pg][ps],FwR[xg][j])); 
//				}
				if ((ps!=-1)&&(!sIsLeftChild))
				{	tmp = MIN(tmp,PLUS(LOSS_COST,R[i][ps])); }
				if ((ps!=-1)&&(sIsLeftChild))
				{	tmp = MIN(tmp,PLUS(LOSS_COST,R[i][ps])); }
				
				R[i][j] = tmp;
			}
		  //cout << "f("<<g->getLabel()<<","<<s->getLabel()<<") -> "<< R[i][j]<<endl;
		}
	}
	return R;
}		


