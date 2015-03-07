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

TYPEDATA** allocateMatricesTASK(Tree * GeneTree, Tree * SpeciesTree)
{
	TYPEDATA** R = new TYPEDATA*[GeneTree->size()];
	for(int i=0;i<GeneTree->size();i++) 
	{ R[i] = new TYPEDATA[SpeciesTree->size()]; }
	return R;
}

void deleteMatricesTASK(Tree * GeneTree, Tree * SpeciesTree, TYPEDATA** R)
{
	for(int i=0;i<GeneTree->size();i++) 
	{
		delete[] R[i];
	}
	delete[] R;
}

void fillMatricesTASK(Tree * GeneTree, Tree * SpeciesTree, TYPEDATA** R)
{
	vector<Tree*> DfoG = computeDepthFirstOrder(GeneTree);
	vector<Tree*> DfoS = computeDepthFirstOrder(SpeciesTree);
	for(int i=0;i<DfoG.size();i++) 
	{
		Tree * g = DfoG[i];
		for(int j=0;j<DfoS.size();j++) 
		{
			Tree * s = DfoS[j];
			// Deplacer le cas Feuille * Noeud Interne vers le cas général
			if (g->isLeaf() && s->isLeaf())
			{				
					if (extentCompatible(g, s))
					{ R[i][j] = ZERO;	}
					else
					{ R[i][j] = INFTY;	}
			}
			else
			{
				TYPEDATA tmp = INFTY;
				int ag = (g->getLeft()?  g->getLeft()->getIndex():-1);
				int bg = (g->getRight()? g->getRight()->getIndex():-1);
				int as = (s->getLeft()?  s->getLeft()->getIndex():-1);
				int bs = (s->getRight()? s->getRight()->getIndex():-1);
				// TODO: Deal with unary nodes in a better way
				if (ag!=-1 && bg!=-1)
				{ tmp = MIN(tmp,PLUS(DUP_COST,PLUS(R[ag][j],R[bg][j]))); }
//				if (ag!=-1 && bg!=-1 && as!=-1) 
//				{ tmp = MIN(tmp,PLUS(R[ag][as],R[bg][as])); }
				if (ag!=-1 && bg!=-1 && as!=-1 && bs!=-1)
				{ tmp = MIN(tmp,PLUS(R[ag][as],R[bg][bs])); }
				if (ag!=-1 && bg!=-1 && as!=-1 && bs!=-1)
				{	tmp = MIN(tmp,PLUS(R[ag][bs],R[bg][as])); }
//				if (ag!=-1 && bg!=-1 && bs!=-1)
//				{	tmp = MIN(tmp,PLUS(R[ag][bs],R[bg][bs])); }
				if (as!=-1)
				{	tmp = MIN(tmp,PLUS(LOSS_COST,R[i][as])); }
				if (bs!=-1)
				{	tmp = MIN(tmp,PLUS(LOSS_COST,R[i][bs])); }
				R[i][j] = tmp;
			}
		  //cout << "f("<<g->getLabel()<<","<<s->getLabel()<<") -> "<< R[i][j]<<endl;
		}
	}
}
