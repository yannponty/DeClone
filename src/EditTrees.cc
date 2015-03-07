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

#include "EditTrees.hh"
#include <iostream>
#include <cstdlib>

string prettyOperationType(NodeType nt){
		switch(nt)
		{
			case (DUP_TYPE): return string("Dup");
			case (LOSS_A_TYPE): return string("Loss_a");
			case (LOSS_B_TYPE): return string("Loss_b");
			//case (SPEC_AA_TYPE): return string("Spec_aa");
			case (SPEC_AB_TYPE): return string("Spec_ab");
			case (SPEC_BA_TYPE): return string("Spec_ba");
			//case (SPEC_BB_TYPE): return string("Spec_bb");
			case (MATCH_TYPE): return string("match");
			default:
				return string("n/a");
		}
}

EditTree::EditTree(string lbl, NodeType t) : Tree(lbl){
	type = t;
}

EditTree::EditTree(string lbl, EditTree * left, EditTree * right, NodeType t) : Tree(lbl,left,right){
	type = t;
}

NodeType EditTree::getNodeType()
{
	return type;
}

EditTree * EditTree::getLeft(){
	return (EditTree *) Tree::getLeft();
}

EditTree * EditTree::getRight(){
	return (EditTree *) Tree::getRight();
}

void EditTree::setLeft(EditTree * t){
	Tree::setLeft(t);
}

void EditTree::setRight(EditTree * t){
	Tree::setRight(t);
}

NodeType**  EditTree::eventsAsMatrix(Tree * geneTree,Tree * speciesTree){
  NodeType ** vals = new NodeType*[geneTree->size()];
  for(int i =0;i<geneTree->size();i++){
		vals[i] = new NodeType[speciesTree->size()];
	}
  for(int i =0;i<geneTree->size();i++){
	  for(int j =0;j<speciesTree->size();j++){
			vals[i][j] = NONE_TYPE;
		}
	}
	eventsAsMatrixRec(geneTree, speciesTree, vals);
	return vals;
}

void EditTree::eventsAsMatrixRec(Tree * geneTree, Tree * speciesTree,  NodeType ** vals){
	vals[geneTree->getIndex()][speciesTree->getIndex()] = this->getNodeType();
	switch(this->getNodeType())
	{
			case (DUP_TYPE):{
				getLeft()->eventsAsMatrixRec(geneTree->getLeft(), speciesTree,  vals);
				getRight()->eventsAsMatrixRec(geneTree->getRight(), speciesTree, vals);
			}
			break;
			case (LOSS_A_TYPE):{
				getRight()->eventsAsMatrixRec(geneTree, speciesTree->getRight(), vals);
			}
			break;
			case (LOSS_B_TYPE):{
				getLeft()->eventsAsMatrixRec(geneTree, speciesTree->getLeft(),  vals);
			}
			break;
//			case (SPEC_AA_TYPE):{
//				getLeft()->eventsAsMatrixRec(geneTree->getLeft(), speciesTree->getLeft(),  vals);
//				getRight()->eventsAsMatrixRec(geneTree->getRight(), speciesTree->getLeft(),  vals);
//			}
//			break;
			case (SPEC_AB_TYPE):{
				getLeft()->eventsAsMatrixRec(geneTree->getLeft(), speciesTree->getLeft(),  vals);
				getRight()->eventsAsMatrixRec(geneTree->getRight(), speciesTree->getRight(), vals);
			}
			break;
			case (SPEC_BA_TYPE):{
				getLeft()->eventsAsMatrixRec(geneTree->getLeft(), speciesTree->getRight(),  vals);
				getRight()->eventsAsMatrixRec(geneTree->getRight(), speciesTree->getLeft(),  vals);
			}
			break;
//			case (SPEC_BB_TYPE):{
//				getLeft()->eventsAsMatrixRec(geneTree->getLeft(), speciesTree->getRight(), vals);
//				getRight()->eventsAsMatrixRec(geneTree->getRight(), speciesTree->getRight(), vals);
//			}
//			break;
			case (MATCH_TYPE):{
			}
			break;
			default:{
				cerr << "Error: Weird edit operation type."<<endl;
				exit(3);
			}
			break;
	}
}


ostream& operator<<(ostream &o , NodeType t)
{
	switch(t)
	{
			case (DUP_TYPE):{
				o<<"D";
			}
			break;
			case (LOSS_A_TYPE):{
				o<<"LA";
			}
			break;
			case (LOSS_B_TYPE):{
				o<<"LB";
			}
			break;
//			case (SPEC_AA_TYPE):{
//				o<<"SAA";
//			}
//			break;
			case (SPEC_AB_TYPE):{
				o<<"SAB";
			}
			break;
			case (SPEC_BA_TYPE):{
				o<<"SBA";
			}
			break;
//			case (SPEC_BB_TYPE):{
//				o<<"SBB";
//			}
//			break;
			case (MATCH_TYPE):{
				o<<"M";
			}
			break;
			case (NONE_TYPE):{
				o<<"-";
			}
			break;
	}
	return o;
}


vector<Tree *> reverse(vector<Tree *> & t)
{
	vector<Tree*> res;
	for (int j=t.size()-1;j>=0;j--)
	{
		res.push_back(t[j]);
	}	
	return res;
}

vector<Tree *>  applyEditTreeRec(Tree * speciesTree, EditTree * ops)
{
	//cerr << "spe:"<<speciesTree<<" ops:"<<ops<<endl;
	vector<Tree*> res;
	if (ops!=NULL)
	{
		vector<Tree*> vl;
		vector<Tree*> vr;
  	//cerr << prettyOperationType(ops->getNodeType())<<endl;
		switch(ops->getNodeType())
		{
			case (DUP_TYPE):{
				vl = applyEditTreeRec(speciesTree, ops->getLeft());
				vr = applyEditTreeRec(speciesTree, ops->getRight());
				res.push_back(new Tree("D",new Tree(speciesTree),new Tree(speciesTree)));
				//Modifications from left edit branch only
				for (int i=0;i<vl.size();i++)
				{
					res.push_back(new Tree("D",new Tree(vl[i]),new Tree(speciesTree)));
				}
				// Modifications from right + all modifs of left edit branch
				for (int j=0;j<vr.size();j++)
				{
					res.push_back(new Tree("D",new Tree(vl[vl.size()-1]),new Tree(vr[j])));
				}
			}
			break;
			case (LOSS_A_TYPE):{
				vr = applyEditTreeRec(speciesTree->getRight(), ops->getRight());
				res.push_back(new Tree("La",new Tree(speciesTree->getRight()),NULL));
				for (int j=0;j<vr.size();j++)
				{
					res.push_back(new Tree(vr[j]));
				}
			}
			break;
			case (LOSS_B_TYPE):{
				vl = applyEditTreeRec(speciesTree->getLeft(), ops->getLeft());
				res.push_back(new Tree("Lb",new Tree(speciesTree->getLeft()),NULL));
				for (int j=0;j<vl.size();j++)
				{
					res.push_back(new Tree(vl[j]));
				}
			}
			break;
//			case (SPEC_AA_TYPE):{
//				vl = applyEditTreeRec(speciesTree->getLeft(), ops->getLeft());
//				vr = applyEditTreeRec(speciesTree->getLeft(), ops->getRight());
//				res.push_back(new Tree("Saa",new Tree(speciesTree->getLeft()),new Tree(speciesTree->getLeft())));
//				// Modifications from left edit branch
//				for (int i=0;i<vl.size();i++)
//				{
//					res.push_back(new Tree("Saa",new Tree(vl[i]),new Tree(speciesTree->getLeft())));
//				}
//				// Modifications from right + all modifs of left edit branch
//				for (int j=0;j<vr.size();j++)
//				{
//					res.push_back(new Tree("Saa",new Tree(vl[vl.size()-1]),new Tree(vr[j])));
//				}
//			}
//			break;
			case (SPEC_AB_TYPE):{
				vector<Tree*> vl = applyEditTreeRec(speciesTree->getLeft(), ops->getLeft());
				vector<Tree*> vr = applyEditTreeRec(speciesTree->getRight(), ops->getRight());
				res.push_back(new Tree("Sab",new Tree(speciesTree->getLeft()),new Tree(speciesTree->getRight())));
				// Modifications from left edit branch
				for (int i=0;i<vl.size();i++)
				{
					res.push_back(new Tree("Sab",new Tree(vl[i]),new Tree(speciesTree->getRight())));
				}
				// Modifications from right + all modifs of left edit branch
				for (int j=0;j<vr.size();j++)
				{
					res.push_back(new Tree("Sab",new Tree(vl[vl.size()-1]),new Tree(vr[j])));
				}
			}
			break;
			case (SPEC_BA_TYPE):{
				vector<Tree*> vl = applyEditTreeRec(speciesTree->getRight(), ops->getLeft());
				vector<Tree*> vr = applyEditTreeRec(speciesTree->getLeft(), ops->getRight());
				res.push_back(new Tree("Sba",new Tree(speciesTree->getRight()),new Tree(speciesTree->getLeft())));
				// Modifications from left edit branch
				for (int i=0;i<vl.size();i++)
				{
					res.push_back(new Tree("Sba",new Tree(vl[i]),new Tree(speciesTree->getLeft())));
				}
				// Modifications from right + all modifs of left edit branch
				for (int j=0;j<vr.size();j++)
				{
					res.push_back(new Tree("Sba",new Tree(vl[vl.size()-1]),new Tree(vr[j])));
				}
			}
			break;
//			case (SPEC_BB_TYPE):{
//				vector<Tree*> vl = applyEditTreeRec(speciesTree->getRight(), ops->getLeft());
//				vector<Tree*> vr = applyEditTreeRec(speciesTree->getRight(), ops->getRight());
//				res.push_back(new Tree("Sbb",new Tree(speciesTree->getLeft()),new Tree(speciesTree->getLeft())));
//				// Modifications from left edit branch
//				for (int i=0;i<vl.size();i++)
//				{
//					res.push_back(new Tree("Sbb",new Tree(vl[i]),new Tree(speciesTree->getRight())));
//				}
//				// Modifications from right + all modifs of left edit branch
//				for (int j=0;j<vr.size();j++)
//				{
//					res.push_back(new Tree("Sbb",new Tree(vl[vl.size()-1]),new Tree(vr[j])));
//				}
//			}
//			break;
			case (MATCH_TYPE):{
						res.push_back(new Tree(speciesTree->getLabel()+string("'")));
			}
			break;
			default:{
				cerr << "Error: Weird edit operation type."<<endl;
				exit(3);
			}
			break;
		}
		for (int i=0;i<vl.size();i++)
		{
			delete vl[i];
		}
		for (int i=0;i<vr.size();i++)
		{
			delete vr[i];
		}
	}
	return res;
}


class StatsAccumulator{
  public:
	long nbDup;
  long nbLoss;
  StatsAccumulator(){
	  nbDup = 0;
	  nbLoss = 0;
	}
};

vector<Tree *>  computeStatsTreeRec(Tree * speciesTree, EditTree * ops, StatsAccumulator & stats)
{
	vector<Tree*> res;
	if (ops!=NULL)
	{
		vector<Tree*> vl;
		vector<Tree*> vr;
		switch(ops->getNodeType())
		{
			case (DUP_TYPE):{
				stats.nbDup++;
				computeStatsTreeRec(speciesTree, ops->getLeft(),stats);
				computeStatsTreeRec(speciesTree, ops->getRight(),stats);
			}
			break;
			case (LOSS_A_TYPE):{
				stats.nbLoss++;
				computeStatsTreeRec(speciesTree->getRight(), ops->getRight(),stats);
			}
			break;
			case (LOSS_B_TYPE):{
				stats.nbLoss++;
				computeStatsTreeRec(speciesTree->getLeft(), ops->getLeft(),stats);
			}
			break;
//			case (SPEC_AA_TYPE):{
//				computeStatsTreeRec(speciesTree->getLeft(), ops->getLeft(),stats);
//				computeStatsTreeRec(speciesTree->getLeft(), ops->getRight(),stats);
//			}
//			break;
			case (SPEC_AB_TYPE):{
				computeStatsTreeRec(speciesTree->getLeft(), ops->getLeft(),stats);
				computeStatsTreeRec(speciesTree->getRight(), ops->getRight(),stats);
			}
			break;
			case (SPEC_BA_TYPE):{
				computeStatsTreeRec(speciesTree->getRight(), ops->getLeft(),stats);
				computeStatsTreeRec(speciesTree->getLeft(), ops->getRight(),stats);
			}
			break;
//			case (SPEC_BB_TYPE):{
//				computeStatsTreeRec(speciesTree->getRight(), ops->getLeft(),stats);
//				computeStatsTreeRec(speciesTree->getRight(), ops->getRight(),stats);
//			}
//			break;
			case (MATCH_TYPE):{
			}
			break;
			default:{
				cerr << "Error: Weird edit operation type."<<endl;
				exit(3);
			}
			break;
		}
		for (int i=0;i<vl.size();i++)
		{
			delete vl[i];
		}
		for (int i=0;i<vr.size();i++)
		{
			delete vr[i];
		}
	}
}



vector<Tree *> applyEditTree(Tree * speciesTree, EditTree * ops)
{
	return applyEditTreeRec(speciesTree, ops);
}


