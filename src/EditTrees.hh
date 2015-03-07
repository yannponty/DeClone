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

#include <vector>
#include <string>
#include <iostream>
#include "Trees.hh"

#ifndef EDIT_TREES_HH
#define EDIT_TREES_HH

using namespace std;

typedef enum {
	DUP_TYPE, 
	LOSS_A_TYPE, 
	LOSS_B_TYPE, 
	//SPEC_AA_TYPE, 
	SPEC_AB_TYPE, 
	SPEC_BA_TYPE, 
	//SPEC_BB_TYPE, 
	MATCH_TYPE, 
	NONE_TYPE} NodeType;

ostream& operator<<(ostream &o , NodeType t);

string prettyOperationType(NodeType nt);

class EditTree: public Tree{
	private:
		NodeType type;

		void eventsAsMatrixRec(Tree * geneTree, Tree * speciesTree,NodeType ** vals);
		
	public:
		EditTree(string lbl, NodeType t);
		EditTree(string lbl, EditTree * left, EditTree * right);
		EditTree(string lbl, EditTree * left, EditTree * right, NodeType nt);
		
		EditTree * getLeft();
		EditTree * getRight();
		void setLeft(EditTree * t);
		void setRight(EditTree * t);

		NodeType ** eventsAsMatrix(Tree * geneTree, Tree * speciesTree);
		
		NodeType getNodeType();
		
};


vector<Tree *> applyEditTree(Tree * speciesTree, EditTree * ops);

#endif
