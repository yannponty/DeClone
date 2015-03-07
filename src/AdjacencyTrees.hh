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
#include <map>
#include "Trees.hh"
#include "RecTrees.hh"

#ifndef ADJ_TREES_HH
#define ADJ_TREES_HH

using namespace std;

class AdjacencyTree: public Tree{
private:
     string gene1;
     string gene2;
     bool isAdj;

public:
     AdjacencyTree(string lbl, string g1, string g2, bool adj);
     AdjacencyTree(string lbl, AdjacencyTree * left, AdjacencyTree * right, string g1, string g2, bool adj);
     
     AdjacencyTree * getLeft();
     AdjacencyTree * getRight();
     AdjacencyTree * getParent();
     void setLeft(AdjacencyTree * t);
     void setRight(AdjacencyTree * t);
     void setParent(AdjacencyTree * t);

     vector< pair<string, string> > listAdjacencies();

     void showB(bool showEOLs=true,ostream & o= cout);
     void showB(bool showEOLs,int level, ostream & o);

     void showMatrix(vector<RecTree*> Dfo1, vector<RecTree*> Dfo2, ostream & o=cout);

     string adjacencyForestString();
     string adjacencyListString();
     string adjacencySpecListString();
};


ostream& operator<<(ostream & o, pair<AdjacencyTree*,double> v);
ostream& operator<<(ostream & o, AdjacencyTree* v);
ostream& operator<<(ostream & o, const vector<pair<AdjacencyTree*,double> > & v);


#endif
