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
#include "SVGDriver.hh"

#ifndef TREES_HH
#define TREES_HH

using namespace std;

class Tree{
private:
     string label;
     string comment;
     Tree * left;
     Tree * right;
     Tree * parent;
     int index;
     string nd;
		
public:
     Tree(string lbl, Tree * left, Tree * right);
     Tree(string lbl);
     Tree(Tree * t);
     ~Tree();
	
     string getLabel();
     void setLabel(string lbl);
	
     string getComment();
     void setComment(string s);
	
     void setLeft(Tree * t);
     Tree * getLeft();
	
     void setRight(Tree * t);
     Tree * getRight();
	
     void setParent(Tree * t);
     Tree * getParent();

     string getND();
     void setND(string n);
	
     bool isLeaf();
     bool isRoot();
		
     int getIndex();	
     void  setIndex(int i );
		
     int size();
     int height();
		
     void show(bool showEOLs=true,ostream & o= cout);
     void show(bool showEOLs,int level, ostream & o);
		
     void asGraphViz(string path);
     void asGraphViz(ostream & o, bool firstCall=true);

     void asPicture(string path, string type);
     void asPDF(string path);
     void asPNG(string path);
     void asJPEG(string path);
		
     void getDupAnnotations(vector<Tree*> & dupNodes);

};


Tree* parseNewickTree(string s);

vector<Tree*> computeDepthFirstOrder(Tree * t);
vector<Tree*> computeInfixOrder(Tree * t);

string formatLabel(string name, Tree * s);

bool extentCompatible(Tree * g, Tree * s);

void reportError(string txt,int nbchar);
void reportWarning(string txt,int nbchar);

void drawProbasSVG(Tree * GeneTree, Tree * SpeciesTree, double** probas, const string & output = "dotplot.svg");

#endif
