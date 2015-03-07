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

#include "Trees.hh"
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "math.h"
#include "utils.hh"

using namespace std;

Tree::Tree(string lbl, Tree * left, Tree * right){
  this->label = lbl;
  setLeft(left);
  setRight(right);
  this->index = -1;
}

Tree::Tree(string lbl){
  this->label = lbl;
  this->left = NULL;
  this->right = NULL;
  this->parent = NULL;
  this->index = -1;
}

Tree::Tree(Tree * t){
  this->label = t->getLabel();
  if (t->getLeft()!=NULL){
  	this->left = new Tree(t->getLeft());
	}
	else{
  	this->left = NULL;
	}
  if (t->getRight()!=NULL){
  	this->right = new Tree(t->getRight());
	}
	else{
  	this->right = NULL;
	}
  this->index = t->getIndex();	
}


Tree::~Tree(){
	if (this->left)
	{ delete this->left; }
	if (this->right)
  { delete this->right; }
}

void Tree::setLabel(string s){
	label = s;
}
string Tree::getLabel(){
	return label;
}

void Tree::setComment(string s){
	comment = s;
}
string Tree::getComment(){
	return comment;
}


void Tree::setParent(Tree * t){
	parent = t;
}
Tree * Tree::getParent(){
	return parent;
}

void Tree::setLeft(Tree * t){
	if (t!=NULL){
		t->setParent(this);
	}
	left = t;
}
Tree * Tree::getLeft(){
	return left;
}

void Tree::setRight(Tree * t){
	if (t!=NULL){
		t->setParent(this);
	}
	right = t;
}	
Tree * Tree::getRight(){
	return right;
}

string Tree::getND(){
     return nd;
}
void Tree::setND(string n){
     nd = n;
}

bool Tree::isLeaf(){
	return (left==NULL)&&(right==NULL);
}

bool Tree::isRoot(){
	return (parent==NULL);
}


int Tree::getIndex(){
	return index;
}

void  Tree::setIndex(int i ){
	index = i;
}

int Tree::size(){
	int count = 1;
	if (left!=NULL)
	{
		count += left->size();
	}
	if (right!=NULL)
	{
		count += right->size();
	}
	return count;
}

int Tree::height(){
	int tmp = 1;
	if (left!=NULL)
	{
		tmp = max(tmp,1+left->height());
	}
	if (right!=NULL)
	{
		tmp = max(tmp,1+right->height());
	}
	return tmp;
}


void Tree::getDupAnnotations(vector<Tree*> & dupNodes){
	if (comment.find("D=Y") != std::string::npos)
	{
		dupNodes.push_back(this);
	}
	if (left!=NULL)
	{
		left->getDupAnnotations(dupNodes);
	}
	if (right!=NULL)
	{
		right->getDupAnnotations(dupNodes);
	}
}

void Tree::show(bool showEOLs, ostream & o){
	this->show(showEOLs,0, o);
}

void Tree::show(bool showEOLs,int level, ostream & o){
     if (showEOLs) 
     {
	  for (int i=0;i<level;i++)
		    o<<" ";
	  o<<"|->["<<label<<","<< comment<<"]";
	  if (showEOLs) o << endl;
	  if (left)
	  {
	       left->show(showEOLs,level+1,o);
	  }
	  if (right)
	  {
	       right->show(showEOLs,level+1,o);
	  }

     }
     else
     {
	
	  for (int i=0;i<level;i++)
	       o<<" ";
	  if (left||right)
	  {
	       o<<"(";
	       if (left)
	       {
		    left->show(showEOLs,0,o);
	       }
	       if (right)
	       {
		    o<<",";	
		    right->show(showEOLs,0,o);
	       }
	       o<<")";	
	  }
	  o<<label;
	  if (comment!="")
	  {
	       o<<"["<<comment<<"]";
	  }	
     }
}

void Tree::asGraphViz(string path)
{
	ofstream out(path.c_str());
	asGraphViz(out, true);
}

void Tree::asGraphViz(ostream & o, bool firstCall)
{
	if (firstCall){
	  o << "digraph {"<<endl<< "node [shape=plaintext];"<<endl;
	}
	o << "\"" <<this<<"\" [label=\""<<getLabel()<<"\"];"<<endl;
	
	if (left!= NULL){
	  o << "\"" <<this<<"\" -> \""<< left<< "\";"<<endl;
	  left->asGraphViz(o,false);
	}
	if (right!= NULL){
	  o << "\"" <<this<<"\" -> \""<< right<< "\";"<<endl;
	  right->asGraphViz(o,false);
	}
	if (firstCall){
	  o << "}"<<endl;
	}
}

void Tree::asPDF(string path)
{
	asGraphViz(path+string(".dot"));
	system((string("dot -Tpdf -o ")+path+string(" ")+path+string(".dot")).c_str());
	system((string("rm ")+path+string(".dot")).c_str());	
}

void Tree::asPNG(string path)
{
	asGraphViz(path+string(".dot"));
	system((string("dot -Tpng -o ")+path+string(" ")+path+string(".dot")).c_str());
	system((string("rm ")+path+string(".dot")).c_str());	
}

void Tree::asJPEG(string path)
{
	asGraphViz(path+string(".dot"));
	system((string("dot -Tjpeg -o ")+path+string(" ")+path+string(".dot")).c_str());
	system((string("rm ")+path+string(".dot")).c_str());	
}

void Tree::asPicture(string path, string type)
{
	if (type=="pdf")
	{
		asPDF(path+string(".")+type);
	}
	else if (type=="png")
	{
		asPNG(path+string(".")+type);
	}
	else if (type=="jpg" || type=="jpeg")
	{
		asJPEG(path+string(".")+type);
	}
}


void reportError(string txt,int nbchar){
		cerr << "Error at char#"<<(nbchar+1)<<":"<<txt<<endl;
}

void reportWarning(string txt,int nbchar){
		cerr << "Warning at char#"<<(nbchar+1)<<": "<<txt<<endl;
}

/*
(,,(,));                               no nodes are named
(A,B,(C,D));                           leaf nodes are named
(A,B,(C,D)E)F;                         all nodes are named
(:0.1,:0.2,(:0.3,:0.4):0.5);           all but root node have a distance to parent
(:0.1,:0.2,(:0.3,:0.4):0.5):0.0;       all have a distance to parent
(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);       distances and leaf names (popular)
(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;     distances and all names
((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;    a tree rooted on a leaf node (rare)
*/

Tree * parseNewickTree(const string & s, int & i);

void parseNewickChildren(const string & s, int & i, Tree * father){
	int ibck = i;
	//cerr<<"  children:"<<i<<endl;
	// Assumes at least two children
	father->setLeft(parseNewickTree(s, i));
	char c = s[i];
	
	if (c!=','){ 
		reportWarning("Single tree node detected. Proceed at your own risk...",ibck+1);
		father->setRight(NULL);
		if (c==')'){
			i++;
		}
	}
	else{
		// Consumes 'comma' char 
		i++;
		father->setRight(parseNewickTree(s, i));
		while(i<s.size())
		{
	  	c = s[i];
			switch(c)
			{
				case ')':
					// Finished parsing
					i++;
					return;
				case ',':
				{
					i++;
					// Parse new subtree and put it as the right of last created node
					Tree * newTree = parseNewickTree(s, i);
					Tree * newFather = new Tree("*");
					Tree * backupRight = father->getRight();
					newFather->setLeft(backupRight);
					newFather->setRight(newTree);
					father->setRight(newFather);
					father = newFather;
					break;
				}
				default:
				{
					reportError("Unexpected char ",i+1);
					break;
				}
			}
		}	
		reportError("Unmatched parenthesis",ibck);
		return ;	
	}
}

string parseNewickComment(const string & s, int & i)
{
	string comment;
	while(i<s.size())
	{
		const char c = s[i];
		switch(c)
		{
			case ']':
				i++;
				return comment;
		}
		comment += c;
		i++;
	}
	return comment;
}
Tree * parseNewickTree(const string & s, int & i)
{
	Tree * result = new Tree("");
	bool inCaption = true;
	while(i<s.size())
	{
		const char c = s[i];
		switch(c)
		{
			case '(':
			{
				i++;
				parseNewickChildren(s,i,result);
				break;
			}
			case ',':
				return result;
			case ')':
				return result;
			case ';':
				return result;
			case ':':
				inCaption = false;
				i++;
				break;
			case '[':
				inCaption = false;
				i++;
				result->setComment(parseNewickComment(s, i));
				break;
			default:
				if (inCaption){
				  result->setLabel(result->getLabel()+s.substr(i,1));
				}
				else{
						// Edge length currently ignored
				}
				i++;
				break;
		}
	}
  return result;
}

Tree* parseNewickTree(string s)
{
	int i=0;
	string content = s;	
	  
	if (existsFile(s))
	{
		std::ifstream ifs(s.c_str());
  	content = string(
		  (std::istreambuf_iterator<char>(ifs) ),
      (std::istreambuf_iterator<char>()    ) );
	}
	Tree * t =parseNewickTree(content,i);
	computeDepthFirstOrder(t);
	return t;
}

void computeInfixOrderRec(Tree * t, vector<Tree*> & result)
{	
	if (t!=NULL)
	{
		computeInfixOrderRec(t->getRight(), result);
		result.push_back(t);
		computeInfixOrderRec(t->getLeft(), result);
	}
}

vector<Tree*> computeInfixOrder(Tree * t)
{
	vector<Tree *> result;
	computeInfixOrderRec(t, result);
	return result;
}



void computeDepthFirstOrderRec(Tree * t, vector<Tree*> & result)
{	
	if (t!=NULL)
	{
		computeDepthFirstOrderRec(t->getRight(), result);
		computeDepthFirstOrderRec(t->getLeft(), result);
		result.push_back(t);
	}
}

vector<Tree*> computeDepthFirstOrder(Tree * t)
{
	vector<Tree *> result;
	computeDepthFirstOrderRec(t, result);
	for(int i=0;i<result.size();i++) 
	{
		Tree * t = result[i];
		t->setIndex(i);
	}
	return result;
}

string formatLabel(string name, Tree * s)
{
	string res = name;
	if (s->getLabel()!="")
	{
		res += string("[")+s->getLabel()+string("]");
	}
	return res;
}


bool extentCompatible(Tree * g, Tree * s)
{
	const string & slbl = s->getLabel();
	//cout <<"  "<<slbl<<" vs "<<g->getComment()<<endl;
	if (slbl == g->getLabel())
	{
		return true;
	}
	else
	{
		const string & gcom = g->getComment();
	  std::size_t found = gcom.find(slbl);
  	if (found!=std::string::npos)
		{
			if (gcom.substr(found-2,2)=="S=")
			return true;
		}
	}
	return false;
};

struct stringbuilder
{
   stringstream ss;
   template<typename T>
   stringbuilder & operator << (const T &data)
   {
        ss << data;
        return *this;
   }
   operator std::string() { return ss.str(); }
};

void drawProbasSVG(Tree * GeneTree, Tree * SpeciesTree, double** probas, const string & output)
{
	SVGFile svg(output);
	double space = 25.;
  string name1 = "Tree 1";
  string name2 = "Tree 2";

	vector<Tree*> InfG = computeInfixOrder(GeneTree);
	vector<Tree*> InfS = computeInfixOrder(SpeciesTree);

	vector<int> heightsG(InfG.size());
	vector<int> posG(InfG.size());
	for(int i=0;i<InfG.size();i++) 
	{
			Tree * g = InfG[i];
			heightsG[g->getIndex()] = g->height();
			posG[g->getIndex()] = i;
	}
	vector<int> heightsS(InfS.size());
	vector<int> posS(InfS.size());
	for(int j=0;j<InfS.size();j++) 
	{
			Tree * s = InfS[j];
			heightsS[s->getIndex()] = s->height();
			posS[s->getIndex()] = j;
	}

  Point delta = space*Point(2.,1.+heightsG[GeneTree->getIndex()]);	

  Point xAxis = space*Point(((double)InfG.size()),0.);
  Point yAxis = space*Point(0.,((double)InfS.size()));

  //EditTree* parsimonious = computeMaxParsimony(GeneTree, SpeciesTree);
  //NodeType** m = parsimonious->eventsAsMatrix(GeneTree,SpeciesTree);
	for(int i=0;i<InfG.size();i++) 
	{
		Tree * g = InfG[i];
		for(int j=0;j<InfS.size();j++) 
		{
			Tree * s = InfS[j];
			double prob = probas[g->getIndex()][s->getIndex()];
			//svg.setColor(Color(1.-prob,1.-prob,1.-prob));
			svg.setColor(Color(1.-prob,.7*prob+1.-prob,1.-prob));
			svg.fillSquare(delta+space*Point(i,j),space);
		}
	}
	// for(int i=0;i<InfG.size();i++) 
	// {
	// 	Tree * g = InfG[i];
	// 	for(int j=0;j<InfS.size();j++) 
	// 	{
	// 		Tree * s = InfS[j];
	// 		bool isParDup = (m[g->getIndex()][s->getIndex()]==DUP_TYPE);
	// 		if(isParDup)
	// 		{
	// 			svg.setColor(Color(1.,.4,.4,.5));
	// 			svg.fillCircle(delta+space*Point(i,j),.7*space);				
	// 			double prob = probas[g->getIndex()][s->getIndex()];
	// 			svg.setColor(Color(1.-prob,1.-prob,1.-prob));
	// 			svg.fillSquare(delta+space*Point(i,j),space);
	// 			svg.setColor(Color(1.,0.,0.));
	// 			svg.drawCircle(delta+space*Point(i,j),.7*space);
	// 		}
	// 	}
	// }

	// freeMatrix(GeneTree, SpeciesTree, m);


	
	svg.setColor(Color(0.,0.,0.));

	for(int i=0;i<InfG.size();i++)
	{
		Tree * g = InfG[i];
		Point p = Point(i, -heightsG[g->getIndex()]);
			
		if (!g->isLeaf())
		{
			Point m = p+Point(0.,.5);
			svg.drawLine(delta+space*p,delta+space*m);
			Point dl = m+Point(posG[g->getLeft()->getIndex()]-posG[g->getIndex()],0.);
			svg.drawLine(delta+space*m,delta+space*dl);
			Point fl = dl+Point(0.,-(heightsG[g->getLeft()->getIndex()]-heightsG[g->getIndex()] +.5));
			svg.drawLine(delta+space*dl,delta+space*fl);
			Point dr = m+Point(posG[g->getRight()->getIndex()]-posG[g->getIndex()],0.);
			svg.drawLine(delta+space*m,delta+space*dr);
			Point fr = dr+Point(0.,-(heightsG[g->getRight()->getIndex()]-heightsG[g->getIndex()] +.5));
			svg.drawLine(delta+space*dr,delta+space*fr);
		}
	}

	for(int i=0;i<=InfG.size();i++) 
	{
		Point orig = space*Point(i-.5,-.5);
		svg.drawLine(delta+orig,delta+orig+yAxis);
		if (i<InfG.size())
		{
			Tree * g = InfG[i];
			string lbl = g->getLabel();
			stringbuilder s;
			s << (i+1);
			lbl = s;
			if (lbl!="")
			{
				Point p = Point(i,-heightsG[g->getIndex()]);
				svg.setColor(Color(1.,1.,1.));
				svg.fillSquare(delta+space*p,space/1.5);
				svg.setColor(Color(0.,0.,0.));
				svg.drawString(delta+space*p,lbl,16.);
			}
		}
	}
	

	for(int j=0;j<InfS.size();j++)
	{
		Tree * s = InfS[j];
		Point p = Point(InfG.size()-1+heightsS[s->getIndex()],j);
		if (!s->isLeaf())
		{
			Point m = p+Point(-.5,0);
			svg.drawLine(delta+space*p,delta+space*m);
			Point dl = m+Point(0.,posS[s->getLeft()->getIndex()]-posS[s->getIndex()]);
			svg.drawLine(delta+space*m,delta+space*dl);
			Point fl = dl+Point(heightsS[s->getLeft()->getIndex()]-heightsS[s->getIndex()] +.5,0);
			svg.drawLine(delta+space*dl,delta+space*fl);
			Point dr = m+Point(0.,posS[s->getRight()->getIndex()]-posS[s->getIndex()]);
			svg.drawLine(delta+space*m,delta+space*dr);
			Point fr = dr+Point(heightsS[s->getRight()->getIndex()]-heightsS[s->getIndex()] +.5,0);
			svg.drawLine(delta+space*dr,delta+space*fr);
		}
	}

	for(int j=0;j<=InfS.size();j++)
	{
		Point orig = space*Point(-.5,j-.5);
		svg.drawLine(delta+orig+xAxis,delta+orig);
		if (j<InfS.size())
		{
			Tree * s = InfS[j];
			Point p = Point(InfG.size()-1+heightsS[s->getIndex()],j);
			svg.setColor(Color(1.,1.,1.));
			svg.fillSquare(delta+space*p,space/1.5);
			svg.setColor(Color(0.,0.,0.));

			string lbl = s->getLabel();
			stringbuilder sb;
			sb << (j+1);
			lbl = sb;
			if (lbl!="")
			{
				svg.drawStringRotated90(delta+space*p,lbl,16.);
			}
		}
	}

  double dx = InfG.size()+heightsS[SpeciesTree->getIndex()];
  double dy = -heightsG[GeneTree->getIndex()];
  dy += 2.;
	//svg.drawStringLeftAlign(delta + space*Point(dx,dy), string("Caption"),25.);
  //dy += 1.5;
	svg.drawStringLeftAlign(delta + space*Point(dx,dy), name1+string(" (horizontal)"),19.);
  dy += 1.2;
	for(int i=0;i<InfG.size();i++) 
	{
		Tree * g = InfG[i];
		if (g->getLabel()!="")
		{
			stringbuilder sb;
			sb << (i+1) << ": " <<g->getLabel();
			svg.drawStringLeftAlign(delta + space*Point(dx,dy), sb,12.);
			dy+=.7;
		}
	}
  dx = InfG.size()+heightsS[SpeciesTree->getIndex()]+10;
  dy = -heightsG[GeneTree->getIndex()]+2;
	svg.drawStringLeftAlign(delta + space*Point(dx,dy), name2+string(" (vertical)"),19.);
  dy += 1.2;
	for(int j=0;j<InfS.size();j++) 
	{
		Tree * s = InfS[j];
		if (s->getLabel()!="")
		{
			stringbuilder sb;
			sb << (j+1) << ": " <<s->getLabel();
			svg.drawStringLeftAlign(delta + space*Point(dx,dy), sb,12.);
			dy+=.7;
		}
	}

}
