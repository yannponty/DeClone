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

#include "AdjacencyTrees.hh"
#include <iostream>
#include <cstdlib>
#include <sstream>
#include "utils.hh"



AdjacencyTree::AdjacencyTree(string lbl, string g1, string g2, bool adj) : Tree(lbl){
     gene1 = g1; 
     gene2 = g2; 
     isAdj = adj;
}

AdjacencyTree::AdjacencyTree(string lbl, AdjacencyTree * left, AdjacencyTree * right, string g1, string g2, bool adj) : Tree(lbl,left,right){
     gene1 = g1; 
     gene2 = g2; 
     isAdj = adj;
}

AdjacencyTree * AdjacencyTree::getLeft(){
     return (AdjacencyTree *) Tree::getLeft();
}

AdjacencyTree * AdjacencyTree::getRight(){
     return (AdjacencyTree *) Tree::getRight();
}

AdjacencyTree * AdjacencyTree::getParent(){
     return (AdjacencyTree *) Tree::getParent();
}

void AdjacencyTree::setLeft(AdjacencyTree * t){
     Tree::setLeft(t);
}

void AdjacencyTree::setRight(AdjacencyTree * t){
     Tree::setRight(t);
}

void AdjacencyTree::setParent(AdjacencyTree * t){
     Tree::setParent(t);
}

vector< pair<string, string> > AdjacencyTree::listAdjacencies(){
     vector< pair<string, string> > resultList;

     vector< pair<string, string> > resultLeft;
     if (getLeft())
     {
	  resultLeft = getLeft()->listAdjacencies();
     }
     
     vector< pair<string, string> > resultRight;
     if (getRight())
     {
	  resultRight = getRight()->listAdjacencies();
     }
     
     if (isAdj)
     {
	  resultList.push_back(pair<string, string>(gene1, gene2));
     }
     resultList.insert(resultList.end(), resultLeft.begin(), resultLeft.end());
     resultList.insert(resultList.end(), resultRight.begin(), resultRight.end());

     return resultList;
}

void AdjacencyTree::showB(bool showEOLs, ostream & o){
     this->showB(showEOLs,0, o);
}

void AdjacencyTree::showB(bool showEOLs,int level, ostream & o){
     o<<"["<<gene1<<"-"<< gene2<<"] Adj:"<<isAdj;
     if (showEOLs) 
     {
	  for (int i=0;i<level;i++)
		    o<<" ";
	  o<<"|->["<<gene1<<"-"<< gene2<<"] Adj:"<<isAdj;
	  if (showEOLs) o << endl;
	  if (getLeft())
	  {
	       getLeft()->showB(showEOLs,level+1,o);
	  }
	  if (getRight())
	  {
	       getRight()->showB(showEOLs,level+1,o);
	  }

     }
     else
     {
	
	  for (int i=0;i<level;i++)
	       o<<" ";
	  if (getLeft()||getRight())
	  {
	       o<<"(";
	       if (getLeft())
	       {
		    getLeft()->showB(showEOLs,0,o);
	       }
	       if (getRight())
	       {
		    o<<",";	
		    getRight()->showB(showEOLs,0,o);
	       }
	       o<<")";	
	  }
	  o<<getLabel();
	  if (getComment()!="")
	  {
	       o<<"["<<getComment()<<"]";
	  }	
     }
}

void AdjacencyTree::showMatrix(vector<RecTree*> Dfo1, vector<RecTree*> Dfo2, ostream & o)
{
     int M[Dfo1.size()][Dfo2.size()];      // Adjacency matrix

     // Initialize values
     for (int i = 0; i < Dfo1.size(); i++)
     {
	  for (int j = 0; j < Dfo2.size(); j++)
	  {
	       M[i][j] = 0;
	  }
     }

     // Get adjacencies
     vector< pair<string, string> > adjList = this->listAdjacencies();
     for (int k = 0; k < adjList.size(); k++)
     {
	  // Get gene labels
	  string g1 = adjList[k].first;
	  string g2 = adjList[k].second;

	  // Look for position
	  int i;
	  int j;
	  for (i = 0; i < Dfo1.size(); i++)
	  {
	       if (g1 == Dfo1[i]->getND())
		    break;
	  }
	  for (j = 0; j < Dfo2.size(); j++)
	  {
	       if (g2 == Dfo2[j]->getND())
		    break;
	  }

	  M[i][j] = 1;
     }


     // Print matrix
     o << "\t";
     for (int j = 0; j < Dfo2.size(); j++)
     {
	  o << Dfo2[j]->getND() << " ";
     }
     o << endl;
     for (int i = 0; i < Dfo1.size(); i++)
     {
	  o << Dfo1[i]->getND() << "\t";
	  for (int j = 0; j < Dfo2.size(); j++)
	  {
	       o << M[i][j] << " "; 
	  }
	  o << endl;
     }
}


// string AdjacencyTree::adjacencyForestString()
// {
//      string result = "";
//      string leftForest = "";
//      string rightForest = "";

//      if (getLeft())
// 	  leftForest = getLeft()->adjacencyForestString();
//      if (getRight())
// 	  rightForest = getRight()->adjacencyForestString();
	 
//      if (isAdj)
//      {
// 	  // Print children
// 	  result += "(";
// 	  result += leftForest;
// 	  if ((leftForest != "") && (rightForest != ""))
// 	       result += ",";
// 	  result += rightForest;
// 	  result += ")";
	       
// 	  // Print node
// 	  result += gene1 + "-" + gene2;
//      }
//      else
//      {
// 	  result += leftForest;
// 	  if ((leftForest != "") && (rightForest != ""))
// 	       result += "\n";
// 	  result += rightForest;
//      }

//      return result;
// }

string AdjacencyTree::adjacencyListString()
{
     string result = "";
     string leftForest = "";
     string rightForest = "";

     if (getLeft())
	  leftForest = getLeft()->adjacencyListString();
     if (getRight())
	  rightForest = getRight()->adjacencyListString();
	 
     if (isAdj)
     {
	  string adjStr = gene1 + "-" + gene2;

	  // // Look for repeated adjacencies
	  // size_t found = (leftForest + " " + rightForest).find(adjStr);
	  // if (found != std::string::npos)
	  //      result += "REPEATED";

	  // Print node
	  result += adjStr + " ";
     }

     // Print children
     result += leftForest;
     result += rightForest;

     return result;
}

string AdjacencyTree::adjacencySpecListString()
{
     string result = "";
     string leftForest = "";
     string rightForest = "";

     //cerr << getLabel() << " " << isAdj << endl;

     if (getLeft())
	  leftForest = getLeft()->adjacencySpecListString();
     if (getRight())
	  rightForest = getRight()->adjacencySpecListString();
	 
     if ((isAdj) && (getLabel().find("SPEC-SPEC") != string::npos))
     {
	  // Print node
	  result += gene1 + "-" + gene2 + " ";
     }

     // Print children
     result += leftForest;
     result += rightForest;

     return result;
}

ostream& operator<<(ostream & o, pair<AdjacencyTree*,double> v)
{
     o <<v.first;
     o <<"Score: "<< v.second<<endl;
     return o;
}

ostream& operator<<(ostream & o, AdjacencyTree* v)
{
     o <<"AdjacencyTree: ";
     v->showB(false,o);
     o <<";"<<endl;
     o <<"AncestralAdjacencies: " << v->adjacencyListString() << endl;
     o <<"SpeciationAdjacencies: " << v->adjacencySpecListString() << endl;
     return o;
}

ostream& operator<<(ostream & o, const vector<pair<AdjacencyTree*,double> > & v)
{
	  for(int i=0;i<v.size();i++)
	  {
	       o << v[i];
	  }
     return o;
}

