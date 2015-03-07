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

#include "RecTrees.hh"
#include <iostream>
#include <cstdlib>
#include <sstream>
#include "utils.hh"

string prettyOperationType(EventType nt){
		switch(nt)
		{
			case (GDup): return string("GDup");
			case (GLos): return string("GLoss");
			case (Spec): return string("Spec");
			case (Extant): return string("Extant");
			default:
				return string("Unknown");
		}
}

RecTree::RecTree(string lbl) : Tree(lbl){
}

RecTree::RecTree(string lbl, EventType t) : Tree(lbl){
	type = t;
}

RecTree::RecTree(string lbl, RecTree * left, RecTree * right, EventType t, string spec) : Tree(lbl,left,right){
	type = t;
	species = spec;
}

EventType RecTree::getEvent()
{
	return type;
}
void RecTree::setEvent(EventType t)
{
	type = t;
}

string RecTree::getSpecies()
{
	return species;
}

string RecTree::getGeneName()
{
    string Gname = split(this->getLabel(),'|')[0];
    return Gname;

}

void RecTree::setSpecies(string spec)
{
	species = spec;
}


RecTree * RecTree::getLeft(){
	return (RecTree *) Tree::getLeft();
}

RecTree * RecTree::getRight(){
	return (RecTree *) Tree::getRight();
}

RecTree * RecTree::getParent(){
	return (RecTree *) Tree::getParent();
}

void RecTree::setLeft(RecTree * t){
	Tree::setLeft(t);
}

void RecTree::setRight(RecTree * t){
	Tree::setRight(t);
}

void RecTree::setParent(RecTree * t){
	Tree::setParent(t);
}

ostream& operator<<(ostream &o , EventType t)
{
	switch(t)
	{
			case (GDup):{
				o<<"D";
			}
			break;
			case (GLos):{
				o<<"L";
			}
			break;
			case (Spec):{
				o<<"S";
			}
			break;
			case (Extant):{
				o<<"E";
			}
			break;
	}
	return o;
}


void computeInfixOrderRec(RecTree * t, vector<RecTree*> & result)
{	
	if (t!=NULL)
	{
		computeInfixOrderRec(t->getRight(), result);
		result.push_back(t);
		computeInfixOrderRec(t->getLeft(), result);
	}
}

vector<RecTree*> computeInfixOrder(RecTree * t)
{
	vector<RecTree *> result;
	computeInfixOrderRec(t, result);
	return result;
}



void computeDepthFirstOrderRec(RecTree * t, vector<RecTree*> & result)
{	
	if (t!=NULL)
	{
		computeDepthFirstOrderRec(t->getRight(), result);
		computeDepthFirstOrderRec(t->getLeft(), result);
		result.push_back(t);
	}
}

vector<RecTree*> computeDepthFirstOrder(RecTree * t)
{
	vector<RecTree *> result;
	computeDepthFirstOrderRec(t, result);
	for(int i=0;i<result.size();i++) 
	{
		RecTree * t = result[i];
		t->setIndex(i);
	}
	return result;
}
 
RecTree * parseNewickRecTree(const string & s, int & i);

void parseNewickChildrenRec(const string & s, int & i, RecTree * father){
	int ibck = i;
	//cerr<<"  children:"<<i<<endl;
	// Assumes at least two children
	father->setLeft(parseNewickRecTree(s, i));
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
		father->setRight(parseNewickRecTree(s, i));
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
					RecTree * newTree = parseNewickRecTree(s, i);
					RecTree * newFather = new RecTree("*");
					RecTree * backupRight = father->getRight();
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


string extractSpecies(string comment)
{
        vector<string> split_comment = split(comment, ':');
	for(std::vector<std::string>::iterator it = split_comment.begin(); it != split_comment.end(); ++it) {
	     string field = *it;
	     if (field.substr(0,2)=="S=")
	     {
		  return field.substr(2,field.size());
	     }
	}
	return "";
}

EventType extractEventType(string comment)
{
        vector<string> split_comment = split(comment, ':');
	for(std::vector<std::string>::iterator it = split_comment.begin(); it != split_comment.end(); ++it)
	{
	     string field = *it;
	     if (field.substr(0,3)=="Ev=")
	     {
		  string event_string = field.substr(3,field.size());
		  if (event_string == "Extant")
		       return Extant;
		  if (event_string == "GLos")
		       return GLos;
		  if (event_string == "GDup") 
		       return GDup;
		  if (event_string == "Spec") 
		       return Spec;
	     }
	}
	return Unknown;
}

string extractND(string comment)
{
        vector<string> split_comment = split(comment, ':');
	for(std::vector<std::string>::iterator it = split_comment.begin(); it != split_comment.end(); ++it) {
	     string field = *it;
	     if (field.substr(0,3)=="ND=")
	     {
		  return field.substr(3,field.size());
	     }
	}
	return "";
}


string parseNewickCommentRec(const string & s, int & i)
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
	reportError("Reached end-of-line while searching for comment-closing bracket",i);
	return comment;
}

RecTree * parseNewickRecTree(const string & s, int & i)
{
	RecTree * result = new RecTree("");
	bool inCaption = true;
	while(i<s.size())
	{
		const char c = s[i];
		switch(c)
		{
			case '(':
			{
				i++;
				parseNewickChildrenRec(s,i,result);
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
			{
				inCaption = false;
				i++;
				string comment = parseNewickCommentRec(s, i);
				result->setComment(comment);
				result->setEvent(extractEventType(comment));
				result->setSpecies(extractSpecies(comment));
				result->setND(extractND(comment));
			}
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

RecTree* parseNewickRecTree(string s)
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
	//cerr << "Read tree: "<< content<<endl;
	RecTree * t = parseNewickRecTree(content,i);
	computeDepthFirstOrder(t);
	return t;
}

bool sameSpecies(RecTree * t1, RecTree * t2)
{
     string s1 = t1->getSpecies();
     string s2 = t2->getSpecies();
     bool b = (s1==s2);
     //cout << "    " << t1->getLabel() << " " << t2->getLabel() << " adjacents? " << b << endl;
     return b;
}

map<string, map<string,string> > loadAdjacencies(string path)
{
  map<string, map<string,string> > adjacencies;
	int i=0;
	string content = "";	
	  
	if (existsFile(path))
	{
	     std::ifstream ifs(path.c_str());
	     content = string(
		  (std::istreambuf_iterator<char>(ifs) ),
		  (std::istreambuf_iterator<char>()    ) );

	     vector<string> lines = split(content, '\n');
	     for(int i=0;i<lines.size();i++)
	     {
		  string line = lines[i];
		  //cout << line << endl;
		  vector<string> pair = split(line,' ');
		  if (pair.size()>=2)
		  {
		       //cout << "Adj: " << pair[0]+pair[1] << endl;
		       // if (adjacencies.count(pair[0])==0)
		       // {
		       // 	    adjacencies[pair[0]] = map<string>();
		       // }
		       adjacencies[pair[0]][pair[1]] = "";
		  }
		  
	     }
	}
	return adjacencies;
}


int RecTree::numNonGDup()
{
     int acc = 0;
     if (getEvent()!=GDup)
     {
	  acc++;
     }
     if (getLeft())
     {
	  acc += getLeft()->numNonGDup();
     }
     if (getRight())
     {
	  acc += getRight()->numNonGDup();
     }
     return acc;
}


bool isAdjacent(RecTree * t1, RecTree * t2,  map<string, map<string,string> > & adjacencies)
{
     string g1 = split(t1->getLabel(),'|')[0];
     string g2 = split(t2->getLabel(),'|')[0];
     bool b = (adjacencies[g1].count(g2)>0)||(adjacencies[g2].count(g1)>0);
     //cout << "    " << t1->getLabel() << " " << t2->getLabel() << " adjacents? " << b << endl;
     return b;
}


