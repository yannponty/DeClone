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

#ifndef REC_TREES_HH
#define REC_TREES_HH

using namespace std;

typedef enum {
	GDup, 
	GLos, 
	Spec, 
	Extant,
        Unknown} EventType;

ostream& operator<<(ostream &o , EventType t);

string prettyOperationType(EventType nt);

class RecTree: public Tree{
	private:
		EventType type;
                string species;

	public:
		RecTree(string lbl);
		RecTree(string lbl, EventType t);
		RecTree(string lbl, RecTree * left, RecTree * right);
                RecTree(string lbl, RecTree * left, RecTree * right, EventType nt, string spec);
		
		RecTree * getLeft();
		RecTree * getRight();
		RecTree * getParent();
		void setLeft(RecTree * t);
		void setRight(RecTree * t);
		void setParent(RecTree * t);
        string getGeneName();

        string getSpecies();
                void setSpecies(string spec);

		EventType getEvent();
		void setEvent(EventType t);

                int numNonGDup();
		
};

vector<RecTree*> computeDepthFirstOrder(RecTree * t);
int findDepthFirstOrder(RecTree * t,string GeneName);

vector<RecTree*> computeInfixOrder(RecTree * t);

RecTree* parseNewickRecTree(string s);

bool sameSpecies(RecTree * t1, RecTree * t2);


//////////// Adjacencies //////////////

map<string, map<string,string> > loadAdjacencies(string path);

bool isAdjacent(RecTree * t1, RecTree * t2,  map<string, map<string,string> > & adjacencies);

#define IsAdj(a,b,adjacencies) (isAdjacent(a,b,adjacencies)? ZERO : INF)
#define IsntAdj(a,b,adjacencies) (isAdjacent(a,b,adjacencies)? INF : ZERO)



#endif
