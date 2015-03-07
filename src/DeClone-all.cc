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

#include "DeClone-coopts.hh"

#include <cfloat>
#include <utility>
#include "RecTrees.hh"
#include "utils.hh"

using namespace std;


vector<pair<Tree *, double> > emptyListAll()
{
     return vector<pair<Tree *, double> >();
}

vector<pair<Tree *, double> > atomicListAll(double score, string label = string(""))
{
     vector<pair<Tree *, double> > result;
     Tree* t = new Tree(label);
     result.push_back(pair<Tree *, double>(t,score));
     return result;
}

vector<pair<Tree *, double> > combineListsAll(vector<pair<Tree *, double> > a, vector<pair<Tree *, double> > b)
{
     vector<pair<Tree *, double> > result;
     for(int i =0;i<a.size();i++)
     {
	  pair<Tree*,double> pi = a[i];
	  for(int j =0;j<b.size();j++)
	  {
	       pair<Tree*,double> pj = b[j];
	       Tree* t = new Tree("C", pi.first,pj.first);
	       result.push_back(pair<Tree *, double>(t,pi.second+pj.second));
	  }
     }
     return result;
}

vector<pair<Tree *, double> > filterListsAll(vector<pair<Tree *, double> > a, vector<pair<Tree *, double> > b, string comment)
{
     vector<pair<Tree *, double> > result;

     for(int i =0;i<a.size();i++)
     {
	  pair<Tree*,double> pi = a[i];
	  pi.first->setLabel(comment);
	  result.push_back(pi);
     }
     for(int j =0;j<b.size();j++)
     {
	  pair<Tree*,double> pj = b[j];
	  pj.first->setLabel(comment);
	  result.push_back(pj);
     }
     return result;
}

#define RESULT_TYPE vector<pair<Tree *, double> >
#define INF emptyListAll()
#define AdjGain atomicListAll(adjacency_gain, "AdjGain")
#define AdjBreak atomicListAll(adjacency_break, "AdjBreak")
#define PLUS(a,b) combineListsAll(a,b)
#define MIN(a,b,c,g1,g2,adj) filterListsAll(a,b,c)
#define ZERO atomicListAll(0.0, "Zero")

#define RESCALING_FACTOR(a) ZERO

#define allocateMatrix allocateMatrixAll
#define deleteMatrix deleteMatrixAll
#define computeMatrix computeMatrixAll

ostream& operator<<(ostream & o, const vector<pair<Tree*,double> > & v)
{
     for(int i=0;i<v.size();i++)
     {
	  v[i].first->show();
	  o << "Score: " << v[i].second << endl;
     }

     return o;
}


#include "DeCoDP.cc"

vector<pair<Tree *, double> > getAllScenarios(RecTree *tree1, RecTree *tree2, map<string, map<string,string> > & adjacencies)
{
     // cout << "Arg1.: " << endl;
     // tree1->show(true,1,cout);
     // cout << endl;
     // cout << "Arg2.: " << endl;
     // tree2->show(true,1,cout);
     // cout << endl;
     return MIN(computeMatrixAll(tree1,tree2,true, adjacencies),
		computeMatrixAll(tree1,tree2,false, adjacencies),
		"Root","","",true);
}
