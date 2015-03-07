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


vector<pair<AdjacencyTree *, double> > emptyList()
{
     return vector<pair<AdjacencyTree *, double> >();
}

vector<pair<AdjacencyTree *, double> > atomicList(double score, string label = string(""))
{
     vector<pair<AdjacencyTree *, double> > result;
     AdjacencyTree* t = new AdjacencyTree(label,"N/A","N/A",false);
     result.push_back(pair<AdjacencyTree *, double>(t,score));
     return result;
}

vector<pair<AdjacencyTree *, double> > combineLists(vector<pair<AdjacencyTree *, double> > a, vector<pair<AdjacencyTree *, double> > b)
{
     vector<pair<AdjacencyTree *, double> > result;
     for(int i =0;i<a.size();i++)
     {
	  pair<AdjacencyTree*,double> pi = a[i];
	  for(int j =0;j<b.size();j++)
	  {
	       pair<AdjacencyTree*,double> pj = b[j];
	       AdjacencyTree* t = new AdjacencyTree("Combined", pi.first,pj.first,"N/A","N/A",false);
	       result.push_back(pair<AdjacencyTree *, double>(t,pi.second+pj.second));
	  }
     }
     return result;
}

#define MAX_COOPT 1000

vector<pair<AdjacencyTree *, double> > filterLists(vector<pair<AdjacencyTree *, double> > a, vector<pair<AdjacencyTree *, double> > b, string comment, bool adj, string g1, string g2)
{
     vector<pair<AdjacencyTree *, double> > result;
     double minVal =  DBL_MAX;
     for(int i =0;i<a.size();i++)
     {
	  pair<AdjacencyTree*,double> pi = a[i];
	  minVal = min(minVal,pi.second);
     }
     for(int j =0;j<b.size();j++)
     {
	  pair<AdjacencyTree*,double> pj = b[j];
	  minVal = min(minVal,pj.second);
     }

     for(int i =0;i<a.size();i++)
     {
	  pair<AdjacencyTree*,double> pi = a[i];
	  //pi.first->setLabel(comment);
	  if (minVal == pi.second)
	  {
	       if (result.size()< MAX_COOPT)
           result.push_back(pi);
	  }
     }
     for(int j =0;j<b.size();j++)
     {
	  pair<AdjacencyTree*,double> pj = b[j];
	  //pj.first->setLabel(comment);
	  AdjacencyTree * vnode = new AdjacencyTree(comment,pj.first,NULL,g1,g2,adj);
	  if (minVal == pj.second)
	  {
	       pair<AdjacencyTree*,double> velem(vnode,pj.second);
	       if (result.size()< MAX_COOPT)
	       result.push_back(velem);
	  }
	  // if (minVal == pj.second)
	  // {
	  //      result.push_back(pj);
	  // }
     }
     return result;
}

// ostream& operator<<(ostream & o, const vector<pair<RecTree*,double> > & v)
// {
//      for(int i=0;i<v.size();i++)
//      {
// 	  //if (i!=0)
// 	  //{o << ',';}
// 	  v[i].first->show(false,0,o);
// 	  o <<":"<<v[i].second;
// 	  o <<";"<<endl;
//      }
//      return o;
// }


#define RESULT_TYPE vector<pair<AdjacencyTree *, double> >
#define INF emptyList()
#define AdjGain atomicList(adjacency_gain, "AdjGain")
#define AdjBreak atomicList(adjacency_break, "AdjBreak")
#define PLUS(a,b) combineLists(a,b)
#define MIN(a,b,c,adj,g1,g2) filterLists(a,b,c,adj,g1,g2)
#define ZERO atomicList(0.0, "Zero")

#define RESCALING_FACTOR(a) ZERO

#define allocateMatrix allocateMatrixCoopts
#define deleteMatrix deleteMatrixCoopts
#define computeMatrix computeMatrixCoopts



#include "DeCoDP.cc"

vector<pair<AdjacencyTree *, double> > getAllOptimalScenarios(RecTree *tree1, RecTree *tree2, map<string, map<string,string> > & adjacencies)
{
     // cout << "Arg1.: " << endl;
     // tree1->show(true,1,cout);
     // cout << endl;
     // cout << "Arg2.: " << endl;
     // tree2->show(true,1,cout);
     // cout << endl;
     vector<pair<AdjacencyTree *, double> > result = MIN(computeMatrixCoopts(tree1,tree2,true, adjacencies),
							 computeMatrixCoopts(tree1,tree2,false, adjacencies),
							 "Root",false,"N/A","N/A");

     // Weird quick fix to show root correctly. Fix later.
     for (int i = 0; i < result.size(); i++)
     {
	  if (result[i].first->getLabel() == "Root")
	       result[i].first->setLabel(result[i].first->getLeft()->getLabel());
     }

     return result;
}
