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

#include "DeClone-stochastic.hh"

#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <utility>
#include "RecTrees.hh"
#include "utils.hh"

using namespace std;



ostream& operator<<(ostream & o, const vector<AdjacencyTree*> & v)
{
     for(int i=0;i<v.size();i++)
     {
	  o<< v[i];
	  o <<";"<<endl;
     }
     return o;
}

int NUM_TREES;

pair< vector<pair<AdjacencyTree *, double> > , double> emptyListSample()
{
     return pair< vector<pair<AdjacencyTree *, double> > , double>( vector<pair<AdjacencyTree *, double> >() , 0.);
}

pair< vector<pair<AdjacencyTree *, double> > , double> atomicListSample(double boltzmannScore, double parsimonyScore, string label = string(""))
{
     pair< vector<pair<AdjacencyTree *, double> > , double> result;
     result.second = boltzmannScore;
     for(int i = 0; i <  NUM_TREES; i++)
     {
	  AdjacencyTree* t = new AdjacencyTree(label,"N/A","N/A",false);
	  result.first.push_back(pair<AdjacencyTree *, double>(t, parsimonyScore));
     }
     return result;
}

//
pair< vector<pair<AdjacencyTree *, double> > , double> combineListsSample(pair< vector<pair<AdjacencyTree *, double> > , double> a, pair< vector<pair<AdjacencyTree *, double> > , double> b)
{
     pair< vector<pair<AdjacencyTree *, double> > , double> result;
     result.second = a.second*b.second;

     for(int i =0;i<min(a.first.size(),b.first.size());i++)
     {
	  AdjacencyTree* pi = a.first[i].first;
	  AdjacencyTree* pj = b.first[i].first;
	  AdjacencyTree* t = new AdjacencyTree("Combined", pi,pj,"N/A","N/A",false);
	  result.first.push_back(pair<AdjacencyTree *, double>(t, a.first[i].second + b.first[i].second));
     }
     return result;
}

pair< vector<pair<AdjacencyTree *, double> > , double> filterListsSample(pair< vector<pair<AdjacencyTree *, double> > , double> a, pair< vector<pair<AdjacencyTree *, double> > , double> b, string comment, bool adj, string g1, string g2)
{
     pair< vector<pair<AdjacencyTree *, double> > , double> result;
     result.second = a.second+b.second;
     for(int i =0;i< NUM_TREES;i++)
     {
	  if (result.second!=0.)
	  {
	       double r  = result.second*((double)rand()/(1.+(double)RAND_MAX));
	       AdjacencyTree* p = NULL;
	       double score = DBL_MAX;
	       if (r-a.second < 0)
	       {
		    //p = a.first[i];
		    p = new AdjacencyTree(comment,a.first[i].first,NULL,g1,g2,adj);
		    score = a.first[i].second;
	       }
	       else
	       {
		    p = new AdjacencyTree(comment,b.first[i].first,NULL,g1,g2,adj);
		    score = b.first[i].second;
	       }
	       result.first.push_back(pair<AdjacencyTree *, double>(p, score));
	  }
     }
     return result;
}


#define RESULT_TYPE pair< vector<pair<AdjacencyTree *, double> > , double>
#define INF emptyListSample()
#define AdjGain atomicListSample( exp(-adjacency_gain/kT), 1.0, "AdjGain")
#define AdjBreak atomicListSample( exp(-adjacency_break/kT), 1.0, "AdjBreak")
#define PLUS(a,b) combineListsSample(a,b)
#define MIN(a,b,c,adj,g1,g2) filterListsSample(a,b,c,adj,g1,g2)
#define ZERO atomicListSample(1.0, 0.0, "Zero")

#define RESCALING_FACTOR(a) atomicListSample( pow(scalingFactor,a), 0.0, "Rescale")

#define allocateMatrix allocateMatrixStochastic
#define deleteMatrix deleteMatrixStochastic
#define computeMatrix computeMatrixStochastic

#include "DeCoDP.cc"

vector<pair<AdjacencyTree *, double> > sample(RecTree *tree1, RecTree *tree2, map<string, map<string,string> > & adjacencies, double kTval, int nbSamples)
{
     NUM_TREES = nbSamples;
     kT = kTval;
     srand((unsigned int)time(NULL));

     pair< vector<pair<AdjacencyTree *, double> > , double> result = MIN(computeMatrixStochastic(tree1,tree2,true, adjacencies),
									 computeMatrixStochastic(tree1,tree2,false, adjacencies),
									 "Root",false,"N/A","N/A");

     return result.first;
}
