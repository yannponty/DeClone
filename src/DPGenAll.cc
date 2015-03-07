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


// Dyn. Prog.

#include "DPGenAll.hh"
#include <iostream>
#include <string>


scenario_t plusScenar(scenario_t a, scenario_t b)
{
	scenario_t result(a.score+b.score,string("(")+a.seq+string(",")+b.seq+string(")"));
	return result;
}

vector<scenario_t> plusVecScenar(vector<scenario_t> a, vector<scenario_t> b)
{
	vector<scenario_t> result;
	for (int i=0;i<a.size();i++)
	{
		for (int j=0;j<b.size();j++)
		{
			scenario_t s = plusScenar(a[i],b[j]);
			if (s.score<DBL_MAX)
			{			
				result.push_back(s);
			}
		}
	}
	return result;
}

vector<scenario_t> minVecScenar(vector<scenario_t> a, vector<scenario_t> b)
{
	vector<scenario_t> result;
	for (int i=0;i<a.size();i++)
	{
		scenario_t s = a[i];
		if (s.score<DBL_MAX)
		{			
			result.push_back(s);
		}
	}
	for (int i=0;i<b.size();i++)
	{
		scenario_t s = b[i];
		if (s.score<DBL_MAX)
		{			
			result.push_back(s);
		}
	}
	return result;
}

vector<scenario_t> singleScenario(scenario_t a)
{
	vector<scenario_t> result;
	result.push_back(a);
	return result;
}

ostream & operator<<(ostream & o, scenario_t s)
{
	o<<"("<<s.seq<<","<<s.score<<")";
	return o;
}


ostream & operator<<(ostream & o, vector<scenario_t> v)
{
	o<<"[";
	for(int i=0;i<v.size();i++)
	{
		o<<v[i]<<",";
	}
	o<<"]";
	return o;
}

#define PLUS(a,b) plusVecScenar(a,b)
#define MIN(a,b) minVecScenar(a,b)
#define ZERO singleScenario(scenario_t(0.,string("M")))

#define TYPEDATA vector<scenario_t>
#define INFTY singleScenario(scenario_t(DBL_MAX,string("")))
#define DUP_COST singleScenario(scenario_t(1.,string("dup")))
#define LOSS_COST singleScenario(scenario_t(1.,string("loss")))

#define allocateMatricesTASK allocateMatricesGenAll
#define fillMatricesTASK fillMatricesGenAll
#define deleteMatricesTASK deleteMatricesGenAll

using namespace std;

#include "DPRaw.cc"

vector<scenario_t> genAllReconciliations(Tree * GeneTree, Tree * SpeciesTree)
{
	TYPEDATA** R = allocateMatricesTASK(GeneTree, SpeciesTree);
  fillMatricesTASK(GeneTree, SpeciesTree, R);
  // Here: HUGE memory leak! 
	return R[GeneTree->getIndex()][SpeciesTree->getIndex()];
}
