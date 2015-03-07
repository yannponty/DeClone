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

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <sstream>

#include "Trees.hh"
#include "DPParsimony.hh"
#include "DPCount.hh"
#include "DPStochasticBacktrack.hh"
#include "DPGenAll.hh"
#include "DPInsideOutside.hh"
#include "SVGDriver.hh"
#include "utils.hh"

#include "math.h"


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


void test1()
{
	string s1 = "(,,(,));";
	string s2 = "(A,B,(C,D));";
	string s3 = "(A,B,(C,D)E)F;";
	string s4 = "(:0.1,:0.2,(:0.3,:0.4):0.5);";
	string s5 = "(:0.1,:0.2,(:0.3,:0.4):0.5):0.0;";
	string s6 = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);";
	string s7 = "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;";
	string s8 = "((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;";
  Tree* t1 = parseNewickTree(s1);
  Tree* t2 = parseNewickTree(s2);
  t1 -> show();
  t2 -> show();
  Tree* t3 = parseNewickTree(s3);
  Tree* t4 = parseNewickTree(s4);
  t3 -> show();
  t4 -> show();
  Tree* t5 = parseNewickTree(s5);
  Tree* t6 = parseNewickTree(s6);
  t5 -> show();
  t6 -> show();
  Tree* t7 = parseNewickTree(s7);
  Tree* t8 = parseNewickTree(s8);
  t7 -> show();
  t8 -> show();
  
	//precompute(gt, st);
	//delete st;
	//delete gt;
}

void test2()
{
	string s1 = "(A,(B,B),(C,C,C),D));";
	string s2 = "(A,B,C,D);";
  Tree* t1 = parseNewickTree(s1);
  Tree* t2 = parseNewickTree(s2);
  //t1 -> show();
  //t2 -> show();
  //computeMaxParsimony(t1, t2);
  //countReconciliations(t1, t2);
  stochasticReconciliations(t1,t2, 100);
	delete t1;
	delete t2;
}



#define GENE_TREE_OPTION_LONG 	    "--gene"
#define GENE_TREE_OPTION_SHORT 	    "-g"
#define SPECIES_TREE_OPTION_LONG    "--species"
#define SPECIES_TREE_OPTION_SHORT   "-s"
#define PARTITION_FUNCTION_OPTION_LONG    "--part-fun"
#define PARTITION_FUNCTION_OPTION_SHORT   "-z"
#define MAX_LIKELIHOOD_OPTION_LONG  "--ML"
#define MAX_LIKELIHOOD_OPTION_SHORT "-m"
#define STOC_BACKTRACK_OPTION_LONG  "--backtrack"
#define STOC_BACKTRACK_OPTION_SHORT "-b"
#define COUNT_OPTION_LONG           "--count"
#define COUNT_OPTION_SHORT          "-c"
#define PRINT_ALL_OPTION_LONG       "--all"
#define PRINT_ALL_OPTION_SHORT      "-a"
#define PRETTY_PRINT_OPTION_LONG    "--pretty"
#define PRETTY_PRINT_OPTION_SHORT   "-p"
#define HELP_OPTION_LONG            "--help"
#define HELP_OPTION_SHORT           "-h"
#define VERBOSE_OPTION_LONG         "--verbose"
#define VERBOSE_OPTION_SHORT        "-v"
#define DUP_ANNOTATIONS_OPTION_LONG             "--dup-matrix"
#define DUP_ANNOTATIONS_OPTION_SHORT             "-j"
#define TREES_OPTION_LONG           "--trees"
#define TREES_OPTION_SHORT          "-t"
#define INSIDE_OUTSIDE_OPTION_LONG  "--in-out"
#define INSIDE_OUTSIDE_OPTION_SHORT "-i"
#define SET_BOLTZMANN_OPTION_SHORT  "-kT"

typedef enum{	PARSIMONY_MODE,
							COUNT_MODE,
							BACKTRACK_MODE, 
							PRINT_ALL_MODE, 
							DUP_ANNOTATIONS_MODE, 
							INSIDE_OUTSIDE_MODE, 
							PARTITION_FUNCTION_MODE
							} runmode;

void usage(string cmd){
	cerr << "Usage: "<<cmd<<" ["<< GENE_TREE_OPTION_SHORT<<"|"<< GENE_TREE_OPTION_LONG<<"] tg ["<< SPECIES_TREE_OPTION_SHORT<<"|"<< SPECIES_TREE_OPTION_LONG<<"] ts [opts]"<<endl;
	cerr << "Where"<<endl;
	cerr << "  tg - Gene Tree (Newick format)"<<endl;
	cerr << "  ts - Species Tree (Newick format)"<<endl;
	cerr << "Options:"<<endl;
	cerr << "  "<<MAX_LIKELIHOOD_OPTION_SHORT<<","<<MAX_LIKELIHOOD_OPTION_LONG<<"          - Maximum parsimony mode (default)"<<endl;
	cerr << "  "<<STOC_BACKTRACK_OPTION_SHORT<<","<<STOC_BACKTRACK_OPTION_LONG<<" k - Stochastic backtrack of k reconciliations (Boltzmann distr.)."<<endl;
	cerr << "  "<<COUNT_OPTION_SHORT<<","<<COUNT_OPTION_LONG<<"       - Count num. reconciliations"<<endl;
	cerr << "  "<<PRINT_ALL_OPTION_SHORT<<","<<PRINT_ALL_OPTION_LONG<<"         - Exhaustive enumeration of reconciliations"<<endl;
	cerr << "  "<<INSIDE_OUTSIDE_OPTION_SHORT<<","<<INSIDE_OUTSIDE_OPTION_LONG<<"      - Computes probability dot-plot for duplications"<<endl;
	cerr << "  "<<PARTITION_FUNCTION_OPTION_SHORT<<","<<PARTITION_FUNCTION_OPTION_LONG<<"      - \"Partition function\" mode"<<endl;
	cerr << "  "<<DUP_ANNOTATIONS_OPTION_SHORT<<","<<DUP_ANNOTATIONS_OPTION_LONG<<"      - Exports gene tree annotated duplications (no species tree needed)"<<endl;

	cerr << "  "<<SET_BOLTZMANN_OPTION_SHORT<<" val  - Sets Boltzmann 'constant' to a given value (def.=1.0)"<<endl;
	cerr << "  "<<VERBOSE_OPTION_SHORT<<","<<VERBOSE_OPTION_LONG<<" - Verbose mode, provides more (possibly unnecessary) information"<<endl;
	cerr << "  "<<PRETTY_PRINT_OPTION_SHORT<<","<<PRETTY_PRINT_OPTION_LONG<<" [pdf|jpeg|png] - Display trees as pictures in files (requires GraphViz)"<<endl;
	cerr << "  "<<TREES_OPTION_SHORT<<","<<TREES_OPTION_LONG<<"      - Output reconciliations as trees (def.= matrices)"<<endl;
	cerr << "  "<<HELP_OPTION_SHORT<<","<<HELP_OPTION_LONG<<"        - Displays help and exits"<<endl;
}


void ensureNextParamAvail(string opt, string param, int curr, int argc, char *argv[])
{
	if (curr+1>=argc)
	{
		cerr << "Error: Missing argument '"<<param<<"' for option "<<opt<<endl;
		usage(argv[0]);
		exit(2);
	}
}

inline bool convertToDouble(std::string const& s, double & d)
{
  std::istringstream i(s);
  double x;
  if (!(i >> x))
    return false;
  d = x;
  return true;
}


template <class T> void freeMatrix(Tree * GeneTree, Tree * SpeciesTree, T** & vals)
{
	int gsize = GeneTree->size();
	for(int i=0;i<gsize;i++)
	{
		delete vals[i];
	}
	delete vals;
	vals = NULL;
}


template <class T> void printMatrix(Tree * GeneTree, Tree * SpeciesTree, T** vals)
{
	vector<Tree*> InfG = computeInfixOrder(GeneTree);
	vector<Tree*> InfS = computeInfixOrder(SpeciesTree);
	cout <<"Columns:";
	for(int j=0;j<InfS.size();j++) 
	{
		Tree * s = InfS[j];
		if(j!=0)
		{ cout <<",";}		
		cout <<s->getLabel();
	}
	cout <<endl;
	cout <<"Rows:";
	for(int j=0;j<InfG.size();j++) 
	{
		Tree * g = InfG[j];
		if(j!=0)
		{ cout <<",";}		
		cout <<g->getLabel();
	}
	cout <<endl;
	cout <<"\t";
	for(int j=0;j<InfS.size();j++) 
	{
		Tree * s = InfS[j];
		cout <<s->getIndex()<<"\t";
	}
	cout <<endl;
	for(int i=0;i<InfG.size();i++) 
	{
		Tree * g = InfG[i];
		cout << g->getIndex() <<"\t";
		for(int j=0;j<InfS.size();j++) 
		{
			Tree * s = InfS[j];
			cout <<vals[g->getIndex()][s->getIndex()]<<"\t";
		}
		cout <<endl;
	}
}

void showReconciliation(EditTree * t, Tree * GeneTree,Tree * SpeciesTree, bool showAsTrees)
{
	if (!showAsTrees)
	{
		NodeType** m = t->eventsAsMatrix(GeneTree,SpeciesTree);
		printMatrix(GeneTree, SpeciesTree, m);
		freeMatrix(GeneTree, SpeciesTree, m);
	}
	else
	{
		t->show(false,2,cout);
		cout << endl;
	}
}


int main(int argc, char *argv[])
{	
	Tree * GeneTree = NULL;
	Tree * SpeciesTree = NULL;
	runmode mode = PARSIMONY_MODE;
	int nbBacktracks = 0;
	bool prettyPrint = false;
	bool verbose = false;
	bool showAsTrees = false;
	string picFormat = "pdf";
	for (int i=1;i<argc;i++)
	{
		string opt(argv[i]);
		if (opt==GENE_TREE_OPTION_SHORT  || opt==GENE_TREE_OPTION_LONG)
		{
			ensureNextParamAvail(opt, "gene tree", i, argc,argv);
			i++;
			GeneTree = parseNewickTree(string(argv[i]));
		}
		else if (opt==SPECIES_TREE_OPTION_SHORT  || opt==SPECIES_TREE_OPTION_LONG)
		{
			ensureNextParamAvail(opt, "species tree", i, argc,argv);
			i++;
			SpeciesTree = parseNewickTree(string(argv[i]));
		}
		else if (opt==MAX_LIKELIHOOD_OPTION_SHORT  || opt==MAX_LIKELIHOOD_OPTION_LONG)
		{
			mode = PARSIMONY_MODE;
		}
		else if (opt==COUNT_OPTION_SHORT  || opt==COUNT_OPTION_LONG)
		{
			mode = COUNT_MODE;
		}
		else if (opt==PARTITION_FUNCTION_OPTION_SHORT  || opt==PARTITION_FUNCTION_OPTION_LONG)
		{
			mode = PARTITION_FUNCTION_MODE;
		}
		else if (opt==PRINT_ALL_OPTION_SHORT  || opt==PRINT_ALL_OPTION_LONG)
		{
			mode = PRINT_ALL_MODE;
		}
		else if (opt==INSIDE_OUTSIDE_OPTION_SHORT  || opt==INSIDE_OUTSIDE_OPTION_LONG)
		{
			mode = INSIDE_OUTSIDE_MODE;
		}
		else if (opt==SET_BOLTZMANN_OPTION_SHORT)
		{
			ensureNextParamAvail(opt, "Boltzmann constant", i, argc,argv);
			i++;
			convertToDouble(string(argv[i]), kT);
		}
		else if (opt==PRETTY_PRINT_OPTION_SHORT  || opt==PRETTY_PRINT_OPTION_LONG)
		{
			ensureNextParamAvail(opt, "file format", i, argc,argv);
			i++;
			prettyPrint = true;
			string format = argv[i];
			std::transform(format.begin(), format.end(), format.begin(), ::tolower);
			picFormat = format;
		}
		else if (opt==HELP_OPTION_SHORT  || opt==HELP_OPTION_LONG)
		{
			usage(argv[0]);
			return EXIT_FAILURE;
		}
		else if (opt==VERBOSE_OPTION_SHORT  || opt==VERBOSE_OPTION_LONG)
		{
			verbose = true;
		}
		else if (opt==TREES_OPTION_SHORT  || opt==TREES_OPTION_LONG)
		{
			showAsTrees = true;
			
		}
		else if (opt==DUP_ANNOTATIONS_OPTION_SHORT  || opt==DUP_ANNOTATIONS_OPTION_LONG)
		{
			mode = DUP_ANNOTATIONS_MODE;
		}
		else if (opt==STOC_BACKTRACK_OPTION_SHORT  || opt==STOC_BACKTRACK_OPTION_LONG)
		{
			ensureNextParamAvail(opt, "#backtracks", i, argc,argv);
			i++;
			mode = BACKTRACK_MODE;
			nbBacktracks = atoi(argv[i]);
		}
		else
		{
			if (GeneTree==NULL)
			{
				GeneTree = parseNewickTree(string(argv[i]));				
			}
			else if (SpeciesTree==NULL)
			{
				SpeciesTree = parseNewickTree(string(argv[i]));				
			}
			
		}	
	}
	if (mode == DUP_ANNOTATIONS_MODE)
	{
		if (GeneTree!=NULL)
		{
			vector<Tree*> dupNodes;
			GeneTree->getDupAnnotations(dupNodes);
			cout << "AnnotatedDup:";
			for (int i=0;i<dupNodes.size();i++)
			{
				if (i!=0)
				{
					cout << ",";
				}
				cout << dupNodes[i]->getIndex();
			}
			return EXIT_SUCCESS;
		}
		else
		{
			cerr << "Error: Missing Gene tree or Species Tree"<<endl;
			return EXIT_FAILURE;
		}
	}
	else 
	{
		if (GeneTree!=NULL && SpeciesTree!=NULL )
		{
			if (verbose)
			{
				cerr << "GeneTree: ";
				GeneTree->show(true,2,cerr);
				cerr << endl;
			}
			if (prettyPrint)
			{	GeneTree->asPicture("GeneTree",picFormat); }
			
			if (verbose)
			{
				cerr << "SpeciesTree: ";
				SpeciesTree->show(true,2,cerr); 
				cerr << endl;
			}
			if (prettyPrint)
			{ SpeciesTree->asPicture("SpeciesTree",picFormat); }
	
			switch(mode)
			{
				case (BACKTRACK_MODE):
					{
						vector<EditTree*> objs = stochasticReconciliations(GeneTree, SpeciesTree, nbBacktracks);
						for (int i=0;i<objs.size();i++)
						{
							EditTree * scenario = objs[i];
							showReconciliation(scenario,GeneTree,SpeciesTree, showAsTrees);
							delete scenario;
						}
					}
					break;
				case (PARSIMONY_MODE):
					{
						EditTree* t = computeMaxParsimony(GeneTree, SpeciesTree);
						if (t)
						{
							showReconciliation(t, GeneTree,SpeciesTree, showAsTrees);
	
							if (prettyPrint)
							{	
								cerr <<endl << "Edit sequence: "<<endl;
								SpeciesTree->show(false,2,cerr); 
								cerr << endl;
								vector<Tree*> v = applyEditTree(SpeciesTree, t);
								SpeciesTree->asPicture("Step-00000-Init",picFormat); 
								for (int i=0;i<v.size();i++)
								{
									v[i]->show(false,2,cerr);
									cout << endl;
									stringbuilder sb;
									sb << "Step-"<< setfill('0') << setw(5)<<(i+1);
									if (i==v.size()-1)
									{ sb<<"-Final"; }
									v[i]->asPicture(sb,picFormat); 
									delete v[i];
								}
							}
						}
						else
						{
							cerr << "Error: No reconciliation possible between these trees."<<endl;
							exit(1);
						}
					}
					break;
				case (COUNT_MODE):
					{
						cout << "#Reconciliations:" << countReconciliations(GeneTree, SpeciesTree)<<endl;
					}
					break;
				case (PARTITION_FUNCTION_MODE):
					{
						cout << "Partition Function:" << getPartitionFunction(GeneTree, SpeciesTree)<<endl;
					}
					break;
				case (INSIDE_OUTSIDE_MODE):
					{
						double** probas = new double*[GeneTree->size()];
						for(int i=0;i<GeneTree->size();i++) 
						{
							probas[i] = new double[SpeciesTree->size()];
						}
						getProbasDupl(GeneTree, SpeciesTree, probas);
						if (prettyPrint)
						{
							drawProbasSVG(GeneTree, SpeciesTree, probas);
						}
						printMatrix(GeneTree, SpeciesTree, probas);
						freeMatrix(GeneTree, SpeciesTree, probas);
					}
					break;
				case (PRINT_ALL_MODE):
					{
						vector<scenario_t> recs = genAllReconciliations(GeneTree, SpeciesTree);
						for(int i=0;i<recs.size();i++) 
						{
							cout << recs[i].score <<'\t'<<recs[i].seq<<endl;
						}
					}
					break;
					
			}
			
			return EXIT_SUCCESS;
		}
		cerr << "Error: Missing Gene tree or Species Tree"<<endl;
		usage(argv[0]);
		return EXIT_FAILURE;
	}
}
