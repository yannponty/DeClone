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

#include <cfloat>
#include <map>
#include <utility>
#include <cstdlib>
#include <sstream>
#include <deque>

#include "RecTrees.hh"
#include "AdjacencyTrees.hh"
#include "OperationsList.hh"

#include "utils.hh"
#include "SVGDriver.hh"

#include "DeClone-all.hh"
#include "DeClone-coopts.hh"
#include "DeClone-count.hh"
#include "DeClone-countcoopts.hh"
#include "DeClone-inside.hh"
#include "DeClone-parsimony.hh"
#include "DeClone-stochastic.hh"
#include "DeClone-outside.hh"

#ifdef USE_POLYTOPE
    #include "DeClone-polytope.hh"
    #include "DeClone-polytope-adj.hh"
    #include "ConvexPolytope.hh"
#endif

#define EXIT_SUCCESS 0


#define TREE_1_OPTION_LONG 	    "--tree1"
#define TREE_1_OPTION_SHORT 	    "-t1"

#define TREE_2_OPTION_LONG    "--tree2"
#define TREE_2_OPTION_SHORT   "-t2"

#define ADJACENCIES_OPTION_LONG    "--adjacencies"
#define ADJACENCIES_OPTION_SHORT   "-a"

#define SCORING_SCHEME_LONG "--score"
#define SCORING_SCHEME_SHORT "-sc"

#define GENE_1_OPTION_LONG    "--gene1"
#define GENE_1_OPTION_SHORT   "-g1"

#define GENE_2_OPTION_LONG    "--gene2"
#define GENE_2_OPTION_SHORT   "-g2"

#define SET_BOLTZMANN_OPTION_SHORT  "-kT"

#define RESCALING_OPTION_SHORT  "-r"
#define RESCALING_OPTION_LONG  "--rescale"

#define OUTPUT_MATRIX_SHORT "-m"
#define OUTPUT_MATRIX_LONG "--matrix"

//// RUN MODES ////
#define PARSIMONY_OPTION_LONG           "--parsimony"
#define PARSIMONY_OPTION_SHORT          "-p"

#define SHOW_COOPTS_OPTION_LONG           "--show-coopts"
#define SHOW_COOPTS_OPTION_SHORT          "-s"

#define COUNT_COOPTS_OPTION_LONG           "--count-coopts"
#define COUNT_COOPTS_OPTION_SHORT          "-c"

#define COUNT_OPTION_LONG           "--count"
#define COUNT_OPTION_SHORT          "-n"

#define PRINT_ALL_OPTION_LONG       "--all"
#define PRINT_ALL_OPTION_SHORT      "-x"

#define STOC_BACKTRACK_OPTION_LONG  "--backtrack"
#define STOC_BACKTRACK_OPTION_SHORT "-b"

#define PARTITION_FUNCTION_OPTION_LONG    "--part-fun"
#define PARTITION_FUNCTION_OPTION_SHORT   "-z"

#define INSIDE_OUTSIDE_OPTION_LONG  "--in-out"
#define INSIDE_OUTSIDE_OPTION_SHORT "-i"

#define POLY_PROP_OPTION_LONG "--polytope"
#define POLY_PROP_OPTION_SHORT "-y"

#define ADJ_POLY_OPTION_LONG "--adjpolytope"
#define ADJ_POLY_OPTION_SHORT "-l"

//#define INTERESTING_ADJ_LONG "--adjy"
//#define INTERESTING_ADJ_SHORT "-j"


//// ADDITIONAL OPTIONS ////

#define DRAW_OPTION_LONG "--draw"
#define DRAW_OPTION_SHORT "-d"


#define HELP_OPTION_LONG            "--help"
#define HELP_OPTION_SHORT           "-h"

#define VERBOSE_OPTION_LONG         "--verbose"
#define VERBOSE_OPTION_SHORT        "-v"

using namespace std;
typedef enum{	PARSIMONY_MODE,
							COUNT_MODE,
							COUNT_COOPTS_MODE,
							SHOW_COOPTS_MODE,
							STOC_BACKTRACK_MODE, 
							PRINT_ALL_MODE, 
							INSIDE_OUTSIDE_MODE, 
							PARTITION_FUNCTION_MODE,
                            POLY_PROP_MODE,
                            ADJ_POLY_MODE
							} RunMode;


void usage(string cmd){
	cerr << "Usage: "<<cmd<<" ["<< TREE_1_OPTION_SHORT<<"|"<< TREE_1_OPTION_LONG<<"] v1 ["<< TREE_2_OPTION_SHORT<<"|"<< TREE_2_OPTION_LONG<<"] v2 ["<< ADJACENCIES_OPTION_SHORT<<"|"<< ADJACENCIES_OPTION_LONG<<"] adj [opts]"<<endl;
	cerr << "Where:"<<endl;
	cerr << "  v1  - (Path to) Gene Tree 1 (Newick format)"<<endl;
	cerr << "  v2  - (Path to) Gene Tree 2 (Newick format)"<<endl;
	cerr << "  adj - Path to a list of adjacent extant genes"<<endl<<endl;
	cerr << "Modes (def.=-p):"<<endl;
	cerr << "  "<<STOC_BACKTRACK_OPTION_SHORT<<","<<STOC_BACKTRACK_OPTION_LONG<<" k   - Stochastic sampling of k adjacency trees"<<endl;
	cerr << "  "<<COUNT_COOPTS_OPTION_SHORT<<","<<COUNT_COOPTS_OPTION_LONG<<"  - Count the number of co-optimal adjacency trees"<<endl;
	cerr << "  "<<HELP_OPTION_SHORT<<","<<HELP_OPTION_LONG<<"          - Displays help and exits"<<endl;
	cerr << "  "<<INSIDE_OUTSIDE_OPTION_SHORT<<","<<INSIDE_OUTSIDE_OPTION_LONG<<"        - Inside-outside mode"<<endl;
  #ifdef USE_POLYTOPE
	  cerr << "  "<<ADJ_POLY_OPTION_SHORT<<","<<ADJ_POLY_OPTION_LONG<<"   - Runs polytope propagation with gain cost, break cost and 2 genes as parameters"<<endl;
    cerr << "                       Requires file with pairs of genes specified by node id in Newick file"<<endl;
  #endif
	cerr << "  "<<COUNT_OPTION_SHORT<<","<<COUNT_OPTION_LONG<<"         - Count the number of valid adjacency trees"<<endl;
	cerr << "  "<<PARSIMONY_OPTION_SHORT<<","<<PARSIMONY_OPTION_LONG<<"     - Maximum parsimony mode (def.)"<<endl;
	cerr << "  "<<SHOW_COOPTS_OPTION_SHORT<<","<<SHOW_COOPTS_OPTION_LONG<<"   - Show all co-optimal adjacency trees"<<endl;
	cerr << "  "<<PRINT_ALL_OPTION_SHORT<<","<<PRINT_ALL_OPTION_LONG<<"           - Exhaustive enumeration of adjacency trees"<<endl;
  #ifdef USE_POLYTOPE
	  cerr << "  "<<POLY_PROP_OPTION_SHORT<<","<<POLY_PROP_OPTION_LONG<<"      - Runs polytope propagation with gain cost and break cost as parameters"<<endl;
  #endif
    cerr << "  " <<PARTITION_FUNCTION_OPTION_SHORT<<","<<PARTITION_FUNCTION_OPTION_LONG<<"      - Computes partition function for instance"<<endl;

    cerr <<endl<< "Parameters:"<<endl;
	cerr << "  "<<DRAW_OPTION_SHORT<<","<<DRAW_OPTION_LONG<<" f        - Draws output to file f (mode-dependent)"<<endl;
	cerr << "  "<<SET_BOLTZMANN_OPTION_SHORT<<" val            - Sets Boltzmann 'constant' (i.e. temperature) to a given value (def.=1.0)"<<endl;
	cerr << "  "<<OUTPUT_MATRIX_SHORT<<","<<OUTPUT_MATRIX_LONG<<"        - Outputs a matrix for the adjacency tree (only for -s and -b modes)"<<endl;
	cerr << "  "<<RESCALING_OPTION_SHORT<<","<<RESCALING_OPTION_LONG<<" val   - Sets rescaling factor (def.=1.0)"<<endl;
	cerr << "  "<<SCORING_SCHEME_SHORT<<","<<SCORING_SCHEME_LONG<<" g b    - Sets costs for adjacency gains (g) and breaks (b) (def.=(1.0,1.0))"<<endl;
	
	cerr << "  "<<VERBOSE_OPTION_SHORT<<","<<VERBOSE_OPTION_LONG<<"       - Verbose mode, provides more (possibly unnecessary) information"<<endl;
}


void ensureNextParamAvail(string opt, string param, int curr, int argc, char *argv[])
{
	if (curr+1>=argc)
	{
		cerr << "Error: Missing argument '"<<param<<"' for option "<<opt<<endl;
		usage(argv[0]);
		exit(EXIT_FAILURE);
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


inline bool convertToInt(std::string const& s, int & d)
{
  std::istringstream i(s);
  int x;
  if (!(i >> x))
    return false;
  d = x;
  return true;
}

void showAdjacencies(std::map<string, map<string,string> > adjacencies)
{
  cerr << "Extant Adjacencies: "<<  endl;
	for(std::map<string, map<string,string> >::iterator iter = adjacencies.begin(); iter != adjacencies.end(); ++iter)
	{
    string a = iter->first;
    for(std::map<string,string>::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); ++iter2)
		{
      string b = iter2->first;
      cerr << "  ("<< a<<","<<b<< ") "<<endl;
    }
  }
}



int main(int argc, char *argv[])
{  
  if (argc>1)
  {
    RecTree * v1 = NULL;
    RecTree * v2 = NULL;
    string gene1 = "";
    string gene2 = "";
    string drawOutput = "";
    map<string, map<string,string> > adjacencies ;
    map<string, map<string,string> > interesting_adjacencies ;
    RunMode mode = PARSIMONY_MODE;
    bool verbose = false;
    bool output_matrix = false;
    int nbSamples;
  	for (int i=1;i<argc;i++)
  	{
  		string opt(argv[i]);
  		if (opt==TREE_1_OPTION_SHORT  || opt==TREE_1_OPTION_LONG)
  		{
  			ensureNextParamAvail(opt, "first tree (v1)", i, argc,argv);
  			i++;
  			v1 = parseNewickRecTree(string(argv[i]));
  		}
  		else if (opt==TREE_2_OPTION_SHORT  || opt==TREE_2_OPTION_LONG)
  		{
  			ensureNextParamAvail(opt, "second tree (v2)", i, argc,argv);
  			i++;
  			v2 = parseNewickRecTree(string(argv[i]));
  		}
  		else if (opt==GENE_1_OPTION_SHORT  || opt==GENE_1_OPTION_LONG)
  		{
//  			ensureNextParamAvail(opt, "first tree (v1)", i, argc,argv);
  			i++;
  			gene1 = string(argv[i]);
  		}
  		else if (opt==GENE_2_OPTION_SHORT  || opt==GENE_2_OPTION_LONG)
  		{
//  			ensureNextParamAvail(opt, "second tree (v2)", i, argc,argv);
  			i++;
  			gene2 = string(argv[i]);
  		}
  		else if (opt==SCORING_SCHEME_LONG || opt==SCORING_SCHEME_SHORT)
  		{
        ensureNextParamAvail(opt, "adj gain cost", i, argc,argv);
  			i++;
        convertToDouble(argv[i],adjacency_gain);
        ensureNextParamAvail(opt, "adj break cost", i, argc,argv);
  			i++;
        convertToDouble(argv[i],adjacency_break);
  		}
  		else if (opt==ADJACENCIES_OPTION_SHORT  || opt==ADJACENCIES_OPTION_LONG)
  		{
  			ensureNextParamAvail(opt, "adjacencies file", i, argc,argv);
  			i++;
  			adjacencies = loadAdjacencies(string(argv[i]));
  			if (verbose)
  			{
          cerr << "Extant Adjacencies: "<<  endl;
          showAdjacencies(adjacencies);
        }
  		}
//      else if (opt==INTERESTING_ADJ_LONG|| opt==INTERESTING_ADJ_SHORT)
//  		{
//  			ensureNextParamAvail(opt, "interesting adjacencies file", i, argc,argv);
//  			i++;
//  			//cout << string(argv[i]) << endl;
//  			interesting_adjacencies = loadAdjacencies(string(argv[i]));
//  			if (verbose)
//  			{
//          cerr << "Interesting Adjacencies: "<<  endl;
//          showAdjacencies(interesting_adjacencies);
//        }
//  		}
  		else if (opt==PARSIMONY_OPTION_SHORT  || opt==PARSIMONY_OPTION_LONG)
  		{
  			mode = PARSIMONY_MODE;
  		}
  		else if (opt==COUNT_OPTION_SHORT  || opt==COUNT_OPTION_LONG)
  		{
  			mode = COUNT_MODE;
  		}
  		else if (opt==COUNT_COOPTS_OPTION_SHORT  || opt==COUNT_COOPTS_OPTION_LONG)
  		{
  			mode = COUNT_COOPTS_MODE;
  		}
        else if (opt==SHOW_COOPTS_OPTION_SHORT  || opt==SHOW_COOPTS_OPTION_LONG)
        {
            mode = SHOW_COOPTS_MODE;
        }
  		else if (opt==PRINT_ALL_OPTION_SHORT  || opt==PRINT_ALL_OPTION_LONG)
  		{
  			mode = PRINT_ALL_MODE;
  		}
  		else if (opt==DRAW_OPTION_SHORT  || opt==DRAW_OPTION_LONG)
  		{
        if (verbose)
        {cerr << "Draw mode set"<<endl;}
  			ensureNextParamAvail(opt, "output file name", i, argc,argv);
        i++;
	  		drawOutput = argv[i];
  		}
  		else if (opt==PARTITION_FUNCTION_OPTION_SHORT  || opt==PARTITION_FUNCTION_OPTION_LONG)
  		{
  			mode = PARTITION_FUNCTION_MODE;
  		}
  		else if (opt==INSIDE_OUTSIDE_OPTION_SHORT  || opt==INSIDE_OUTSIDE_OPTION_LONG)
  		{
  			mode = INSIDE_OUTSIDE_MODE;
  		}
      else if (opt==POLY_PROP_OPTION_SHORT  || opt==POLY_PROP_OPTION_LONG)
  		{
  			mode = POLY_PROP_MODE;
  		}
  		else if (opt==ADJ_POLY_OPTION_SHORT  || opt==ADJ_POLY_OPTION_LONG)
  		{
  			mode = ADJ_POLY_MODE;
  			ensureNextParamAvail(opt, "interesting adjacencies file", i, argc,argv);
  			i++;
  			//cout << string(argv[i]) << endl;
  			interesting_adjacencies = loadAdjacencies(string(argv[i]));
  			if (verbose)
  			{
            cerr << "Interesting Adjacencies: "<<  endl;
            showAdjacencies(interesting_adjacencies);
            }
  		}
      else if (opt==STOC_BACKTRACK_OPTION_SHORT  || opt==STOC_BACKTRACK_OPTION_LONG)
  		{
  			mode = STOC_BACKTRACK_MODE;
  			ensureNextParamAvail(opt, "number of samples", i, argc,argv);
	  		i++;
			convertToInt(string(argv[i]), nbSamples);
  		}
  		else if ((opt==HELP_OPTION_SHORT)  || (opt==HELP_OPTION_LONG))
  		{
  			usage(argv[0]);
  			return EXIT_SUCCESS;
  		}
  		else if (opt==SET_BOLTZMANN_OPTION_SHORT)
		  {
			  ensureNextParamAvail(opt, "Boltzmann constant", i, argc,argv);
			  i++;
			  convertToDouble(string(argv[i]), kT);
		  }
  		else if ((opt==RESCALING_OPTION_SHORT)||(opt==RESCALING_OPTION_LONG))
		  {
			  ensureNextParamAvail(opt, "rescaling factor", i, argc,argv);
			  i++;
			  convertToDouble(string(argv[i]), scalingFactor);
		  }
 		  else if (opt==VERBOSE_OPTION_SHORT  || opt==VERBOSE_OPTION_LONG)
  		{
  			verbose = true;
  		}
      else if (opt==OUTPUT_MATRIX_SHORT  || opt==OUTPUT_MATRIX_LONG)
  		{
  			output_matrix = true;
  		}
  		else
  		{
  			if (v1==NULL)
  			{
  				v1 = parseNewickRecTree(string(argv[i]));				
  			}
  			else if (v2==NULL)
  			{
          //cerr << "Here: "<<string(argv[i])<<endl;
  				v2 = parseNewickRecTree(string(argv[i]));				
  			}
  			else if (adjacencies.size()==0)
  			{
  				adjacencies = loadAdjacencies(string(argv[i]));	
  			}
  		}	 
    }
    if (v1!=NULL && v2!=NULL )
    {
    	if (verbose)
    	{
    		cerr << "Tree 1: ";
    		v1->show(true,2,cerr);
    		cerr << endl;
        cerr.flush();
    	}
    	
    	if (verbose)
    	{
    		cerr << "Tree 2: ";
    		v2->show(true,2,cerr); 
    		cerr << endl;
        cerr.flush();
    	}
    
	cout.precision(10);
    	switch(mode)
	{
        case PARSIMONY_MODE:
        {
      	  cout << computeMaxParsimony(v1,v2,adjacencies) << endl;            
        }
        break;
        case COUNT_COOPTS_MODE:
        {
	     cout << countCooptimalAdjacencyTrees(v1,v2,adjacencies) << endl;
        }
        break;
        case SHOW_COOPTS_MODE:
        {
	     vector<pair<AdjacencyTree *, double> > result = getAllOptimalScenarios(v1,v2,adjacencies);

	     if (output_matrix)
	     {
		  vector<RecTree*> Dfo1 = computeDepthFirstOrder(v1);
		  vector<RecTree*> Dfo2 = computeDepthFirstOrder(v2);

		       for(int i=0; i<result.size(); i++)
		       {
			    result[i].first->showMatrix(Dfo1, Dfo2);
			    cout << endl;
		       }
	     }
	     else
	     {
		  cout << result << endl;
	     }
        }
        break;
        case COUNT_MODE:
        {
      	  cout << countValidAdjacencyTrees(v1,v2,adjacencies) << endl;
        }
        break;
        case PARTITION_FUNCTION_MODE:
        {
      	  cout << computeInside(v1,v2,adjacencies,kT) << endl;
        }
        break;
        case STOC_BACKTRACK_MODE:
        {
	     vector<pair<AdjacencyTree *, double> > result = sample(v1,v2,adjacencies,kT,nbSamples);

	     if (output_matrix)
	     {
		  vector<RecTree*> Dfo1 = computeDepthFirstOrder(v1);
		  vector<RecTree*> Dfo2 = computeDepthFirstOrder(v2);

		       for(int i=0; i<result.size(); i++)
		       {
			    result[i].first->showMatrix(Dfo1, Dfo2);
			    cout << endl;
		       }
	     }
	     else
	     {
		  cout << result << endl;
	     }
	     
        }
        break;
        case PRINT_ALL_MODE:
        {
	     cout << getAllScenarios(v1,v2,adjacencies) << endl;
        }
        break;
        case INSIDE_OUTSIDE_MODE:
        {
	     double ***W = computeOutside(v1,v2,adjacencies);

	     double Z = computeInside(v1,v2,adjacencies,kT);

	     vector<RecTree*> Dfo1 = computeDepthFirstOrder(v1);
	     vector<RecTree*> Dfo2 = computeDepthFirstOrder(v2);
///////////////////////////////////////////////////////////////////
//   
         for(int i=0;i<Dfo1.size();i++)
         {
            if(Dfo1[i]->isLeaf()){
               string firstgene = Dfo1[i]->getGeneName();
               std::map<string,map<string,string> >::iterator it1 = adjacencies.find(firstgene);
               for(std::map<string,string>::iterator it2=it1->second.begin(); it2 != it1->second.end();it2++)
               {
                    for(int j =0;j<Dfo2.size();j++){
                        if(Dfo2[j]->isLeaf() and it2->first == Dfo2[j]->getGeneName()){
                            cout<< "> "<<it1->first<<" "<<it2->first<<" "<<Dfo1[i]->getND()<<" "<<Dfo2[j]->getND()<<endl;
                        }
                    }
               }
            }
         }
         cout<< endl;
//
//////////////////////////////////////////////////////////////////

	     double** probas = new double*[Dfo1.size()];
	     for(int i =0;i<Dfo1.size();i++)
	     {
		  probas[i] = new double[Dfo2.size()];
	     }
	     vector<CaseLabel> c1Lbls =  C1Labels();
	     cout << "\t";
	     for (int j = 0; j < Dfo2.size(); j++)
	     {
		  cout << Dfo2[j]->getND() << " ";
	     }
		       cout << endl;
	     for (int i = 0; i < Dfo1.size(); i++)
	     {
		  cout << Dfo1[i]->getND() << "\t";
	     	  for (int j = 0; j < Dfo2.size(); j++)
	     	  {
		       double c1Weight = 0.0;
	     	       //cout << i << " " << j << " --- ";
	     	       for (int k = 0; k < c1Lbls.size(); k++)
	     	       {
	     		    c1Weight +=  W[Dfo1[i]->getIndex()][Dfo2[j]->getIndex()][c1Lbls[k]];
	     	       }
		       cout << c1Weight / Z << " ";
		       probas[Dfo1[i]->getIndex()][Dfo2[j]->getIndex()] = c1Weight/Z;
	     	  }
			    cout << endl;
	     }
	     if (drawOutput.length()!=0)
	     {  
          if (verbose)
          {
            cerr << "Drawing dot plot to '"<<drawOutput<<"'"<<endl;
          }
          drawProbasSVG(v1,v2,probas,drawOutput);
       }
	     // vector<RecTree*> Dfo1 = computeDepthFirstOrder(v1);
	     // vector<RecTree*> Dfo2 = computeDepthFirstOrder(v2);
	     // double total_weight = 0.;
	     // for (int i=0;i<Dfo1.size();i++) 
	     // {
	     // 	  RecTree* v1_spec = Dfo1[i];
	     // 	  for (int j=0;j<Dfo2.size();j++)
	     // 	  {
	     // 	       RecTree* v2_spec = Dfo2[j];
	     // 	       if ((v2_spec->getEvent() == Spec)&&(v1_spec->getEvent() == Spec))
	     // 	       {
	     // 		    for (int k = 0; k < c1Lbls.size(); k++)
	     // 		    {
	     // 			 int v1_idx = v1_spec->getIndex();
	     // 			 int v2_idx = v2_spec->getIndex();
	     // 			 total_weight += W[v1_idx][v2_idx][c1Lbls[k]];
	     // 		    }
	     // 	       }
	     // 	  }
	     // }
	     // cout << v1->size() + v2->size() << " ";
	     // cout << total_weight << " ";
	     // cout << computeInside(v1,v2,adjacencies,kT) << endl;
        }
        break;
        case POLY_PROP_MODE:
        { 
        #ifdef USE_POLYTOPE
          Polytope p = polycomputeValidAdjacencyTrees(v1, v2, adjacencies);
          cout << "Polygon: "<< p << endl;
          vector<NormalVector> normals = p.normalVectors();
          cout << "Normals (+Signatures): "<<endl<<"{"<<endl;
          for (int i=0;i<normals.size();i++)
          {
            cout << "  ";
            cout << normals[i];
            if (i<normals.size()-1) cout << "," ;
            cout << endl;
          } 
          cout << "}";
        #endif
        #ifndef USE_POLYTOPE
          cerr << "Error : Option ["<<POLY_PROP_OPTION_SHORT<<"|" << POLY_PROP_OPTION_LONG<< "] not-available with current compilation mode."<<endl<<"Please recompile using one of the 'Polytope-aware' compilation targets."<<endl;
    			usage(argv[0]);
          return EXIT_FAILURE;
        #endif
	      }
        break;
        case ADJ_POLY_MODE:
        {
        #ifdef USE_POLYTOPE
        for(std::map<string, map<string,string> >::iterator iter = interesting_adjacencies.begin(); iter != interesting_adjacencies.end(); ++iter)
    		{
          string a = iter->first;
          for(std::map<string,string>::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); ++iter2)
      		{
            string b = iter2->first;
            cout << "Adjacency: "<< a<<","<<b<< endl;
            Polytope p = adjpolycomputeValidAdjacencyTrees(v1, v2,  adjacencies, a, b);
            cout << "Polygon: "<< p << endl;
            vector<NormalVector> normals = p.normalVectors();
            cout << "Normals (+Signatures): "<<endl<<"{"<<endl;
            for (int i=0;i<normals.size();i++)
            {
              cout << "  ";
              cout << normals[i];
              if (i<normals.size()-1) cout << "," ;
              cout << endl;
            } 
            cout << "}" << endl;
            cout << "----------------"<< endl;
          }
        }
        #endif
        #ifndef USE_POLYTOPE
          cerr << "Error : Option ["<<ADJ_POLY_OPTION_SHORT<<"|" << ADJ_POLY_OPTION_LONG<< "] not-available with current compilation mode."<<endl<<"Please recompile using one of the 'Polytope-aware' compilation targets."<<endl;
    			usage(argv[0]);
          return EXIT_FAILURE;
        #endif
  	    }
        break;
    }
    return EXIT_SUCCESS;
  }
  }
  cerr << "Error: Missing arguments"<<endl; 
	usage(argv[0]);
  return EXIT_FAILURE;
}
