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

#include "DeClone-polytope-adj.hh"

#include <cfloat>
#include <map>
#include <utility>
#include "RecTrees.hh"
#include "utils.hh"
#include "ConvexPolytope.hh"
#include <fstream>
#include <iostream>
#include <cstdlib>

#define RESULT_TYPE Polytope
#define DIMENSION 3

#define allocateMatrix allocateMatrixAdjPolytope
#define deleteMatrix deleteMatrixAdjPolytope
#define computeMatrix computeMatrixAdjPolytope

#define AdjGain Polytope(DIMENSION,adjgain)
#define AdjBreak Polytope(DIMENSION,adjbreak)
#define ZERO Polytope(DIMENSION,origin)
#define INF Polytope(DIMENSION)
#define RESCALING_FACTOR(a) ZERO 

#define PLUS(a,b) minkowski_addition_adj(a,b) 
#define MIN(a,b,c,adj,g1,g2) convex_hull_adj(a,b,adj,g1,g2)


std::pair<std::string,std::string> adj_pair;
const double adjgain[] = {1.,0.,0.};
const double adjbreak[] = {0.,1.,0.};
const double origin[] = {0.,0.,0.};


Polytope minkowski_addition_adj(Polytope p1, Polytope p2)
{
  return p1.minkovskiSum(p2);
}

Polytope adjacency_parameter_shift(Polytope p, const bool adj, const std::string& g1, const std::string& g2){
  
    //if (((g1 == adj_pair.first)) ){
    //  cout << g1 << " ("<< g1.size() << ") "<<g2<< " ("<< g2.size() << ") "<<" "<<adj_pair.first<<" ("<<adj_pair.first.size()<<") "<<adj_pair.second<<" ("<<adj_pair.second.size() << ") -> YES "  << endl;
    //}
    
    if (adj && ((g1 == adj_pair.first && g2 == adj_pair.second)) ){
       const double adjacency_param[] = {0.,0.,1.};
       return p.minkovskiSum(Polytope(DIMENSION,adjacency_param));     
    }
    else{
        return p;
    }
}

Polytope convex_hull_adj(Polytope p1, \
                            Polytope p2, \
                            const bool adj, \
                            const std::string& g1, \
                            const std::string& g2){
    //cout << g1 << " "<<g2<<" "<<adj_pair.first<<" "<<adj_pair.second << " -> " << adjacency_parameter_shift(p2,adj,g1,g2) << endl;
    //cout << "Convx: "<< p1 << " "<<p2<<" ["<< adjacency_parameter_shift(p2,adj,g1,g2)<<"]=> "<< p1.convexSum(adjacency_parameter_shift(p2,adj,g1,g2)) << endl;
    Polytope result = p1.convexSum(adjacency_parameter_shift(p2,adj,g1,g2));
    return result;
}


#include "DeCoDP.cc"


Polytope adjpolycomputeValidAdjacencyTrees(RecTree *tree1, RecTree *tree2,  map<string, map<string,string> > & adjacencies, const std::string& gene1, const std::string& gene2)
{
    adj_pair.first = gene1;
    adj_pair.second = gene2;

    return MIN(MIN(INF,computeMatrixAdjPolytope(tree1,tree2,true, adjacencies),"Root",true,"N/A","N/A"),
               computeMatrixAdjPolytope(tree1,tree2,false, adjacencies),
  	           "Root",false,"N/A","N/A");

}




/*int
main(int argc, char * argv[])
{
 
    size_t dim = 3;
    points_type points_(5);                      // Actual list of points
    
    // Assigning points
    for (size_type i=0;i<points_.size();i++)
    {
      point_type & point_ = points_[i];
      point_.resize(dim);
    }
    {
    point_type & point_ = points_[0];
    point_[0] = 0.0;
    point_[1] = 0.0;
    point_[2] = 0.0;
    }
    {
    point_type & point_ = points_[1];
    point_[0] = 0.0;
    point_[1] = 0.0;
    point_[2] = 1.0;
    }
    {
    point_type & point_ = points_[2];
    point_[0] = 0.0;
    point_[1] = 1.0;
    point_[2] = 0.0;
    }
    {
    point_type & point_ = points_[3];
    point_[0] = 1.0;
    point_[1] = 0.0;
    point_[2] = 0.0;
    }
    {
    point_type & point_ = points_[4];
    point_[0] = 0.1;
    point_[1] = 0.1;
    point_[2] = 0.1;
    }
    
    Polytope p1(dim,points_);
    std::cout << "Polytope p1:"<< std::endl;
    std::cout << p1<< std::endl;
    std::cout << "Convex Hull (p1):"<< std::endl;
    std::cout << p1.convexHull()<<std::endl;

    const double np[] = {1.3,1.6,1.7};
    Polytope p2(3,np);

    std::cout << "Polytope p2:"<< std::endl;
    std::cout << p2<< std::endl;
    
    std::cout << "Minkovski sum (p1+p2):"<< std::endl;
    std::cout << p1.minkovskiSum(p2)<<std::endl;
    std::cout << "Minkovski sum (p2+p1):"<< std::endl;
    std::cout << p2.minkovskiSum(p1)<<std::endl;

    std::cout << "p3 := Minkovski sum ((p1+p2)+(p1+p2)):"<< std::endl;
    Polytope p3 = (p1.minkovskiSum(p2)).minkovskiSum(p1.minkovskiSum(p2));
    std::cout << p3 <<std::endl;
    
    std::cout << "Convex Hull (p3):"<< std::endl;
    std::cout << p3.convexHull()<<std::endl;

    std::cout << "Convex Hull (p3):"<< std::endl;
    std::cout << p3.convexHull()<<std::endl;
    
    // Debugging DP
    
    adj_pair.first = "";
    adj_pair.second = "";

    points_.resize(3);
    // Assigning points
    for (size_type i=0;i<points_.size();i++)
    {
      point_type & point_ = points_[i];
      point_.resize(dim);
    }
    {
    point_type & point_ = points_[0];
    point_[0] = 0.0;
    point_[1] = 0.0;
    point_[2] = 0.0;
    }
    {
    point_type & point_ = points_[1];
    point_[0] = 0.0;
    point_[1] = 1.0;
    point_[2] = 0.0;
    }
    {
    point_type & point_ = points_[2];
    point_[0] = 2.0;
    point_[1] = 1.0;
    point_[2] = 0.0;
    }
    Polytope p4(DIMENSION,points_);

    points_.resize(1);
    // Assigning points
    for (size_type i=0;i<points_.size();i++)
    {
      point_type & point_ = points_[i];
      point_.resize(dim);
    }
    {
    point_type & point_ = points_[0];
    point_[0] = 2.0;
    point_[1] = 2.0;
    point_[2] = 0.0;
    }
    Polytope p5(DIMENSION,points_);

    cout << p4 << " " << p5 << "=>" << p4.convexSum(adjacency_parameter_shift(p5,true,"x","")) << endl;
//    cout << p4 << " " << p5 <<"("<<adjacency_parameter_shift(p5,true,"","")<<")" << "=>" << p4.convexSum(adjacency_parameter_shift(p5,true,"","")) >> endl;
    {
    point_type point1(3);
    point1[0] = 2.0;
    point1[1] = 2.0;
    point1[2] = 0.0;
    point_type point2(3);
    point2[0] = 2.0;
    point2[1] = 1.0;
    point2[2] = 0.0;
    }


    //Batch mode
    RecTree * v1 = NULL;
    RecTree * v2 = NULL;
    map<string, map<string,string> > adjacencies ;
    map<string, map<string,string> > interesting_adjacencies ;

    if (argc>4)
    {
  		v1 = parseNewickRecTree(string(argv[1]));				
  		v2 = parseNewickRecTree(string(argv[2]));				
  		adjacencies = loadAdjacencies(string(argv[3]));	
  		interesting_adjacencies = loadAdjacencies(string(argv[4]));	
      for(std::map<string, map<string,string> >::iterator iter = interesting_adjacencies.begin(); iter != interesting_adjacencies.end(); ++iter)
  		{
        string a = iter->first;
        for(std::map<string,string>::iterator iter2 = iter->second.begin(); iter2 != iter->second.end(); ++iter2)
    		{
          string b = iter2->first;
          cout << "Adjacency: "<< a<<","<<b<< endl;
          cout << adjpolycomputeValidAdjacencyTrees(v1, v2,  adjacencies, a, b) << endl;
          cout << "----------------"<< endl;
        }
      }
    }
    else
    {
  		cerr << "Usage: "<<argv[0]<<" tree1 tree2 adjacencies interesting_adjancencies"<<endl; 
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}*/

