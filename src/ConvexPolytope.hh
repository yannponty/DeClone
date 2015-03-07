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

#ifndef YPOLYTOPE_HH
#define YPOLYTOPE_HH

#include <iostream>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <valarray>
#include <set>
#include <chrono>

#include <cstdlib>
#include <cstdio>
#include <cassert>

using size_type = std::size_t;
using G = double;                                  // Type for individual coordinate
using point_type = std::valarray< G >;            // Type for point
using points_type = std::vector< point_type >;  // Type for list of points

class NormalVector{
  public:
    point_type support;
    point_type vec;
    points_type sigs;

    NormalVector(size_t dimension);
    NormalVector(std::vector<G> vals);
//    NormalVector(point_type v);
//    NormalVector(point_type v,point_type supp);
    NormalVector(const NormalVector & v);
    void normalize();
    friend std::ostream & operator<<(std::ostream & o, const NormalVector & p);
    
    void setSignatures(const points_type p);
};

class Polytope{
  private:
    size_type dimension;
    points_type points;                      // Actual list of points
    NormalVector computeNormal(std::vector<size_t> vertices_);
   
   public: 
    Polytope();
    Polytope(size_type dim);
    Polytope(const size_t dim, const G* point);
    Polytope(size_type dim,points_type p);
    Polytope(const Polytope & p);
    
    std::vector<std::vector<size_t> > convexHullFacets();
    std::vector<NormalVector> normalVectors();
    
    Polytope convexHull();

    Polytope unionSum(Polytope p2);

    Polytope minkovskiSum(Polytope p2);

    Polytope convexSum(Polytope p2);

    friend std::ostream & operator<<(std::ostream & o, const Polytope & p);
};

std::ostream & operator<<(std::ostream & o, const point_type & p);
std::ostream & operator<<(std::ostream & o, const points_type & p);

G l2Norm(std::valarray<G> vector);

#endif
