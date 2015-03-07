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

#include "ConvexPolytope.hh"
#include "quickhull.hpp"


#include <algorithm>
#include <vector>
#include <iterator>     // std::distance


using namespace std;

std::ostream & operator<<(std::ostream & o, const point_type & p)
{
  o<<"{";
  bool first = true;
  for (G const coord : p) {
    if(!first)
    { o <<","; }
    o << coord;
    first = false;
  }
  o<<"}";
  return o;
}

std::ostream & operator<<(std::ostream & o, const points_type & p)
{
  o <<"{";
  bool first = true;
  for (point_type const pt : p) {
    if(!first)
    { o <<","; }
    o << pt;
    first = false;
  }
  o<< "}"; 
  return o;
}


std::ostream & operator<<(std::ostream & o, const Polytope & p)
{
//  o<<"Polygon["<<p.points<<"]";
  o<<""<<p.points<<"";
  return o;
}

std::ostream & operator<<(std::ostream & o, const NormalVector & v)
{
  // Uncomment below to place the normals on the polytope
  // o<<"{"<<v.support<<","<<(v.vec+v.support)<<"}";
  o<<v.vec << " -> "<< v.sigs;
  return o;
}


point_type operator+(const point_type & x, const point_type & y)
{
  point_type res(x.size());
  for(size_t i=0;i<x.size();i++)
  {
    res[i]=x[i]+y[i];
  }
  return res;
}

point_type operator/(const point_type & x, G k)
{
  point_type res(x.size());
  for(size_t i=0;i<x.size();i++)
  {
    res[i]=x[i]/k;
  }
  return res;
}


NormalVector::NormalVector(size_t dimension)
{
  vec.resize(dimension);
  support.resize(dimension);
  sigs.resize(0);
  for(size_t i=0;i<dimension;i++)
  {
    vec[i]=0.;
    support[i]=0.;
  }
}

NormalVector::NormalVector(vector<G> vals)
{
  vec.resize(vals.size());
  support.resize(vals.size());
  sigs.resize(0);
  for(size_t i=0;i<vals.size();i++)
  {
    vec[i]=vals[i];
    support[i]=0.;
  }
  normalize();
}


/*NormalVector::NormalVector(point_type v)
{
  vec.resize(v.size());
  support.resize(v.size());
  for(size_t i=0;i<v.size();i++)
  {
    vec[i]=v[i];
    support[i]=0.;
  }
}

NormalVector::NormalVector(point_type v,point_type supp)
{
  vec.resize(v.size());
  support.resize(supp.size());
  for(size_t i=0;i<vec.size();i++)
  {
    vec[i]=v[i];
    support[i]=supp[i];
  }
}*/

NormalVector::NormalVector(const NormalVector & v)
{

  vec.resize(v.vec.size());
  support.resize(v.support.size());
  sigs.resize(v.sigs.size());
  for(size_t i=0;i<vec.size();i++)
  {
    vec[i]=v.vec[i];
    support[i]=v.support[i];
  }
  for(size_t i=0;i<sigs.size();i++)
  {
    sigs[i].resize(v.sigs[i].size());
    for(size_t j=0;j<v.sigs[i].size();j++)
    { 
      sigs[i][j]=v.sigs[i][j];
    }
  }
}


void NormalVector::setSignatures(const points_type p)
{
  sigs.resize(p.size()); 
  for(size_t i=0;i<sigs.size();i++)
  {
    sigs[i].resize(p[i].size());
    for(size_t j=0;j<sigs[i].size();j++)
    { 
      sigs[i][j]=p[i][j];
    }
  }
}


bool pointsCompRec (const point_type & p1, const point_type & p2, int d=0)
{
  if (p1.size()!=p2.size())
  {
    return p1.size()<p2.size();
  }
  if (d<p1.size())
  {
    if (p1[d]<p2[d])
    {
      return true;
    }
    else if (p1[d]>p2[d])
    {
      return false;
    }
    else

    {
      return pointsCompRec(p1,p2,d+1);
    }
  }
  return false;
}


struct pointsComp {
  bool operator() (const point_type& lhs, const point_type& rhs) const
  {return pointsCompRec(lhs,rhs);}
};



Polytope::Polytope()
{
  dimension = -1;
}


Polytope::Polytope(size_type dim)
{
  dimension = dim;
  points.resize(0);
}

Polytope::Polytope(const size_t dim, const G* point)
{
  dimension = dim;
  points.resize(1);
  points[0].resize(dimension);
  for (int i=0;i<dim;i++)
  {
    points[0][i] = point[i];
  } 
}

Polytope::Polytope(const Polytope & p)
{
  dimension = p.dimension;
  points.resize(p.points.size());
  for (int i=0;i<p.points.size();i++)
  {
    points[i].resize(dimension);
    for (int j=0;j<dimension;j++)
    {
      points[i][j] = p.points[i][j];
    } 
  } 
}



Polytope::Polytope(size_type dim,points_type p)
{
  dimension = dim;
  std::set<point_type,pointsComp> tmppoints;
    
  for (point_type const pt : p) {
    tmppoints.insert(pt);
  }

  points.resize(tmppoints.size());
  int i=0;
  for (point_type const p : tmppoints) {
    points[i] = p;
    i++;
  }
}


bool pointsEqual (const point_type & p1, const point_type & p2,  int d=0)
{
  if (p1.size()!=p2.size())
  {
    return false;
  }

  if (d<p1.size())
  {
    if (p1[d]!=p2[d])
    {
      return false;
    }
    else
    {
      return pointsEqual(p1,p2,d+1);
    }
  }
  return true;
}


points_type dropLastDim(const points_type & p)
{
  //std::cerr << "DropDim:"<< p;
  
  points_type res(p.size());  
  for (int i=0;i<p.size();i++) {
    res[i].resize(p[i].size()-1);
    for (int j=0;j<res[i].size();j++) {
      res[i][j] = p[i][j];
    }
  }
  //std::cerr << "=>"<< res<<std::endl;
  return res;
}

bool lastDimHomogenous(const points_type & p)
{
  G val;
  bool firstval = true;
  for (int i=0;i<p.size();i++) {
    if (firstval){
      val = p[i][p[i].size()-1];
      firstval = true;
    }
    else
    {
      if (val != p[i][p[i].size()-1])
      return false;
    }
  }
  return true;
}


Polytope Polytope::convexHull()
{
  if (points.size()<= dimension )
  {    
    return Polytope(dimension,points);
  }
  else
  {
    std::vector<std::vector<size_t> > facets_ = convexHullFacets();
    std::set<size_t> indices;
    for (size_type i = 0; i < facets_.size(); ++i) 
    {
      auto const & vertices_ = facets_[i];          
      for (size_type const vertex_ : vertices_) 
      {
        indices.insert(vertex_);
      }
    }

    if (indices.size()==0)
    {
      return Polytope(dimension,points);
    }    
    points_type respoints(indices.size());
    size_type i=0;
    for (size_type const index : indices) 
    {
      respoints[i] = points[index];
      i++;
    }
    return Polytope(dimension,respoints);
  }
}

void printMatrix(std::valarray<std::valarray<G> > & a ){
    for (int k = 0; k < a.size(); ++k){
        for (int l = 0; l < a.size() + 1; ++l){
            cout << a[k][l] <<" ";
        }
        cout << endl;
    }
    cout << endl;   
}

void forwardSubstitution(std::valarray<std::valarray<G> > & a) {
    int i, j, k, max;
    float t;
    for (i = 0; i < a.size(); ++i) {
        max = i;
        // Find largest at coordinate i among system of equations.         
        for (j = i + 1; j < a.size(); ++j) {
            if (abs(a[j][i]) > abs(a[max][i])){
                max = j;
            }
        }
        // Swap rows         
        for (j = 0; j < a.size() + 1; ++j) {
            t = a[max][j];
            a[max][j] = a[i][j];
            a[i][j] = t;
        }
        for (j = a.size(); j >= i; --j){
            for (k = i + 1; k < a.size(); ++k){
                a[k][j] -= a[k][i]/a[i][i] * a[i][j];
            }
        }
    }
}
 
std::vector<G> reverseElimination(std::valarray<std::valarray<G> > & a) {
    int i, j;
    int temp_dimension = a[0].size();
    std::vector<G> result;
    // Subtract last row from rest so that (i,j)^th terms for j>i are 0.
    for (i = a.size() - 1; i >= 0; --i) {
        // Check if (i,i)^th term is 0.
        if(a[i][i] != 0.){
            // If not, use the row to reduce the rows occurring before. 
//            cout<<"Here\n";
//            printMatrix(a);
            a[i][a.size()] = a[i][a.size()] / a[i][i];
            a[i][i] = 1;
            for (j = i - 1; j >= 0; --j) {
                a[j][a.size()] -= a[j][i] * a[i][a.size()];
                a[j][i] = 0.;
            }
//            cout<<"Test elimination\n";
//            printMatrix(a);
        }
        else{
//            cout<<"Here\n";
//            printMatrix(a);
            // Otherwise, i^th entry in normal vector will be 0.
                std::valarray<std::valarray<G> > system(a.size()-1);
                for (size_t k=0; k< a.size()-1; k++){
                    system[k].resize(a[k].size()-1);
                    for (size_t l = 0; l<a[k].size()-1; l++){
///                       if (l<a.size()-1){
                          system[k][l] = a[k][l];
///                       }
///                       else{
///                           system[k][l] = -a[k][l];
///                       }
                   }
                }
                // Run reverseElimination on the submatrix on k<i.
                result = reverseElimination(system);
//            cout<<"Test elimination\n";
//            printMatrix(a);
                result.push_back(0.);
//                cout<< result<<endl;
                break;

        }
    }
    if (result.size() > 0){
        return result;
    }
    else{
        for (size_t i=0;i<a.size();i++)
        {
          G comp = a[i][temp_dimension-1];
          result.push_back(comp);
        } 
        result.push_back(-1.);
//        cout<< result<< endl;
        return result;
    }
}

std::vector<G> crossProduct(std::valarray<std::valarray<G> > & a){
    std::vector<G> result;
    if(a.size() > 2){
        cout << "System is too big"<<endl;
        // Raise exfception here.
        result.push_back(0.);
        result.push_back(0.);
        result.push_back(0.);
        return result; 
    }
    
    G temp_coord = a[0][1]*a[1][2]-a[1][1]*a[0][2];
    result.push_back(temp_coord);

    temp_coord = a[0][2]*a[1][0]-a[0][0]*a[1][2];
    result.push_back(temp_coord);

    temp_coord = a[0][0]*a[1][1]-a[1][0]*a[0][1];
    result.push_back(temp_coord);

    return result; 

}


#define EPS 0.00001

int testDirectionality(const point_type & planep, const vector<G> & normalvec, const point_type & anotherp)
{
  G sum = 0.;
//  vector<G> difference;
  for (size_t i=0;i<planep.size();i++)
  {
    sum += (anotherp[i]-planep[i])*normalvec[i];
//    difference.push_back(anotherp[i]-planep[i]);
  }
  if (abs(sum)<EPS)
  {return 0;}
  else if (sum>0)
  { return 1; }
  else if (sum<0)
  { return -1; }
}

bool badDirectionality(const point_type & planep, const vector<G> & normalvec, const point_type & anotherp, const int d)
{
  G sum = 0.;
  G direction = 0.;
  if (d>0){
      direction = 1.;
  }
  else{
      direction = -1.;
  }
//  vector<G> difference;
  for (size_t i=0;i<planep.size();i++)
  {
    sum += (anotherp[i]-planep[i])*normalvec[i];
//    difference.push_back(anotherp[i]-planep[i]);
  }
  if (abs(sum)<EPS)
  {
    //cerr << "Colinear stuff" <<endl;
    return true;
  }
  else if (sum*d > 0)
  { return true; }
  else if (sum*d < 0)
  { 
    cerr << "Bad facets" <<endl;
    return false; 
  }

}


G l2Norm(valarray<G> vector){
  G sum = 0.;
  for (size_t i=0;i<vector.size();i++)
  {
    sum += vector[i]*vector[i];
  }
  sum = sqrt(sum);
  return sum;
}


void NormalVector::normalize()
{
//  G sum = 0.;
//  for (size_t i=0;i<normalvec.size();i++)
//  {
//    sum += normalvec[i]*normalvec[i];
//  }
  G sum  = l2Norm(vec);
  if(sum != 0.){
  for (size_t i=0;i<vec.size();i++)
  {
    vec[i] /= sum;
  }
  }
}



size_t findNullDim(std::valarray<std::valarray<G> > & a)
{
  if (a.size()>0)
  {
    for (size_t d=0;d<a[0].size();d++)
    {
      bool nul = true;
      for (size_t i=0;i<a.size();i++)
      {
        if (a[i][d]!=0.)
        {
          nul = false;
          break;
        }
      }
      if (nul) return d;
    }
  }
  return -1;
}

NormalVector Polytope::computeNormal(std::vector<size_t> vertices_)
{
  std::vector<G> result;
  //cerr << "[A]" << flush;
  std::valarray<std::valarray<G> > system(vertices_.size()-1);
  auto const & plast = points[vertices_[vertices_.size()-1]];
  /*cerr << "  Points: " ;
  for (size_t i=0;i<vertices_.size();i++)
  {
   cerr << points[vertices_[i]]<< "  ";  
  }
  cerr << endl;
  cerr << "[B]" << flush;*/
  for (size_t i=0;i<vertices_.size()-1;i++)
  {
    size_t pindex = vertices_[i];
    auto const & p = points[pindex];
    system[i].resize(dimension);
    for (size_t j=0;j<dimension;j++)
    {
      system[i][j] = p[j]-plast[j];
    }
  }
  int intDim = findNullDim(system);
  if (intDim!=-1)
  {
    for (size_t i=0;i<dimension;i++)
    {
      if (i==intDim)      
        result.push_back(1.0);
      else
        result.push_back(0.);
    } 
  }
  else
  {
//    if(dimension == 3){
//        result = crossProduct(system);
//
//    }
//    else{
    forwardSubstitution(system);
    result = reverseElimination(system);
//    }

//    for (size_t i=0;i<system.size();i++)
//    {
//      G comp = system[i][dimension-1];
//      result.push_back(comp);
//    } 
//    result.push_back(-1.);
//    normalize(result);    
  }
  
  
  int d = 0;
  size_t k = 0;
  while((d==0)&&(k<points.size()))
  {
    point_type & p = points[k];
    d = testDirectionality(plast,result,p);
    k++;
  }
  
  k = 0;
  for(k=0; k<points.size();k++)
  {
    point_type & p = points[k];
    if (not(badDirectionality(plast,result,p,d))){
        cout << "BadPolyComp\n";
        for (size_t current_dim = 0; current_dim < result.size(); current_dim++){
        result[current_dim] = 0.;
        }
        return result; 
    }
  }
  


  if (d<0)
  {
    for (size_t i=0;i<result.size();i++)
    {
      result[i] *= -1.;  
    }    
  }
  
  /*cerr << "  Normal vec: " ;
  for (size_t i=0;i<result.size();i++)
  {
    cerr << result[i]<< "  ";  
  }
  cerr << endl;*/
  NormalVector finRes(result);
  point_type acc(dimension);
  for (size_t i=0;i<vertices_.size();i++)
  {
    size_t pindex = vertices_[i];
    acc = acc + points[pindex];
  }
  finRes.support = acc/((G)dimension);
  
  points_type ref(vertices_.size());
  for (size_t i=0;i<vertices_.size();i++)
  {
    size_t pindex = vertices_[i];
    point_type & p = points[pindex];
    ref[i].resize(p.size());
    for (size_t j=0;j<p.size();j++)
    {
      ref[i][j] = p[j];
    }
  }
  finRes.setSignatures(ref);
  
  return finRes;
}


std::vector<NormalVector> Polytope::normalVectors()
{
  std::vector<NormalVector > normals;
  std::vector<std::vector<size_t> > facets_ = convexHullFacets();
  if (facets_.size()>0)
  {
    for (size_type i = 0; i < facets_.size(); ++i) 
    {
//      cout<< "Computed facet"<<endl;
      auto const & vertices_ = facets_[i];
//      for(size_t k=0; k< vertices_.size(); k++){
//        cout << points[vertices_[k]]<< endl;
//      }
      NormalVector current_normal = NormalVector(computeNormal(vertices_));      
      if(l2Norm(current_normal.vec)!= 0.){
        normals.push_back(current_normal);
      }
    }
  }    
  return normals;
}


std::vector<std::vector<size_t> > getConvexHullHullFacets(points_type points, int dimension)
{
  std::vector<std::vector<size_t> > result;
  using quick_hull_type = quick_hull< points_type::const_iterator >;
  using std::sqrt;
  quick_hull_type quick_hull_(dimension, sqrt(std::numeric_limits< G >::epsilon()));
  typename quick_hull_type::point_list initial_simplex_;
  initial_simplex_ = quick_hull_.create_simplex(std::begin(points), std::end(points));
  size_type const basis_size_ = initial_simplex_.size();
  if (basis_size_ == dimension+1)
  {
    quick_hull_.create_convex_hull();
    auto const & facets_ = quick_hull_.facets_;
  
    for (size_type i = 0; i < facets_.size(); ++i) 
    {
      auto const & vertices_ = facets_[i].vertices_;
      std::vector<size_t> v;
      for (auto const vertex_ : vertices_) 
      {
        points_type::const_iterator c = points.begin();
        v.push_back(std::distance(c,vertex_));
      }
      result.push_back(v);
    }
  }
  return result;
}

std::vector<std::vector<size_t> > Polytope::convexHullFacets()
{
  std::vector<std::vector<size_t> > result;
  if (points.size()> dimension )
  {    
    result = getConvexHullHullFacets(points,dimension);
    /*if ((result.size()==0) && lastDimHomogenous(points) && dimension > 2)
    {
      result = getConvexHullHullFacets(dropLastDim(points),dimension-1); 
    }*/
  }
  return result;
}


Polytope Polytope::unionSum(Polytope p2)
{
  points_type pres(points.size()+p2.points.size());
  for (int i=0;i<points.size();i++)
  {
    pres[i].resize(dimension);
    for (int j=0;j<dimension;j++)
    {
      pres[i][j] = points[i][j];
    } 
  } 
  for (int i=0;i<p2.points.size();i++)
  {
    pres[points.size()+i].resize(dimension);
    for (int j=0;j<dimension;j++)
    {
      pres[points.size()+i][j] = p2.points[i][j];
    } 
  } 
  return Polytope(dimension,pres);
}

Polytope Polytope::minkovskiSum(Polytope p2)
{
  size_t n1 = points.size();
  size_t n2 = p2.points.size();
  
  points_type pres(n1*n2);
  for (int i=0;i<n1;i++)
  {
    const point_type & pt1 = points[i];
    for (int j=0;j<n2;j++)
    {
      const point_type & pt2 = p2.points[j];
      pres[i*n2+j].resize(dimension);
      for (int k=0;k<dimension;k++)
      {
        pres[i*n2+j][k] = pt1[k]+pt2[k];
      }
    } 
  } 
  return Polytope(dimension,pres);
}

Polytope Polytope::convexSum(Polytope p2)
{
  return unionSum(p2).convexHull();
}

