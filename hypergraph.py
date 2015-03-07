#   DeClone: A software for computing and analyzing ancestral adjacency scenarios.
#   Copyright (C) 2015 Cedric Chauve, Yann Ponty, Ashok Rajaraman, Joao P.P. Zanetti
# 
#   This file is part of DeClone.
#   
#   DeClone is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
# 
#   DeClone is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
# 
#   You should have received a copy of the GNU General Public License
#   along with DeClone.  If not, see <http://www.gnu.org/licenses/>.
# 
#   Contact: <yann.ponty@lix.polytechnique.fr>.
# 
# 
#   DeClone uses the Quickhull algorithm implementation programmed by 
#   Anatoly V. Tomilov. The code is available on <https://bitbucket.org/tomilov/quickhull/src/585267abb3a63794c04fc8325aa9ec9f726112ed/include/quickhull.hpp?at=master>.
# 
#   Contact: <tomilovanatoliy@gmail.com>
# 


import sys,os


class Vertex:
    """Represents a Vertex, identified by a unique name."""

    def __init__(self, name, hypergraph):
        """Creates a new Vertex.

        Keyword arguments:
        name       -- the unique name of the vertex. Must be hashable!
        hypergraph -- the parent hypergraph
        """
        self.name = name
        self.hypergraph = hypergraph
        
    def getName(self):
        """Returns the unique name of this Vertex."""
        return self.name

    def getOutList(self):
        """Returns a list of FEdges leaving this Vertex."""
        return self.hypergraph.getOutList(self.name)

    def __str__(self):
        """Returns a string representation of this Vertex."""
        return str(self.name)

    def __repr__(self):
        """Returns a string representation of this Vertex."""
        return repr(self.name)


class FArc:
    """
    Represents a forward-arc, ie an hyperarc linking
    a single origin vertex to multiple destination vertices.
    An FArc is weighted, and can be assigned a type for the users'
    convenience.
    """

    def __init__(self,hypergraph,origin,destinations,weight=0.0,type=None):
        """
        Creates a new FArc from the name of its single origin vertex
        and the names of its multiple destinations vertices.
        Vertices are identified by their unique names and must be hashable.
        Vertices that are absent from the graph will be created.

        Keyword arguments:
        hypergraph   -- the parent hypergraph
        origin       -- the unique name of the origin vertex. Must be hashable!
        destinations -- the unique names of the destination vertices. Must be hashable!
        weight       -- the real-valued weight of this arc (default: 0.0)
        type         -- user-provided type for this arc (default: None)
        """
        self.hypergraph = hypergraph
        self.weight = weight
        self.type = type
        vo = self.hypergraph.addVertex(origin)
        vdest = []
        for name in destinations:
            vdest.append(self.hypergraph.addVertex(name))
        self.origin = vo
        self.destinations = vdest

    def getOrigin(self):
        """Returns the unique name of the origin Vertex."""
        return self.origin.getName()

    def getDestinations(self):
        """Returns as a list the unique names of the destination vertices."""
        return [v.getName() for v in self.destinations]

    def getType(self):
        """Returns the type of this Arc."""
        return self.type

    def getWeight(self):
        """Returns the weight of this Arc."""
        return self.weight
        
    def __repr__(self):
        """Returns a string representation of this Arc."""
        return repr(self.origin)+"->"+repr(self.destinations)+"(w:%.1f)"%(self.weight)+" t:%s"%(repr(self.type))
        #return repr(self.origin)+" -> "+repr(self.destinations)
            
class FGraph:
    """ Represents a Forward Hypergraph, ie a graph whose arcs join
    a single origin vertex to multiple destination vertices.
    """
    def __init__(self):
        """Creates a new FGraph."""
        self.vertices = {}
        self.arcs = {}

    def getVertices(self):
        """Returns the list of vertices."""
        return self.vertices.keys()

    def addVertex(self, name):
        """Adds a new vertex to the FGraph.
        Keyword arguments:
        name   -- the vertex name
        """
        if name not in self.vertices:
            self.vertices[name] = Vertex(name,self)
        return self.vertices[name]
    
    
    def addFArc(self, origin, destinations, weight=0.0, type=None):
        """
        Creates a new FArce from the names of its single origin vertex
        and its multiple destinations vertices. Vertices are identified by
        their names and must be hashable. Vertices that are absent from the
        graph will be created on demand.

        Keyword arguments:
        origin       -- the unique name of the origin vertex. Must be hashable!
        destinations -- the unique names of the destination vertices. Must be hashable!
        weight       -- the real-valued weight of this arc (default: 0.0)
        type         -- user-provided type for this arc (default: None)
        """
        if origin not in self.arcs:
            self.arcs[origin] = []
        narc = FArc(self,origin,destinations,weight,type)
        self.arcs[origin].append(narc)
        return narc


    def getOutList(self,name):
        """Returns a list of FEdges leaving a given Vertex identified by its unique name.

          Keyword arguments:
          name       -- the unique name of the origin vertex. Must be hashable!
        """
        if name in self.arcs:
            return self.arcs[name]
        return []
    

    def getParents(self,name):
        """Returns a list of FEdges leading to a given Vertex identified by its unique name.

          Keyword arguments:
          name       -- the unique name of the origin vertex. Must be hashable!
        """
        result = []
        names = self.arcs.keys()
        names.sort()
        for v in names:
            for e in self.arcs[v]:
              if name in e.getDestinations():
                  result.append(e)
        return result
    
    def __str__(self):
        """Returns a string representation of this HyperGraph"""
        tmp = ""
        names = self.arcs.keys()
        names.sort()
        for v in names:
            for e in self.arcs[v]:
              tmp += "  "+str(e)+"\n"
        return "Vertices:\n"+str(self.vertices.keys())+"\nArcs:\n"+tmp

def memoize(function):
    """Generically memoizes a function results."""
    cache = {}
    def decorated_function(*args):
        if args in cache:
            return cache[args]
        else:
            val = function(*args)
            cache[args] = val
            return val
    return decorated_function

def testHypergraph():
    """Creates and return a dummy hypergraph"""
    g = FGraph()
    g.addFArc(0,[1])
    g.addFArc(0,[2])
    g.addFArc(0,[3,4])
    g.addFArc(1,[5,6])
    g.addFArc(2,[5,7],2)
    g.addFArc(3,[5])
    g.addFArc(4,[7])
    g.addFArc(4,[8])
    g.addFArc(5,[],1)
    g.addFArc(5,[9])
    g.addFArc(5,[10])
    g.addFArc(6,[])
    g.addFArc(7,[])    
    g.addFArc(8,[])
    g.addFArc(9,[])
    g.addFArc(10,[])
    return g

