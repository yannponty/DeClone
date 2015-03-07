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


from hypergraph import *


def addCaseFromString(h, caseString):
    def costTranslation(inStr):
        transDict = {
            'G' : 'AdjGain', 
            'B' : 'AdjBreak', 
            '0' : 'ZERO',
            'IsAdj' : 'IsAdj(v1,v2,adjacencies)',
            'IsntAdj' : 'IsntAdj(v1,v2,adjacencies)'}
        if inStr in transDict:
            return transDict[inStr]
        else:
            return inStr

    splitCaseAtEqual = caseString.split('=')

    LHSStr = splitCaseAtEqual[0].strip()
    if 'c0' in LHSStr:
        LHS = 0
    elif 'c1' in LHSStr:
        LHS = 1

    splitInfo = splitCaseAtEqual[1].split(';')

    RHSList = splitInfo[0].split(':')
    type1List = splitInfo[1].split()
    type2List = splitInfo[2].split()

    # adjacency = None
    # if len(splitInfo) > 3:
    #     adjacency = splitInfo[3]
    # # TODO: deal with the adjacency

    for RHS in RHSList:
        destinations = []
        arguments = []
        costs = []

        for addend in RHS.split('+'):
            addend = addend.strip()
            if 'c0' in addend:
                destinations.append(0)
                arguments.append(tuple(addend[2:].strip('()').split(',')))
            elif 'c1' in addend:
                destinations.append(1)
                arguments.append(tuple(addend[2:].strip('()').split(',')))
            else:
                costs.append(costTranslation(addend.strip()))

        h.addFArc(LHS,
                  destinations,
                  0.,
                  ((type1List,type2List),
                   costs,
                   arguments))

def buildHypergraphFromFile(inputFilename):
    inFile = open(inputFilename, 'r')

    h = FGraph()

    fullCase = ''
    for raw_line in inFile:
        line = raw_line.strip()

        # Ignore empty lines and lines starting with #
        if len(line) > 0 and line[0] != '#':

            fullCase += line

            if line[-1] != ':' and line[-1] != ';':
                # End of case, parse it
                addCaseFromString(h, fullCase)

                # Initialize next case
                fullCase = ''

    inFile.close()
    return h
                
    

def buildBackwardDP(hg,node):
    def transformArgs(j,rhs1,rhs2,(X1,X2),trans=None):
        if rhs1=="v1":
            nX1 = X1
        elif rhs1=="v1a":
            if X1 == "v1a":
                nX1 = "v1"
            elif X1 == "v1b":
                nX1 = "v1s"
            elif X1 == "v1":
                nX1 = "v1p"
        elif rhs1=="v1b":
            if X1 == "v1b":
                nX1 = "v1"
            elif X1 == "v1a":
                nX1 = "v1s"
            elif X1 == "v1":
                nX1 = "v1p"
        if rhs2=="v2":
            nX2 = X2
        elif rhs2=="v2a":
            if X2 == "v2a":
                nX2 = "v2"
            elif X2 == "v2b":
                nX2 = "v2s"
            elif X2 == "v2":
                nX2 = "v2p"
        elif rhs2=="v2b":
            if X2 == "v2b":
                nX2 = "v2"
            elif X2 == "v2a":
                nX2 = "v2s"
            elif X2 == "v2":
                nX2 = "v2p"
        if trans is not None:
            return (j,trans[nX1],trans[nX2])
        else:
            return (j,nX1,nX2)

    trans = {
        "v1" : "i",
        "v2" : "j",
        "v1a" : "a1",
        "v1b" : "b1",
        "v1p" : "p1",
        "v1s" : "s1",
        "v2a" : "a2",
        "v2b" : "b2",
        "v2p" : "p2",
        "v2s" : "s2",
             }

            
    prec = g.getParents(node)
    print "              tmp = INF;"
    print """              if ((p1 == -1) && (p2 == -1))
	      {
		   tmp = ZERO;
	      }"""
    for e in prec:
        (conditions, events, relationship) = e.getType()
        cond1, cond2 = conditions    
        lbl = formatLabel(e.getOrigin(),cond1,cond2, e.getDestinations(),relationship)

        for i,n in enumerate (e.getDestinations()):
            if n==node:
                parent = e.getOrigin()
                nconds = []
                (rhs1, rhs2) = relationship[i]
                if rhs1=='v1':
                    b1 = 'v1' 
                elif rhs1=='v1a':
                    b1 = 'v1p'
                    nconds.append("(p1 != -1)")
                    nconds.append("(v1IsLeftChild)")
                elif rhs1=='v1b':
                    b1 = 'v1p' 
                    nconds.append("(p1 != -1)")
                    nconds.append("(v1IsRightChild)")
                if rhs2=='v2':
                    b2 = 'v2' 
                elif rhs2=='v2a':
                    b2 = 'v2p'
                    nconds.append("(p2 != -1)")
                    nconds.append("(v2IsLeftChild)")
                elif  rhs2=='v2b':
                    b2 = 'v2p' 
                    nconds.append("(p2 != -1)")
                    nconds.append("(v2IsRightChild)")
                nconds.append("(" + "||".join(["(event1[%s] == %s)" % (trans[b1], cond) for cond in cond1]) + ")")
                nconds.append("(" + "||".join(["(event2[%s] == %s)" % (trans[b2], cond) for cond in cond2]) + ")")

                print (("              // b%s[v1][v2] =  " % node) 
                       + " + ".join(["c%s(%s,%s)" % transformArgs(e.getDestinations()[j],rhs1,rhs2,e.getType()[2][j]) for j,m in enumerate (e.getDestinations()) if i!=j])
                       + " +  b%s(%s,%s)" % (parent,b1,b2)
                       + "".join([" + "+str(ev) for ev in events]))  
                
                rhs = reduce(lambda a,b : 'PLUS(%s,%s)' % (a, b),
                             ["C%s[%s][%s]" % transformArgs(e.getDestinations()[j],rhs1,rhs2,e.getType()[2][j],trans) for j,m in enumerate (e.getDestinations()) if i!=j] + 
                             ["B%s[%s][%s]" % (parent, trans[b1], trans[b2])] + 
                             [str(ev) for ev in events] +
                             ["RESCALING_FACTOR((((event1[%s]!=GDup)&&(event2[%s]!=GDup))?1:0))" % (trans[b1], trans[b2])])
                
                print "              if (" + " && ".join(nconds) + ")"
                print "              {"
                print "                  tmp = MIN(tmp, %s, %s, %s, v1->getND(), v2->getND());" % (rhs, '"COMMENT"', 'false')
                print "              }"
    print "              B%s[i][j] = tmp;\n" % node

def formatLabel(case,v1cond,v2cond, dest,args):
    def shorthand(r):
        return r
        if r[-1] in ["a","b"]:
            return r[-1]
        else:
            return ""
    if len(v1cond)!=1:
        v1c = "".join([s[0] for s in v1cond])
    else:
        v1c = v1cond[0]
    if len(v2cond)!=1:
        v2c = "".join([s[0] for s in v2cond])
    else:
        v2c = v2cond[0]
    recCalls = ["_C%s%s%s"%(d,shorthand(x1),shorthand(x2)) for (d,(x1,x2)) in zip(dest, args)]
    recCalls.sort()
    result = ("C%s_%s_%s%s"%(case,v1c,v2c,"".join(recCalls))).upper()
    #sys.stderr.write(result+"\n")
    return result

def printEnum(hg):
    labelSet = set()
    for case in hg.getVertices():
        for e in hg.getOutList(case):
            ((v1cond,v2cond),costs,args) = e.getType()
            dest = e.getDestinations()
            labelSet.add(formatLabel(case,v1cond,v2cond, dest,args))
    print "#ifndef OPERATIONS_HH" 
    print "#define OPERATIONS_HH" 
    print "#include <string>"
    print "#include <vector>"
    print 'using namespace std;'
    print "#define NUM_CASES %d" % len(labelSet)
    print "typedef enum {"
    for lbl in sorted(labelSet):
        print "   %s," % (lbl)
    print "} CaseLabel;"
    print "string label2String(CaseLabel c);"
    print "vector<CaseLabel> C1Labels();"
    print "#endif" 

def printPrintingLabels(hg):
    print '#include "OperationsList.hh"'
    print 'using namespace std;'
    print "string label2String(CaseLabel c) {"
    print "    switch(c) {"

    for case in hg.getVertices():
        for e in hg.getOutList(case):
            ((v1cond,v2cond),costs,args) = e.getType()
            dest = e.getDestinations()
            lbl = formatLabel(case,v1cond,v2cond, dest,args)
            print "     case (%s):{" % (lbl)
            print '       return "%s";' % (lbl)
            print "     }"
            print "     break;"

    print "    }"
    print "}"
    print "vector<CaseLabel> C1Labels() {"
    print "  vector<CaseLabel> result;"
    for case in hg.getVertices():
        for e in hg.getOutList(case):
            ((v1cond,v2cond),costs,args) = e.getType()
            dest = e.getDestinations()
            lbl = formatLabel(case,v1cond,v2cond, dest,args)
            if lbl.startswith("C1"):
                print "  result.push_back(%s);" % (lbl)
    print "  return result;"
    print "}"
    

def GenCPP(hg, backwardsSwitch):
    trans = {
        "v1" : "i",
        "v2" : "j",
        "v1a" : "a1",
        "v1b" : "b1",
        "v1p" : "p1",
        "v1s" : "s1",
        "v2a" : "a2",
        "v2b" : "b2",
        "v2p" : "p2",
        "v2s" : "s2",
             }
    outType = "RESULT_TYPE"
    if (backwardsSwitch):
        outType = "RESULT_TYPE***"
    print """

RESULT_TYPE** allocateMatrix(RecTree * t1, RecTree * t2)
{
	RESULT_TYPE** R = new RESULT_TYPE*[t1->size()];
	for(int i=0;i<t1->size();i++) 
	{ R[i] = new RESULT_TYPE[t2->size()]; }
	return R;
}

void deleteMatrix(RESULT_TYPE** R, RecTree * t1, RecTree * t2)
{
	for(int i=0;i<t1->size();i++) 
	{
		delete[] R[i];
	}
	delete[] R;
}

%s computeMatrix(RecTree * t1, RecTree * t2, bool adjacent, map<string, map<string,string> > & adjacencies){
    RESULT_TYPE** C0 = allocateMatrix(t1,t2);
    RESULT_TYPE** C1 = allocateMatrix(t1,t2);"""%(outType)

    if (backwardsSwitch):
        print """
    RESULT_TYPE** B0 = allocateMatrix(t1,t2);
    RESULT_TYPE** B1 = allocateMatrix(t1,t2);"""

    print """
    vector<RecTree*> Dfo1 = computeDepthFirstOrder(t1);
    vector<RecTree*> Dfo2 = computeDepthFirstOrder(t2);

    vector<EventType> event1;
    vector<EventType> event2;
    for(int i=0;i<Dfo1.size();i++) 
    {
        event1.push_back(Dfo1[i]->getEvent());
    }
    for(int j=0;j<Dfo2.size();j++) 
    {
        event2.push_back(Dfo2[j]->getEvent());
    }

    for(int i=0;i<Dfo1.size();i++) 
    {
        RecTree * v1 = Dfo1[i];
        for(int j=0;j<Dfo2.size();j++) 
        {
            RESULT_TYPE tmp;
	    RecTree * v2 = Dfo2[j];
            if (!sameSpecies(v1,v2))
            { C1[i][j] = INF; C0[i][j] = INF; }
            else
            {
              EventType t1 = event1[i];
              EventType t2 = event2[j];
              int a1 = (v1->getLeft()?  v1->getLeft()->getIndex():-1);
              int b1 = (v1->getRight()? v1->getRight()->getIndex():-1);
              int a2 = (v2->getLeft()?  v2->getLeft()->getIndex():-1);
              int b2 = (v2->getRight()? v2->getRight()->getIndex():-1);
"""
    nbop = 1
    for case in hg.getVertices():
        print "              tmp = INF;"
        for e in hg.getOutList(case):
            ((v1cond,v2cond),costs,args) = e.getType()
            dest = e.getDestinations()
            precond  = "("+" || ".join(["(t1==%s)"%(et) for et in v1cond ])+")"
            precond += " && ("+" || ".join(["(t2==%s)"%(et) for et in v2cond ])+")"
            if "v1a" in [n1 for n1,n2 in args]:
                precond += " && (a1 != -1)"
            if "v1b" in [n1 for n1,n2 in args]:
                precond += " && (b1 != -1)"
            if "v2a" in [n2 for n1,n2 in args]:
                precond += " && (a2 != -1)"
            if "v2b" in [n2 for n1,n2 in args]:
                precond += " && (b2 != -1)"
            rescaling_factor = "(((t1!=GDup)&&(t2!=GDup))?1:0)"
            if (("GLos" in v1cond) or ("GLos" in v2cond)):
                rescaling_factor = "(v1->numNonGDup() + v2->numNonGDup() - 1)"
                
            rhs = reduce(lambda a,b : 'PLUS(%s,%s)' % (a, b),
                         ["C%s[%s][%s]" % (d, trans[x1], trans[x2]) for (d,(x1,x2)) in zip(dest, args)] + costs + ["RESCALING_FACTOR(%s)"%(rescaling_factor)]  )
            print "              // Op#%s: c%s[v1,v2] = %s"%(nbop,case,(" + ".join(["c%s[%s,%s]"%(d, x1, x2) for (d,(x1,x2)) in zip(dest,args)]+costs)))
            print "              if (%s)"%(precond)
            print "              {"
            print "                 //cout << \"     \" << \"%s; Old val:\"<< tmp << \"; Cand: \" << %s << endl;"%(precond,rhs)
            print "                 tmp = MIN(tmp, %s, (\"%s|\"+v1->getSpecies()), %s, v1->getND(), v2->getND());"%(rhs, formatLabel(case,v1cond,v2cond, dest,args), (repr(case==1)).lower())
            print "              }"
            nbop += 1
        print "              C%s[%s][%s] = tmp;\n"%(case,trans["v1"], trans["v2"])

    print """
            }
	    //cerr << "Node: " << v1->getLabel() << " (" << v1->getSpecies() << ") " << v2->getLabel() << " (" << v2->getSpecies() << ")" << endl;
        }
    }"""

    if (backwardsSwitch):
        print """
    for(int i=Dfo1.size()-1; i>=0; i--) 
    {
        RecTree * v1 = Dfo1[i];
        for(int j=Dfo2.size()-1; j>=0; j--) 
        {
            RESULT_TYPE tmp;
	    RecTree * v2 = Dfo2[j];
            if (!sameSpecies(v1,v2))
            { B1[i][j] = INF; B0[i][j] = INF; }
            else
            {
              EventType t1 = event1[i];
              EventType t2 = event2[j];
              int a1 = (v1->getLeft()?  v1->getLeft()->getIndex():-1);
              int b1 = (v1->getRight()? v1->getRight()->getIndex():-1);
              int a2 = (v2->getLeft()?  v2->getLeft()->getIndex():-1);
              int b2 = (v2->getRight()? v2->getRight()->getIndex():-1);

              RecTree *  v1p = v1->getParent();
              RecTree *  v2p = v2->getParent();
              int p1 = -1;
              bool v1IsLeftChild = false;
              bool v1IsRightChild = false;
              int s1 = -1;
	      if (v1p)
	      {
		   p1 = v1p->getIndex();
		   v1IsLeftChild = (v1p->getLeft()==v1);
		   v1IsRightChild = !v1IsLeftChild;
		   s1 = (v1IsLeftChild?v1p->getRight()->getIndex():v1p->getLeft()->getIndex());
	      }
	      int p2 = -1;
	      bool v2IsLeftChild = false;
	      bool v2IsRightChild = false;
	      int s2 = -1;
	      if (v2p)
	      {
		   p2 = v2p->getIndex();
		   v2IsLeftChild = (v2p->getLeft()==v2);
		   v2IsRightChild = !v2IsLeftChild;
		   s2 = (v2IsLeftChild?v2p->getRight()->getIndex():v2p->getLeft()->getIndex());
	      }
"""
        buildBackwardDP(hg,1)
        buildBackwardDP(hg,0)
        print """
            }
	}
    }"""

        print """
    RESULT_TYPE*** W = new RESULT_TYPE**[t1->size()];
    for(int i=0;i<t1->size();i++) 
    { 
      W[i] = new RESULT_TYPE*[t2->size()]; 
      for(int j=0;j<t2->size();j++) 
      { 
        W[i][j] = new RESULT_TYPE[%s]; 
        for(int k=0;k<%s;k++)
        {
            W[i][j][k] = INF;
        }
      }
    }
    for(int i=0;i<Dfo1.size();i++) 
    {
	RecTree * v1 = Dfo1[i];
	for(int j=0;j<Dfo2.size();j++) 
	{
            RESULT_TYPE tmp;
	    RecTree * v2 = Dfo2[j];
            if (!sameSpecies(v1,v2))
            {
                for(int k=0;k<%s;k++)
                {
                    W[i][j][k] = INF;
                }
            }
            else
            {
              EventType t1 = event1[i];
              EventType t2 = event2[j];
              int a1 = (v1->getLeft()?  v1->getLeft()->getIndex():-1);
              int b1 = (v1->getRight()? v1->getRight()->getIndex():-1);
              int a2 = (v2->getLeft()?  v2->getLeft()->getIndex():-1);
              int b2 = (v2->getRight()? v2->getRight()->getIndex():-1);
"""%(nbop-1,nbop-1,nbop-1)

        nbop = 1
        for case in hg.getVertices():
            for e in hg.getOutList(case):
                ((v1cond,v2cond),costs,args) = e.getType()
                dest = e.getDestinations()
                precond  = "("+" || ".join(["(t1==%s)"%(et) for et in v1cond ])+")"
                precond += " && ("+" || ".join(["(t2==%s)"%(et) for et in v2cond ])+")"
                if "v1a" in [n1 for n1,n2 in args]:
                    precond += " && (a1 != -1)"
                if "v1b" in [n1 for n1,n2 in args]:
                    precond += " && (b1 != -1)"
                if "v2a" in [n2 for n1,n2 in args]:
                    precond += " && (a2 != -1)"
                if "v2b" in [n2 for n1,n2 in args]:
                    precond += " && (b2 != -1)"
                rescaling_factor = "(((t1!=GDup)&&(t2!=GDup))?1:0)"
                if (("GLos" in v1cond) or ("GLos" in v2cond)):
                   rescaling_factor = "(v1->numNonGDup() + v2->numNonGDup() - 1)"

                rhsW = reduce(lambda a,b : 'PLUS(%s,%s)' % (a, b),
                              ["C%s[%s][%s]" % (d, trans[x1], trans[x2]) for (d,(x1,x2)) in zip(dest, args)] 
                              + costs 
                              + ["RESCALING_FACTOR(%s)"%(rescaling_factor)]
                              + ["B%s[%s][%s]" % (case, trans["v1"], trans["v2"])])

                lbl = formatLabel(case,v1cond,v2cond, dest,args)

                print "              if (%s)"%(precond)
                print "              {"
                print "                 W[%s][%s][%s]"%(trans["v1"], trans["v2"], lbl);
                print "                     = %s;" % (rhsW)
                print "              }"
                nbop += 1
        print """
            }
        }
    }"""



        print"""
    RESULT_TYPE*** finalResult = W;
    // We don't want to do that: this matrix is the output!
    /*for(int i=0;i<t1->size();i++) 
    { 
      for(int j=0;j<t2->size();j++) 
      { 
        delete W[i][j]; 
      }
      delete W[i]; 
    }
    delete W;*/ 
    deleteMatrix(B0,t1,t2);
    deleteMatrix(B1,t1,t2);"""

    else:
        print """
    RESULT_TYPE finalResult = (adjacent?C1[Dfo1.size()-1][Dfo2.size()-1]:C0[Dfo1.size()-1][Dfo2.size()-1]);"""

    print """
    deleteMatrix(C0,t1,t2);
    deleteMatrix(C1,t1,t2);
    return finalResult;
}
"""


if __name__ == "__main__":
    if len(sys.argv)>1:
        g = buildHypergraphFromFile('DeCo-DP.txt')
        if sys.argv[1] == "-c":
            GenCPP(g, False)
        elif sys.argv[1] == "-b":
            GenCPP(g, True)
        elif sys.argv[1] == "-e":
            printEnum(g)
        elif sys.argv[1] == "-l":
            printPrintingLabels(g)
        else:
            print g
