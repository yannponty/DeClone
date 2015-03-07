
MAIN_SOURCE = src/ProbaReconciliations.cc
SOURCE = src/Trees.cc src/DPParsimony.cc src/DPCount.cc src/DPStochasticBacktrack.cc src/DPGenAll.cc src/EditTrees.cc src/DPInsideOutside.cc src/SVGDriver.cc src/utils.cc
OBJS = $(SOURCE:.cc=.o)
EXEC = $(MAIN_SOURCE:.cc=)

DECO_SOURCES = src/RecTrees.cc src/DeClone-parsimony.cc src/DeClone-coopts.cc src/DeClone-all.cc src/DeClone-count.cc src/DeClone-countcoopts.cc src/DeClone-inside.cc src/DeClone-outside.cc src/DeClone-stochastic.cc src/AdjacencyTrees.cc src/OperationsList.cc
ALL_DECO_SOURCES = $(DECO_SOURCES) $(SOURCE)
DECO_OBJS = $(ALL_DECO_SOURCES:.cc=.o)

POLYTOPE_SOURCES = src/ConvexPolytope.cc src/DeClone-polytope-adj.cc src/DeClone-polytope.cc
POLYTOPE_OBJS = $(POLYTOPE_SOURCES:.cc=.o)

GENERATED_HH_ENUM = src/OperationsList.hh
GENERATED_CC_ENUM = src/OperationsList.cc
GENERATED_DP = src/DeCoDP.cc
GENERATED_DP_OUTSIDE = src/DeCoDP-outside.cc
PYTHON_HG_GEN = DPDesign.py
DP_SOURCE = DeCo-DP.txt

PRODUCED = $(OBJS) $(EXEC) $(DECO_OBJS)

COMPILER = g++ -g -static-libgcc -static-libstdc++
PYTHON = python

all: UsePolytope 

ProbaReconciliations: $(OBJS) $(MAIN_SOURCE)
	$(COMPILER) $(OBJS) $(MAIN_SOURCE) -o $(EXEC)

%.o: %.cc
	$(COMPILER) -c $< -o $@

$(GENERATED_DP):  $(PYTHON_HG_GEN) $(DP_SOURCE) 
	$(PYTHON) $(PYTHON_HG_GEN) -c >  $(GENERATED_DP)

$(GENERATED_DP_OUTSIDE):  $(PYTHON_HG_GEN) $(DP_SOURCE) 
	$(PYTHON) $(PYTHON_HG_GEN) -b >  $(GENERATED_DP_OUTSIDE)

$(GENERATED_HH_ENUM):   $(PYTHON_HG_GEN) $(DP_SOURCE)
	$(PYTHON) $(PYTHON_HG_GEN) -e >  $(GENERATED_HH_ENUM)

$(GENERATED_CC_ENUM):   $(PYTHON_HG_GEN) $(DP_SOURCE)
	$(PYTHON) $(PYTHON_HG_GEN) -l >  $(GENERATED_CC_ENUM)

clean: DeClone-clean
	rm -f $(PRODUCED) 

DeClone: $(GENERATED_DP) $(GENERATED_DP_OUTSIDE) $(GENERATED_HH_ENUM) $(GENERATED_CC_ENUM) ProbaReconciliations $(DECO_OBJS) src/DeClone.cc
	$(COMPILER)  $(DECO_OBJS) src/DeClone.cc -o DeClone

UsePolytope: COMPILER += -std=c++0x -DUSE_POLYTOPE
UsePolytope: $(GENERATED_DP) $(GENERATED_DP_OUTSIDE) $(GENERATED_HH_ENUM) $(GENERATED_CC_ENUM) ProbaReconciliations $(DECO_OBJS) $(POLYTOPE_OBJS) src/DeClone.cc 
	$(COMPILER)  $(DECO_OBJS) $(POLYTOPE_OBJS) src/DeClone.cc -o DeClone

DeClone-clean:
	rm -f src/RecTrees.o src/DeClone src/DeClone.o $(GENERATED_DP) $(GENERATED_DP_OUTSIDE) $(GENERATED_HH_ENUM) $(GENERATED_CC_ENUM) $(POLYTOPE_OBJS)
	rm -f $(PRODUCED)