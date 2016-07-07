# CXX=g++
CXX=$(shell root-config --cxx)
CXXFLAGS=-Wall -O2 -g -std=c++0x #-std=c++11 # no C++11 possible for gcc version 4.4.7

ROOTCXXFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs)
ROOTINC=$(shell root-config --incdir)
ROOFITLIBS=-lRooFit -lRooFitCore

# conditionally define the preprocessor flag. CHIC_BINNING is not defined a default value
# is defined in commonVar.h
PREPFLAGS=
ifneq ($(strip ${CHIC_BINNING}),)
	PREPFLAGS=-D"USE_CHIC_BINNING=${CHIC_BINNING}"
endif

# %.o : %.C
# 	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) $(ROOTINC) -c $<

all: compPol compWS simpleFit MCFitComp printWS testStatEff plotEff plotJpsiDists

compPol: compPol.C
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) -o $@ -I$(ROOTINC) $(ROOTLIBS) $^ $(PREPFLAGS)

compWS: compWorkspace.C
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) -o $@ -I$(ROOTINC) $(ROOTLIBS) $(ROOFITLIBS) $^ $(PREPFLAGS)

simpleFit: simpleChicMassFit.cc
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) -o $@ -I$(ROOTINC) -I../interface $(ROOTLIBS) $(ROOFITLIBS) $^ $(PREPFLAGS)

MCFitComp: MCfitComp.cc
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) -o $@ -I$(ROOTINC) -I../interface $(ROOTLIBS) $(ROOFITLIBS) $^ $(PREPFLAGS)

printWS: printWorkspace.cc
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) -o $@ -I$(ROOTINC) -I../interface -I$(HOME) $(ROOTLIBS) $(ROOFITLIBS) $^ $(PREPFLAGS)

testStatEff: testStatVarEff.cc ../interface/effsAndCuts_chi.h
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) -o $@ -I$(ROOTINC) -I../interface -I$(HOME) $(ROOTLIBS) testStatVarEff.cc $(PREPFLAGS)

plotEff: plotEff.C
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) -o $@ -I$(ROOTINC) $(ROOTLIBS) $^ $(PREPFLAGS)

plotJpsiDists: plotJpsiDistributions.cc
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) -o $@ -I$(ROOTINC) -I../interface -I$(HOME)/snippets $(ROOTLIBS) $(ROOFITLIBS) $^ $(PREPFLAGS)

clean:
	rm compPol compWS simpleFit MCFitComp printWS testStatEff plotEff plotJpsiDists
