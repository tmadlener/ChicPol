# CXX=g++
CXX=$(shell root-config --cxx)
CXXFLAGS=-Wall -O2 -g -std=c++0x #-std=c++11 # no C++11 possible for gcc version 4.4.7

ROOTCXXFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs)
ROOTINC=$(shell root-config --incdir)

# %.o : %.C
# 	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) $(ROOTINC) -c $<

all: compPol

compPol: compPol.C
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) -o $@ -I$(ROOTINC) $(ROOTLIBS) compPol.C

clean:
	rm compPol
