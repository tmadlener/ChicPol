CXX=$(shell root-config --cxx --cflags) -Wall
LIBS=$(shell root-config --libs)

all: runMuonTests

runMuonTests: runMuonTests.cc
	$(CXX) $^ -o $@ $(LIBS)
