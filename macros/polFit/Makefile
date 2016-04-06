CXX=$(shell root-config --cxx --cflags)
LIBS=$(shell root-config --libs)
INCDIR=-I../../interface

%.o : %.cc
	$(CXX) -c $<

#all: polRapPtPlot polGenRecFitPlot EvaluateSyst PlotFinalResults MattsResults AverageSystematics ChangeTGraph AlterPPD
all: polRapPtPlot polGenRecFitPlot PlotFinalResults

polRapPtPlot: polRapPtPlot.cc
	$(CXX) $^ $(INCDIR) -o $@ $(LIBS) -lFoam -lMinuit

polGenRecFitPlot: polGenRecFitPlot.cc
	$(CXX) $^ $(INCDIR) -o $@ $(LIBS) -lFoam -lMinuit

EvaluateSyst: EvaluateSyst.cc
	$(CXX) $^ $(INCDIR-) -o $@ $(LIBS) -lFoam -lMinuit

PlotFinalResults: PlotFinalResults.cc
	$(CXX) $^ $(INCDIR) -o $@ $(LIBS) -lFoam -lMinuit

MattsResults: MattsResults.cc
	$(CXX) $^ $(INCDIR) -o $@ $(LIBS) -lFoam -lMinuit

AverageSystematics: AverageSystematics.cc
	$(CXX) $^ $(INCDIR) -o $@ $(LIBS) -lFoam -lMinuit

ChangeTGraph: ChangeTGraph.cc
	$(CXX) $^ $(INCDIR) -o $@ $(LIBS) -lFoam -lMinuit

AlterPPD: AlterPPD.cc
	$(CXX) $^ $(INCDIR) -o $@ $(LIBS) -lFoam -lMinuit

clean:
	rm polRapPtPlot polGenRecFitPlot EvaluateSyst PlotFinalResults MattsResults AverageSystematics ChangeTGraph AlterPPD*.o
