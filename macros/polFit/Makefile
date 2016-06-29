CXX=$(shell root-config --cxx --cflags)
LIBS=$(shell root-config --libs)
INCDIR=-I../../interface

%.o : %.cc
	$(CXX) -c $<

#all: polRapPtPlot polGenRecFitPlot EvaluateSyst PlotFinalResults MattsResults AverageSystematics ChangeTGraph AlterPPD
all: polRapPtPlot polGenRecFitPlot PlotFinalResults

polRapPtPlot: polRapPtPlot.cc
	$(CXX) $^ $(INCDIR) -o $@ $(LIBS) -lFoam -lMinuit -D"USE_CHIC_BINNING=${CHIC_BINNING}"

polGenRecFitPlot: polGenRecFitPlot.cc
	$(CXX) $^ $(INCDIR) -o $@ $(LIBS) -lFoam -lMinuit -D"USE_CHIC_BINNING=${CHIC_BINNING}"

EvaluateSyst: EvaluateSyst.cc
	$(CXX) $^ $(INCDIR-) -o $@ $(LIBS) -lFoam -lMinuit -D"USE_CHIC_BINNING=${CHIC_BINNING}"

PlotFinalResults: PlotFinalResults.cc
	$(CXX) $^ $(INCDIR) -o $@ $(LIBS) -lFoam -lMinuit -D"USE_CHIC_BINNING=${CHIC_BINNING}"

MattsResults: MattsResults.cc
	$(CXX) $^ $(INCDIR) -o $@ $(LIBS) -lFoam -lMinuit -D"USE_CHIC_BINNING=${CHIC_BINNING}"

AverageSystematics: AverageSystematics.cc
	$(CXX) $^ $(INCDIR) -o $@ $(LIBS) -lFoam -lMinuit -D"USE_CHIC_BINNING=${CHIC_BINNING}"

ChangeTGraph: ChangeTGraph.cc
	$(CXX) $^ $(INCDIR) -o $@ $(LIBS) -lFoam -lMinuit -D"USE_CHIC_BINNING=${CHIC_BINNING}"

AlterPPD: AlterPPD.cc
	$(CXX) $^ $(INCDIR) -o $@ $(LIBS) -lFoam -lMinuit -D"USE_CHIC_BINNING=${CHIC_BINNING}"

clean:
	rm polRapPtPlot polGenRecFitPlot EvaluateSyst PlotFinalResults MattsResults AverageSystematics ChangeTGraph AlterPPD*.o
