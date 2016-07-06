CXX=$(shell root-config --cxx --cflags)
LIBS=$(shell root-config --libs)
INCDIR=-I../../interface

# conditionally define the preprocessor flag. CHIC_BINNING is not defined a default value
# is defined in commonVar.h
PREPFLAGS=
ifdef ${CHIC_BINNING}
	PREFLAGS=-D"USE_CHIC_BINNING=${CHIC_BINNING}"
endif

%.o : %.cc
	$(CXX) -c $<

#all: polRapPtPlot polGenRecFitPlot EvaluateSyst PlotFinalResults MattsResults AverageSystematics ChangeTGraph AlterPPD
all: polRapPtPlot polGenRecFitPlot PlotFinalResults

polRapPtPlot: polRapPtPlot.cc
	$(CXX) $^ $(INCDIR) -o $@ $(LIBS) -lFoam -lMinuit $(PREPFLAGS)

polGenRecFitPlot: polGenRecFitPlot.cc
	$(CXX) $^ $(INCDIR) -o $@ $(LIBS) -lFoam -lMinuit $(PREPFLAGS)

EvaluateSyst: EvaluateSyst.cc
	$(CXX) $^ $(INCDIR-) -o $@ $(LIBS) -lFoam -lMinuit $(PREPFLAGS)

PlotFinalResults: PlotFinalResults.cc
	$(CXX) $^ $(INCDIR) -o $@ $(LIBS) -lFoam -lMinuit $(PREPFLAGS)

MattsResults: MattsResults.cc
	$(CXX) $^ $(INCDIR) -o $@ $(LIBS) -lFoam -lMinuit $(PREPFLAGS)

AverageSystematics: AverageSystematics.cc
	$(CXX) $^ $(INCDIR) -o $@ $(LIBS) -lFoam -lMinuit $(PREPFLAGS)

ChangeTGraph: ChangeTGraph.cc
	$(CXX) $^ $(INCDIR) -o $@ $(LIBS) -lFoam -lMinuit $(PREPFLAGS)

AlterPPD: AlterPPD.cc
	$(CXX) $^ $(INCDIR) -o $@ $(LIBS) -lFoam -lMinuit $(PREPFLAGS)

clean:
	rm polRapPtPlot polGenRecFitPlot EvaluateSyst PlotFinalResults MattsResults AverageSystematics ChangeTGraph AlterPPD*.o
