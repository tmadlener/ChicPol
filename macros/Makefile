CXX=$(shell root-config --cxx --cflags) -Wall
LIBS=$(shell root-config --libs)
INCS=-I../interface

%.o : %.cc
	$(CXX) -c $<

all: runData runChiData runWorkspace runMassFit runLifetimeFit runPlotMassLifetime runBkgHistos PlotFitPar PlotCosThetaPhiBG PlotCosThetaPhiDistribution runChiMassLifetimeFit runDefineRegionsAndFractions runPlotDataDistributions runPlotJpsiMassLifetime PlotJpsiFitPar PlotMassRapPtBG runMC runBkgHistos_new PlotPsiPtDist

PlotJpsiFitPar: PlotJpsiFitPar.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore -D"USE_CHIC_BINNING=${CHIC_BINNING}"

runPlotDataDistributions: runPlotDataDistributions.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore -D"USE_CHIC_BINNING=${CHIC_BINNING}"

runPlotMassLifetime: runPlotMassLifetime.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore -D"USE_CHIC_BINNING=${CHIC_BINNING}"

runPlotJpsiMassLifetime: runPlotJpsiMassLifetime.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore -D"USE_CHIC_BINNING=${CHIC_BINNING}"

runData: runData.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -D"USE_CHIC_BINNING=${CHIC_BINNING}"

runChiData: runChiData.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -D"USE_CHIC_BINNING=${CHIC_BINNING}"

runWorkspace: runWorkspace.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore -D"USE_CHIC_BINNING=${CHIC_BINNING}"

runMassFit: runMassFit.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore -D"USE_CHIC_BINNING=${CHIC_BINNING}"

runLifetimeFit: runLifetimeFit.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore -D"USE_CHIC_BINNING=${CHIC_BINNING}"

runChiMassLifetimeFit: runChiMassLifetimeFit.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore -D"USE_CHIC_BINNING=${CHIC_BINNING}"

runBkgHistos: runBkgHistos.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore -D"USE_CHIC_BINNING=${CHIC_BINNING}"

PlotFitPar: PlotFitPar.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore -D"USE_CHIC_BINNING=${CHIC_BINNING}"

PlotCosThetaPhiBG: PlotCosThetaPhiBG.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore -D"USE_CHIC_BINNING=${CHIC_BINNING}"

PlotCosThetaPhiDistribution: PlotCosThetaPhiDistribution.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore -D"USE_CHIC_BINNING=${CHIC_BINNING}"

runDefineRegionsAndFractions: runDefineRegionsAndFractions.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore -D"USE_CHIC_BINNING=${CHIC_BINNING}"

PlotMassRapPtBG: PlotMassRapPtBG.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore -D"USE_CHIC_BINNING=${CHIC_BINNING}"

runMC: runMC.cc runMC_helperStructs.h
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -D"USE_CHIC_BINNING=${CHIC_BINNING}"

runBkgHistos_new: runBkgHistos_new.cc BkgHistoProducerFactory.h BkgHistoProducer.h BkgHistoProducer_helper.h bkgHistos_helper.h
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lRooFit -lRooFitCore -D"USE_CHIC_BINNING=${CHIC_BINNING}"

PlotPsiPtDist: PlotPsiPtDistribution.cc
	$(CXX) $^ $(INCS) -I$(HOME) -o $@ $(LIBS) -D"USE_CHIC_BINNING=${CHIC_BINNING}"

clean:
	rm runData runChiData runWorkspace runMassFit runLifetimeFit runPlotMassLifetime runBkgHistos PlotFitPar PlotCosThetaPhiBG PlotCosThetaPhiDistribution runChiMassLifetimeFit runDefineRegionsAndFractions runPlotDataDistributions runPlotJpsiMassLifetime PlotJpsiFitPar PlotMassRapPtBG runMC runBkgHistos_new PlotPsiPtDist
