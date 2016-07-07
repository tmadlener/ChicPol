CXX=$(shell root-config --cxx --cflags) -Wall
LIBS=$(shell root-config --libs)
INCS=-I../interface

# conditionally define the preprocessor flag. CHIC_BINNING is not defined a default value
# is defined in commonVar.h
PREPFLAGS=
ifneq ($(strip ${CHIC_BINNING}),)
	PREPFLAGS=-D"USE_CHIC_BINNING=${CHIC_BINNING}"
endif

%.o : %.cc
	$(CXX) -c $<

all: runData runChiData runWorkspace runMassFit runLifetimeFit runPlotMassLifetime runBkgHistos PlotFitPar PlotCosThetaPhiBG PlotCosThetaPhiDistribution runChiMassLifetimeFit runDefineRegionsAndFractions runPlotDataDistributions runPlotJpsiMassLifetime PlotJpsiFitPar PlotMassRapPtBG runMC runBkgHistos_new PlotPsiPtDist

PlotJpsiFitPar: PlotJpsiFitPar.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore $(PREPFLAGS)

runPlotDataDistributions: runPlotDataDistributions.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore $(PREPFLAGS)

runPlotMassLifetime: runPlotMassLifetime.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore $(PREPFLAGS)

runPlotJpsiMassLifetime: runPlotJpsiMassLifetime.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore $(PREPFLAGS)

runData: runData.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit $(PREPFLAGS)

runChiData: runChiData.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit $(PREPFLAGS)

runWorkspace: runWorkspace.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore $(PREPFLAGS)

runMassFit: runMassFit.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore $(PREPFLAGS)

runLifetimeFit: runLifetimeFit.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore $(PREPFLAGS)

runChiMassLifetimeFit: runChiMassLifetimeFit.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore $(PREPFLAGS)

runBkgHistos: runBkgHistos.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore $(PREPFLAGS)

PlotFitPar: PlotFitPar.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore $(PREPFLAGS)

PlotCosThetaPhiBG: PlotCosThetaPhiBG.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore $(PREPFLAGS)

PlotCosThetaPhiDistribution: PlotCosThetaPhiDistribution.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore $(PREPFLAGS)

runDefineRegionsAndFractions: runDefineRegionsAndFractions.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore $(PREPFLAGS)

PlotMassRapPtBG: PlotMassRapPtBG.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore $(PREPFLAGS)

runMC: runMC.cc runMC_helperStructs.h
	$(CXX) $^ $(INCS) -o $@ $(LIBS) $(PREPFLAGS)

runBkgHistos_new: runBkgHistos_new.cc BkgHistoProducerFactory.h BkgHistoProducer.h BkgHistoProducer_helper.h bkgHistos_helper.h
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lRooFit -lRooFitCore $(PREPFLAGS)

# PlotPsiPtDist: PlotPsiPtDistribution.cc
# 	$(CXX) $^ $(INCS) -I$(HOME) -o $@ $(LIBS) $(PREPFLAGS)

clean:
	rm runData runChiData runWorkspace runMassFit runLifetimeFit runPlotMassLifetime runBkgHistos PlotFitPar PlotCosThetaPhiBG PlotCosThetaPhiDistribution runChiMassLifetimeFit runDefineRegionsAndFractions runPlotDataDistributions runPlotJpsiMassLifetime PlotJpsiFitPar PlotMassRapPtBG runMC runBkgHistos_new PlotPsiPtDist
