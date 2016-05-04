CXX=$(shell root-config --cxx --cflags) -Wall
LIBS=$(shell root-config --libs)
INCS=-I../interface

%.o : %.cc
	$(CXX) -c $<

all: runData runChiData runWorkspace runMassFit runLifetimeFit runPlotMassLifetime runBkgHistos PlotFitPar PlotCosThetaPhiBG PlotCosThetaPhiDistribution runChiMassLifetimeFit runDefineRegionsAndFractions runPlotDataDistributions runPlotJpsiMassLifetime PlotJpsiFitPar PlotMassRapPtBG runMC runBkgHistos_new

PlotJpsiFitPar: PlotJpsiFitPar.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore

runPlotDataDistributions: runPlotDataDistributions.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore

runPlotMassLifetime: runPlotMassLifetime.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore

runPlotJpsiMassLifetime: runPlotJpsiMassLifetime.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore

runData: runData.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit

runChiData: runChiData.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit

runWorkspace: runWorkspace.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore

runMassFit: runMassFit.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore

runLifetimeFit: runLifetimeFit.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore

runChiMassLifetimeFit: runChiMassLifetimeFit.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore

runBkgHistos: runBkgHistos.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore

PlotFitPar: PlotFitPar.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore

PlotCosThetaPhiBG: PlotCosThetaPhiBG.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore

PlotCosThetaPhiDistribution: PlotCosThetaPhiDistribution.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore

runDefineRegionsAndFractions: runDefineRegionsAndFractions.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore

PlotMassRapPtBG: PlotMassRapPtBG.cc
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lFoam -lMinuit -lRooFit -lRooFitCore

runMC: runMC.cc runMC_helperStructs.h
	$(CXX) $^ $(INCS) -o $@ $(LIBS)

runBkgHistos_new: runBkgHistos_new.cc BkgHistoProducerFactory.h BkgHistoProducer.h BkgHistoProducer_helper.h bkgHistos_helper.h
	$(CXX) $^ $(INCS) -o $@ $(LIBS) -lRooFit -lRooFitCore

clean:
	rm runData runChiData runWorkspace runMassFit runLifetimeFit runPlotMassLifetime runBkgHistos PlotFitPar PlotCosThetaPhiBG PlotCosThetaPhiDistribution runChiMassLifetimeFit runDefineRegionsAndFractions runPlotDataDistributions runPlotJpsiMassLifetime PlotJpsiFitPar PlotMassRapPtBG runMC runBkgHistos_new
