CXX=$(shell root-config --cxx --cflags)
LIBS=$(shell root-config --libs) 

%.o : %.cc
	$(CXX) -c $<

all: runData runChiData runWorkspace runMassFit runLifetimeFit runPlotMassLifetime runBkgHistos PlotFitPar PlotCosThetaPhiBG PlotCosThetaPhiDistribution runChiMassLifetimeFit runDefineRegionsAndFractions runPlotDataDistributions runPlotJpsiMassLifetime PlotJpsiFitPar

PlotJpsiFitPar: PlotJpsiFitPar.cc
	$(CXX) $^ -o $@ $(LIBS) `root-config --cflags` -lFoam -lMinuit -lRooFit -lRooFitCore

runPlotDataDistributions: runPlotDataDistributions.cc
	$(CXX) $^ -o $@ $(LIBS) `root-config --cflags` -lFoam -lMinuit -lRooFit -lRooFitCore

runPlotMassLifetime: runPlotMassLifetime.cc
	$(CXX) $^ -o $@ $(LIBS) `root-config --cflags` -lFoam -lMinuit -lRooFit -lRooFitCore

runPlotJpsiMassLifetime: runPlotJpsiMassLifetime.cc
	$(CXX) $^ -o $@ $(LIBS) `root-config --cflags` -lFoam -lMinuit -lRooFit -lRooFitCore

runData: runData.cc 
	$(CXX) $^ -o $@ $(LIBS) `root-config --cflags` -lFoam -lMinuit

runChiData: runChiData.cc 
	$(CXX) $^ -o $@ $(LIBS) `root-config --cflags` -lFoam -lMinuit

runWorkspace: runWorkspace.cc
	$(CXX) $^ -o $@ $(LIBS) `root-config --cflags` -lFoam -lMinuit -lRooFit -lRooFitCore

runMassFit: runMassFit.cc
	$(CXX) $^ -o $@ $(LIBS) `root-config --cflags` -lFoam -lMinuit -lRooFit -lRooFitCore

runLifetimeFit: runLifetimeFit.cc
	$(CXX) $^ -o $@ $(LIBS) `root-config --cflags` -lFoam -lMinuit -lRooFit -lRooFitCore

runChiMassLifetimeFit: runChiMassLifetimeFit.cc
	$(CXX) $^ -o $@ $(LIBS) `root-config --cflags` -lFoam -lMinuit -lRooFit -lRooFitCore

runBkgHistos: runBkgHistos.cc 
	$(CXX) $^ -o $@ $(LIBS) `root-config --cflags` -lFoam -lMinuit -lRooFit -lRooFitCore

PlotFitPar: PlotFitPar.cc
	$(CXX) $^ -o $@ $(LIBS) `root-config --cflags` -lFoam -lMinuit -lRooFit -lRooFitCore

PlotCosThetaPhiBG: PlotCosThetaPhiBG.cc
	$(CXX) $^ -o $@ $(LIBS) `root-config --cflags` -lFoam -lMinuit -lRooFit -lRooFitCore

PlotCosThetaPhiDistribution: PlotCosThetaPhiDistribution.cc
	$(CXX) $^ -o $@ $(LIBS) `root-config --cflags` -lFoam -lMinuit -lRooFit -lRooFitCore

runDefineRegionsAndFractions: runDefineRegionsAndFractions.cc
	$(CXX) $^ -o $@ $(LIBS) `root-config --cflags` -lFoam -lMinuit -lRooFit -lRooFitCore




clean: 
	runData runChiData runWorkspace runMassFit runLifetimeFit runPlotMassLifetime runBkgHistos PlotFitPar PlotCosThetaPhiBG PlotCosThetaPhiDistribution runChiMassLifetimeFit runDefineRegionsAndFractions runPlotDataDistributions runPlotJpsiMassLifetime PlotJpsiFitPar

