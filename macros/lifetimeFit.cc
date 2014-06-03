#include "calculatePar.cc"

void buildLifetimePDF(RooWorkspace *ws);
void doFit(RooWorkspace *ws, int nState, double BkgRatio3Sig, double fracBkgInLSB, double fracBkgInRSB, int rapBin, int ptBin);

//=================================================================================
void lifetimeFit(const std::string &infilename, int rapBin, int ptBin, int nState){

	TFile* infile = new TFile(infilename.c_str(), "UPDATE");
	if(infile->IsZombie()){
		std::cout << "Error: failed to open mass root file" << std::endl;
		return;
	}
	std::string workspacename = "ws_masslifetime";
	RooWorkspace *ws=dynamic_cast<RooWorkspace*>(infile->Get(workspacename.c_str()));
	if(ws == 0){
		std::cout << "Error: failed to get workspace " << workspacename << " from file " << infilename << std::endl;
		infile->ls();
		return;
	}

	RooRealVar *JpsiMass = ws->var("JpsiMass");
	assert( 0 != JpsiMass );
	RooRealVar *JpsiRap = ws->var("JpsiRap");
	assert( 0 != JpsiRap );

	std::stringstream dataname;
	dataname << "data_rap" << rapBin << "_pt" << ptBin;
	RooAbsData* data = ws->data(dataname.str().c_str());

	std::stringstream binNameSR, binNameLSB, binNameRSB;
	binNameSR  << "data_rap" << rapBin << "_pt" << ptBin << "_SR";
	binNameLSB << "data_rap" << rapBin << "_pt" << ptBin << "_LSB";
	binNameRSB << "data_rap" << rapBin << "_pt" << ptBin << "_RSB";

	RooAbsData* dataSR = ws->data(binNameSR.str().c_str());
	RooAbsData* dataLSB = ws->data(binNameLSB.str().c_str());
	RooAbsData* dataRSB = ws->data(binNameRSB.str().c_str());


	std:: cout << "----------------------------" << "\n"
		<< "events in SR: " << dataSR->numEntries() << "\n"
		<< "events in LSB: " << dataLSB->numEntries() << "\n"
		<< "events in RSB: " << dataRSB->numEntries() << "\n"
		<< "----------------------------" << std::endl;

	// caculating median of different regions by filling events into histogram and getting the mean
	TH1* histSR =  dataSR->createHistogram("histSR", *JpsiMass,  Binning(120));
	TH1* histLSB = dataLSB->createHistogram("histLSB", *JpsiMass, Binning(120));
	TH1* histRSB = dataRSB->createHistogram("histRSB", *JpsiMass, Binning(120));

	double meanSR = histSR->GetMean();
	double meanLSB = histLSB->GetMean();
	double meanRSB = histRSB->GetMean();

	RooRealVar* MeanSR = new RooRealVar("MeanJpsiSR","MeanJpsiSR",meanSR);
	RooRealVar* MeanLSB = new RooRealVar("MeanJpsiLSB","MeanJpsiLSB",meanLSB);
	RooRealVar* MeanRSB = new RooRealVar("MeanJpsiRSB","MeanJpsiRSB",meanRSB);
	ws->import(RooArgList(*MeanSR, *MeanLSB, *MeanRSB));

	std::cout << "----------------------------" << "\n"
		<< "meanSR: " << meanSR << "\n"
		<< "meanLSB: " << meanLSB << "\n"
		<< "meanRSB: " << meanRSB
		<< "----------------------------" << std::endl;


	// building lifetime pdf
	std::cout << ">>>Building Mass and LifeTime PDF" << std::endl;
	buildLifetimePDF(ws);

	// fitting
	std::cout << ">>>Fitting" << std::endl;
	doFit(ws, nState, ws->var("var_frac_jpsi_BGInSR")->getVal(), ws->var("var_frac_jpsi_BGInLSB")->getVal(), ws->var("var_frac_jpsi_BGInRSB")->getVal(), rapBin, ptBin);

	std::cout << ">>>Writing results to root file" << std::endl;
	ws->Write();
	infile->Close();
}

//=================================================================================

void buildLifetimePDF(RooWorkspace *ws){


	RooRealVar ResDeformation("ResDeformation","ResDeformation",1.);
	RooFormulaVar JpsictErrDeformed("JpsictErrDeformed","@0*@1",RooArgList(*ws->var("JpsictErr"),ResDeformation));
	ws->import(JpsictErrDeformed);
	ws->import(ResDeformation);

	RooRealVar jpsi_ctResolution_p0("jpsi_ctResolution_p0","jpsi_ctResolution_p0",0.8,0.5,1.15);
	RooRealVar jpsi_ctResolution_p1("jpsi_ctResolution_p1","jpsi_ctResolution_p1",-5.);
	ws->import(RooArgList(jpsi_ctResolution_p0,jpsi_ctResolution_p1));
	RooFormulaVar jpsi_ctResolution_param("jpsi_ctResolution_param","abs(@0+(sign(-@2)+1)/2*@1*@2)",RooArgList(*ws->var("jpsi_ctResolution_p0"), *ws->var("jpsi_ctResolution_p1"), *ws->var("JpsictErr")));
	//RooFormulaVar jpsi_ctResolution_param("jpsi_ctResolution_param","@0",RooArgList(*ws->var("jpsi_ctResolution_p0")));
	ws->import(jpsi_ctResolution_param);

//	RooRealVar jpsi_ctResolution("jpsi_ctResolution","jpsi_ctResolution",0.8,0.5,1.15);
//	ws->import(jpsi_ctResolution);

	//----prompt
	//resolution function
	ws->factory("RooGaussModel::jpsi_promptLifetime(Jpsict,jpsi_promptMean[0.,-0.1, 0.1],jpsi_ctResolution[0.8,0.5,1.15],JpsictErr)");
	//ws->factory("RooGaussModel::jpsi_promptLifetime(Jpsict,jpsi_promptMean[0.],jpsi_ctResolution_param,JpsictErr)");
	((RooGaussModel*)ws->pdf("jpsi_promptLifetime"))->advertiseFlatScaleFactorIntegral(true);

	ws->factory("RooGaussModel::jpsi_promptLifetime2(Jpsict,jpsi_promptMean, jpsi_ctResolution2[1.3,1.15,2.5],JpsictErr)");
	((RooGaussModel*)ws->pdf("jpsi_promptLifetime2"))->advertiseFlatScaleFactorIntegral(true);

	ws->factory("RooGaussModel::jpsi_promptLifetime3(Jpsict,jpsi_promptMean, jpsi_ctResolution3[1.8,1.5,5.],JpsictErr)");
	((RooGaussModel*)ws->pdf("jpsi_promptLifetime3"))->advertiseFlatScaleFactorIntegral(true);

	RooGaussModel* jpsi_promptLifetime = (RooGaussModel*)ws->pdf("jpsi_promptLifetime");
	RooGaussModel* jpsi_promptLifetime2 = (RooGaussModel*)ws->pdf("jpsi_promptLifetime2");
	RooGaussModel* jpsi_promptLifetime3 = (RooGaussModel*)ws->pdf("jpsi_promptLifetime3");
	RooRealVar jpsi_fracGauss2("jpsi_fracGauss2","jpsi_fracGauss2",0.219,0.,0.45);
	RooRealVar jpsi_fracGauss3("jpsi_fracGauss3","jpsi_fracGauss3",0.1,0.,0.45);
	//RooAddModel jpsi_TotalPromptLifetime("jpsi_TotalPromptLifetime","jpsi_TotalPromptLifetime", RooArgList(*jpsi_promptLifetime3,*jpsi_promptLifetime2,*jpsi_promptLifetime),RooArgList(jpsi_fracGauss3,jpsi_fracGauss2));
	RooAddModel jpsi_TotalPromptLifetime("jpsi_TotalPromptLifetime","jpsi_TotalPromptLifetime", RooArgList(*jpsi_promptLifetime2,*jpsi_promptLifetime),RooArgList(jpsi_fracGauss2));
	ws->import(jpsi_fracGauss2);
	ws->import(jpsi_fracGauss3);
	ws->import(jpsi_TotalPromptLifetime);
	jpsi_TotalPromptLifetime.Print();

	//ws->factory("RooDecay::jpsi_promptSignalDSD(Jpsict,jpsi_signalTauDSD[.016,0.001,0.2],jpsi_TotalPromptLifetime,RooDecay::DoubleSided)");

	//RooAddModel jpsi_TotalPromptLifetime("jpsi_TotalPromptLifetime","jpsi_TotalPromptLifetime", RooArgList(*jpsi_promptGauss1,*jpsi_promptGauss2),RooArgList(jpsi_fracGauss2));

	//ws->factory("RooGaussModel::jpsi_promptGauss1(Jpsict,mean[0],width1[0.005,0.001,0.03])");
	//ws->factory("RooGaussModel::jpsi_promptGauss2(Jpsict,mean,width2[0.01,0.005,0.05])");
	//ws->factory("RooGaussModel::jpsi_promptGauss3(Jpsict,mean,width3[0.03,0.02,0.08])");
	//ws->factory("RooGaussModel::jpsi_promptGauss4(Jpsict,mean,width4[0.06,0.06,0.5])");
	//RooGaussModel* jpsi_promptGauss1 = (RooGaussModel*)ws->pdf("jpsi_promptGauss1");
	//RooGaussModel* jpsi_promptGauss2 = (RooGaussModel*)ws->pdf("jpsi_promptGauss2");
	//RooGaussModel* jpsi_promptGauss3 = (RooGaussModel*)ws->pdf("jpsi_promptGauss3");
	//RooGaussModel* jpsi_promptGauss4 = (RooGaussModel*)ws->pdf("jpsi_promptGauss4");
	//RooRealVar jpsi_fracPromptGauss1("jpsi_fracPromptGauss1","jpsi_fracPromptGauss1",0.1,0.01,0.7);
	//RooRealVar jpsi_fracPromptGauss2("jpsi_fracPromptGauss2","jpsi_fracPromptGauss2",0.25,0.01,0.7);
	//RooRealVar jpsi_fracPromptGauss3("jpsi_fracPromptGauss3","jpsi_fracPromptGauss3",0.6,0.01,0.7);
	//ws->import(jpsi_fracPromptGauss1);
	//ws->import(jpsi_fracPromptGauss2);
	//ws->import(jpsi_fracPromptGauss3);
	//RooAddModel jpsi_TotalPromptLifetime("jpsi_TotalPromptLifetime","jpsi_TotalPromptLifetime", RooArgList(*jpsi_promptGauss1,*jpsi_promptGauss2,*jpsi_promptGauss3,*jpsi_promptGauss4),RooArgList(jpsi_fracPromptGauss1,jpsi_fracPromptGauss2,jpsi_fracPromptGauss3));
	//ws->import(jpsi_TotalPromptLifetime);


	//----noprompt
	ws->factory("RooDecay::jpsi_nonPromptSSD(Jpsict,jpsi_nonPromptTau[.377,.15,0.6],jpsi_TotalPromptLifetime,RooDecay::SingleSided)");

	//----background
	//LSB
	ws->factory("RooDecay::jpsi_backgroundSSD(Jpsict,jpsi_bkgTauSSD[.355,0.1,0.6],jpsi_TotalPromptLifetime,RooDecay::SingleSided)");
	ws->factory("RooDecay::jpsi_backgroundFD(Jpsict,jpsi_bkgTauFD[.11,.01,0.3],jpsi_TotalPromptLifetime,RooDecay::Flipped)");
	ws->factory("RooDecay::jpsi_backgroundDSD(Jpsict,jpsi_bkgTauDSD[.016,0.001,0.05],jpsi_TotalPromptLifetime,RooDecay::DoubleSided)");
	//ws->factory("SUM::jpsi_backgroundlifetimeLpre(jpsi_fBkgSSDR_LSB[.4,0,1.]*jpsi_backgroundSSD,jpsi_fBkgDSD_LSB[.2,0,1.]*jpsi_backgroundDSD,jpsi_backgroundFD)");
	ws->factory("SUM::jpsi_backgroundlifetimeLpre(jpsi_fBkgSSDR[.7,0.45,0.95]*jpsi_backgroundSSD,jpsi_fBkgDSD[.25,0.05,0.35]*jpsi_backgroundDSD,jpsi_backgroundFD)");

	RooFormulaVar jpsi_fBkgSSDL("jpsi_fBkgSSDL","1-@0-@1",RooArgList(*ws->var("jpsi_fBkgSSDR"),*ws->var("jpsi_fBkgDSD")));
	ws->import(jpsi_fBkgSSDL);

	////RSB
	//ws->factory("RooDecay::jpsi_backgroundSSD_RSB(Jpsict,jpsi_bkgTauSSD_RSB[.4,0,3],jpsi_TotalPromptLifetime,RooDecay::SingleSided)");
	//ws->factory("SUM::jpsi_backgroundlifetimeRpre(jpsi_fBkgSSDR_RSB[.4,0,1.]*jpsi_backgroundSSD,jpsi_fBkgDSD_RSB[.2,0,1.]*jpsi_backgroundDSD,jpsi_backgroundFD)");
	ws->factory("SUM::jpsi_backgroundlifetimeRpre(jpsi_fBkgSSDR*jpsi_backgroundSSD,jpsi_fBkgDSD*jpsi_backgroundDSD,jpsi_backgroundFD)");

	//Signal region
	//interpolation
	//ws->factory("expr::jpsi_fBkgSSDR('@0+(@2-@3)*(@1-@0)/(@4-@3)',jpsi_fBkgSSDR_LSB,jpsi_fBkgSSDR_RSB,MeanSR,MeanLSB,MeanRSB)");
	//ws->factory("expr::jpsi_fBkgDSD('@0+(@2-@3)*(@1-@0)/(@4-@3)',jpsi_fBkgDSD_LSB,jpsi_fBkgDSD_RSB,MeanSR,MeanLSB,MeanRSB)");
	//ws->factory("expr::bkgTauSSD('@0+(@2-@3)*(@1-@0)/(@4-@3)',bkgTauSSD_LSB,bkgTauSSD_RSB,MeanSR,MeanLSB,MeanRSB)");

	//ws->factory("RooDecay::jpsi_backgroundSSD(Jpsict,jpsi_bkgTauSSD,jpsi_TotalPromptLifetime,RooDecay::SingleSided)");
	ws->factory("SUM::jpsi_backgroundlifetime(jpsi_fBkgSSDR*jpsi_backgroundSSD,jpsi_fBkgDSD*jpsi_backgroundDSD,jpsi_backgroundFD)");

	//---final pdf
	//signal region
	ws->factory("SUM::jpsi_fulllifetimeSRnoC(jpsi_fBkg[0.05,0.,1.]*jpsi_backgroundlifetime,jpsi_fPrompt[0.55,0.,1.]*jpsi_TotalPromptLifetime,jpsi_nonPromptSSD)");
	//ws->factory("SUM::jpsi_fulllifetimeSRnoC(jpsi_fBkg[0.05,0.,1.]*jpsi_backgroundlifetime,jpsi_fPrompt[0.55,0.,1.]*jpsi_promptSignalDSD,jpsi_nonPromptSSD)");

	RooFormulaVar jpsi_fNonPrompt("jpsi_fNonPrompt","1-@0-@1",RooArgList(*ws->var("jpsi_fBkg"),*ws->var("jpsi_fPrompt")));
	ws->import(jpsi_fNonPrompt);

	RooRealVar jpsi_FracBkg("jpsi_FracBkg","jpsi_FracBkg",ws->var("var_frac_jpsi_BGInSR")->getVal());
	ws->import(jpsi_FracBkg);
	RooRealVar jpsi_FracBkgErr("jpsi_FracBkgErr","jpsi_FracBkgErr",ws->var("var_frac_jpsi_BGInSR")->getError());
	ws->import(jpsi_FracBkgErr);
	ws->factory("RooGaussian::jpsi_fBkgConstraint(jpsi_fBkg, jpsi_FracBkg, jpsi_FracBkgErr)");
	ws->factory("PROD::jpsi_fulllifetimeSR(jpsi_fulllifetimeSRnoC, jpsi_fBkgConstraint)");

	RooRealVar* jpsi_fBkgLSB = new RooRealVar("jpsi_fBkgLSB","jpsi_fBkgLSB",0.8,0.,1.);
	RooRealVar* jpsi_fBkgRSB = new RooRealVar("jpsi_fBkgRSB","jpsi_fBkgRSB",0.8,0.,1.);
	ws->import(*jpsi_fBkgLSB); ws->import(*jpsi_fBkgRSB);

	//f_P / f_NP should be same in LSB,SR,RSB, to interpolate f_P in L(R)SB from SR
	ws->factory("expr::jpsi_fPromptLSB('@0*(1.-@1)/(1.-@2)',jpsi_fPrompt,jpsi_fBkgLSB,jpsi_fBkg)");
	ws->factory("expr::jpsi_fPromptRSB('@0*(1.-@1)/(1.-@2)',jpsi_fPrompt,jpsi_fBkgRSB,jpsi_fBkg)");

	RooFormulaVar jpsi_fNonPromptLSB("jpsi_fNonPromptLSB","1-@0-@1",RooArgList(*ws->function("jpsi_fPromptLSB"),*ws->var("jpsi_fBkgLSB")));
	ws->import(jpsi_fNonPromptLSB);
	RooFormulaVar jpsi_fNonPromptRSB("jpsi_fNonPromptRSB","1-@0-@1",RooArgList(*ws->function("jpsi_fPromptRSB"),*ws->var("jpsi_fBkgRSB")));
	ws->import(jpsi_fNonPromptRSB);

	//LSB
	ws->factory("SUM::jpsi_backgroundlifetimeLnoC(jpsi_fBkgLSB*jpsi_backgroundlifetimeLpre,jpsi_fPromptLSB*jpsi_TotalPromptLifetime,jpsi_nonPromptSSD)");
	//ws->factory("SUM::jpsi_backgroundlifetimeLnoC(jpsi_fBkgLSB*jpsi_backgroundlifetimeLpre,jpsi_fPromptLSB*jpsi_promptSignalDSD,jpsi_nonPromptSSD)");
	RooRealVar jpsi_FracBkgLSB("jpsi_FracBkgLSB","jpsi_FracBkgLSB",ws->var("var_frac_jpsi_BGInLSB")->getVal());
	ws->import(jpsi_FracBkgLSB);
	RooRealVar jpsi_FracBkgLSBErr("jpsi_FracBkgLSBErr","jpsi_FracBkgLSBErr",ws->var("var_frac_jpsi_BGInLSB")->getError());
	ws->import(jpsi_FracBkgLSBErr);
	ws->factory("RooGaussian::jpsi_fBkgLSBConstraint(jpsi_fBkgLSB, jpsi_FracBkgLSB, jpsi_FracBkgLSBErr)");
	ws->factory("PROD::jpsi_backgroundlifetimeL(jpsi_backgroundlifetimeLnoC, jpsi_fBkgLSBConstraint)");

	//RSB
	ws->factory("SUM::jpsi_backgroundlifetimeRnoC(jpsi_fBkgRSB*jpsi_backgroundlifetimeRpre,jpsi_fPromptRSB*jpsi_TotalPromptLifetime,jpsi_nonPromptSSD)");
	//ws->factory("SUM::jpsi_backgroundlifetimeRnoC(jpsi_fBkgRSB*jpsi_backgroundlifetimeRpre,jpsi_fPromptRSB*jpsi_promptSignalDSD,jpsi_nonPromptSSD)");
	RooRealVar jpsi_FracBkgRSB("jpsi_FracBkgRSB","jpsi_FracBkgRSB",ws->var("var_frac_jpsi_BGInRSB")->getVal());
	ws->import(jpsi_FracBkgRSB);
	RooRealVar jpsi_FracBkgRSBErr("jpsi_FracBkgRSBErr","jpsi_FracBkgRSBErr",ws->var("var_frac_jpsi_BGInRSB")->getError());
	ws->import(jpsi_FracBkgRSBErr);
	ws->factory("RooGaussian::jpsi_fBkgRSBConstraint(jpsi_fBkgRSB, jpsi_FracBkgRSB, jpsi_FracBkgRSBErr)");
	ws->factory("PROD::jpsi_backgroundlifetimeR(jpsi_backgroundlifetimeRnoC, jpsi_fBkgRSBConstraint)");

	ws->Print("v");
}


//=================================================================================
void doFit(RooWorkspace *ws, int nState, double BkgRatio3Sig, double fracBkgInLSB, double fracBkgInRSB, int rapBin, int ptBin){

	RooRealVar JpsiMass(*ws->var("JpsiMass"));
	RooRealVar Jpsict(*ws->var("Jpsict"));

	stringstream binNameSR, binNameLSB, binNameRSB;
	binNameSR  << "data_rap" << rapBin << "_pt" << ptBin <<"_SR";
	binNameLSB << "data_rap" << rapBin << "_pt" << ptBin <<"_LSB";
	binNameRSB << "data_rap" << rapBin << "_pt" << ptBin <<"_RSB";
	RooDataSet *dataSR = (RooDataSet*)ws->data(binNameSR.str().c_str());
	RooDataSet *dataLSB = (RooDataSet*)ws->data(binNameLSB.str().c_str());
	RooDataSet *dataRSB = (RooDataSet*)ws->data(binNameRSB.str().c_str());

	//starting parameter of constrained parameters
	ws->var("jpsi_fBkg")->setVal(BkgRatio3Sig);
	ws->var("jpsi_fBkgLSB")->setVal(fracBkgInLSB);
	ws->var("jpsi_fBkgRSB")->setVal(fracBkgInRSB);
	ws->var("jpsi_fBkgLSB")->setConstant(kTRUE);
	ws->var("jpsi_fBkgRSB")->setConstant(kTRUE);

	//parameter juggling of free parameters

	//ws->var("jpsi_fBkg")->setVal(0.);
	//ws->var("jpsi_fPrompt")->setVal(0.999999);
	//ws->var("jpsi_fBkgLSB")->setVal(0.);
	//ws->var("jpsi_fBkgRSB")->setVal(0.);
	//ws->var("jpsi_ctResolution2")->setMax(2.5);

	//ws->var("jpsi_fBkg")->setConstant(kTRUE);
	//ws->var("jpsi_fPrompt")->setConstant(kTRUE);
	//ws->var("jpsi_fBkgSSDR")->setConstant(kTRUE);
	//ws->var("jpsi_fBkgDSD")->setConstant(kTRUE);
	//ws->var("jpsi_bkgTauSSD")->setConstant(kTRUE);
	//ws->var("jpsi_bkgTauFD")->setConstant(kTRUE);
	//ws->var("jpsi_bkgTauDSD")->setConstant(kTRUE);
	//ws->var("jpsi_nonPromptTau")->setConstant(kTRUE);

	//ws->var("jpsi_ctResolution")->setVal(0.848*1.0525);
	//ws->var("jpsi_ctResolution")->setConstant(kTRUE);
	//ws->var("jpsi_fracGauss2")->setConstant(kTRUE);

	//ws->var("jpsi_ctResolution2")->setVal(1.55);
	//ws->var("jpsi_ctResolution2")->setConstant(kTRUE);

	//RooAbsPdf *ModelLifeSR = (RooAbsPdf*)ws->pdf("jpsi_fulllifetimeSR");
	//RooAbsPdf *ModelLifeLSB = (RooAbsPdf*)ws->pdf("jpsi_backgroundlifetimeL");
	//RooAbsPdf *ModelLifeRSB = (RooAbsPdf*)ws->pdf("jpsi_backgroundlifetimeR");

	RooAbsPdf *ModelLifeSR = (RooAbsPdf*)ws->pdf("jpsi_fulllifetimeSR");
	RooAbsPdf *ModelLifeLSB = (RooAbsPdf*)ws->pdf("jpsi_backgroundlifetimeLnoC");
	RooAbsPdf *ModelLifeRSB = (RooAbsPdf*)ws->pdf("jpsi_backgroundlifetimeRnoC");



	RooArgSet *NLLs = new RooArgSet();
	RooAbsReal *MLNLLSR = NULL, *MLNLLLSB = NULL, *MLNLLRSB = NULL;

	MLNLLSR = (RooAbsReal *)ModelLifeSR->createNLL(*dataSR,
			//ConditionalObservables(RooArgSet(*ws->function("JpsictErrDeformed"))),
			ConditionalObservables(RooArgSet(*ws->var("JpsictErr"))),
			Constrain(RooArgSet(*ws->var("jpsi_fBkg"))),
			Extended(kFALSE),
			NumCPU(6));
	MLNLLLSB = (RooAbsReal *)ModelLifeLSB->createNLL(*dataLSB,
			//ConditionalObservables(RooArgSet(*ws->function("JpsictErrDeformed"))),
			ConditionalObservables(RooArgSet(*ws->var("JpsictErr"))),
			//Constrain(RooArgSet(*ws->var("jpsi_fBkgLSB"))),
			Extended(kFALSE),
			NumCPU(6));
	MLNLLRSB = (RooAbsReal *)ModelLifeRSB->createNLL(*dataRSB,
			//ConditionalObservables(RooArgSet(*ws->function("JpsictErrDeformed"))),
			ConditionalObservables(RooArgSet(*ws->var("JpsictErr"))),
			//Constrain(RooArgSet(*ws->var("jpsi_fBkgRSB"))),
			Extended(kFALSE),
			NumCPU(6));

	cout<<"SR NLL before fit = "<<MLNLLSR->getVal()<<endl;

	NLLs->add(*MLNLLSR);
	NLLs->add(*MLNLLLSB);
	NLLs->add(*MLNLLRSB);

	RooAddition *simNLL = new RooAddition("add","add",*NLLs);
	RooMinuit *lMinuit = new RooMinuit(*simNLL);
	//RooMinuit *lMinuit = new RooMinuit(*MLNLLSR);

	lMinuit->setStrategy(1); 
	//if(nState==4 && ptBin > 7) lMinuit->setStrategy(2);
	lMinuit->setPrintEvalErrors(-1);
	lMinuit->setEvalErrorWall(false);
	lMinuit->setVerbose(false);
	lMinuit->setPrintLevel(-1);

	lMinuit->migrad();

    std::stringstream rltNameMigrad;
    rltNameMigrad<< "fitresultaftermigrad_rap"<<rapBin<<"_pt"<<ptBin;
    RooFitResult *fitresultaftermigrad = (RooFitResult*)lMinuit->save(rltNameMigrad.str().c_str());
    double covQualMigrad=fitresultaftermigrad->covQual();
	cout<<"covQualMigrad = "<<covQualMigrad<<endl;

    if(covQualMigrad<3){
    	cout<<"-> repeat migrad() one more time..."<<endl;
    	lMinuit->migrad();
        fitresultaftermigrad = (RooFitResult*)lMinuit->save(rltNameMigrad.str().c_str());
        double covQualMigradTen=10*covQualMigrad;
        covQualMigrad=fitresultaftermigrad->covQual();
    	cout<<"second turn migrad -> covQualMigrad = "<<covQualMigrad<<endl;
        covQualMigrad+=covQualMigradTen;
    	cout<<"-> MigradCovQualHistory = "<<covQualMigrad<<endl;
    }

    lMinuit->hesse();


    std::cout<< "End of fit" <<std::endl;
    std::stringstream rltName, snapshotName;
    rltName<< "jpsi_l_fitresult_rap"<<rapBin<<"_pt"<<ptBin;
    snapshotName<< "jpsi_l_snapshot_rap"<<rapBin<<"_pt"<<ptBin;
    RooFitResult *fitresult = (RooFitResult*)lMinuit->save(rltName.str().c_str());
    ws->import(*fitresult,kTRUE);
    fitresult->Print();
    double covQualHesse=fitresult->covQual();
	cout<<"covQualHesse = "<<covQualHesse<<endl;


	cout<<"SR NLL after fit = "<<MLNLLSR->getVal()<<endl;
	cout<<"LSB NLL after fit = "<<MLNLLLSB->getVal()<<endl;
	cout<<"RSB NLL after fit = "<<MLNLLRSB->getVal()<<endl;

	//ws->var("jpsi_ctResolution")->setVal(0.89);
	//ws->var("jpsi_ctResolution2")->setVal(1.5);

	cout<<"SR NLL after manually changing resolutions = "<<MLNLLSR->getVal()<<endl;
	cout<<"LSB NLL after manually changing resolutions = "<<MLNLLLSB->getVal()<<endl;
	cout<<"RSB NLL after manually changing resolutions = "<<MLNLLRSB->getVal()<<endl;


	//double BkgTauFD = ((RooRealVar*)ws->var("bkgTauFD"))->getVal();
	//double BkgTauDSD = ((RooRealVar*)ws->var("bkgTauDSD"))->getVal();
	//double CtResolution = ((RooRealVar*)ws->var("ctResolution"))->getVal();
	//double CtResolution2 = ((RooRealVar*)ws->var("ctResolution2"))->getVal();
    //
	//int refitCount = 0;
	//while(BkgTauFD < BkgTauDSD){
	//	if(refitCount>3) break;
	//	cout<<"bkgTauFD < bkgTauDSD"<<endl;
	//	ws->var("bkgTauFD")->setVal(BkgTauDSD);
	//	ws->var("bkgTauDSD")->setVal(BkgTauFD);
    //
	//	lMinuit->migrad();
	//	fitresult = (RooFitResult*)lMinuit->save(rltName.str().c_str());
	//	cout<<"fitresult->covQual(): "<<fitresult->covQual()<<endl;
	//	lMinuit->migrad();
	//	fitresult = (RooFitResult*)lMinuit->save(rltName.str().c_str());
	//	cout<<"fitresult->covQual(): "<<fitresult->covQual()<<endl;
	//	lMinuit->hesse();
	//	fitresult = (RooFitResult*)lMinuit->save(rltName.str().c_str());
	//	cout<<"fitresult->covQual(): "<<fitresult->covQual()<<endl;
    //
	//	BkgTauFD = ((RooRealVar*)ws->var("bkgTauFD"))->getVal();
	//	BkgTauDSD = ((RooRealVar*)ws->var("bkgTauDSD"))->getVal();
    //
	//	refitCount++;
	//}
    //
	//while(CtResolution > CtResolution2){
	//	if(refitCount>3) break;
	//	cout<<"CtResolution > CtResolution2"<<endl;
	//	ws->var("ctResolution")->setVal(CtResolution2);
	//	ws->var("ctResolution2")->setVal(CtResolution);
    //
	//	lMinuit->migrad();
	//	fitresult = (RooFitResult*)lMinuit->save(rltName.str().c_str());
	//	cout<<"fitresult->covQual(): "<<fitresult->covQual()<<endl;
	//	lMinuit->migrad();
	//	fitresult = (RooFitResult*)lMinuit->save(rltName.str().c_str());
	//	cout<<"fitresult->covQual(): "<<fitresult->covQual()<<endl;
	//	lMinuit->hesse();
	//	fitresult = (RooFitResult*)lMinuit->save(rltName.str().c_str());
	//	cout<<"fitresult->covQual(): "<<fitresult->covQual()<<endl;
    //
	//	CtResolution = ((RooRealVar*)ws->var("ctResolution"))->getVal();
	//	CtResolution2 = ((RooRealVar*)ws->var("ctResolution2"))->getVal();
    //
	//	refitCount++;
	//}
    //
	//cout<<">>>>>>Refitted "<<refitCount<<" times"<<endl;
	//if(refitCount>3) cout<<">>>>>>Fit not converged"<<endl;

	fitresult->Print();
	ws->import(*fitresult,kTRUE);

	ws->saveSnapshot(snapshotName.str().c_str(), ws->allVars());
}

