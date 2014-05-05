#include "rootIncludes.inc"
#include "commonVar.h"
#include "RooUtils.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"

using namespace RooFit;

void massFit(const std::string &infilename, int rapBin, int ptBin, int nState, bool fitMassPR, bool fitMassNP, bool MC){

    TFile *infile = new TFile(infilename.c_str(), "UPDATE");
    if(!infile){
        std::cout << "Error: failed to open file with dataset" << std::endl;
    }
    RooWorkspace *ws=(RooWorkspace *)infile->Get("ws_masslifetime");
    if(!ws){
        std::cout << "Error: failed to open workspace " << std::endl;
    }

	//------------------------------------------------------------------------------------------------------------------
	// mass pdf
	std::cout << "Building mass pdf" << std::endl;

	//define mass shape
	//signal (prompt+noprompt) 

	RooRealVar CBsigma_p0_jpsi("CBsigma_p0_jpsi","CBsigma_p0_jpsi",0.02,0.015,0.025);
	RooRealVar CBsigma_p1_jpsi("CBsigma_p1_jpsi","CBsigma_p1_jpsi",0.,-0.01,0.01);
	RooRealVar CBsigma_p2_jpsi("CBsigma_p2_jpsi","CBsigma_p2_jpsi",0.01,0.005,0.025);
	ws->import(RooArgList(CBsigma_p0_jpsi,CBsigma_p1_jpsi,CBsigma_p2_jpsi));
	RooFormulaVar CBsigma_jpsi("CBsigma_jpsi","@0+@1*abs(@3)+@2*abs(@3)*abs(@3)",RooArgList(*ws->var("CBsigma_p0_jpsi"), *ws->var("CBsigma_p1_jpsi"), *ws->var("CBsigma_p2_jpsi"), *ws->var("JpsiRap")));
	ws->import(CBsigma_jpsi);

	RooRealVar CBalpha_p0_jpsi("CBalpha_p0_jpsi","CBalpha_p0_jpsi",1.729,1.2,2.5);
	RooRealVar CBalpha_p1_jpsi("CBalpha_p1_jpsi","CBalpha_p1_jpsi",0.191,0.,0.5);
	ws->import(RooArgList(CBalpha_p0_jpsi,CBalpha_p1_jpsi));
	RooFormulaVar CBalpha_jpsi("CBalpha_jpsi","@0+@1*abs(@2)",RooArgList(*ws->var("CBalpha_p0_jpsi"), *ws->var("CBalpha_p1_jpsi"), *ws->var("JpsiRap")));
	ws->import(CBalpha_jpsi);

	RooRealVar CBmass_p0_jpsi("CBmass_p0_jpsi","CBmass_p0_jpsi",3.094,3.086,3.098);
	RooRealVar CBmass_p1_jpsi("CBmass_p1_jpsi","CBmass_p1_jpsi",0.001,-0.002,0.002);
	RooRealVar CBmass_p2_jpsi("CBmass_p2_jpsi","CBmass_p2_jpsi",-0.003,-0.005,0.001);
	ws->import(RooArgList(CBmass_p0_jpsi,CBmass_p1_jpsi,CBmass_p2_jpsi));
	RooFormulaVar CBmass_jpsi("CBmass_jpsi","@0+@1*abs(@3)+@2*abs(@3)*abs(@3)",RooArgList(*ws->var("CBmass_p0_jpsi"), *ws->var("CBmass_p1_jpsi"), *ws->var("CBmass_p2_jpsi"), *ws->var("JpsiRap")));
	ws->import(CBmass_jpsi);


	ws->factory("RooCBShape::sigMassShape_jpsi(JpsiMass,CBmass_jpsi,CBsigma_jpsi,CBalpha_jpsi,CBn_jpsi[2.5])");
	ws->factory("Gaussian::gaussMassShape_jpsi(JpsiMass,CBmass_jpsi,CBsigma_jpsi)");

	//ws->factory("RooCBShape::massCBShape2_jpsi(JpsiMass,CBmass_jpsi,CBsigma2_jpsi[0.02,0.0001,0.1],CBalpha_jpsi,CBn_jpsi)");
	//ws->factory("SUM::sigMassShape_jpsi(fracCB1_jpsi[0.5,0.,1.]*massCBShape_jpsi, massCBShape2_jpsi)");

    RooRealVar fracCB1_jpsi("fracCB1_jpsi","fracCB1_jpsi",0.5,0.,1.);
    ws->import(fracCB1_jpsi);
    RooRealVar CBsigma2_jpsi("CBsigma2_jpsi","CBsigma2_jpsi",0.5,0.,1.);
    ws->import(CBsigma2_jpsi);

	//background
	ws->factory("RooExponential::bkgMassShape_jpsi(JpsiMass,bkgLambda_jpsi[0,-5,5])");
	//full Mass shape
	ws->factory("SUM::massModel_jpsi(fracBkg_jpsi[0.085,0.,1.]*bkgMassShape_jpsi, sigMassShape_jpsi)");

	ws->Print("v");

	cout<<"test1"<<endl;
    RooRealVar test1("test1","test1",1.);
    ws->import(test1);

	//------------------------------------------------------------------------------------------------------------------
	// do fit
	std::cout << "Fitting" << std::endl;

	//RooRealVar JpsiMass(*JpsiMass);
	//RooRealVar JpsiRap(*JpsiRap);

   // if(MC){
   //     ws->var("fracBkg_jpsi")->setVal(1.e-10);
   //     ws->var("fracBkg_jpsi")->setConstant(kTRUE);
   //     ws->var("bkgLambda_jpsi")->setVal(1.e-10);
   //     ws->var("bkgLambda_jpsi")->setConstant(kTRUE);
   // }





		//ws->var("CBsigma_p0_jpsi")->setMax(0.06);
        //ws->var("CBsigma_p1_jpsi")->setVal(0.);
        //ws->var("CBsigma_p1_jpsi")->setConstant(kTRUE);
        //ws->var("CBsigma_p2_jpsi")->setVal(0.);
        //ws->var("CBsigma_p2_jpsi")->setConstant(kTRUE);

		ws->var("CBmass_p1_jpsi")->setVal(0.);
		ws->var("CBmass_p1_jpsi")->setConstant(kTRUE);
		ws->var("CBmass_p2_jpsi")->setVal(0.);
		ws->var("CBmass_p2_jpsi")->setConstant(kTRUE);

		//ws->var("CBalpha_p1_jpsi")->setVal(0.);
		//ws->var("CBalpha_p1_jpsi")->setConstant(kTRUE);



		//ws->var("CBsigma")->setMax(0.07);
		//ws->var("CBsigma2")->setMax(0.07);
		//ws->var("CBsigma")->setVal(.02);
		//ws->var("CBsigma2")->setVal(.02);
		//ws->var("fracCB1")->setVal(0.5);

	//if(nState==4){
	//}
	//else if(nState==5){
	//	ws->var("CBmass")->setMax(3.75);
	//	ws->var("CBmass")->setMin(3.65);
	//	ws->var("CBmass")->setVal(3.69);
	//	ws->var("CBalpha")->setVal(2.);
	//	ws->var("CBn")->setVal(2.5);
	//	ws->var("CBn")->setConstant(kTRUE);
    //
	//	//if(rapBin == 1 && ptBin == 1){
	//	//  ws->var("CBsigma")->setMax(0.035);
	//	//  ws->var("CBsigma2")->setMax(0.035);
	//	//}
    //
	//	ws->var("CBsigma")->setMax(0.07);
	//	ws->var("CBsigma2")->setMax(0.07);
	//	if(rapBin == 3){
	//		ws->var("CBsigma")->setMin(0.03);
	//		ws->var("CBsigma2")->setMin(0.03);
	//	}
	//	if(rapBin == 3 && ptBin > 3)
	//		ws->var("fracCB1")->setMin(0.1);
    //
	//	ws->var("CBsigma")->setVal(.03);
	//	ws->var("CBsigma2")->setVal(.04);
	//	if(rapBin == 1 && ptBin == 1){
	//		ws->var("fracCB1")->setVal(0.5);
	//		ws->var("bkgLambda")->setVal(-1.5);
	//		ws->var("fracBkg")->setVal(0.1);
	//	}
	//}

		std::cout << "..." << std::endl;
	//RooAddPdf *MPdf = (RooAddPdf*)ws->pdf("massModel_jpsi");
	//MPdf->SetNameTitle("MPdf","MPdf");

	std::stringstream binName;
	binName << "jpsi_data_rap" << rapBin << "_pt" << ptBin;
	RooDataSet *data = (RooDataSet*)ws->data(binName.str().c_str());

	std::cout << "..." << std::endl;
	data->Print("v");

	RooArgSet *NLLs = new RooArgSet();
	RooAbsReal *MassNLL = NULL; 

	MassNLL = (RooAbsReal *)ws->pdf("massModel_jpsi")->createNLL(*data, ConditionalObservables(*ws->var("JpsiRap")), NumCPU(4));

	NLLs->add(*MassNLL);

	RooAddition *simNLL = new RooAddition("add","add",*NLLs);
	RooMinuit *mMinuit = new RooMinuit(*simNLL);

	mMinuit->setStrategy(1);
	mMinuit->setPrintEvalErrors(-1);
	mMinuit->setEvalErrorWall(false);
	mMinuit->setVerbose(false);

    mMinuit->migrad();

    std::stringstream rltNameMigrad;
    rltNameMigrad<< "fitresultaftermigrad_rap"<<rapBin<<"_pt"<<ptBin;
    RooFitResult *fitresultaftermigrad = (RooFitResult*)mMinuit->save(rltNameMigrad.str().c_str());
    double covQualMigrad=fitresultaftermigrad->covQual();
	cout<<"covQualMigrad = "<<covQualMigrad<<endl;

    if(covQualMigrad<3){
    	cout<<"-> repeat migrad() one more time..."<<endl;
        mMinuit->migrad();
        fitresultaftermigrad = (RooFitResult*)mMinuit->save(rltNameMigrad.str().c_str());
        double covQualMigradTen=10*covQualMigrad;
        covQualMigrad=fitresultaftermigrad->covQual();
    	cout<<"second turn migrad -> covQualMigrad = "<<covQualMigrad<<endl;
        covQualMigrad+=covQualMigradTen;
    	cout<<"-> MigradCovQualHistory = "<<covQualMigrad<<endl;
    }

    mMinuit->hesse();


    std::cout<< "End of fit" <<std::endl;
    std::stringstream rltName, snapshotName;
    rltName<< "m_fitresult_rap"<<rapBin<<"_pt"<<ptBin;
    snapshotName<< "m_snapshot_rap"<<rapBin<<"_pt"<<ptBin;
    RooFitResult *fitresult = (RooFitResult*)mMinuit->save(rltName.str().c_str());
    ws->import(*fitresult,kTRUE);
    fitresult->Print();
    double covQualHesse=fitresult->covQual();
	cout<<"covQualHesse = "<<covQualHesse<<endl;

    std::cout<< "save snapshot" <<std::endl;
    ws->saveSnapshot(snapshotName.str().c_str(), ws->allVars());

    RooRealVar var_covQualMigrad("var_jpsi_mass_covQualMigrad","var_jpsi_mass_covQualMigrad",covQualMigrad); if(!ws->var("var_jpsi_mass_covQualMigrad")) ws->import(var_covQualMigrad); else ws->var("var_jpsi_mass_covQualMigrad")->setVal(covQualMigrad);
    RooRealVar var_covQualHesse("var_jpsi_mass_covQualHesse","var_jpsi_mass_covQualHesse",covQualHesse); if(!ws->var("var_jpsi_mass_covQualHesse")) ws->import(var_covQualHesse); else ws->var("var_jpsi_mass_covQualHesse")->setVal(covQualHesse);




    //Define regions and calculate fractions

	RooRealVar *JpsiMass = ws->var("JpsiMass");
	assert( 0 != JpsiMass );
	RooRealVar *JpsiRap = ws->var("JpsiRap");
	assert( 0 != JpsiRap );

	RooAbsPdf *massPdf = ws->pdf("massModel_jpsi");
	assert ( 0 != massPdf );
	RooAbsPdf *bkgMassShape = ws->pdf("bkgMassShape_jpsi");
	assert ( 0 != bkgMassShape );
	RooAbsPdf *sigMassShape_jpsi = ws->pdf("sigMassShape_jpsi");
	assert ( 0 != sigMassShape_jpsi );
	RooAbsPdf *gaussMassShape_jpsi = ws->pdf("gaussMassShape_jpsi");
	assert ( 0 != gaussMassShape_jpsi );

	int nEntries = data->numEntries();


	RooDataSet *dataJpsiRap = (RooDataSet*)data->reduce(SelectVars(RooArgSet(*JpsiRap)),Name("dataJpsiRap"));

	int nevt=1000000;
	cout<<"generating datasets"<<endl;
	RooDataSet *GaussPseudoData = gaussMassShape_jpsi->generate(*JpsiMass,ProtoData(*dataJpsiRap),NumEvents(nevt/100.));
	RooDataSet *BackgroundPseudoData = bkgMassShape->generate(*JpsiMass,ProtoData(*dataJpsiRap),NumEvents(nevt/100.));
	RooDataSet *SignalPseudoData = sigMassShape_jpsi->generate(*JpsiMass,ProtoData(*dataJpsiRap),NumEvents(nevt));
	cout<<"finished generating datasets"<<endl;

	int nbinsHists=100;
	TH1F* sigMassShape_jpsi_asHist = new TH1F("sigMassShape_jpsi_asHist","sigMassShape_jpsi_asHist", nbinsHists, onia::massMin, onia::massMax);
	SignalPseudoData->fillHistogram(sigMassShape_jpsi_asHist,RooArgList(*JpsiMass));
	RooDataHist* signalMassShape_asRooDataHist = new RooDataHist("signalMassShape_asRooDataHist","signalMassShape_asRooDataHist", RooArgList(*JpsiMass), sigMassShape_jpsi_asHist);
	RooHistPdf* signalMassShape_asHistPdf = new RooHistPdf("signalMassShape_asHistPdf","signalMassShape_asHistPdf", RooArgSet(*JpsiMass), *signalMassShape_asRooDataHist, 3);

    ws->import(*signalMassShape_asHistPdf);

	ws->factory("SUM::FullmassPdf(fracBkg_jpsi*bkgMassShape_jpsi, signalMassShape_asHistPdf)");


    std::stringstream cutSR;
    std::stringstream cutLSB;
    std::stringstream cutRSB;
    cutSR << "TMath::Abs(JpsiMass-" << ws->var("CBmass_p0_jpsi")->getVal() << ") < " << onia::nSigMass << "*(" << ws->var("CBsigma_p0_jpsi")->getVal() << "+"<< ws->var("CBsigma_p1_jpsi")->getVal() << "*TMath::Abs(JpsiRap) + "<< ws->var("CBsigma_p2_jpsi")->getVal() << "*TMath::Abs(JpsiRap)*TMath::Abs(JpsiRap))";
    cutLSB << "JpsiMass < " << ws->var("CBmass_p0_jpsi")->getVal() << " - " << onia::nSigBkgLow << "*(" << ws->var("CBsigma_p0_jpsi")->getVal() << "+"<< ws->var("CBsigma_p1_jpsi")->getVal() << "*TMath::Abs(JpsiRap) + "<< ws->var("CBsigma_p2_jpsi")->getVal() << "*TMath::Abs(JpsiRap)*TMath::Abs(JpsiRap))";
    cutRSB << "JpsiMass > " << ws->var("CBmass_p0_jpsi")->getVal() << " + " << onia::nSigBkgHigh << "*(" << ws->var("CBsigma_p0_jpsi")->getVal() << "+"<< ws->var("CBsigma_p1_jpsi")->getVal() << "*TMath::Abs(JpsiRap) + "<< ws->var("CBsigma_p2_jpsi")->getVal() << "*TMath::Abs(JpsiRap)*TMath::Abs(JpsiRap))";

    cout<<cutSR.str().c_str()<<endl;
    cout<<cutLSB.str().c_str()<<endl;
    cout<<cutRSB.str().c_str()<<endl;

	std::stringstream binNameSR, binNameLSB, binNameRSB;
	binNameSR  << "data_rap" << rapBin << "_pt" << ptBin << "_SR";
	binNameLSB << "data_rap" << rapBin << "_pt" << ptBin << "_LSB";
	binNameRSB << "data_rap" << rapBin << "_pt" << ptBin << "_RSB";

	RooAbsData* dataSR, *dataLSB, *dataRSB;
	int events=0;

	dataSR	= data->reduce(Cut(cutSR.str().c_str()));
	dataLSB = data->reduce(Cut(cutLSB.str().c_str()));
	dataRSB = data->reduce(Cut(cutRSB.str().c_str()));

	dataSR->SetNameTitle(binNameSR.str().c_str(), "data in signal region");
	dataLSB->SetNameTitle(binNameLSB.str().c_str(), "data in LSB");
	dataRSB->SetNameTitle(binNameRSB.str().c_str(), "data in RSB");

	ws->import(*dataSR);
	ws->import(*dataLSB);
	ws->import(*dataRSB);


    RooAbsData* SignalPseudoDataSR = SignalPseudoData->reduce(Cut(cutSR.str().c_str()));
    RooAbsData* SignalPseudoDataLSB = SignalPseudoData->reduce(Cut(cutLSB.str().c_str()));
    RooAbsData* SignalPseudoDataRSB = SignalPseudoData->reduce(Cut(cutRSB.str().c_str()));

    SignalPseudoDataSR->Print();

	double fracSigEventsInSR= double(SignalPseudoDataSR->numEntries()) / double(SignalPseudoData->numEntries());
	cout<<"fracSigEventsInSR = "<<fracSigEventsInSR<<endl;
	double fracSigEventsInLSB= double(SignalPseudoDataLSB->numEntries()) / double(SignalPseudoData->numEntries());
	cout<<"fracSigEventsInLSB = "<<fracSigEventsInLSB<<endl;
	double fracSigEventsInRSB= double(SignalPseudoDataRSB->numEntries()) / double(SignalPseudoData->numEntries());
	cout<<"fracSigEventsInRSB = "<<fracSigEventsInRSB<<endl;

    RooAbsData* BackgroundPseudoDataSR = BackgroundPseudoData->reduce(Cut(cutSR.str().c_str()));
    RooAbsData* BackgroundPseudoDataLSB = BackgroundPseudoData->reduce(Cut(cutLSB.str().c_str()));
    RooAbsData* BackgroundPseudoDataRSB = BackgroundPseudoData->reduce(Cut(cutRSB.str().c_str()));

	double fracBGEventsInSR= double(BackgroundPseudoDataSR->numEntries()) / double(BackgroundPseudoData->numEntries());
	cout<<"fracBGEventsInSR = "<<fracBGEventsInSR<<endl;
	double fracBGEventsInLSB= double(BackgroundPseudoDataLSB->numEntries()) / double(BackgroundPseudoData->numEntries());
	cout<<"fracBGEventsInLSB = "<<fracBGEventsInLSB<<endl;
	double fracBGEventsInRSB= double(BackgroundPseudoDataRSB->numEntries()) / double(BackgroundPseudoData->numEntries());
	cout<<"fracBGEventsInRSB = "<<fracBGEventsInRSB<<endl;


	int nbinsSigmaDef=200;
	TH1F* hist1D = (TH1F*)GaussPseudoData->createHistogram("hist1D",*JpsiMass,Binning(nbinsSigmaDef));
	hist1D->Scale(1./hist1D->Integral());
	hist1D->SetLineColor(kRed);
	hist1D->SetMarkerColor(kRed);
	hist1D->GetXaxis()->SetLimits(-.2,.2);

	double massres = hist1D->GetRMS();
	double err_massres = hist1D->GetRMSError();



	double fracBG = ws->var("fracBkg_jpsi")->getVal();
	double fracBGErr = ws->var("fracBkg_jpsi")->getError();

	double relativeErr_fracBG=fracBGErr/fracBG;
	double relativeErr_fracSig=fracBGErr/(1.-fracBG);

	double n_SigInSR = nEntries*(1.-fracBG)*fracSigEventsInSR;
	double n_BGInSR = nEntries*(fracBG)*fracBGEventsInSR;

	double n_SigInLSB = nEntries*(1.-fracBG)*fracSigEventsInLSB;
	double n_BGInLSB = nEntries*(fracBG)*fracBGEventsInLSB;

	double n_SigInRSB = nEntries*(1.-fracBG)*fracSigEventsInRSB;
	double n_BGInRSB = nEntries*(fracBG)*fracBGEventsInRSB;

	double frac_jpsi_SigInSR = n_SigInSR / (n_SigInSR+n_BGInSR);
	double frac_jpsi_BGInSR = n_BGInSR / (n_SigInSR+n_BGInSR);

	double frac_jpsi_SigInLSB = n_SigInLSB / (n_SigInLSB+n_BGInLSB);
	double frac_jpsi_BGInLSB = n_BGInLSB / (n_SigInLSB+n_BGInLSB);

	double frac_jpsi_SigInRSB = n_SigInRSB / (n_SigInRSB+n_BGInRSB);
	double frac_jpsi_BGInRSB = n_BGInRSB / (n_SigInRSB+n_BGInRSB);

	cout<<"frac_jpsi_SigInSR = "<<frac_jpsi_SigInSR<<endl;
	cout<<"frac_jpsi_BGInSR = "<<frac_jpsi_BGInSR<<endl;

	cout<<"frac_jpsi_SigInLSB = "<<frac_jpsi_SigInLSB<<endl;
	cout<<"frac_jpsi_BGInLSB = "<<frac_jpsi_BGInLSB<<endl;

	cout<<"frac_jpsi_SigInRSB = "<<frac_jpsi_SigInRSB<<endl;
	cout<<"frac_jpsi_BGInRSB = "<<frac_jpsi_BGInRSB<<endl;

	cout<<"massres "<<massres<<endl;
	cout<<"err_massres "<<err_massres<<endl;

    RooRealVar var_frac_jpsi_SigInSR("var_frac_jpsi_SigInSR","var_frac_jpsi_SigInSR",frac_jpsi_SigInSR);  var_frac_jpsi_SigInSR.setError(frac_jpsi_SigInSR*relativeErr_fracSig); if(!ws->var("var_frac_jpsi_SigInSR")) ws->import(var_frac_jpsi_SigInSR); else {ws->var("var_frac_jpsi_SigInSR")->setVal(frac_jpsi_SigInSR); ws->var("var_frac_jpsi_SigInSR")->setError(frac_jpsi_SigInSR*relativeErr_fracSig);}
    RooRealVar var_frac_jpsi_BGInSR("var_frac_jpsi_BGInSR","var_frac_jpsi_BGInSR",frac_jpsi_BGInSR);  var_frac_jpsi_BGInSR.setError(frac_jpsi_BGInSR*relativeErr_fracBG); if(!ws->var("var_frac_jpsi_BGInSR")) ws->import(var_frac_jpsi_BGInSR); else {ws->var("var_frac_jpsi_BGInSR")->setVal(frac_jpsi_BGInSR); ws->var("var_frac_jpsi_BGInSR")->setError(frac_jpsi_BGInSR*relativeErr_fracBG);}
    RooRealVar var_frac_jpsi_SigInLSB("var_frac_jpsi_SigInLSB","var_frac_jpsi_SigInLSB",frac_jpsi_SigInLSB);  var_frac_jpsi_SigInLSB.setError(frac_jpsi_SigInLSB*relativeErr_fracSig); if(!ws->var("var_frac_jpsi_SigInLSB")) ws->import(var_frac_jpsi_SigInLSB); else {ws->var("var_frac_jpsi_SigInLSB")->setVal(frac_jpsi_SigInLSB); ws->var("var_frac_jpsi_SigInLSB")->setError(frac_jpsi_SigInLSB*relativeErr_fracSig);}
    RooRealVar var_frac_jpsi_BGInLSB("var_frac_jpsi_BGInLSB","var_frac_jpsi_BGInLSB",frac_jpsi_BGInLSB);  var_frac_jpsi_SigInLSB.setError(frac_jpsi_BGInLSB*relativeErr_fracBG); if(!ws->var("var_frac_jpsi_BGInLSB")) ws->import(var_frac_jpsi_BGInLSB); else {ws->var("var_frac_jpsi_BGInLSB")->setVal(frac_jpsi_BGInLSB); ws->var("var_frac_jpsi_BGInLSB")->setError(frac_jpsi_BGInLSB*relativeErr_fracBG);}
    RooRealVar var_frac_jpsi_SigInRSB("var_frac_jpsi_SigInRSB","var_frac_jpsi_SigInRSB",frac_jpsi_SigInRSB);  var_frac_jpsi_SigInRSB.setError(frac_jpsi_SigInRSB*relativeErr_fracSig); if(!ws->var("var_frac_jpsi_SigInRSB")) ws->import(var_frac_jpsi_SigInRSB); else {ws->var("var_frac_jpsi_SigInRSB")->setVal(frac_jpsi_SigInRSB); ws->var("var_frac_jpsi_SigInRSB")->setError(frac_jpsi_SigInRSB*relativeErr_fracSig);}
    RooRealVar var_frac_jpsi_BGInRSB("var_frac_jpsi_BGInRSB","var_frac_jpsi_BGInRSB",frac_jpsi_BGInRSB);  var_frac_jpsi_BGInRSB.setError(frac_jpsi_BGInRSB*relativeErr_fracBG); if(!ws->var("var_frac_jpsi_BGInRSB")) ws->import(var_frac_jpsi_BGInRSB); else {ws->var("var_frac_jpsi_BGInRSB")->setVal(frac_jpsi_BGInRSB); ws->var("var_frac_jpsi_BGInRSB")->setError(frac_jpsi_BGInRSB*relativeErr_fracBG);}

    RooRealVar var_massres("var_massres","var_massres",massres);  var_massres.setError(err_massres); if(!ws->var("var_massres")) ws->import(var_massres); else {ws->var("var_massres")->setVal(massres); ws->var("var_massres")->setError(err_massres);}


	ws->Write();
	infile->Close();
}
