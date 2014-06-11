/*
 * chiMassLifetimeFit.cc
 *
 *  Created on: Jan 23, 2014
 *      Author: valentinknuenz
 */

#include "rootIncludes.inc"
#include "commonVar.h"
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooAbsPdf.h"
#include "RooFFTConvPdf.h"
#include "RooBreitWigner.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooCBShape.h"


using namespace RooFit;

void chiMassLifetimeFit(const std::string &infilename, int rapBin, int ptBin, int nState, bool runChiMassFitOnly, bool MC){
	cout<<"chiMassLifetimeFit"<<endl;

    TFile *infile = new TFile(infilename.c_str(), "UPDATE");
    if(!infile){
        std::cout << "Error: failed to open file with dataset" << std::endl;
    }
    RooWorkspace *ws=(RooWorkspace *)infile->Get("ws_masslifetime");
    if(!ws){
        std::cout << "Error: failed to open workspace " << std::endl;
    }

    std::stringstream binName;
    binName << "data_rap" << rapBin << "_pt" << ptBin<< "_SR";
    RooDataSet *data = (RooDataSet*)ws->data(binName.str().c_str());
    double ev = data->numEntries();
    std::cout<<"Number of Events in dataset: "<<ev<<endl;

    double ExtensionFactor=2.; //Allow total number of events to float from ev/ExtensionFactor to ev*ExtensionFactor

    bool ConstrainMassDifferenceWithPESandSigmaWithScale=true;

    if(MC) ConstrainMassDifferenceWithPESandSigmaWithScale=false;

    //------------------------------------------------------------------------------------------------------------------
    // mass pdf
    std::cout << "Building mass pdf" << std::endl;

    //define mass shape

    if(ConstrainMassDifferenceWithPESandSigmaWithScale){

    	ws->factory("RooCBShape::M_chic1(chicMass,CBmass1[3.5,3.45,3.54],CBsigma1[0.008,0.003,0.02],CBalpha1[0.6,.2,1.1],CBn[2.5,1.8,3.2])");

		RooFormulaVar PES("PES",Form("(@0-%f)/%f",onia::MpsiPDG, onia::Mchi1PDG-onia::MpsiPDG),RooArgList(*ws->var("CBmass1")));
		ws->import(PES);

		RooFormulaVar CBmass0("CBmass0",Form("@0*%f+%f",onia::Mchi0PDG-onia::MpsiPDG, onia::MpsiPDG),RooArgList(*ws->function("PES")));
		RooFormulaVar CBsigma0("CBsigma0",Form("@0*%f",(onia::Mchi0PDG-onia::MpsiPDG)/(onia::Mchi1PDG-onia::MpsiPDG)),RooArgList(*ws->var("CBsigma1")));
		ws->import(CBmass0);
		ws->import(CBsigma0);
		RooFormulaVar CBmass2("CBmass2",Form("@0*%f+%f",onia::Mchi2PDG-onia::MpsiPDG, onia::MpsiPDG),RooArgList(*ws->function("PES")));
		RooFormulaVar CBsigma2("CBsigma2",Form("@0*%f",(onia::Mchi2PDG-onia::MpsiPDG)/(onia::Mchi1PDG-onia::MpsiPDG)),RooArgList(*ws->var("CBsigma1")));
		ws->import(CBmass2);
		ws->import(CBsigma2);

		RooFormulaVar CBn2("CBn2",Form("@0"),RooArgList(*ws->var("CBn")));
		ws->import(CBn2);
		ws->factory("RooCBShape::M_chic2(chicMass,CBmass2,CBsigma2,CBalpha2[0.6,.2,1.1],CBn2)");

		RooFormulaVar CBalpha0("CBalpha0","(@0+@1)/2.",RooArgList(*ws->var("CBalpha1"),*ws->var("CBalpha2")));
		ws->import(CBalpha0);

		ws->factory("RooVoigtian::M_chic0(chicMass,CBmass0,CBsigma0, CBwidth0[0.0104])");

    }
    else{

    	ws->factory("RooCBShape::M_chic1(chicMass,CBmass1[3.5,3.45,3.54],CBsigma1[0.008,0.003,0.02],CBalpha1[0.6,.2,1.1],CBn[2.5,1.8,3.2])");

		RooFormulaVar PES("PES",Form("(@0-%f)/%f",onia::MpsiPDG, onia::Mchi1PDG-onia::MpsiPDG),RooArgList(*ws->var("CBmass1")));
		ws->import(PES);

		RooFormulaVar CBmass0("CBmass0",Form("@0*%f+%f",onia::Mchi0PDG-onia::MpsiPDG, onia::MpsiPDG),RooArgList(*ws->function("PES")));
		RooFormulaVar CBsigma0("CBsigma0",Form("@0*%f",(onia::Mchi0PDG-onia::MpsiPDG)/(onia::Mchi1PDG-onia::MpsiPDG)),RooArgList(*ws->var("CBsigma1")));
		ws->import(CBmass0);
		ws->import(CBsigma0);

		ws->factory("RooCBShape::M_chic2(chicMass,CBmass2[3.54,3.49,3.58],CBsigma2[0.008,0.003,0.02],CBalpha2[0.6,.2,1.1],CBn2[2.5,1.8,3.2])");

		RooFormulaVar CBalpha0("CBalpha0","(@0+@1)/2.",RooArgList(*ws->var("CBalpha1"),*ws->var("CBalpha2")));
		ws->import(CBalpha0);

		ws->factory("RooVoigtian::M_chic0(chicMass,CBmass0,CBsigma0, CBwidth0[0.0104])");

    }


    //background

    RooAbsPdf* M_background;
    double q01S_Start=3.1;
    RooRealVar q01S("q01S","q01S",q01S_Start,3.0,3.2);

    RooRealVar alpha1("alpha1","alpha1",0.6,0.,5.);
    RooRealVar beta1("beta1","beta1",-2.5,-10.,0.);
    RooFormulaVar a1("a1","TMath::Abs(@0-@1)",RooArgList(*ws->var("chicMass"),q01S));
    RooFormulaVar b1("b1","@0*(@1-@2)",RooArgList(beta1,*ws->var("chicMass"),q01S));
    RooFormulaVar signum1("signum1","(TMath::Sign(-1.,@0-@1)+1)/2.",RooArgList(*ws->var("chicMass"),q01S));
    ws->import(alpha1);
    ws->import(beta1);
    ws->import(q01S);
    ws->var("alpha1")->setConstant(kTRUE);
    ws->var("beta1")->setConstant(kTRUE);
    ws->var("q01S")->setConstant(kTRUE);

    //M_background = new RooGenericPdf("M_background","signum1*pow(a1,alpha1)*exp(b1)",RooArgSet(signum1,a1,alpha1,b1));

    RooRealVar BK_p1("BK_p1","BK_p1",0.,-1.,1.);
    RooRealVar BK_p2("BK_p2","BK_p2",0.,-.1,.1);

    M_background = new RooPolynomial("M_background","M_background",*ws->var("chicMass"),RooArgList(BK_p1,BK_p2));

    ws->import(BK_p1);
    ws->import(BK_p2);
    ws->import(*M_background);

    //ws->factory("RooPolynomial::M_background(chicMass,poly1[0., -10., 10.])");



    //lifetime pdf
    std::cout << "Building lifetime pdf" << std::endl;

    //define lifetime shape
    //prompt
    //resolution function
    ws->factory("RooGaussModel::promptLifetime(Jpsict,promptMean[0,-.01,.01],ctResolution[0.8,.5,1.], 1, JpsictErr)");
    ((RooGaussModel*)ws->pdf("promptLifetime"))->advertiseFlatScaleFactorIntegral(true);
    ws->factory("RooGaussModel::promptLifetime2(Jpsict,promptMean, ctResolution2[1.5,1.,2.5], 1, JpsictErr)");
    ((RooGaussModel*)ws->pdf("promptLifetime2"))->advertiseFlatScaleFactorIntegral(true);
    RooGaussModel* promptLifetime = (RooGaussModel*)ws->pdf("promptLifetime");
    RooGaussModel* promptLifetime2 = (RooGaussModel*)ws->pdf("promptLifetime2");
    RooRealVar fracGauss2("fracGauss2","fracGauss2",0.2,0.1,0.45);
    RooAddModel L_TotalPromptLifetime_punzi("L_TotalPromptLifetime_punzi","L_TotalPromptLifetime_punzi",
                                    RooArgSet(*promptLifetime2,*promptLifetime),fracGauss2);
    ws->import(fracGauss2);
    ws->import(L_TotalPromptLifetime_punzi);
    L_TotalPromptLifetime_punzi.Print();
	ws->factory("PROD::L_TotalPromptLifetime(L_TotalPromptLifetime_punzi|JpsictErr, pdf_ctauerrModelPR)");

    //nonprompt signal
    ws->factory("RooDecay::L_chic1_NP_punzi(Jpsict,NP_TauChic[.3,.25,0.5],L_TotalPromptLifetime_punzi,RooDecay::SingleSided)");
    ws->factory("RooDecay::L_chic2_NP_punzi(Jpsict,NP_TauChic,L_TotalPromptLifetime_punzi,RooDecay::SingleSided)");
    ws->factory("RooDecay::L_chic0_NP_punzi(Jpsict,NP_TauChic,L_TotalPromptLifetime_punzi,RooDecay::SingleSided)");
	ws->factory("PROD::L_chic1_NP(L_chic1_NP_punzi|JpsictErr, pdf_ctauerrModelNP)");
	ws->factory("PROD::L_chic2_NP(L_chic1_NP_punzi|JpsictErr, pdf_ctauerrModelNP)");
	ws->factory("PROD::L_chic0_NP(L_chic1_NP_punzi|JpsictErr, pdf_ctauerrModelNP)");
    //Jpsi background
    ws->factory("RooDecay::L_background_NP_punzi(Jpsict,NP_TauBkg[.37,.35,0.425],L_TotalPromptLifetime_punzi,RooDecay::SingleSided)");
	ws->factory("PROD::L_background_NP(L_background_NP_punzi|JpsictErr, pdf_ctauerrModelNP)");
    ws->factory("SUM::L_background(fBkgNP[.5,0,1.]*L_background_NP, L_TotalPromptLifetime)");
    //Combinatorial mumugamma background
    ws->factory("RooDecay::L_comb_backgroundSSD(Jpsict,jpsi_bkgTauSSD,L_TotalPromptLifetime_punzi,RooDecay::SingleSided)");
	ws->factory("RooDecay::L_comb_backgroundFD(Jpsict,jpsi_bkgTauFD,L_TotalPromptLifetime_punzi,RooDecay::Flipped)");
	ws->factory("RooDecay::L_comb_backgroundDSD(Jpsict,jpsi_bkgTauDSD,L_TotalPromptLifetime_punzi,RooDecay::DoubleSided)");
	ws->factory("SUM::L_comb_background_punzi(jpsi_fBkgSSDR*L_comb_backgroundSSD,jpsi_fBkgSSDL*L_comb_backgroundFD,L_comb_backgroundDSD)");
	ws->factory("PROD::L_comb_background(L_comb_background_punzi|JpsictErr, pdf_ctauerrModelBG)");


    //combine mass and lifetime
    ws->factory("PROD::ML_chic1_PR(M_chic1,L_TotalPromptLifetime)");
    ws->factory("PROD::ML_chic1_NP(M_chic1,L_chic1_NP)");
    ws->factory("SUM::ML_chic1(fracNP_chic1[0.25,0.1,0.4]*ML_chic1_NP, ML_chic1_PR)");

    ws->factory("PROD::ML_chic2_PR(M_chic2,L_TotalPromptLifetime)");
    ws->factory("PROD::ML_chic2_NP(M_chic2,L_chic2_NP)");
    ws->factory("SUM::ML_chic2(fracNP_chic2[0.15,0.,0.3]*ML_chic2_NP, ML_chic2_PR)");

    ws->factory("PROD::ML_chic0_PR(M_chic0,L_TotalPromptLifetime)");
    ws->factory("PROD::ML_chic0_NP(M_chic0,L_chic0_NP)");
    ws->factory("SUM::ML_chic0(fracNP_chic0[0.5,0.2,0.7]*ML_chic0_NP, ML_chic0_PR)");

    ws->factory("PROD::ML_background(M_background,L_background)");
    ws->factory("PROD::ML_comb_background(M_background,L_comb_background)");
    ws->factory("SUM::ML_signal(fracSignal_chic1[0.7,0.6,0.8]*ML_chic1, fracSignal_chic0[0.03,0.,0.1]*ML_chic0, ML_chic2)");

	RooFormulaVar fracSignal_chic2("fracSignal_chic2","1-@0-@1",RooArgList(*ws->var("fracSignal_chic0"),*ws->var("fracSignal_chic1")));
	ws->import(fracSignal_chic2);

    //full mass lifetime shape
    ws->factory("SUM::ML_fullModel(fracBackground[0.6,0.,1.]*ML_background, jpsi_fBkg*ML_comb_background, ML_signal)");
    //ws->factory(Form("ExtendPdf::ML_fullModel(ML_fullModelNonE, NumEvE[%f,%f,%f])",ev, ev/ExtensionFactor, ev*ExtensionFactor));

    //full mass shape
    ws->factory("SUM::M_signal(fracSignal_chic1*M_chic1, fracSignal_chic0*M_chic0, M_chic2)");
    ws->factory("SUM::M_fullModel(fracBackground*M_background, jpsi_fBkg*M_background, M_signal)");
    //ws->factory("ExtendPdf::M_fullModel(M_fullModelNonE, NumEvE)");


    RooRealVar NumEvE("NumEvE","NumEvE",0.);
	ws->import(NumEvE);



    ws->Print("v");

    //------------------------------------------------------------------------------------------------------------------
    // do fit
    std::cout << "Fitting" << std::endl;

    RooRealVar chicMass(*ws->var("chicMass"));


    //ws->var("q01S")->setVal(3.2);
    //ws->var("q01S")->setConstant(kTRUE);

    ws->var("promptMean")->setVal(0.);
    ws->var("promptMean")->setConstant(kTRUE);

    ws->var("BK_p2")->setVal(1.e-10);
    ws->var("BK_p2")->setConstant(kTRUE);

    if(MC){
        ws->var("fracBackground")->setVal(1.e-10);
        ws->var("fracBackground")->setConstant(kTRUE);

        ws->var("fracNP_chic0")->setVal(1.e-10);
        ws->var("fracNP_chic0")->setConstant(kTRUE);

        ws->var("fracSignal_chic0")->setVal(1.e-10);
        ws->var("fracSignal_chic0")->setConstant(kTRUE);

        ws->var("BK_p1")->setVal(1.e-10);
        ws->var("BK_p1")->setConstant(kTRUE);

        ws->var("BK_p2")->setVal(1.e-10);
        ws->var("BK_p2")->setConstant(kTRUE);



    }

    //ws->var("DSD_TauBkg")->setVal(0.04);
    //ws->var("DSD_TauBkg")->setConstant(kTRUE);



	if(ptBin==5){
		ws->var("CBalpha1")->setVal(0.55);
		ws->var("CBalpha1")->setConstant(kTRUE);
		ws->var("CBalpha2")->setVal(0.475);
		ws->var("CBalpha2")->setConstant(kTRUE);

		//ws->var("fBkgNP")->setVal(0.65);
		//ws->var("fBkgNP")->setConstant(kTRUE);

        //ws->var("fracNP_chic1")->setVal(0.34);
        //ws->var("fracNP_chic1")->setConstant(kTRUE);
        //ws->var("fracNP_chic2")->setVal(0.16);
        //ws->var("fracNP_chic2")->setConstant(kTRUE);
        //ws->var("fracBackground")->setVal(0.687);
        //ws->var("fracBackground")->setConstant(kTRUE);

	}

	if(ptBin>3){
		ws->var("CBmass1")->setVal(3.5062);
		ws->var("CBmass1")->setConstant(kTRUE);
	}


    //ws->var("FD_TauBkg")->setVal(0.1);
    //ws->var("FD_TauBkg")->setConstant(kTRUE);


    //ws->var("ctResolution")->setVal(0.9);
    //ws->var("ctResolution")->setConstant(kTRUE);
    //if(rapBin == 1){
    //    ws->var("ctResolution2")->setVal(1.1);
    //}
    //else if(rapBin == 2){
    //    ws->var("ctResolution2")->setVal(1.5);
    //}


    ws->var("CBn")->setVal(2.75);
    ws->var("CBn")->setConstant(kTRUE);
    ws->var("fracSignal_chic0")->setVal(0.03);
    ws->var("fracSignal_chic0")->setConstant(kTRUE);



    //Fix comb background from jpsi fit

    ws->var("ctResolution2")->setVal(ws->var("jpsi_ctResolution2")->getVal());
    ws->var("ctResolution2")->setConstant(kTRUE);
    ws->var("ctResolution")->setVal(ws->var("jpsi_ctResolution")->getVal());
    ws->var("ctResolution")->setConstant(kTRUE);
    ws->var("fracGauss2")->setVal(ws->var("jpsi_fracGauss2")->getVal());
    ws->var("fracGauss2")->setConstant(kTRUE);

    //Fix comb background from jpsi fit

	ws->var("jpsi_fBkg")->setConstant(kTRUE);
	ws->var("jpsi_fBkgSSDR")->setConstant(kTRUE);
	ws->var("jpsi_fBkgSSDL")->setConstant(kTRUE);
	ws->var("jpsi_bkgTauSSD")->setConstant(kTRUE);
	ws->var("jpsi_bkgTauFD")->setConstant(kTRUE);
	ws->var("jpsi_bkgTauDSD")->setConstant(kTRUE);


    //ws->var("fracGauss2")->setVal(0.);
    //ws->var("fracGauss2")->setConstant(kTRUE);

    //ws->var("fracSignal_chic0")->setVal(0.);
    //ws->var("fracSignal_chic0")->setConstant(kTRUE);

    //ws->var("fracNP_chic0")->setConstant(kTRUE);



    /*
    if(nState==4){
        //    ws->var("CBalpha")->setVal(2.);
    //    ws->var("bkgLambda")->setVal(-1.5);
        ws->var("CBmass")->setMax(3.15);
        ws->var("CBmass")->setMin(3.05);
        ws->var("CBmass")->setVal(3.1);
        ws->var("CBn")->setVal(2.5);
        ws->var("CBn")->setConstant(kTRUE);
        ws->var("CBsigma")->setMax(0.07);
        ws->var("CBsigma")->setVal(.02);
        ws->var("CBsigma2")->setVal(.02);
        ws->var("CBsigma2")->setMax(.07);
        ws->var("fracCB1")->setVal(0.5);
        ws->var("fracBkg")->setVal(0.1);

        ws->var("fBkgSSD")->setVal(.7); //0.4
        ws->var("fBkgDSD")->setVal(.2);
        ws->var("bkgTauSSD")->setVal(.4);
        ws->var("bkgTauFD")->setVal(.1);
        ws->var("bkgTauDSD")->setVal(.01);
        ws->var("fracNP")->setVal(.3);
        ws->var("nonPromptTau")->setVal(.4);
        //** set limit
        ws->var("bkgTauSSD")->setMax(1.);
        ws->var("bkgTauFD")->setMax(0.3);
        ws->var("nonPromptTau")->setMax(1.);
        ws->var("promptMean")->setVal(0.);
        ws->var("ctResolution")->setVal(.9);
        ws->var("promptMean")->setConstant(kTRUE);
        ws->var("ctResolution")->setConstant(kTRUE);
        if(rapBin == 1){
            ws->var("ctResolution2")->setVal(1.1);
            ws->var("CBalpha")->setVal(1.85);
            ws->var("bkgLambda")->setVal(2.0);
        }
        else if(rapBin == 2){
            ws->var("ctResolution2")->setVal(1.5);
            ws->var("CBalpha")->setVal(1.95);
            ws->var("bkgLambda")->setVal(1.75);
        }
        ws->var("ctResolution2")->setConstant(kTRUE);
        ws->var("CBalpha")->setConstant(kTRUE);
        ws->var("bkgLambda")->setConstant(kTRUE);
        ws->var("fracGauss2")->setVal(0.2);
        ws->var("fracGauss2")->setMin(0.04);
    }
    else if(nState==5){
        if(rapBin == 1){
            ws->var("CBalpha")->setVal(1.85);
            ws->var("bkgLambda")->setVal(2.0);
        }
        else if(rapBin == 2){
            ws->var("CBalpha")->setVal(1.95);
            ws->var("bkgLambda")->setVal(1.75);
        }
        ws->var("CBalpha")->setConstant(kTRUE);
        ws->var("bkgLambda")->setConstant(kTRUE);
        //ws->var("CBalpha")->setVal(2.);
        //ws->var("bkgLambda")->setVal(-1.5);
        ws->var("CBmass")->setMax(3.75);
        ws->var("CBmass")->setMin(3.65);
        ws->var("CBmass")->setVal(3.69);
        ws->var("CBalpha")->setVal(2.);
        ws->var("CBn")->setVal(2.5);
        ws->var("CBn")->setConstant(kTRUE);
        ws->var("CBsigma")->setMax(0.07);
        ws->var("CBsigma2")->setMax(0.07);
        ws->var("CBsigma")->setVal(.03);
        ws->var("CBsigma2")->setVal(.04);

        ws->var("fBkgSSD")->setVal(.6);
        ws->var("fBkgDSD")->setVal(.4);
        ws->var("bkgTauSSD")->setVal(.4);
        ws->var("bkgTauFD")->setVal(.1);
        ws->var("bkgTauDSD")->setVal(.01);
        ws->var("nonPromptTau")->setVal(.4);
        //** set limit
        ws->var("bkgTauSSD")->setMax(1.);
        ws->var("bkgTauFD")->setMax(0.3);
        ws->var("bkgTauDSD")->setMax(0.05);
        ws->var("nonPromptTau")->setMax(1.);
        ws->var("fBkgDSD")->setMin(0.2);
        ws->var("promptMean")->setVal(0.);
        ws->var("ctResolution")->setVal(.9);
        ws->var("promptMean")->setConstant(kTRUE);
        ws->var("ctResolution")->setConstant(kTRUE);
        ws->var("ctResolution2")->setVal(3.);
        ws->var("ctResolution2")->setConstant(kTRUE);
        ws->var("fracGauss2")->setVal(0.01);

    }

*/


    RooAbsPdf *fullPdf;

    if(!runChiMassFitOnly)
    	fullPdf = (RooAbsPdf*)ws->pdf("ML_fullModel");
    else
    	fullPdf = (RooAbsPdf*)ws->pdf("M_fullModel");


    data->Print("v");

    RooArgSet *NLLs = new RooArgSet();
    RooAbsReal *MassNLL = NULL;

    MassNLL = (RooAbsReal *)fullPdf->createNLL(*data, /*Extended(true), */NumCPU(4));

    NLLs->add(*MassNLL);

    RooAddition *simNLL = new RooAddition("add","add",*NLLs);
    RooMinuit *mMinuit = new RooMinuit(*simNLL);

    //mMinuit->setStrategy(2);
    mMinuit->setPrintEvalErrors(-1);
    mMinuit->setEvalErrorWall(false);
    mMinuit->setVerbose(false);

    //mMinuit->simplex();
    //mMinuit->migrad();
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
    rltName<< "fitresult_rap"<<rapBin<<"_pt"<<ptBin;
    snapshotName<< "snapshot_rap"<<rapBin<<"_pt"<<ptBin;
    RooFitResult *fitresult = (RooFitResult*)mMinuit->save(rltName.str().c_str());
    ws->import(*fitresult,kTRUE);
    fitresult->Print();
    double covQualHesse=fitresult->covQual();
	cout<<"covQualHesse = "<<covQualHesse<<endl;

    std::cout<< "save snapshot" <<std::endl;
    ws->saveSnapshot(snapshotName.str().c_str(), ws->allVars());

    RooRealVar var_covQualMigrad("var_covQualMigrad","var_covQualMigrad",covQualMigrad); if(!ws->var("var_covQualMigrad")) ws->import(var_covQualMigrad); else ws->var("var_covQualMigrad")->setVal(covQualMigrad);
    RooRealVar var_covQualHesse("var_covQualHesse","var_covQualHesse",covQualHesse); if(!ws->var("var_covQualHesse")) ws->import(var_covQualHesse); else ws->var("var_covQualHesse")->setVal(covQualHesse);


    //cout<<"Fitted Chic0 mass = "<<CBmass0.getVal()<<" GeV --> "<<(onia::Mchi0PDG-CBmass0.getVal())*1000<<" MeV smaller than PDG mass"<<endl;
    //cout<<"Fitted Chic1 mass = "<<ws->var("CBmass1")->getVal()<<" GeV --> "<<(onia::Mchi1PDG-ws->var("CBmass1")->getVal())*1000<<" MeV smaller than PDG mass"<<endl;
    //cout<<"Fitted Chic2 mass = "<<CBmass2.getVal()<<" GeV --> "<<(onia::Mchi2PDG-CBmass2.getVal())*1000<<" MeV smaller than PDG mass"<<endl;




    std::cout<< "write ws" <<std::endl;
    ws->Write();
    std::cout<< "print ws" <<std::endl;
    ws->Print("v");

    infile->Close();

}



