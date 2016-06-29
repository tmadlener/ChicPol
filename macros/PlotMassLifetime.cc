#include "rootIncludes.inc"
#include "commonVar.h"
#include "RooUtils.h"
#include "THStack.h"
#include "TCutG.h"

using namespace RooFit;

double plotMass(RooWorkspace *ws, int rapBin, int ptBin, int nState);
void plot2DMassLifetime(RooWorkspace *ws, int rapBin, int ptBin, int nState);
void plot2DMassLifetimePedagogical(RooWorkspace *ws, int rapBin, int ptBin, int nState);
double plotLifetime(RooWorkspace *ws, int rapBin, int ptBin, int nState, bool SpeedPlotting, int nSpeedPlotting, bool zoom);
double plotLifetimeSR1(RooWorkspace *ws, int rapBin, int ptBin, int nState, bool SpeedPlotting, int nSpeedPlotting, bool zoom);
double plotLifetimeSR2(RooWorkspace *ws, int rapBin, int ptBin, int nState, bool SpeedPlotting, int nSpeedPlotting, bool zoom);
double plotLifetimeLSB(RooWorkspace *ws, int rapBin, int ptBin, int nState, bool SpeedPlotting, int nSpeedPlotting, bool zoom);
double plotLifetimeRSB(RooWorkspace *ws, int rapBin, int ptBin, int nState, bool SpeedPlotting, int nSpeedPlotting, bool zoom);
//==============================================
void latexFloatingLifetimePars(RooWorkspace *ws, TLatex* latex);

void PlotMassLifetime(const std::string &infilename, int rapBin, int ptBin, int nState, int Plotting){

    TFile *infile = new TFile(infilename.c_str(), "UPDATE");
    if(!infile){
        std::cout << "Error: failed to open file with dataset" << std::endl;
    }
    RooWorkspace *ws=(RooWorkspace *)infile->Get("ws_masslifetime");
    if(!ws){
        std::cout << "Error: failed to open workspace " << std::endl;
    }


	bool SpeedPlotting=true;
	int nSpeedPlotting=250;
	bool plotZoom=true;

	double chi2ndfMass=0.;
	double chi2ndfLifetime=0.;
	double chi2ndfLifetimeSR1=0.;
	double chi2ndfLifetimeSR2=0.;
	double chi2ndfLifetimeLSB=0.;
	double chi2ndfLifetimeRSB=0.;
	double buff_chi2ndf=0.;

	switch (Plotting) {
		case 1:
			//std::cout << ">>>>Plotting 2D mass-lifetime all data" << std::endl;
			//plot2DMassLifetime(ws, rapBin, ptBin, nState);
			std::cout << ">>>>Plotting 2D mass-lifetime pedagogical" << std::endl;
			if(ptBin==0) plot2DMassLifetimePedagogical(ws, rapBin, ptBin, nState);
			std::cout << ">>>>Plotting mass all data" << std::endl;
			chi2ndfMass = plotMass(ws, rapBin, ptBin, nState);
			std::cout << ">>>>chi2ndfMass = " << chi2ndfMass << std::endl;
			std::cout << ">>>>Plotting lifetime: SR2" << std::endl;
			chi2ndfLifetimeSR2 = plotLifetimeSR2(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false);
			std::cout << ">>>>Plotting lifetime: SR1" << std::endl;
			chi2ndfLifetimeSR1 = plotLifetimeSR1(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false);
			std::cout << ">>>>Plotting lifetime: LSB" << std::endl;
			chi2ndfLifetimeLSB = plotLifetimeLSB(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false);
			std::cout << ">>>>Plotting lifetime: RSB" << std::endl;
			chi2ndfLifetimeRSB = plotLifetimeRSB(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false);
			std::cout << ">>>>Plotting lifetime all data" << std::endl;
			chi2ndfLifetime = plotLifetime(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false);


			if(plotZoom){
				std::cout << ">>>>Plotting lifetime: SR2 zoom" << std::endl;
				buff_chi2ndf = plotLifetimeSR2(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true);
				std::cout << ">>>>Plotting lifetime: SR1 zoom" << std::endl;
				buff_chi2ndf = plotLifetimeSR1(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true);
				std::cout << ">>>>Plotting lifetime: LSB zoom" << std::endl;
				buff_chi2ndf = plotLifetimeLSB(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true);
				std::cout << ">>>>Plotting lifetime: RSB zoom" << std::endl;
				buff_chi2ndf = plotLifetimeRSB(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true);
				std::cout << ">>>>Plotting lifetime all data zoom" << std::endl;
				buff_chi2ndf = plotLifetime(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true);
			}
			break;
		case 2:
			std::cout << ">>>>Plotting 2D mass-lifetime pedagogical" << std::endl;
			if(ptBin==0) plot2DMassLifetimePedagogical(ws, rapBin, ptBin, nState);
			std::cout << ">>>>Plotting mass all data" << std::endl;
			chi2ndfMass = plotMass(ws, rapBin, ptBin, nState);
			break;
		case 3:
			std::cout << ">>>>Plotting lifetime: SR2" << std::endl;
			chi2ndfLifetimeSR2 = plotLifetimeSR2(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false);
			std::cout << ">>>>Plotting lifetime: SR1" << std::endl;
			chi2ndfLifetimeSR1 = plotLifetimeSR1(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false);
			std::cout << ">>>>Plotting lifetime: LSB" << std::endl;
			chi2ndfLifetimeLSB = plotLifetimeLSB(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false);
			std::cout << ">>>>Plotting lifetime: RSB" << std::endl;
			chi2ndfLifetimeRSB = plotLifetimeRSB(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false);
			std::cout << ">>>>Plotting lifetime all data" << std::endl;
			chi2ndfLifetime = plotLifetime(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false);
			if(plotZoom){
				std::cout << ">>>>Plotting lifetime: SR2 zoom" << std::endl;
				buff_chi2ndf = plotLifetimeSR2(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true);
				std::cout << ">>>>Plotting lifetime: SR1 zoom" << std::endl;
				buff_chi2ndf = plotLifetimeSR1(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true);
				std::cout << ">>>>Plotting lifetime: LSB zoom" << std::endl;
				buff_chi2ndf = plotLifetimeLSB(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true);
				std::cout << ">>>>Plotting lifetime: RSB zoom" << std::endl;
				buff_chi2ndf = plotLifetimeRSB(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true);
				std::cout << ">>>>Plotting lifetime all data zoom" << std::endl;
				buff_chi2ndf = plotLifetime(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true);
			}
			break;
		case 4:
			std::cout << ">>>>Plotting mass all data" << std::endl;
			chi2ndfMass = plotMass(ws, rapBin, ptBin, nState);
			std::cout << ">>>>Plotting lifetime: SR1" << std::endl;
			chi2ndfLifetimeSR1 = plotLifetimeSR1(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false);
			if(plotZoom){
				std::cout << ">>>>Plotting lifetime: SR1 zoom" << std::endl;
				buff_chi2ndf = plotLifetimeSR1(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true);
			}
			break;
		case 5:
			std::cout << ">>>>Plotting mass all data" << std::endl;
			chi2ndfMass = plotMass(ws, rapBin, ptBin, nState);
			std::cout << ">>>>Plotting lifetime: SR2" << std::endl;
			chi2ndfLifetimeSR2 = plotLifetimeSR2(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false);
			if(plotZoom){
				std::cout << ">>>>Plotting lifetime: SR2 zoom" << std::endl;
				buff_chi2ndf = plotLifetimeSR2(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true);
			}
			break;
		case 6:
			std::cout << ">>>>Plotting mass all data" << std::endl;
			chi2ndfMass = plotMass(ws, rapBin, ptBin, nState);
			std::cout << ">>>>Plotting lifetime: LSB" << std::endl;
			chi2ndfLifetimeLSB = plotLifetimeLSB(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false);
			if(plotZoom){
				std::cout << ">>>>Plotting lifetime: LSB zoom" << std::endl;
				buff_chi2ndf = plotLifetimeLSB(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true);
			}
			break;
		case 7:
			std::cout << ">>>>Plotting mass all data" << std::endl;
			chi2ndfMass = plotMass(ws, rapBin, ptBin, nState);
			std::cout << ">>>>Plotting lifetime: RSB" << std::endl;
			chi2ndfLifetimeRSB = plotLifetimeRSB(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false);
			if(plotZoom){
				std::cout << ">>>>Plotting lifetime: RSB zoom" << std::endl;
				buff_chi2ndf = plotLifetimeRSB(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true);
			}
			break;
		case 8:
			std::cout << ">>>>Plotting mass all data" << std::endl;
			chi2ndfMass = plotMass(ws, rapBin, ptBin, nState);
			std::cout << ">>>>Plotting lifetime all data" << std::endl;
			chi2ndfLifetime = plotLifetime(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false);
			if(plotZoom){
				std::cout << ">>>>Plotting lifetime all data zoom" << std::endl;
				buff_chi2ndf = plotLifetime(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true);
			}
			break;
		default:
			std::cerr << "I do not know what do do with this value of Plotting" << std::endl;
	}

    infile->Close();

    TFile *infileOut = new TFile(infilename.c_str(), "UPDATE");
    if(!infile){
        std::cout << "Error: failed to open file with dataset" << std::endl;
    }
    RooWorkspace *wsOut=(RooWorkspace *)infileOut->Get("ws_masslifetime");
    if(!ws){
        std::cout << "Error: failed to open workspace " << std::endl;
    }

	RooRealVar var_chi2ndf_Mass("var_chi2ndf_Mass","var_chi2ndf_Mass",chi2ndfMass); var_chi2ndf_Mass.setVal(chi2ndfMass); if(!wsOut->var("var_chi2ndf_Mass")) wsOut->import(var_chi2ndf_Mass); else wsOut->var("var_chi2ndf_Mass")->setVal(chi2ndfMass);
	if(Plotting==1 || Plotting==3 || Plotting==4){ RooRealVar var_chi2ndf_Lifetime_SR1("var_chi2ndf_Lifetime_SR1","var_chi2ndf_Lifetime_SR1",chi2ndfLifetimeSR1); var_chi2ndf_Lifetime_SR1.setVal(chi2ndfLifetimeSR1); if(!wsOut->var("var_chi2ndf_Lifetime_SR1")) wsOut->import(var_chi2ndf_Lifetime_SR1); else wsOut->var("var_chi2ndf_Lifetime_SR1")->setVal(chi2ndfLifetimeSR1);}
	if(Plotting==1 || Plotting==3 || Plotting==5){ RooRealVar var_chi2ndf_Lifetime_SR2("var_chi2ndf_Lifetime_SR2","var_chi2ndf_Lifetime_SR2",chi2ndfLifetimeSR2); var_chi2ndf_Lifetime_SR2.setVal(chi2ndfLifetimeSR2); if(!wsOut->var("var_chi2ndf_Lifetime_SR2")) wsOut->import(var_chi2ndf_Lifetime_SR2); else wsOut->var("var_chi2ndf_Lifetime_SR2")->setVal(chi2ndfLifetimeSR2);}
	if(Plotting==1 || Plotting==3 || Plotting==6){ RooRealVar var_chi2ndf_Lifetime_LSB("var_chi2ndf_Lifetime_LSB","var_chi2ndf_Lifetime_LSB",chi2ndfLifetimeLSB); var_chi2ndf_Lifetime_LSB.setVal(chi2ndfLifetimeLSB); if(!wsOut->var("var_chi2ndf_Lifetime_LSB")) wsOut->import(var_chi2ndf_Lifetime_LSB); else wsOut->var("var_chi2ndf_Lifetime_LSB")->setVal(chi2ndfLifetimeLSB);}
	if(Plotting==1 || Plotting==3 || Plotting==7){ RooRealVar var_chi2ndf_Lifetime_RSB("var_chi2ndf_Lifetime_RSB","var_chi2ndf_Lifetime_RSB",chi2ndfLifetimeRSB); var_chi2ndf_Lifetime_RSB.setVal(chi2ndfLifetimeRSB); if(!wsOut->var("var_chi2ndf_Lifetime_RSB")) wsOut->import(var_chi2ndf_Lifetime_RSB); else wsOut->var("var_chi2ndf_Lifetime_RSB")->setVal(chi2ndfLifetimeRSB);}
	if(Plotting==1 || Plotting==3 || Plotting==8){ RooRealVar var_chi2ndf_Lifetime("var_chi2ndf_Lifetime","var_chi2ndf_Lifetime",chi2ndfLifetime); var_chi2ndf_Lifetime.setVal(chi2ndfLifetime); if(!wsOut->var("var_chi2ndf_Lifetime")) wsOut->import(var_chi2ndf_Lifetime); else wsOut->var("var_chi2ndf_Lifetime")->setVal(chi2ndfLifetime);}

    std::cout<< "write wsOut" <<std::endl;
    wsOut->Write();
    std::cout<< "print wsOut" <<std::endl;
    wsOut->Print("v");

    infileOut->Close();


	delete ws;
}



void plot2DMassLifetime(RooWorkspace *ws, int rapBin, int ptBin, int nState){
	TGaxis::SetMaxDigits(3);

	RooRealVar *chicMass = ws->var("chicMass");
	assert( 0 != chicMass );

	RooAbsData *data= ws->data(Form("data_rap%d_pt%d_SR",rapBin,ptBin));
	assert ( 0 != data );
	RooFitResult* fitRlt = dynamic_cast<RooFitResult*>(ws->obj(Form("fitresult_rap%d_pt%d",rapBin,ptBin)));
	assert ( 0 != fitRlt);

    // get variables
    RooRealVar *ct = ws->var("Jpsict");
    RooRealVar *ctErr = ws->var("JpsictErr");


	RooAbsPdf *fullMassPdf = ws->pdf("M_fullModel");
	assert ( 0 != fullMassPdf );
	RooAbsPdf *fullPdf = ws->pdf("ML_fullModel");
	assert ( 0 != fullPdf );

	int nEntries = data->numEntries();

	//calc global chi2 of 2D fit
	int parsFit=(fitRlt->floatParsFinal()).getSize(); //this used the full p.d.f.

	int n2DbinsX=50;
	int n2DbinsY=50;
	Double_t scaleFactor;

	scaleFactor=ws->var("var_ev")->getVal();
	TH2* model2Dhist = chicMass->createHistogram("mass vs ct pdf",*ct);
	model2Dhist->SetBins(n2DbinsX, onia::chimassMin, onia::chimassMax, n2DbinsY, onia::ctPlotMin, onia::ctPlotMax);
	fullPdf->fillHistogram(model2Dhist, RooArgList(*chicMass, *ct), scaleFactor);//, &RooArgSet(*ctErr));

	//model1Dhist = new TH1(n2DbinsX, onia::chimassMin, onia::chimassMax, n2DbinsY, -0.02, 0.18);
	//int genNevents=100000;
	//RooDataSet *genDataFromModel = fullPdf->generate(RooArgSet(*chicMass, *ct), genNevents);
	//TH2* model2Dhist = (TH2*)genDataFromModel->createHistogram(*chicMass, *ct, n2DbinsX, n2DbinsY,Form("chicMass>%f&&chicMass<%f&&Jpsict>%f&&Jpsict<%f",onia::chimassMin, onia::chimassMax, n2DbinsY, -0.02, 0.18));
	//model2Dhist->Scale(scaleFactor/double(genNevents));

	TH2* data2Dhist = chicMass->createHistogram("mass vs ct pdf",*ct);
	data2Dhist->SetBins(n2DbinsX, onia::chimassMin, onia::chimassMax, n2DbinsY, onia::ctPlotMin, onia::ctPlotMax);
	data->fillHistogram(data2Dhist, RooArgList(*chicMass, *ct));
	data2Dhist->SetMarkerStyle(20);
	//data2Dhist->Print("all");

	double chi2Full=0;
	int nBinsForNdf=0;
	for(int ix=1; ix<n2DbinsX+1; ix++){
		for(int iy=1; iy<n2DbinsY+1; iy++){

			if(data2Dhist->GetBinContent(ix,iy)>0){
				double chicFullElement=(model2Dhist->GetBinContent(ix,iy)-data2Dhist->GetBinContent(ix,iy))/data2Dhist->GetBinError(ix,iy);
				chicFullElement*=chicFullElement;
				chi2Full+=chicFullElement;
				nBinsForNdf++;
			}

		}
	}

	double ndfFull=nBinsForNdf-parsFit;

	cout<<"chi2Full = "<<chi2Full<<endl;
	cout<<"ndfFull = "<<ndfFull<<endl;
	cout<<"chi2Full/ndfFull = "<<chi2Full/ndfFull<<endl;

	TCanvas *c2;
	bool plot2D=true;
	if(plot2D){

		double minZ = 5e-1;
		double maxZ = data2Dhist->GetMaximum();


		for(int LinLog=0; LinLog<2; LinLog++){

			c2=new TCanvas("c2","",1000,900);

			model2Dhist->SetStats(0);
			model2Dhist->SetTitle("");
			model2Dhist->SetMinimum(minZ);
			model2Dhist->SetMaximum(maxZ);
			model2Dhist->GetXaxis()->SetTitle(chicMass->GetTitle());
			model2Dhist->GetXaxis()->SetTitleOffset(1.8);
			model2Dhist->GetYaxis()->SetTitle(ct->GetTitle());
			model2Dhist->GetYaxis()->SetTitleOffset(1.8);
			model2Dhist->Draw("SURF1");
			data2Dhist->Draw("esame");

			c2->SetTheta(25.);
			c2->SetPhi(130.);

			double left=0.7, top=0.885, textSize=0.03;
			TLatex *latex=new TLatex();
			latex->SetTextFont(42);
			latex->SetNDC(kTRUE);
			latex->SetTextSize(textSize);
			double step=textSize*1.6;


			//gStyle->SetPadBottomMargin(0.08); //0.12
			//gStyle->SetPadLeftMargin(0.09); //0.12
			//gStyle->SetPadRightMargin(0.035); //0.05
			//gStyle->SetPadTopMargin(0.05); //0.05


			textSize=0.03; latex->SetTextSize(textSize);
			if(rapBin<1) left=0.1125;
			else left=0.100;
			top=0.885;
			if(rapBin==0) latex->DrawLatex(left,top,Form("|y%s| < %.1f",onia::KinParticleChar,onia::rapForPTRange[onia::kNbRapForPTBins]));
			else if(rapBin==1) latex->DrawLatex(left,top,Form("|y%s| < %.1f",onia::KinParticleChar,onia::rapForPTRange[rapBin]));
			else latex->DrawLatex(left,top,Form("%.1f < |y%s| < %.1f",onia::rapForPTRange[rapBin-1],onia::KinParticleChar,onia::rapForPTRange[rapBin]));

			left=0.0755; top-=textSize*1.7;
			if(ptBin==0)
				latex->DrawLatex(left,top,Form("%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin],onia::KinParticleChar,onia::pTRange[rapBin][onia::kNbPTMaxBins]));
			else
				latex->DrawLatex(left,top,Form("%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin-1],onia::KinParticleChar,onia::pTRange[rapBin][ptBin]));

			top+=textSize*1.7/2.;
			left=0.65;
			textSize=0.03; latex->SetTextSize(textSize);
			latex->DrawLatex(left,top,Form("#chi^{2}/ndf = %.1f / %d", chi2Full, ndfFull));



			if(LinLog==0) c2->SetLogz(0);
			if(LinLog==1) c2->SetLogz(1);


			std::stringstream saveMassLifetime;
			std::stringstream saveMassLifetimeRoot;
			if(LinLog==0){
				saveMassLifetime << "Fit/chicFit/2D_mass_ctau_lin_rap" << rapBin << "_pt" << ptBin << ".pdf";
				//saveMassLifetimeRoot << "Fit/root/chic_2D_mass_ctau_lin_rap" << rapBin << "_pt" << ptBin << ".root";
			}
			if(LinLog==1){
				saveMassLifetime << "Fit/chicFit/2D_mass_ctau_log_rap" << rapBin << "_pt" << ptBin << ".pdf";
				saveMassLifetimeRoot << "Fit/root/chic_2D_mass_ctau_log_rap" << rapBin << "_pt" << ptBin << ".root";
			}
			c2->SaveAs(saveMassLifetime.str().c_str());
			c2->SaveAs(saveMassLifetimeRoot.str().c_str());

		}


	}

delete c2;
return;

}



void plot2DMassLifetimePedagogical(RooWorkspace *ws, int rapBin, int ptBin, int nState){
	TGaxis::SetMaxDigits(3);

	RooRealVar *chicMass = ws->var("chicMass");
	assert( 0 != chicMass );

	RooAbsData *data= ws->data(Form("data_rap%d_pt%d_SR",rapBin,ptBin));
	assert ( 0 != data );
	RooFitResult* fitRlt = dynamic_cast<RooFitResult*>(ws->obj(Form("fitresult_rap%d_pt%d",rapBin,ptBin)));
	assert ( 0 != fitRlt);

    // get variables
    RooRealVar *ct = ws->var("Jpsict");
    RooRealVar *ctErr = ws->var("JpsictErr");


	RooAbsPdf *fullMassPdf = ws->pdf("M_fullModel");
	assert ( 0 != fullMassPdf );
	RooAbsPdf *fullPdf = ws->pdf("ML_fullModel");
	assert ( 0 != fullPdf );

	int nEntries = data->numEntries();

	//calc global chi2 of 2D fit
	int parsFit=(fitRlt->floatParsFinal()).getSize(); //this used the full p.d.f.

	int n2DbinsX=100;
	int n2DbinsY=100;
	Double_t scaleFactor;

	double ctauMinPedag=-0.1999;
	double ctauMaxPedag=0.4999;
	double massMinPedag=onia::chimassMin;
	double massMaxPedag=onia::chimassMax;

    double sig1MaxMass = ws->var("var_sig1MaxMass")->getVal();
    double sig1MinMass = ws->var("var_sig1MinMass")->getVal();
    double sig2MaxMass = ws->var("var_sig2MaxMass")->getVal();
    double sig2MinMass = ws->var("var_sig2MinMass")->getVal();
    double lsbMaxMass = ws->var("var_lsbMaxMass")->getVal();
    double lsbMinMass = ws->var("var_lsbMinMass")->getVal();
    double rsbMaxMass = ws->var("var_rsbMaxMass")->getVal();
    double rsbMinMass = ws->var("var_rsbMinMass")->getVal();

    double PRMin = ws->var("var_PRMin")->getVal();
    double PRMax = ws->var("var_PRMax")->getVal();
    double NPMin = ws->var("var_NPMin")->getVal();
    double NPMax = ws->var("var_NPMax")->getVal();



    int colorPedag_Basis=18;
    int colorPedag_PRLSB=432-9;
    int colorPedag_NPLSB=432-6;
    int colorPedag_PRRSB=432-9;
    int colorPedag_NPRSB=432-6;
    int colorPedag_PRSR1=416-7;
    int colorPedag_NPSR1=800+0;
    int colorPedag_PRSR2=632-4;
    int colorPedag_NPSR2=900+5;

    //colorPedag_NPSR2=800;
    //colorPedag_NPSR1=600;

	scaleFactor=1.;//ws->var("var_ev")->getVal();
	TH2* model2Dhist = chicMass->createHistogram("mass vs ct pdf",*ct);
	model2Dhist->SetBins(n2DbinsX, massMinPedag, massMaxPedag, n2DbinsY, ctauMinPedag, ctauMaxPedag);
	fullPdf->fillHistogram(model2Dhist, RooArgList(*chicMass, *ct), scaleFactor);//, &RooArgSet(*ctErr));
	model2Dhist->SetFillColor(colorPedag_Basis);


	TH2* model2Dhist_NPLSB = chicMass->createHistogram("mass vs ct pdf NPLSB",*ct);
	model2Dhist_NPLSB->SetBins(n2DbinsX, massMinPedag, massMaxPedag, n2DbinsY, ctauMinPedag, ctauMaxPedag);
	fullPdf->fillHistogram(model2Dhist_NPLSB, RooArgList(*chicMass, *ct), scaleFactor);
	model2Dhist_NPLSB->SetFillColor(colorPedag_NPLSB);

	TH2* model2Dhist_PRLSB = chicMass->createHistogram("mass vs ct pdf PRLSB",*ct);
	model2Dhist_PRLSB->SetBins(n2DbinsX, massMinPedag, massMaxPedag, n2DbinsY, ctauMinPedag, ctauMaxPedag);
	fullPdf->fillHistogram(model2Dhist_PRLSB, RooArgList(*chicMass, *ct), scaleFactor);
	model2Dhist_PRLSB->SetFillColor(colorPedag_PRLSB);

	TH2* model2Dhist_NPRSB = chicMass->createHistogram("mass vs ct pdf NPRSB",*ct);
	model2Dhist_NPRSB->SetBins(n2DbinsX, massMinPedag, massMaxPedag, n2DbinsY, ctauMinPedag, ctauMaxPedag);
	fullPdf->fillHistogram(model2Dhist_NPRSB, RooArgList(*chicMass, *ct), scaleFactor);
	model2Dhist_NPRSB->SetFillColor(colorPedag_NPRSB);

	TH2* model2Dhist_PRRSB = chicMass->createHistogram("mass vs ct pdf PRRSB",*ct);
	model2Dhist_PRRSB->SetBins(n2DbinsX, massMinPedag, massMaxPedag, n2DbinsY, ctauMinPedag, ctauMaxPedag);
	fullPdf->fillHistogram(model2Dhist_PRRSB, RooArgList(*chicMass, *ct), scaleFactor);
	model2Dhist_PRRSB->SetFillColor(colorPedag_PRRSB);

	TH2* model2Dhist_NPSR1 = chicMass->createHistogram("mass vs ct pdf NPSR1",*ct);
	model2Dhist_NPSR1->SetBins(n2DbinsX, massMinPedag, massMaxPedag, n2DbinsY, ctauMinPedag, ctauMaxPedag);
	fullPdf->fillHistogram(model2Dhist_NPSR1, RooArgList(*chicMass, *ct), scaleFactor);
	model2Dhist_NPSR1->SetFillColor(colorPedag_NPSR1);

	TH2* model2Dhist_PRSR1 = chicMass->createHistogram("mass vs ct pdf PRSR1",*ct);
	model2Dhist_PRSR1->SetBins(n2DbinsX, massMinPedag, massMaxPedag, n2DbinsY, ctauMinPedag, ctauMaxPedag);
	fullPdf->fillHistogram(model2Dhist_PRSR1, RooArgList(*chicMass, *ct), scaleFactor);
	model2Dhist_PRSR1->SetFillColor(colorPedag_PRSR1);

	TH2* model2Dhist_NPSR2 = chicMass->createHistogram("mass vs ct pdf NPSR2",*ct);
	model2Dhist_NPSR2->SetBins(n2DbinsX, massMinPedag, massMaxPedag, n2DbinsY, ctauMinPedag, ctauMaxPedag);
	fullPdf->fillHistogram(model2Dhist_NPSR2, RooArgList(*chicMass, *ct), scaleFactor);
	model2Dhist_NPSR2->SetFillColor(colorPedag_NPSR2);

	TH2* model2Dhist_PRSR2 = chicMass->createHistogram("mass vs ct pdf PRSR2",*ct);
	model2Dhist_PRSR2->SetBins(n2DbinsX, massMinPedag, massMaxPedag, n2DbinsY, ctauMinPedag, ctauMaxPedag);
	fullPdf->fillHistogram(model2Dhist_PRSR2, RooArgList(*chicMass, *ct), scaleFactor);
	model2Dhist_PRSR2->SetFillColor(colorPedag_PRSR2);

	TCutG *cutg_NPLSB = new TCutG("cutg_NPLSB",5);
	cutg_NPLSB->SetPoint(0,lsbMinMass,NPMin);
	cutg_NPLSB->SetPoint(1,lsbMaxMass,NPMin);
	cutg_NPLSB->SetPoint(2,lsbMaxMass,ctauMaxPedag);
	cutg_NPLSB->SetPoint(3,lsbMinMass,ctauMaxPedag);
	cutg_NPLSB->SetPoint(4,lsbMinMass,NPMin);

	TCutG *cutg_PRLSB = new TCutG("cutg_PRLSB",5);
	cutg_PRLSB->SetPoint(0,lsbMinMass,PRMin);
	cutg_PRLSB->SetPoint(1,lsbMaxMass,PRMin);
	cutg_PRLSB->SetPoint(2,lsbMaxMass,PRMax);
	cutg_PRLSB->SetPoint(3,lsbMinMass,PRMax);
	cutg_PRLSB->SetPoint(4,lsbMinMass,PRMin);

	TCutG *cutg_NPRSB = new TCutG("cutg_NPRSB",5);
	cutg_NPRSB->SetPoint(0,rsbMinMass,NPMin);
	cutg_NPRSB->SetPoint(1,rsbMaxMass,NPMin);
	cutg_NPRSB->SetPoint(2,rsbMaxMass,ctauMaxPedag);
	cutg_NPRSB->SetPoint(3,rsbMinMass,ctauMaxPedag);
	cutg_NPRSB->SetPoint(4,rsbMinMass,NPMin);

	TCutG *cutg_PRRSB = new TCutG("cutg_PRRSB",5);
	cutg_PRRSB->SetPoint(0,rsbMinMass,PRMin);
	cutg_PRRSB->SetPoint(1,rsbMaxMass,PRMin);
	cutg_PRRSB->SetPoint(2,rsbMaxMass,PRMax);
	cutg_PRRSB->SetPoint(3,rsbMinMass,PRMax);
	cutg_PRRSB->SetPoint(4,rsbMinMass,PRMin);

	TCutG *cutg_NPSR1 = new TCutG("cutg_NPSR1",5);
	cutg_NPSR1->SetPoint(0,sig1MinMass,NPMin);
	cutg_NPSR1->SetPoint(1,sig1MaxMass,NPMin);
	cutg_NPSR1->SetPoint(2,sig1MaxMass,ctauMaxPedag);
	cutg_NPSR1->SetPoint(3,sig1MinMass,ctauMaxPedag);
	cutg_NPSR1->SetPoint(4,sig1MinMass,NPMin);

	TCutG *cutg_PRSR1 = new TCutG("cutg_PRSR1",5);
	cutg_PRSR1->SetPoint(0,sig1MinMass,PRMin);
	cutg_PRSR1->SetPoint(1,sig1MaxMass,PRMin);
	cutg_PRSR1->SetPoint(2,sig1MaxMass,PRMax);
	cutg_PRSR1->SetPoint(3,sig1MinMass,PRMax);
	cutg_PRSR1->SetPoint(4,sig1MinMass,PRMin);

	TCutG *cutg_NPSR2 = new TCutG("cutg_NPSR2",5);
	cutg_NPSR2->SetPoint(0,sig2MinMass,NPMin);
	cutg_NPSR2->SetPoint(1,sig2MaxMass,NPMin);
	cutg_NPSR2->SetPoint(2,sig2MaxMass,ctauMaxPedag);
	cutg_NPSR2->SetPoint(3,sig2MinMass,ctauMaxPedag);
	cutg_NPSR2->SetPoint(4,sig2MinMass,NPMin);

	TCutG *cutg_PRSR2 = new TCutG("cutg_PRSR2",5);
	cutg_PRSR2->SetPoint(0,sig2MinMass,PRMin);
	cutg_PRSR2->SetPoint(1,sig2MaxMass,PRMin);
	cutg_PRSR2->SetPoint(2,sig2MaxMass,PRMax);
	cutg_PRSR2->SetPoint(3,sig2MinMass,PRMax);
	cutg_PRSR2->SetPoint(4,sig2MinMass,PRMin);

	gStyle->SetPadBottomMargin(0.115); //0.12
	gStyle->SetPadLeftMargin(0.1125); //0.12
	gStyle->SetPadRightMargin(0.025); //0.05
	gStyle->SetPadTopMargin(0.025); //0.05


	TCanvas *c2;
	bool plot2D=true;
	if(plot2D){

		double minZ = 3.501e-6;
		double maxZ = model2Dhist->GetMaximum()*1.05;

		model2Dhist->SetMinimum(minZ);
		model2Dhist->SetMaximum(maxZ);

		model2Dhist_PRLSB->SetMinimum(minZ);
		model2Dhist_PRLSB->SetMaximum(maxZ);
		model2Dhist_NPLSB->SetMinimum(minZ);
		model2Dhist_NPLSB->SetMaximum(maxZ);

		model2Dhist_PRSR1->SetMinimum(minZ);
		model2Dhist_PRSR1->SetMaximum(maxZ);
		model2Dhist_NPSR1->SetMinimum(minZ);
		model2Dhist_NPSR1->SetMaximum(maxZ);

		model2Dhist_PRSR2->SetMinimum(minZ);
		model2Dhist_PRSR2->SetMaximum(maxZ);
		model2Dhist_NPSR2->SetMinimum(minZ);
		model2Dhist_NPSR2->SetMaximum(maxZ);

		model2Dhist_PRRSB->SetMinimum(minZ);
		model2Dhist_PRRSB->SetMaximum(maxZ);
		model2Dhist_NPRSB->SetMinimum(minZ);
		model2Dhist_NPRSB->SetMaximum(maxZ);

		for(int LinLog=0; LinLog<2; LinLog++){

			//c2=new TCanvas("c2","",1000,900);
			c2=new TCanvas("c2","",4000,3600);

			model2Dhist->SetStats(0);
			model2Dhist->SetTitle("");
			model2Dhist->GetXaxis()->SetTitle(chicMass->GetTitle());
			model2Dhist->GetXaxis()->SetTitleOffset(1.8);
			model2Dhist->GetYaxis()->SetTitle(ct->GetTitle());
			model2Dhist->GetYaxis()->SetTitleOffset(2.1);
			model2Dhist->GetZaxis()->SetTitle("arb. units");
			model2Dhist->GetZaxis()->SetTitleOffset(1.7);
			model2Dhist->Draw("SURF4 FB");



			model2Dhist_PRLSB->Draw("SURF4 FB BB A same [cutg_PRLSB]");
			model2Dhist_NPLSB->Draw("SURF4 FB BB A same [cutg_NPLSB]");
			model2Dhist_PRSR1->Draw("SURF4 FB BB A same [cutg_PRSR1]");
			model2Dhist_NPSR1->Draw("SURF4 FB BB A same [cutg_NPSR1]");
			model2Dhist_PRSR2->Draw("SURF4 FB BB A same [cutg_PRSR2]");
			model2Dhist_NPSR2->Draw("SURF4 FB BB A same [cutg_NPSR2]");
			model2Dhist_PRRSB->Draw("SURF4 FB BB A same [cutg_PRRSB]");
			model2Dhist_NPRSB->Draw("SURF4 FB BB A same [cutg_NPRSB]");




			c2->SetTheta(25.);
			//c2->SetPhi(130.);
			c2->SetPhi(200.);

			double left=0.7, top=0.885, textSize=0.03;
			TLatex *latex=new TLatex();
			latex->SetTextFont(42);
			latex->SetNDC(kTRUE);
			latex->SetTextSize(textSize);
			double step=textSize*1.6;




			if(LinLog==0) c2->SetLogz(0);
			if(LinLog==1) c2->SetLogz(1);


			std::stringstream saveMassLifetime;
			std::stringstream saveMassLifetimePng;
			std::stringstream saveMassLifetimeRoot;
			if(LinLog==0){
				saveMassLifetime << "Fit/chicFit/2D_mass_ctau_pedagogical_lin_rap" << rapBin << "_pt" << ptBin << ".pdf";
				saveMassLifetimePng << "Fit/chicFit/2D_mass_ctau_pedagogical_lin_rap" << rapBin << "_pt" << ptBin << ".gif";
			}
			if(LinLog==1){
				saveMassLifetime << "Fit/chicFit/2D_mass_ctau_pedagogical_log_rap" << rapBin << "_pt" << ptBin << ".pdf";
				saveMassLifetimePng << "Fit/chicFit/2D_mass_ctau_pedagogical_log_rap" << rapBin << "_pt" << ptBin << ".gif";
				saveMassLifetimeRoot << "Fit/root/chic_2D_mass_ctau_pedagogical_log_rap" << rapBin << "_pt" << ptBin << ".root";
				//c2->SaveAs(saveMassLifetimeRoot.str().c_str());
			}
			c2->SaveAs(saveMassLifetimePng.str().c_str());
			c2->SaveAs(saveMassLifetime.str().c_str());

		}


	}

delete c2;
return;

}


//==============================================
double plotMass(RooWorkspace *ws, int rapBin, int ptBin, int nState){
	int nbins=onia::ChicMassPlotBins;
	TGaxis::SetMaxDigits(3);

	double binWidth=(onia::chimassMax-onia::chimassMin)/double(nbins)*1000;

	RooRealVar *chicMass = ws->var("chicMass");
	assert( 0 != chicMass );

	RooPlot *massFrame = chicMass->frame(Range(onia::chimassMin, onia::chimassMax),Bins(nbins));
	assert ( 0 != massFrame );
	massFrame->SetName(Form("mass_plot_rap%d_pt%d",rapBin,ptBin));
	massFrame->SetTitle("");
	massFrame->GetYaxis()->SetTitle(Form("Events / %1.0f MeV",binWidth));
	massFrame->GetYaxis()->SetTitleOffset(1.3);

	RooPlot *massFramePull = chicMass->frame(Range(onia::chimassMin, onia::chimassMax),Bins(nbins));
	assert ( 0 != massFramePull );
	massFramePull->SetName(Form("pullmass_plot_rap%d_pt%d",rapBin,ptBin));
	massFramePull->SetTitle("");
	massFramePull->GetYaxis()->SetTitle("pull");
	massFramePull->GetXaxis()->SetTitleSize(0.08);
	massFramePull->GetYaxis()->SetTitleSize(0.08);
	massFramePull->GetXaxis()->SetLabelSize(0.08);
	massFramePull->GetYaxis()->SetLabelSize(0.08);
	massFramePull->GetYaxis()->SetTitleOffset(0.4);
	massFramePull->GetYaxis()->SetRangeUser(-5.99,5.99);


	RooAbsData *data= ws->data(Form("data_rap%d_pt%d_SR",rapBin,ptBin));
	assert ( 0 != data );
	RooFitResult* fitRlt = dynamic_cast<RooFitResult*>(ws->obj(Form("fitresult_rap%d_pt%d",rapBin,ptBin)));
	assert ( 0 != fitRlt);

    // get variables
    RooRealVar *ct = ws->var("Jpsict");
    RooRealVar *ctErr = ws->var("JpsictErr");

    //calculate mass ranges
    double sig1MaxMass = ws->var("var_sig1MaxMass")->getVal();
    double sig1MinMass = ws->var("var_sig1MinMass")->getVal();
    double sig2MaxMass = ws->var("var_sig2MaxMass")->getVal();
    double sig2MinMass = ws->var("var_sig2MinMass")->getVal();
    double lsbMaxMass = ws->var("var_lsbMaxMass")->getVal();
    double lsbMinMass = ws->var("var_lsbMinMass")->getVal();
    double rsbMaxMass = ws->var("var_rsbMaxMass")->getVal();
    double rsbMinMass = ws->var("var_rsbMinMass")->getVal();

	RooAbsPdf *fullMassPdf = ws->pdf("M_fullModel");
	assert ( 0 != fullMassPdf );
	// RooAbsPdf *fullPdf = ws->pdf("ML_fullModel");
	// assert ( 0 != fullPdf );

	int nEntries = data->numEntries();

	double minY = 0.;
	double maxY = 0.;
	double enlargeYby=onia::enlargeYby_ML;






	data->plotOn(massFrame,MarkerSize(onia::markerSize_ML), Name("myHist"));
	maxY = massFrame->GetMaximum()*enlargeYby;
	minY = 5e-1;
	massFrame->GetYaxis()->SetRangeUser(minY,maxY);

	fullMassPdf->plotOn(massFrame,
			LineWidth(onia::lineWidth_ML),
			ProjWData(*data), Name("myCurve"));

	int parsFit=(fitRlt->floatParsFinal()).getSize(); //this used the full p.d.f.
	int nBins_Mass=massFrame->GetNbinsX();
	double chi2Pre_Mass=massFrame->chiSquare(parsFit);  //reduced chi-squared = chi2/ndof
	int ndof_Mass=nBins_Mass-parsFit;  //num of degree of freedom
	double chi2_Mass=chi2Pre_Mass*ndof_Mass;

	TH1F* pull = new TH1F("pull","pull distribution", 100,-10.,10.);
	gSystem->mkdir("Fit/root",kTRUE);gSystem->mkdir("Fit/chicFit",kTRUE);
	TFile *pullFile = new TFile(Form("Fit/root/pull_mass_rap%d_pt%d.root",rapBin,ptBin),"RECREATE");

	RooHist* hpull_mass = massFrame->pullHist("myHist","myCurve",kTRUE);

	hpull_mass->SetMarkerSize(onia::markerSize_ML);
	for(int i=0;i<hpull_mass->GetN();i++){
		hpull_mass->SetPointEYlow(i,0.);
		hpull_mass->SetPointEYhigh(i,0.);
		double x,y;
		hpull_mass->GetPoint(i,x,y);
		pull->Fill(y);
	}
	pullFile->cd();
	pull->Write();
	pullFile->Close();

	double pullMean=pull->GetMean();
	double pullRMS=pull->GetRMS();
	double err_pullMean=pull->GetMeanError();
	double err_pullRMS=pull->GetRMSError();

	massFramePull->addPlotable(hpull_mass,"P");


	fullMassPdf->plotOn(massFrame,
			Components("M_background"),
			LineStyle(onia::lineStyle_subComps_ML),
			LineColor(onia::colorBackground),
			LineWidth(onia::lineWidth_ML),
			ProjWData(*data));

	fullMassPdf->plotOn(massFrame,
			Components("M_chic0"),
			LineStyle(onia::lineStyle_subComps_ML),
			LineColor(onia::colorChic0),
			LineWidth(onia::lineWidth_ML),
			ProjWData(*data));
	fullMassPdf->plotOn(massFrame,
			Components("M_chic1"),
			LineStyle(onia::lineStyle_subComps_ML),
			LineColor(onia::colorChic1),
			LineWidth(onia::lineWidth_ML),
			ProjWData(*data));

	fullMassPdf->plotOn(massFrame,
			Components("M_chic2"),
			LineStyle(onia::lineStyle_subComps_ML),
			LineColor(onia::colorChic2),
			LineWidth(onia::lineWidth_ML),
			ProjWData(*data));

	//fullMassPdf->paramOn(massFrame, Layout(0.6,0.95,0.95), Format("NE",AutoPrecision(2)));


	double lineWidthRegions = 1.5;

	gStyle->SetLineStyleString(11,"12 50");

	int lineStyleRegions=11;

	TLine *lineSBLow = new TLine(lsbMaxMass, minY, lsbMaxMass, maxY);
	TLine *lineSBHigh = new TLine(rsbMinMass, minY, rsbMinMass, maxY);
	TLine *lineSig1Low = new TLine(sig1MinMass, minY, sig1MinMass, maxY);
	TLine *lineSig1High = new TLine(sig1MaxMass, minY, sig1MaxMass, maxY);
	TLine *lineSig2Low = new TLine(sig2MinMass, minY, sig2MinMass, maxY);
	TLine *lineSig2High = new TLine(sig2MaxMass, minY, sig2MaxMass, maxY);
	lineSBLow->SetLineWidth(lineWidthRegions);
	lineSBLow->SetLineColor(onia::colorBackground);
	lineSBLow->SetLineStyle(lineStyleRegions);
	lineSBHigh->SetLineWidth(lineWidthRegions);
	lineSBHigh->SetLineColor(onia::colorBackground);
	lineSBHigh->SetLineStyle(lineStyleRegions);
	lineSig1Low->SetLineWidth(lineWidthRegions);
	lineSig1Low->SetLineColor(onia::colorChic1);
	lineSig1Low->SetLineStyle(lineStyleRegions);
	lineSig1High->SetLineWidth(lineWidthRegions);
	lineSig1High->SetLineColor(onia::colorChic1);
	lineSig1High->SetLineStyle(lineStyleRegions);
	lineSig2Low->SetLineWidth(lineWidthRegions);
	lineSig2Low->SetLineColor(onia::colorChic2);
	lineSig2Low->SetLineStyle(lineStyleRegions);
	lineSig2High->SetLineWidth(lineWidthRegions);
	lineSig2High->SetLineColor(onia::colorChic2);
	lineSig2High->SetLineStyle(lineStyleRegions);

	if(sig1MaxMass==sig2MinMass)
		lineSig1High->SetLineWidth(lineWidthRegions*2);

	TH1* legend_Tot = data->createHistogram("legend_Tot",*chicMass,Binning(50)) ; legend_Tot->SetLineColor(kBlue) ; legend_Tot->SetLineStyle(1) ; legend_Tot->SetLineWidth(2.) ;
	TH1* legend_Background = data->createHistogram("legend_Background",*chicMass,Binning(50)) ; legend_Background->SetLineColor(onia::colorBackground) ; legend_Background->SetLineStyle(onia::lineStyle_subComps_ML) ; legend_Background->SetLineWidth(onia::lineWidth_ML) ;
	TH1* legend_Chic0 = data->createHistogram("legend_Chic0",*chicMass,Binning(50)) ; legend_Chic0->SetLineColor(onia::colorChic0) ; legend_Chic0->SetLineStyle(onia::lineStyle_subComps_ML) ; legend_Chic0->SetLineWidth(onia::lineWidth_ML) ;
	TH1* legend_Chic1 = data->createHistogram("legend_Chic1",*chicMass,Binning(50)) ; legend_Chic1->SetLineColor(onia::colorChic1) ; legend_Chic1->SetLineStyle(onia::lineStyle_subComps_ML) ; legend_Chic1->SetLineWidth(onia::lineWidth_ML) ;
	TH1* legend_Chic2 = data->createHistogram("legend_Chic2",*chicMass,Binning(50)) ; legend_Chic2->SetLineColor(onia::colorChic2) ; legend_Chic2->SetLineStyle(onia::lineStyle_subComps_ML) ; legend_Chic2->SetLineWidth(onia::lineWidth_ML) ;

	TLegend* MassLegend=new TLegend(0.18,0.5,0.38,0.7);
	MassLegend->SetFillColor(kWhite);
	MassLegend->SetFillStyle(0);
	MassLegend->SetTextFont(42);
	MassLegend->SetTextSize(0.035);
	MassLegend->SetBorderSize(0.);
	MassLegend->AddEntry(legend_Tot,"sum","l");
	MassLegend->AddEntry(legend_Background,"BG","l");
	MassLegend->AddEntry(legend_Chic0,"#chi_{c0}","l");
	MassLegend->AddEntry(legend_Chic1,"#chi_{c1}","l");
	MassLegend->AddEntry(legend_Chic2,"#chi_{c2}","l");


	double left=0.7, top=0.885, textSize=0.03;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	double step=textSize*1.6;


	double leftmargin=0.09;
	double topmargin=0.05;
	double bottommargin=0.08;
	double rightmargin=0.02;

	gStyle->SetPadBottomMargin(0.08); //0.12
	gStyle->SetPadLeftMargin(0.09); //0.12
	gStyle->SetPadRightMargin(0.035); //0.05
	gStyle->SetPadTopMargin(0.05); //0.05
	
	TCanvas *c1;
	TPad *pad1;
	TPad *pad2;

	massFrame->SetMinimum(minY);
	massFrame->SetMaximum(maxY);


	for(int LinLog=0; LinLog<1; LinLog++){

		c1=new TCanvas("c1","",1000,900);

		c1->cd();
		pad1 = new TPad("pad1","pad1",0.,0.,1.,0.3);
		pad1->SetGridy();
		pad1->SetBottomMargin(0.2);
		pad1->Draw();
		c1->cd();
		pad2 = new TPad("pad2","pad2",0.,0.3,1.,1.);
		pad2->Draw();


		//Sig
		pad2->cd(0);
		if(LinLog==0) pad2->SetLogy(false);
		if(LinLog==1) pad2->SetLogy(true);


		massFrame->Draw(); MassLegend->Draw("same");
		lineSBLow->Draw("same");
		lineSBHigh->Draw("same");
		lineSig1Low->Draw("same");
		lineSig1High->Draw("same");
		lineSig2Low->Draw("same");
		lineSig2High->Draw("same");


		left=0.15; top=0.885; textSize=0.030; latex->SetTextSize(textSize);
		if(rapBin==0) latex->DrawLatex(left,top,Form("|y%s| < %.1f",onia::KinParticleChar,onia::rapForPTRange[onia::kNbRapForPTBins]));
		else if(rapBin==1) latex->DrawLatex(left,top,Form("|y%s| < %.1f",onia::KinParticleChar,onia::rapForPTRange[rapBin]));
		else latex->DrawLatex(left,top,Form("%.1f < |y%s| < %.1f",onia::rapForPTRange[rapBin-1],onia::KinParticleChar,onia::rapForPTRange[rapBin]));
		top-=step;
		if(ptBin==0)
			latex->DrawLatex(left,top,Form("%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin],onia::KinParticleChar,onia::pTRange[rapBin][onia::kNbPTBins[rapBin]]));
		else
			latex->DrawLatex(left,top,Form("%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin-1],onia::KinParticleChar,onia::pTRange[rapBin][ptBin]));

		top-=step;
		latex->SetTextColor(kRed);
		latex->DrawLatex(left,top,Form("J/#psi SR, #chi all"));
		latex->SetTextColor(kBlack);
		top-=step;
		latex->DrawLatex(left,top,Form("#chi^{2}/ndf = %.1f / %d", chi2_Mass, ndof_Mass));
		top-=step;


		double fullMass=(onia::chimassMax-onia::chimassMin)/(1-(rightmargin+leftmargin));
		double fullMassMin=onia::chimassMin-fullMass*leftmargin;

		textSize=0.02875; latex->SetTextSize(textSize);
		left=(sig1MinMass-fullMassMin)/fullMass+0.0125; top=0.885;
		latex->DrawLatex(left,top,Form("N^{SR1}_{#chi_{c1}} = %1.0f", ws->var("var_nChic1")->getVal()*ws->var("var_fChic1InSR1")->getVal()));
		top-=textSize*2.25;
		latex->DrawLatex(left,top,Form("f^{SR1}_{BG} = %1.2f", ws->var("var_fracBackgroundInSR1")->getVal()));

		left=(sig2MinMass-fullMassMin)/fullMass+0.0125; top=0.885;
		latex->DrawLatex(left,top,Form("N^{SR2}_{#chi_{c2}} = %1.0f", ws->var("var_nChic2")->getVal()*ws->var("var_fChic2InSR2")->getVal()));
		top-=textSize*2.25;
		latex->DrawLatex(left,top,Form("f^{SR2}_{BG} = %1.2f", ws->var("var_fracBackgroundInSR2")->getVal()));


		textSize=0.03; latex->SetTextSize(textSize);
		left=0.7575; top=0.885;
		double stepsizeTimes=1.9;
		latex->SetTextSize(textSize);

		latex->DrawLatex(left,top,Form("n^{Tot}_{data}  =  %.0f",ws->var("var_data_ev")->getVal()));
		top-=textSize*stepsizeTimes;
		//if(!ws->var("NumEvE")->isConstant()){
		//latex->DrawLatex(left,top,Form("n^{Tot}_{fit}  =  %.1f #pm %.1f",ws->var("NumEvE")->getVal(), ws->var("NumEvE")->getError()));
		//top-=textSize*stepsizeTimes;
		//}
		if(!ws->var("fracBackground")->isConstant()){
		latex->DrawLatex(left,top,Form("f^{#chi}_{#psiBG}  =  %.3f #pm %.3f",ws->var("fracBackground")->getVal(), ws->var("fracBackground")->getError()));
		top-=textSize*stepsizeTimes;
		}
		if(!ws->var("fracSignal_chic0")->isConstant()){
		latex->DrawLatex(left,top,Form("f^{Sig}_{#chi_{c0}}  =  %.3f #pm %.3f",ws->var("fracSignal_chic0")->getVal(), ws->var("fracSignal_chic0")->getError()));
		top-=textSize*stepsizeTimes;
		}
		if(!ws->var("fracSignal_chic1")->isConstant()){
		latex->DrawLatex(left,top,Form("f^{Sig}_{#chi_{c1}}  =  %.3f #pm %.3f",ws->var("fracSignal_chic1")->getVal(), ws->var("fracSignal_chic1")->getError()));
		top-=textSize*stepsizeTimes;
		}
		if(!ws->var("CBmass1")->isConstant()){
		latex->DrawLatex(left,top,Form("m_{#chi_{c1}}  =  %.1f #pm %.1f MeV",ws->var("CBmass1")->getVal()*1000, ws->var("CBmass1")->getError()*1000));
		top-=textSize*stepsizeTimes;
		}
		if(ws->var("CBmass2")!=NULL)
		if(!ws->var("CBmass2")->isConstant()){
		latex->DrawLatex(left,top,Form("m_{#chi_{c2}}  =  %.1f #pm %.1f MeV",ws->var("CBmass2")->getVal()*1000, ws->var("CBmass2")->getError()*1000));
		top-=textSize*stepsizeTimes;
		}
		if(!ws->var("CBsigma1")->isConstant()){
		latex->DrawLatex(left,top,Form("#sigma_{#chi_{c1}}  =  %.2f #pm %.2f MeV",ws->var("CBsigma1")->getVal()*1000, ws->var("CBsigma1")->getError()*1000));
		top-=textSize*stepsizeTimes;
		}
		if(ws->var("CBsigma2")!=NULL)
		if(!ws->var("CBsigma2")->isConstant()){
		latex->DrawLatex(left,top,Form("#sigma_{#chi_{c2}}  =  %.2f #pm %.2f MeV",ws->var("CBsigma2")->getVal()*1000, ws->var("CBsigma2")->getError()*1000));
		top-=textSize*stepsizeTimes;
		}
		if(ws->var("CBn")!=NULL)
		if(!ws->var("CBn")->isConstant()){
		latex->DrawLatex(left,top,Form("n^{CB}_{#chi_{c1}}  =  %.4f #pm %.4f",ws->var("CBn")->getVal(), ws->var("CBn")->getError()));
		top-=textSize*stepsizeTimes;
		}
		if(ws->var("CBn2")!=NULL)
		if(!ws->var("CBn2")->isConstant()){
		latex->DrawLatex(left,top,Form("n^{CB}_{#chi_{c2}}  =  %.4f #pm %.4f",ws->var("CBn2")->getVal(), ws->var("CBn2")->getError()));
		top-=textSize*stepsizeTimes;
		}
		if(!ws->var("CBalpha1")->isConstant()){
		latex->DrawLatex(left,top,Form("#alpha^{CB}_{#chi_{c1}}  =  %.3f #pm %.3f",ws->var("CBalpha1")->getVal(), ws->var("CBalpha1")->getError()));
		top-=textSize*stepsizeTimes;
		}
		if(!ws->var("CBalpha2")->isConstant()){
		latex->DrawLatex(left,top,Form("#alpha^{CB}_{#chi_{c2}}  =  %.3f #pm %.3f",ws->var("CBalpha2")->getVal(), ws->var("CBalpha2")->getError()));
		top-=textSize*stepsizeTimes;
		}
		if(!ws->var("alpha1")->isConstant()){
		latex->DrawLatex(left,top,Form("#alpha^{BG}  =  %.3f #pm %.3f",ws->var("alpha1")->getVal(), ws->var("alpha1")->getError()));
		top-=textSize*stepsizeTimes;
		}
		if(!ws->var("beta1")->isConstant()){
		latex->DrawLatex(left,top,Form("#beta^{BG}  =  %.3f #pm %.3f",ws->var("beta1")->getVal(), ws->var("beta1")->getError()));
		top-=textSize*stepsizeTimes;
		}
		if(!ws->var("q01S")->isConstant()){
		latex->DrawLatex(left,top,Form("Q_{0}^{BG}  =  %.4f #pm %.4f",ws->var("q01S")->getVal(), ws->var("q01S")->getError()));
		top-=textSize*stepsizeTimes;
		}
		if(!ws->var("BK_p1")->isConstant()){
		latex->DrawLatex(left,top,Form("p_{1}^{BG}  =  %.4f #pm %.4f",ws->var("BK_p1")->getVal(), ws->var("BK_p1")->getError()));
		top-=textSize*stepsizeTimes;
		}
		if(!ws->var("BK_p2")->isConstant()){
		latex->DrawLatex(left,top,Form("p_{2}^{BG}  =  %.4f #pm %.4f",ws->var("BK_p2")->getVal(), ws->var("BK_p2")->getError()));
		top-=textSize*stepsizeTimes;
		}


		textSize=0.015; latex->SetTextSize(textSize);
		left=0.925; top=0.925;
		latex->SetTextSize(textSize);
		latex->DrawLatex(left,top,Form("M%1.0fH%1.0f",ws->var("var_covQualMigrad")->getVal(), ws->var("var_covQualHesse")->getVal()));

		pad1->cd(0); pad1->SetLogy(0);
		massFramePull->Draw();

		textSize=0.06; latex->SetTextSize(textSize);
		left=0.385; top=0.875;
		latex->SetTextSize(textSize);
		//latex->DrawLatex(left,top,Form("#mu_{pull}  =  %.3f #pm %.3f,  #sigma_{pull}  =  %.3f #pm %.3f",pullMean, err_pullMean, pullRMS, err_pullRMS));

		c1->cd();

		std::stringstream saveMass;
		if(LinLog==0) saveMass << "Fit/chicFit/mass_lin_rap" << rapBin << "_pt" << ptBin << ".pdf";
		if(LinLog==1) saveMass << "Fit/chicFit/mass_log_rap" << rapBin << "_pt" << ptBin << ".pdf";
		c1->SaveAs(saveMass.str().c_str());


	}

	delete c1;
	delete legend_Tot;
	delete legend_Background;
	delete legend_Chic0;
	delete legend_Chic1;
	delete legend_Chic2;

	double returnChi2Ndf=double(chi2_Mass/ndof_Mass);

	return returnChi2Ndf;
}






//==============================================
double plotLifetime(RooWorkspace *ws, int rapBin, int ptBin, int nState, bool SpeedPlotting, int nSpeedPlotting, bool zoom){
	int nbins=onia::LifetimePlotBins;
	TGaxis::SetMaxDigits(3);

	double PlotMin, PlotMax;
	if(zoom){
		PlotMin=onia::ctPlotMinZoom;
		PlotMax=onia::ctPlotMaxZoom;
	}
	else{
		PlotMin=onia::ctPlotMin;
		PlotMax=onia::ctPlotMax;
	}

	double binWidth=(PlotMax-PlotMin)/double(nbins)*1000;


	bool correctResolutionForPlotting=false;
	double resCorrFactor=1.0525;
	if(ptBin==0) resCorrFactor=1.0825;
	if(ptBin==1) resCorrFactor=1.0725;
	if(correctResolutionForPlotting){
		ws->var("ctResolution")->setVal(ws->var("ctResolution")->getVal()*resCorrFactor);
		ws->var("ctResolution2")->setVal(ws->var("ctResolution2")->getVal()*resCorrFactor);
	}

	RooRealVar *chicMass = ws->var("chicMass");
	RooRealVar *Jpsict = ws->var("Jpsict");
	RooRealVar *JpsictErr = ws->var("JpsictErr");
	assert( 0 != Jpsict );

	RooPlot *ctauFrame = ((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins),Range(PlotMin, PlotMax));
	assert ( 0 != ctauFrame );
	ctauFrame->SetName(Form("ctau_plot_rap%d_pt%d",rapBin,ptBin));
	ctauFrame->SetTitle("");
	ctauFrame->GetYaxis()->SetTitle(Form("Events / %1.0f micron",binWidth));
	ctauFrame->GetYaxis()->SetTitleOffset(1.3);
	ctauFrame->GetXaxis()->SetRangeUser(PlotMin, PlotMax);

	RooPlot *ctauFramePull = ((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins),Range(PlotMin, PlotMax));
	assert ( 0 != ctauFramePull );
	ctauFramePull->SetName(Form("pullctau_plot_rap%d_pt%d",rapBin,ptBin));
	ctauFramePull->SetTitle("");
	ctauFramePull->GetYaxis()->SetTitle("pull");
	ctauFramePull->GetXaxis()->SetTitleSize(0.08);
	ctauFramePull->GetYaxis()->SetTitleSize(0.08);
	ctauFramePull->GetXaxis()->SetLabelSize(0.08);
	ctauFramePull->GetYaxis()->SetLabelSize(0.08);
	ctauFramePull->GetYaxis()->SetTitleOffset(0.4);
	ctauFramePull->GetYaxis()->SetRangeUser(-5.99,5.99);
	ctauFramePull->GetXaxis()->SetRangeUser(PlotMin, PlotMax);


	RooAbsData *data= ws->data(Form("data_rap%d_pt%d_SR",rapBin,ptBin));
	assert ( 0 != data );
	RooFitResult* fitRlt = dynamic_cast<RooFitResult*>(ws->obj(Form("fitresult_rap%d_pt%d",rapBin,ptBin)));
	assert ( 0 != fitRlt);

	fitRlt->Print();

    RooAbsData* dataFullRegionProj;
    if(SpeedPlotting){
    	dataFullRegionProj = data->reduce(EventRange(0,nSpeedPlotting));
    }
    else{
    	dataFullRegionProj = data->reduce(EventRange(0,1e10));
    }
    cout<<"number of events in data = "<<data->numEntries()<<endl;
    cout<<"number of events in dataFullRegionProj = "<<dataFullRegionProj->numEntries()<<endl;


	RooAbsPdf *fullPdf = ws->pdf("ML_fullModel");
	assert ( 0 != fullPdf );

	int nEntries = data->numEntries();

	double minY = 0.;
	double maxY = 0.;
	double enlargeYby=onia::enlargeYby_ML;

	data->plotOn(ctauFrame,MarkerSize(onia::markerSize_ML), Name("myHist"));
	maxY = ctauFrame->GetMaximum()*enlargeYby;
	minY = 5e-1;
	ctauFrame->GetYaxis()->SetRangeUser(minY,maxY);

	cout<<"Plotting sum"<<endl;
	fullPdf->plotOn(ctauFrame,
			//Normalization(ws->var("var_ev")->getVal(),2),
			LineWidth(onia::lineWidth_ML),
			ProjWData(*JpsictErr, *dataFullRegionProj),
			NumCPU(1),
			Name("myCurve"));
	cout<<"Plotting sum finished"<<endl;

	//------get chi2------------
	int parsFit=(fitRlt->floatParsFinal()).getSize(); //this used the full p.d.f.
	int nBins_Ctau=ctauFrame->GetNbinsX();
	double chi2Pre_Ctau=ctauFrame->chiSquare(parsFit);  //reduced chi-squared = chi2/ndof
	int ndof_Ctau=nBins_Ctau-parsFit;  //num of degree of freedom
	double chi2_Ctau=chi2Pre_Ctau*ndof_Ctau;

	TH1F* pull = new TH1F("pull","pull distribution", 100,-10.,10.);
	gSystem->mkdir("Fit/root",kTRUE);gSystem->mkdir("Fit/chicFit",kTRUE);
	TFile *pullFile = new TFile(Form("Fit/root/pull_ctau_all_rap%d_pt%d.root",rapBin,ptBin),"RECREATE");

	RooHist* hpull_ctau = ctauFrame->pullHist("myHist","myCurve",kTRUE);

	hpull_ctau->SetMarkerSize(onia::markerSize_ML);
	for(int i=0;i<hpull_ctau->GetN();i++){
		hpull_ctau->SetPointEYlow(i,0.);
		hpull_ctau->SetPointEYhigh(i,0.);
		double x,y;
		hpull_ctau->GetPoint(i,x,y);
		pull->Fill(y);
	}
	pullFile->cd();
	pull->Write();
	pullFile->Close();

	double pullMean=pull->GetMean();
	double pullRMS=pull->GetRMS();
	double err_pullMean=pull->GetMeanError();
	double err_pullRMS=pull->GetRMSError();

	ctauFramePull->addPlotable(hpull_ctau,"P");


	bool plotBackground=true;
	bool plotChic0=true;
	bool plotChic1=true;
	bool plotChic2=true;

	if(plotBackground){
		cout<<"Plotting background"<<endl;
		fullPdf->plotOn(ctauFrame,
				Components("ML_background,ML_comb_background"),
				//Normalization(ws->var("var_ev")->getVal(),2),
				LineStyle(onia::lineStyle_subComps_ML),
				LineColor(onia::colorBackground),
				LineWidth(onia::lineWidth_ML),
				ProjWData(*JpsictErr, *dataFullRegionProj), NumCPU(1));
		cout<<"Plotting background finished"<<endl;
	}

	if(plotChic0){
		cout<<"Plotting Chic0"<<endl;
		fullPdf->plotOn(ctauFrame,
				Components("ML_chic0"),
				//Normalization(ws->var("var_ev")->getVal(),2),
				LineStyle(onia::lineStyle_subComps_ML),
				LineColor(onia::colorChic0),
				LineWidth(onia::lineWidth_ML),
				ProjWData(*JpsictErr, *dataFullRegionProj), NumCPU(1));
		cout<<"Plotting Chic0 finished"<<endl;
	}

	if(plotChic1){
		cout<<"Plotting Chic1"<<endl;
		fullPdf->plotOn(ctauFrame,
				Components("ML_chic1"),
				//Normalization(ws->var("var_ev")->getVal(),2),
				LineStyle(onia::lineStyle_subComps_ML),
				LineColor(onia::colorChic1),
				LineWidth(onia::lineWidth_ML),
				ProjWData(*JpsictErr, *dataFullRegionProj), NumCPU(1));
		cout<<"Plotting Chic1 finished"<<endl;
	}

	if(plotChic2){
		cout<<"Plotting Chic2"<<endl;
		fullPdf->plotOn(ctauFrame,
				Components("ML_chic2"),
				//Normalization(ws->var("var_ev")->getVal(),2),
				LineStyle(onia::lineStyle_subComps_ML),
				LineColor(onia::colorChic2),
				LineWidth(onia::lineWidth_ML),
				ProjWData(*JpsictErr, *dataFullRegionProj), NumCPU(1));
		cout<<"Plotting Chic2 finished"<<endl;
	}




	//fullPdf->plotOn(ctauFrame,
	//		Components("ML_chic0_PR, ML_chic1_PR, ML_chic2_PR, ML_background"),
	//		//Normalization(ws->var("var_ev")->getVal(),2),
	//		LineStyle(onia::lineStyle_subComps_ML),
	//		LineColor(kOrange),
	//		LineWidth(onia::lineWidth_ML),
	//		ProjWData(*JpsictErr, *dataFullRegionProj), NumCPU(1));
	//fullPdf->plotOn(ctauFrame,
	//		Components("ML_chic0_NP, ML_chic1_NP, ML_chic2_NP, ML_background"),
	//		//Normalization(ws->var("var_ev")->getVal(),2),
	//		LineStyle(onia::lineStyle_subComps_ML),
	//		LineColor(kMagenta),
	//		LineWidth(onia::lineWidth_ML),
	//		ProjWData(*JpsictErr, *dataFullRegionProj), NumCPU(1));



	if(correctResolutionForPlotting){
		ws->var("ctResolution")->setVal(ws->var("ctResolution")->getVal()/resCorrFactor);
		ws->var("ctResolution2")->setVal(ws->var("ctResolution2")->getVal()/resCorrFactor);
	}

	double lineWidthRegions = 1.5;

	gStyle->SetLineStyleString(11,"12 50");

	int lineStyleRegions=11;

	TLine *linePRLow = new TLine(ws->var("var_PRMin")->getVal(), minY, ws->var("var_PRMin")->getVal(), maxY);
	TLine *linePRHigh = new TLine(ws->var("var_PRMax")->getVal(), minY, ws->var("var_PRMax")->getVal(), maxY);
	TLine *lineNPLow = new TLine(ws->var("var_NPMin")->getVal(), minY, ws->var("var_NPMin")->getVal(), maxY);
	TLine *lineNPHigh = new TLine(ws->var("var_NPMax")->getVal(), minY, ws->var("var_NPMax")->getVal(), maxY);

	linePRLow->SetLineWidth(lineWidthRegions);
	linePRLow->SetLineColor(onia::colorPR);
	linePRLow->SetLineStyle(lineStyleRegions);
	linePRHigh->SetLineWidth(lineWidthRegions);
	linePRHigh->SetLineColor(onia::colorPR);
	linePRHigh->SetLineStyle(lineStyleRegions);
	lineNPLow->SetLineWidth(lineWidthRegions);
	lineNPLow->SetLineColor(onia::colorNP);
	lineNPLow->SetLineStyle(lineStyleRegions);
	lineNPHigh->SetLineWidth(lineWidthRegions);
	lineNPHigh->SetLineColor(onia::colorNP);
	lineNPHigh->SetLineStyle(lineStyleRegions);


	if(ws->var("var_NPMin")->getVal()==ws->var("var_PRMax")->getVal())
		lineNPLow->SetLineWidth(lineWidthRegions*2);

	TH1* legend_Tot = data->createHistogram("legend_Tot",*Jpsict,Binning(50)) ; legend_Tot->SetLineColor(kBlue) ; legend_Tot->SetLineStyle(1) ; legend_Tot->SetLineWidth(2.) ;
	TH1* legend_Background = data->createHistogram("legend_Background",*Jpsict,Binning(50)) ; legend_Background->SetLineColor(onia::colorBackground) ; legend_Background->SetLineStyle(onia::lineStyle_subComps_ML) ; legend_Background->SetLineWidth(onia::lineWidth_ML) ;
	TH1* legend_Chic0 = data->createHistogram("legend_Chic0",*Jpsict,Binning(50)) ; legend_Chic0->SetLineColor(onia::colorChic0) ; legend_Chic0->SetLineStyle(onia::lineStyle_subComps_ML) ; legend_Chic0->SetLineWidth(onia::lineWidth_ML) ;
	TH1* legend_Chic1 = data->createHistogram("legend_Chic1",*Jpsict,Binning(50)) ; legend_Chic1->SetLineColor(onia::colorChic1) ; legend_Chic1->SetLineStyle(onia::lineStyle_subComps_ML) ; legend_Chic1->SetLineWidth(onia::lineWidth_ML) ;
	TH1* legend_Chic2 = data->createHistogram("legend_Chic2",*Jpsict,Binning(50)) ; legend_Chic2->SetLineColor(onia::colorChic2) ; legend_Chic2->SetLineStyle(onia::lineStyle_subComps_ML) ; legend_Chic2->SetLineWidth(onia::lineWidth_ML) ;

	TLegend* CtauLegend=new TLegend(0.13,0.71,0.23,0.91);
	CtauLegend->SetFillColor(kWhite);
	CtauLegend->SetFillStyle(0);
	CtauLegend->SetTextFont(42);
	CtauLegend->SetTextSize(0.035);
	CtauLegend->SetBorderSize(0.);
	CtauLegend->AddEntry(legend_Tot,"sum","l");
	if(plotBackground) CtauLegend->AddEntry(legend_Background,"BG","l");
	if(plotChic0) CtauLegend->AddEntry(legend_Chic0,"#chi_{c0}","l");
	if(plotChic1) CtauLegend->AddEntry(legend_Chic1,"#chi_{c1}","l");
	if(plotChic2) CtauLegend->AddEntry(legend_Chic2,"#chi_{c2}","l");


	double left=0.7, top=0.885, textSize=0.03;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	double step=textSize*1.6;


	gStyle->SetPadBottomMargin(0.08); //0.12
	gStyle->SetPadLeftMargin(0.09); //0.12
	gStyle->SetPadRightMargin(0.035); //0.05
	gStyle->SetPadTopMargin(0.05); //0.05

	TCanvas *c1;
	TPad *pad1;
	TPad *pad2;

	ctauFrame->SetMinimum(minY);
	ctauFrame->SetMaximum(maxY);

	for(int LinLog=0; LinLog<2; LinLog++){
		cout<<"LinLog "<<LinLog<<endl;

		c1=new TCanvas("c1","",1000,900);

		c1->cd();
		pad1 = new TPad("pad1","pad1",0.,0.,1.,0.3);
		pad1->SetGridy();
		pad1->SetBottomMargin(0.2);
		pad1->Draw();
		c1->cd();
		pad2 = new TPad("pad2","pad2",0.,0.3,1.,1.);
		pad2->Draw();


		//Sig
		pad2->cd(0);
		if(LinLog==0) pad2->SetLogy(false);
		if(LinLog==1) pad2->SetLogy(true);


		ctauFrame->Draw(); CtauLegend->Draw("same");
		lineNPLow->Draw("same");
		linePRLow->Draw("same");
		linePRHigh->Draw("same");



		left=0.54; top=0.885; textSize=0.030; latex->SetTextSize(textSize);
		if(rapBin==0) latex->DrawLatex(left,top,Form("|y%s| < %.1f",onia::KinParticleChar,onia::rapForPTRange[onia::kNbRapForPTBins]));
		else if(rapBin==1) latex->DrawLatex(left,top,Form("|y%s| < %.1f",onia::KinParticleChar,onia::rapForPTRange[rapBin]));
		else latex->DrawLatex(left,top,Form("%.1f < |y%s| < %.1f",onia::rapForPTRange[rapBin-1],onia::KinParticleChar,onia::rapForPTRange[rapBin]));
		top-=step;
		if(ptBin==0)
			latex->DrawLatex(left,top,Form("%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin],onia::KinParticleChar,onia::pTRange[rapBin][onia::kNbPTBins[rapBin]]));
		else
			latex->DrawLatex(left,top,Form("%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin-1],onia::KinParticleChar,onia::pTRange[rapBin][ptBin]));

		top-=step;
		latex->SetTextColor(kRed);
		latex->DrawLatex(left,top,Form("J/#psi SR, #chi all"));
		latex->SetTextColor(kBlack);
		top-=step;
		latex->DrawLatex(left,top,Form("#chi^{2}/ndf = %.1f / %d", chi2_Ctau, ndof_Ctau));
		top-=step;







	    RooRealVar var_chi2ndf_Lifetime("var_chi2ndf_Lifetime","var_chi2ndf_Lifetime",double(chi2_Ctau/ndof_Ctau)); if(!ws->var("var_chi2ndf_Lifetime")) ws->import(var_chi2ndf_Lifetime); else ws->var("var_chi2ndf_Lifetime")->setVal(double(chi2_Ctau/ndof_Ctau));

		latexFloatingLifetimePars(ws, latex);


		pad1->cd(0); pad1->SetLogy(0);
		ctauFramePull->Draw();

		textSize=0.06; latex->SetTextSize(textSize);
		left=0.385; top=0.875;
		latex->SetTextSize(textSize);
		//latex->DrawLatex(left,top,Form("#mu_{pull}  =  %.3f #pm %.3f,  #sigma_{pull}  =  %.3f #pm %.3f",pullMean, err_pullMean, pullRMS, err_pullRMS));

		c1->cd();

		std::stringstream saveCtau;
		if(!zoom){
			if(LinLog==0) saveCtau << "Fit/chicFit/ctau_all_lin_rap" << rapBin << "_pt" << ptBin << ".pdf";
			if(LinLog==1) saveCtau << "Fit/chicFit/ctau_all_log_rap" << rapBin << "_pt" << ptBin << ".pdf";
		}
		else{
			if(LinLog==0) saveCtau << "Fit/chicFit/ctau_all_zoom_lin_rap" << rapBin << "_pt" << ptBin << ".pdf";
			if(LinLog==1) saveCtau << "Fit/chicFit/ctau_all_zoom_log_rap" << rapBin << "_pt" << ptBin << ".pdf";
		}
		c1->SaveAs(saveCtau.str().c_str());


	}

	delete c1;
	delete legend_Tot;
	delete legend_Background;
	delete legend_Chic0;
	delete legend_Chic1;
	delete legend_Chic2;

	double returnChi2Ndf=double(chi2_Ctau/ndof_Ctau);

	return returnChi2Ndf;

}










//==============================================
double plotLifetimeSR1(RooWorkspace *ws, int rapBin, int ptBin, int nState, bool SpeedPlotting, int nSpeedPlotting, bool zoom){
	int nbins=onia::LifetimePlotBins;
	TGaxis::SetMaxDigits(3);

	double PlotMin, PlotMax;
	if(zoom){
		PlotMin=onia::ctPlotMinZoom;
		PlotMax=onia::ctPlotMaxZoom;
	}
	else{
		PlotMin=onia::ctPlotMin;
		PlotMax=onia::ctPlotMax;
	}

	double binWidth=(PlotMax-PlotMin)/double(nbins)*1000;

	bool correctResolutionForPlotting=false;
	double resCorrFactor=1.0525;
	if(ptBin==0) resCorrFactor=1.0825;
	if(ptBin==1) resCorrFactor=1.0725;
	if(correctResolutionForPlotting){
		ws->var("ctResolution")->setVal(ws->var("ctResolution")->getVal()*resCorrFactor);
		ws->var("ctResolution2")->setVal(ws->var("ctResolution2")->getVal()*resCorrFactor);
	}


	RooRealVar *chicMass = ws->var("chicMass");
	RooRealVar *Jpsict = ws->var("Jpsict");
	RooRealVar *JpsictErr = ws->var("JpsictErr");
	assert( 0 != Jpsict );

	RooPlot *ctauFrame = ((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins),Range(PlotMin, PlotMax));
	assert ( 0 != ctauFrame );
	ctauFrame->SetName(Form("ctau_plot_rap%d_pt%d",rapBin,ptBin));
	ctauFrame->SetTitle("");
	ctauFrame->GetYaxis()->SetTitle(Form("Events / %1.0f micron",binWidth));
	ctauFrame->GetYaxis()->SetTitleOffset(1.3);
	ctauFrame->GetXaxis()->SetRangeUser(PlotMin, PlotMax);

	RooPlot *ctauFramePull = ((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins),Range(PlotMin, PlotMax));
	assert ( 0 != ctauFramePull );
	ctauFramePull->SetName(Form("pullctau_plot_rap%d_pt%d",rapBin,ptBin));
	ctauFramePull->SetTitle("");
	ctauFramePull->GetYaxis()->SetTitle("pull");
	ctauFramePull->GetXaxis()->SetTitleSize(0.08);
	ctauFramePull->GetYaxis()->SetTitleSize(0.08);
	ctauFramePull->GetXaxis()->SetLabelSize(0.08);
	ctauFramePull->GetYaxis()->SetLabelSize(0.08);
	ctauFramePull->GetYaxis()->SetTitleOffset(0.4);
	ctauFramePull->GetYaxis()->SetRangeUser(-5.99,5.99);
	ctauFramePull->GetXaxis()->SetRangeUser(PlotMin, PlotMax);


	RooAbsData *data= ws->data(Form("data_rap%d_pt%d_SR",rapBin,ptBin));
	assert ( 0 != data );
	RooFitResult* fitRlt = dynamic_cast<RooFitResult*>(ws->obj(Form("fitresult_rap%d_pt%d",rapBin,ptBin)));
	assert ( 0 != fitRlt);

	fitRlt->Print();


    //calculate ctau ranges
    double sig1MaxMass = ws->var("var_sig1MaxMass")->getVal();
    double sig1MinMass = ws->var("var_sig1MinMass")->getVal();
    double sig2MaxMass = ws->var("var_sig2MaxMass")->getVal();
    double sig2MinMass = ws->var("var_sig2MinMass")->getVal();
    double lsbMaxMass = ws->var("var_lsbMaxMass")->getVal();
    double lsbMinMass = ws->var("var_lsbMinMass")->getVal();
    double rsbMaxMass = ws->var("var_rsbMaxMass")->getVal();
    double rsbMinMass = ws->var("var_rsbMinMass")->getVal();


    // define data in different regions
    std::stringstream cutSR1;
    cutSR1 << "chicMass > " << sig1MinMass << " && chicMass < " << sig1MaxMass;

    std::stringstream binNameSR1;
    binNameSR1  << "data_rap" << rapBin << "_pt" << ptBin << "_SR1";

    RooAbsData* dataSR1 = data->reduce(Cut(cutSR1.str().c_str()));

    RooAbsData* dataSR1Proj;
    if(SpeedPlotting){
    	//dataSR1Proj = dataSR1->reduce(EventRange(0,nSpeedPlotting));
    	dataSR1Proj = data->reduce(EventRange(0,nSpeedPlotting));
    }
    else{
    	//dataSR1Proj = dataSR1->reduce(EventRange(0,1e10));
    	dataSR1Proj = data->reduce(EventRange(0,1e10));
    }
    cout<<"number of events in dataSR1 = "<<dataSR1->numEntries()<<endl;
    cout<<"number of events in dataSR1Proj = "<<dataSR1Proj->numEntries()<<endl;

	RooAbsPdf *fullPdf = ws->pdf("ML_fullModel_SR1");
	assert ( 0 != fullPdf );

	int nEntries = dataSR1->numEntries();

	double minY = 0.;
	double maxY = 0.;
	double enlargeYby=onia::enlargeYby_ML;

	dataSR1->plotOn(ctauFrame,MarkerSize(onia::markerSize_ML), Name("myHist"));
	maxY = ctauFrame->GetMaximum()*enlargeYby;
	minY = 5e-1;
	ctauFrame->GetYaxis()->SetRangeUser(minY,maxY);

	cout<<"Plotting sum"<<endl;
	fullPdf->plotOn(ctauFrame,
			Normalization(ws->var("var_ev")->getVal()*ws->var("var_fTotInSR1")->getVal(),2),
			LineWidth(onia::lineWidth_ML),
			ProjWData(*JpsictErr, *dataSR1Proj),
			NumCPU(1),
			Name("myCurve"));
	cout<<"Plotting sum finished"<<endl;

	//------get chi2------------
	int parsFit=(fitRlt->floatParsFinal()).getSize(); //this used the full p.d.f.
	int nBins_Ctau=ctauFrame->GetNbinsX();
	double chi2Pre_Ctau=ctauFrame->chiSquare(parsFit);  //reduced chi-squared = chi2/ndof
	int ndof_Ctau=nBins_Ctau-parsFit;  //num of degree of freedom
	double chi2_Ctau=chi2Pre_Ctau*ndof_Ctau;

	TH1F* pull = new TH1F("pull","pull distribution", 100,-10.,10.);
	gSystem->mkdir("Fit/root",kTRUE);gSystem->mkdir("Fit/chicFit",kTRUE);
	TFile *pullFile = new TFile(Form("Fit/root/pull_ctau_SR1_rap%d_pt%d.root",rapBin,ptBin),"RECREATE");

	RooHist* hpull_ctau = ctauFrame->pullHist("myHist","myCurve",kTRUE);
	hpull_ctau->SetMarkerSize(onia::markerSize_ML);
	for(int i=0;i<hpull_ctau->GetN();i++){
		hpull_ctau->SetPointEYlow(i,0.);
		hpull_ctau->SetPointEYhigh(i,0.);
		double x,y;
		hpull_ctau->GetPoint(i,x,y);
		pull->Fill(y);
	}
	pullFile->cd();
	pull->Write();
	pullFile->Close();

	double pullMean=pull->GetMean();
	double pullRMS=pull->GetRMS();
	double err_pullMean=pull->GetMeanError();
	double err_pullRMS=pull->GetRMSError();

	ctauFramePull->addPlotable(hpull_ctau,"P");


	cout<<"nBackgroundInSR1: "<<ws->var("var_ev")->getVal()*ws->var("var_fTotInSR1")->getVal()*ws->var("var_fracBackgroundInSR1")->getVal()<<endl;
	cout<<"nChic0InSR1: "<<ws->var("var_ev")->getVal()*ws->var("var_fTotInSR1")->getVal()*ws->var("var_fracChic0InSR1")->getVal()<<endl;
	cout<<"nChic1InSR1: "<<ws->var("var_ev")->getVal()*ws->var("var_fTotInSR1")->getVal()*ws->var("var_fracChic1InSR1")->getVal()<<endl;
	cout<<"nChic2InSR1: "<<ws->var("var_ev")->getVal()*ws->var("var_fTotInSR1")->getVal()*ws->var("var_fracChic2InSR1")->getVal()<<endl;

	double minCompFrac=0.03;

	bool plotBackground=false;
	bool plotChic0=false;
	bool plotChic1=false;
	bool plotChic2=false;
	if(ws->var("var_fracBackgroundInSR1")->getVal()>minCompFrac) plotBackground=true;
	if(ws->var("var_fracChic0InSR1")->getVal()>minCompFrac) plotChic0=true;
	if(ws->var("var_fracChic1InSR1")->getVal()>minCompFrac) plotChic1=true;
	if(ws->var("var_fracChic2InSR1")->getVal()>minCompFrac) plotChic2=true;

	if(plotBackground){
		cout<<"Plotting background"<<endl;
		fullPdf->plotOn(ctauFrame,
				Components("ML_background,ML_comb_background"),
				//Normalization(ws->var("var_ev")->getVal()*ws->var("var_fTotInSR1")->getVal(),2),
				LineStyle(onia::lineStyle_subComps_ML),
				LineColor(onia::colorBackground),
				LineWidth(onia::lineWidth_ML),
				ProjWData(*JpsictErr, *dataSR1Proj), NumCPU(1));
		cout<<"Plotting background finished"<<endl;
	}

	if(plotChic0){
		cout<<"Plotting Chic0"<<endl;
		fullPdf->plotOn(ctauFrame,
				Components("ML_chic0"),
				//Normalization(ws->var("var_ev")->getVal()*ws->var("var_fTotInSR1")->getVal(),2),
				LineStyle(onia::lineStyle_subComps_ML),
				LineColor(onia::colorChic0),
				LineWidth(onia::lineWidth_ML),
				ProjWData(*JpsictErr, *dataSR1Proj), NumCPU(1));
		cout<<"Plotting Chic0 finished"<<endl;
	}

	if(plotChic1){
		cout<<"Plotting Chic1"<<endl;
		fullPdf->plotOn(ctauFrame,
				Components("ML_chic1"),
				//Normalization(ws->var("var_ev")->getVal()*ws->var("var_fTotInSR1")->getVal(),2),
				LineStyle(onia::lineStyle_subComps_ML),
				LineColor(onia::colorChic1),
				LineWidth(onia::lineWidth_ML),
				ProjWData(*JpsictErr, *dataSR1Proj), NumCPU(1));
		cout<<"Plotting PR Chic1"<<endl;
		fullPdf->plotOn(ctauFrame,
				Components("ML_chic1_PR"),
				//Normalization(ws->var("var_ev")->getVal()*ws->var("var_fTotInSR1")->getVal(),2),
				LineStyle(onia::lineStyle_subCompPRsignal_ML),
				LineColor(onia::colorChic1),
				LineWidth(onia::lineWidth_PRsignal_ML),
				ProjWData(*JpsictErr, *dataSR1Proj), NumCPU(1));
		cout<<"Plotting NP Chic1"<<endl;
		fullPdf->plotOn(ctauFrame,
				Components("ML_chic1_NP"),
				//Normalization(ws->var("var_ev")->getVal()*ws->var("var_fTotInSR1")->getVal(),2),
				LineStyle(onia::lineStyle_subCompPRsignal_ML),
				LineColor(onia::colorChic1),
				LineWidth(onia::lineWidth_PRsignal_ML),
				ProjWData(*JpsictErr, *dataSR1Proj), NumCPU(1));
		cout<<"Plotting Chic1 finished"<<endl;
	}

	if(plotChic2){
		cout<<"Plotting Chic2"<<endl;
		fullPdf->plotOn(ctauFrame,
				Components("ML_chic2"),
				//Normalization(ws->var("var_ev")->getVal()*ws->var("var_fTotInSR1")->getVal(),2),
				LineStyle(onia::lineStyle_subComps_ML),
				LineColor(onia::colorChic2),
				LineWidth(onia::lineWidth_ML),
				ProjWData(*JpsictErr, *dataSR1Proj), NumCPU(1));
		cout<<"Plotting Chic2 finished"<<endl;
	}

	if(correctResolutionForPlotting){
		ws->var("ctResolution")->setVal(ws->var("ctResolution")->getVal()/resCorrFactor);
		ws->var("ctResolution2")->setVal(ws->var("ctResolution2")->getVal()/resCorrFactor);
	}


	double lineWidthRegions = 1.5;

	gStyle->SetLineStyleString(11,"12 50");

	int lineStyleRegions=11;

	TLine *linePRLow = new TLine(ws->var("var_PRMin")->getVal(), minY, ws->var("var_PRMin")->getVal(), maxY);
	TLine *linePRHigh = new TLine(ws->var("var_PRMax")->getVal(), minY, ws->var("var_PRMax")->getVal(), maxY);
	TLine *lineNPLow = new TLine(ws->var("var_NPMin")->getVal(), minY, ws->var("var_NPMin")->getVal(), maxY);
	TLine *lineNPHigh = new TLine(ws->var("var_NPMax")->getVal(), minY, ws->var("var_NPMax")->getVal(), maxY);

	linePRLow->SetLineWidth(lineWidthRegions);
	linePRLow->SetLineColor(onia::colorPR);
	linePRLow->SetLineStyle(lineStyleRegions);
	linePRHigh->SetLineWidth(lineWidthRegions);
	linePRHigh->SetLineColor(onia::colorPR);
	linePRHigh->SetLineStyle(lineStyleRegions);
	lineNPLow->SetLineWidth(lineWidthRegions);
	lineNPLow->SetLineColor(onia::colorNP);
	lineNPLow->SetLineStyle(lineStyleRegions);
	lineNPHigh->SetLineWidth(lineWidthRegions);
	lineNPHigh->SetLineColor(onia::colorNP);
	lineNPHigh->SetLineStyle(lineStyleRegions);


	if(ws->var("var_NPMin")->getVal()==ws->var("var_PRMax")->getVal())
		lineNPLow->SetLineWidth(lineWidthRegions*2);

	TH1* legend_Tot = data->createHistogram("legend_Tot",*Jpsict,Binning(50)) ; legend_Tot->SetLineColor(kBlue) ; legend_Tot->SetLineStyle(1) ; legend_Tot->SetLineWidth(2.) ;
	TH1* legend_Background = data->createHistogram("legend_Background",*Jpsict,Binning(50)) ; legend_Background->SetLineColor(onia::colorBackground) ; legend_Background->SetLineStyle(onia::lineStyle_subComps_ML) ; legend_Background->SetLineWidth(onia::lineWidth_ML) ;
	TH1* legend_Chic0 = data->createHistogram("legend_Chic0",*Jpsict,Binning(50)) ; legend_Chic0->SetLineColor(onia::colorChic0) ; legend_Chic0->SetLineStyle(onia::lineStyle_subComps_ML) ; legend_Chic0->SetLineWidth(onia::lineWidth_ML) ;
	TH1* legend_Chic1 = data->createHistogram("legend_Chic1",*Jpsict,Binning(50)) ; legend_Chic1->SetLineColor(onia::colorChic1) ; legend_Chic1->SetLineStyle(onia::lineStyle_subComps_ML) ; legend_Chic1->SetLineWidth(onia::lineWidth_ML) ;
	TH1* legend_Chic2 = data->createHistogram("legend_Chic2",*Jpsict,Binning(50)) ; legend_Chic2->SetLineColor(onia::colorChic2) ; legend_Chic2->SetLineStyle(onia::lineStyle_subComps_ML) ; legend_Chic2->SetLineWidth(onia::lineWidth_ML) ;

	TH1* legend_ChicPR = data->createHistogram("legend_ChicPR",*Jpsict,Binning(50)) ; legend_ChicPR->SetLineColor(onia::colorChic1) ; legend_ChicPR->SetLineStyle(onia::lineStyle_subCompPRsignal_ML) ; legend_ChicPR->SetLineWidth(onia::lineWidth_PRsignal_ML) ;

	TLegend* CtauLegend=new TLegend(0.13,0.71,0.23,0.91);
	CtauLegend->SetFillColor(kWhite);
	CtauLegend->SetFillStyle(0);
	CtauLegend->SetTextFont(42);
	CtauLegend->SetTextSize(0.035);
	CtauLegend->SetBorderSize(0.);
	CtauLegend->AddEntry(legend_Tot,"sum","l");
	if(plotBackground) CtauLegend->AddEntry(legend_Background,"BG","l");
	if(plotChic0) CtauLegend->AddEntry(legend_Chic0,"#chi_{c0}","l");
	if(plotChic1){
		CtauLegend->AddEntry(legend_Chic1,"#chi_{c1}","l");
		CtauLegend->AddEntry(legend_ChicPR,"#chi^{PR}_{c1}, #chi^{NP}_{c1}","l");
	}
	if(plotChic2) CtauLegend->AddEntry(legend_Chic2,"#chi_{c2}","l");


	double left=0.7, top=0.885, textSize=0.03;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	double step=textSize*1.6;


	gStyle->SetPadBottomMargin(0.08); //0.12
	gStyle->SetPadLeftMargin(0.09); //0.12
	gStyle->SetPadRightMargin(0.035); //0.05
	gStyle->SetPadTopMargin(0.05); //0.05

	TCanvas *c1;
	TPad *pad1;
	TPad *pad2;

	ctauFrame->SetMinimum(minY);
	ctauFrame->SetMaximum(maxY);

	for(int LinLog=0; LinLog<2; LinLog++){
		cout<<"LinLog "<<LinLog<<endl;

		c1=new TCanvas("c1","",1000,900);

		c1->cd();
		pad1 = new TPad("pad1","pad1",0.,0.,1.,0.3);
		pad1->SetGridy();
		pad1->SetBottomMargin(0.2);
		pad1->Draw();
		c1->cd();
		pad2 = new TPad("pad2","pad2",0.,0.3,1.,1.);
		pad2->Draw();


		//Sig
		pad2->cd(0);
		if(LinLog==0) pad2->SetLogy(false);
		if(LinLog==1) pad2->SetLogy(true);


		ctauFrame->Draw(); CtauLegend->Draw("same");
		lineNPLow->Draw("same");
		linePRLow->Draw("same");
		linePRHigh->Draw("same");


		left=0.54; top=0.885; textSize=0.030; latex->SetTextSize(textSize);
		if(rapBin==0) latex->DrawLatex(left,top,Form("|y%s| < %.1f",onia::KinParticleChar,onia::rapForPTRange[onia::kNbRapForPTBins]));
		else if(rapBin==1) latex->DrawLatex(left,top,Form("|y%s| < %.1f",onia::KinParticleChar,onia::rapForPTRange[rapBin]));
		else latex->DrawLatex(left,top,Form("%.1f < |y%s| < %.1f",onia::rapForPTRange[rapBin-1],onia::KinParticleChar,onia::rapForPTRange[rapBin]));
		top-=step;
		if(ptBin==0)
			latex->DrawLatex(left,top,Form("%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin],onia::KinParticleChar,onia::pTRange[rapBin][onia::kNbPTBins[rapBin]]));
		else
			latex->DrawLatex(left,top,Form("%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin-1],onia::KinParticleChar,onia::pTRange[rapBin][ptBin]));

		top-=step;
		latex->SetTextColor(kRed);
		latex->DrawLatex(left,top,Form("J/#psi SR, #chi SR1"));
		latex->SetTextColor(kBlack);
		top-=step;
		latex->DrawLatex(left,top,Form("#chi^{2}/ndf = %.1f / %d", chi2_Ctau, ndof_Ctau));
		top-=step;


	    RooRealVar var_chi2ndf_Lifetime_SR1("var_chi2ndf_Lifetime_SR1","var_chi2ndf_Lifetime_SR1",double(chi2_Ctau/ndof_Ctau)); if(!ws->var("var_chi2ndf_Lifetime_SR1")) ws->import(var_chi2ndf_Lifetime_SR1); else ws->var("var_chi2ndf_Lifetime_SR1")->setVal(double(chi2_Ctau/ndof_Ctau));


		textSize=0.03; latex->SetTextSize(textSize);
		left=0.805; top=0.885;
		double stepsizeTimes=1.9;
		//latex->DrawLatex(left,top,Form("n^{SR1}_{tot} = %.0f", ws->var("var_ev")->getVal()*ws->var("var_fTotInSR1")->getVal()));
		//top-=textSize*stepsizeTimes;
		//latex->DrawLatex(left,top,Form("n^{PRSR1}_{tot} = %.0f", ws->var("var_ev")->getVal()*ws->var("var_fTotInSR1")->getVal()*ws->var("var_fPRSR1InSR1")->getVal()));
		//top-=textSize*stepsizeTimes;
		//latex->DrawLatex(left,top,Form("n^{NPSR1}_{tot} = %.0f", ws->var("var_ev")->getVal()*ws->var("var_fTotInSR1")->getVal()*ws->var("var_fNPSR1InSR1")->getVal()));
		//top-=textSize*stepsizeTimes;
		latex->SetTextColor(onia::colorPR);
		latex->SetTextSize(textSize);
		latex->DrawLatex(left,top,Form("f^{PRSR1}_{PR#chi_{c1}}  =  %.3f",ws->var("var_fracPRChic1InPRSR1")->getVal()));
		top-=textSize*stepsizeTimes;
		latex->DrawLatex(left,top,Form("f^{PRSR1}_{NP#chi_{c1}}  =  %.3f",ws->var("var_fracNPChic1InPRSR1")->getVal()));
		top-=textSize*stepsizeTimes;
		latex->DrawLatex(left,top,Form("f^{PRSR1}_{BG}  =  %.3f",ws->var("var_fracBackgroundInPRSR1")->getVal()));
		top-=textSize*stepsizeTimes;
		if(plotChic2){
			latex->DrawLatex(left,top,Form("f^{PRSR1}_{PR#chi_{c2}}  =  %.3f",ws->var("var_fracPRChic2InPRSR1")->getVal()));
			top-=textSize*stepsizeTimes;
			//latex->DrawLatex(left,top,Form("f^{PRSR1}_{NP#chi_{c2}}  =  %.3f",ws->var("var_fracNPChic2InPRSR1")->getVal()));
			//top-=textSize*stepsizeTimes;
		}
		latex->SetTextColor(onia::colorNP);
		latex->SetTextSize(textSize);
		latex->DrawLatex(left,top,Form("f^{NPSR1}_{NP#chi_{c1}}  =  %.3f",ws->var("var_fracNPChic1InNPSR1")->getVal()));
		top-=textSize*stepsizeTimes;
		latex->DrawLatex(left,top,Form("f^{NPSR1}_{PR#chi_{c1}}  =  %.3f",ws->var("var_fracPRChic1InNPSR1")->getVal()));
		top-=textSize*stepsizeTimes;
		latex->DrawLatex(left,top,Form("f^{NPSR1}_{BG}  =  %.3f",ws->var("var_fracBackgroundInNPSR1")->getVal()));
		top-=textSize*stepsizeTimes;
		//if(plotChic2){
		//	latex->DrawLatex(left,top,Form("f^{NPSR1}_{NP#chi_{c2}}  =  %.3f",ws->var("var_fracNPChic2InNPSR1")->getVal()));
		//	top-=textSize*stepsizeTimes;
		//	latex->DrawLatex(left,top,Form("f^{NPSR1}_{PR#chi_{c2}}  =  %.3f",ws->var("var_fracPRChic2InNPSR1")->getVal()));
		//	top-=textSize*stepsizeTimes;
		//}
		latex->SetTextColor(1);


		pad1->cd(0); pad1->SetLogy(0);
		ctauFramePull->Draw();

		textSize=0.06; latex->SetTextSize(textSize);
		left=0.385; top=0.875;
		latex->SetTextSize(textSize);
		//latex->DrawLatex(left,top,Form("#mu_{pull}  =  %.3f #pm %.3f,  #sigma_{pull}  =  %.3f #pm %.3f",pullMean, err_pullMean, pullRMS, err_pullRMS));

		c1->cd();

		std::stringstream saveCtau;
		if(!zoom){
			if(LinLog==0) saveCtau << "Fit/chicFit/ctau_SR1_lin_rap" << rapBin << "_pt" << ptBin << ".pdf";
			if(LinLog==1) saveCtau << "Fit/chicFit/ctau_SR1_log_rap" << rapBin << "_pt" << ptBin << ".pdf";
		}
		else{
			if(LinLog==0) saveCtau << "Fit/chicFit/ctau_SR1_zoom_lin_rap" << rapBin << "_pt" << ptBin << ".pdf";
			if(LinLog==1) saveCtau << "Fit/chicFit/ctau_SR1_zoom_log_rap" << rapBin << "_pt" << ptBin << ".pdf";
		}
		c1->SaveAs(saveCtau.str().c_str());


	}

	delete c1;
	delete legend_Tot;
	delete legend_Background;
	delete legend_Chic0;
	delete legend_Chic1;
	delete legend_Chic2;

	double returnChi2Ndf=double(chi2_Ctau/ndof_Ctau);

	return returnChi2Ndf;

}




//==============================================
double plotLifetimeSR2(RooWorkspace *ws, int rapBin, int ptBin, int nState, bool SpeedPlotting, int nSpeedPlotting, bool zoom){
	int nbins=onia::LifetimePlotBins;
	TGaxis::SetMaxDigits(3);

	double PlotMin, PlotMax;
	if(zoom){
		PlotMin=onia::ctPlotMinZoom;
		PlotMax=onia::ctPlotMaxZoom;
	}
	else{
		PlotMin=onia::ctPlotMin;
		PlotMax=onia::ctPlotMax;
	}

	double binWidth=(PlotMax-PlotMin)/double(nbins)*1000;

	bool correctResolutionForPlotting=false;
	double resCorrFactor=1.0525;
	if(ptBin==0) resCorrFactor=1.0825;
	if(ptBin==1) resCorrFactor=1.0725;
	if(correctResolutionForPlotting){
		ws->var("ctResolution")->setVal(ws->var("ctResolution")->getVal()*resCorrFactor);
		ws->var("ctResolution2")->setVal(ws->var("ctResolution2")->getVal()*resCorrFactor);
	}


	RooRealVar *chicMass = ws->var("chicMass");
	RooRealVar *Jpsict = ws->var("Jpsict");
	RooRealVar *JpsictErr = ws->var("JpsictErr");
	assert( 0 != Jpsict );

	RooPlot *ctauFrame = ((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins),Range(PlotMin, PlotMax));
	assert ( 0 != ctauFrame );
	ctauFrame->SetName(Form("ctau_plot_rap%d_pt%d",rapBin,ptBin));
	ctauFrame->SetTitle("");
	ctauFrame->GetYaxis()->SetTitle(Form("Events / %1.0f micron",binWidth));
	ctauFrame->GetYaxis()->SetTitleOffset(1.3);
	ctauFrame->GetXaxis()->SetRangeUser(PlotMin, PlotMax);

	RooPlot *ctauFramePull = ((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins),Range(PlotMin, PlotMax));
	assert ( 0 != ctauFramePull );
	ctauFramePull->SetName(Form("pullctau_plot_rap%d_pt%d",rapBin,ptBin));
	ctauFramePull->SetTitle("");
	ctauFramePull->GetYaxis()->SetTitle("pull");
	ctauFramePull->GetXaxis()->SetTitleSize(0.08);
	ctauFramePull->GetYaxis()->SetTitleSize(0.08);
	ctauFramePull->GetXaxis()->SetLabelSize(0.08);
	ctauFramePull->GetYaxis()->SetLabelSize(0.08);
	ctauFramePull->GetYaxis()->SetTitleOffset(0.4);
	ctauFramePull->GetYaxis()->SetRangeUser(-5.99,5.99);
	ctauFramePull->GetXaxis()->SetRangeUser(PlotMin, PlotMax);


	RooAbsData *data= ws->data(Form("data_rap%d_pt%d_SR",rapBin,ptBin));
	assert ( 0 != data );
	RooFitResult* fitRlt = dynamic_cast<RooFitResult*>(ws->obj(Form("fitresult_rap%d_pt%d",rapBin,ptBin)));
	assert ( 0 != fitRlt);

	fitRlt->Print();


    //calculate ctau ranges
    double sig1MaxMass = ws->var("var_sig1MaxMass")->getVal();
    double sig1MinMass = ws->var("var_sig1MinMass")->getVal();
    double sig2MaxMass = ws->var("var_sig2MaxMass")->getVal();
    double sig2MinMass = ws->var("var_sig2MinMass")->getVal();
    double lsbMaxMass = ws->var("var_lsbMaxMass")->getVal();
    double lsbMinMass = ws->var("var_lsbMinMass")->getVal();
    double rsbMaxMass = ws->var("var_rsbMaxMass")->getVal();
    double rsbMinMass = ws->var("var_rsbMinMass")->getVal();


    // define data in different regions
    std::stringstream cutSR2;
    cutSR2 << "chicMass > " << sig2MinMass << " && chicMass < " << sig2MaxMass;

    std::stringstream binNameSR2;
    binNameSR2  << "data_rap" << rapBin << "_pt" << ptBin << "_SR2";

    RooAbsData* dataSR2 = data->reduce(Cut(cutSR2.str().c_str()));

    RooAbsData* dataSR2Proj;
    if(SpeedPlotting){
    	//dataSR2Proj = dataSR2->reduce(EventRange(0,nSpeedPlotting));
    	dataSR2Proj = data->reduce(EventRange(0,nSpeedPlotting));
    }
    else{
    	//dataSR2Proj = dataSR2->reduce(EventRange(0,1e10));
    	dataSR2Proj = data->reduce(EventRange(0,1e10));
    }
    cout<<"number of events in dataSR2 = "<<dataSR2->numEntries()<<endl;
    cout<<"number of events in dataSR2Proj = "<<dataSR2Proj->numEntries()<<endl;

	RooAbsPdf *fullPdf = ws->pdf("ML_fullModel_SR2");
	assert ( 0 != fullPdf );

	int nEntries = dataSR2->numEntries();

	double minY = 0.;
	double maxY = 0.;
	double enlargeYby=onia::enlargeYby_ML;

	dataSR2->plotOn(ctauFrame,MarkerSize(onia::markerSize_ML), Name("myHist"));
	maxY = ctauFrame->GetMaximum()*enlargeYby;
	minY = 5e-1;
	ctauFrame->GetYaxis()->SetRangeUser(minY,maxY);

	cout<<"Plotting sum"<<endl;
	fullPdf->plotOn(ctauFrame,
			//Normalization(ws->var("var_ev")->getVal()*ws->var("var_fTotInSR2")->getVal(),2),
			LineWidth(onia::lineWidth_ML),
			ProjWData(*JpsictErr, *dataSR2Proj),
			NumCPU(1),
			Name("myCurve"));
	cout<<"Plotting sum finished"<<endl;

	//------get chi2------------
	int parsFit=(fitRlt->floatParsFinal()).getSize(); //this used the full p.d.f.
	int nBins_Ctau=ctauFrame->GetNbinsX();
	double chi2Pre_Ctau=ctauFrame->chiSquare(parsFit);  //reduced chi-squared = chi2/ndof
	int ndof_Ctau=nBins_Ctau-parsFit;  //num of degree of freedom
	double chi2_Ctau=chi2Pre_Ctau*ndof_Ctau;

	TH1F* pull = new TH1F("pull","pull distribution", 100,-10.,10.);
	gSystem->mkdir("Fit/root",kTRUE);gSystem->mkdir("Fit/chicFit",kTRUE);
	TFile *pullFile = new TFile(Form("Fit/root/pull_ctau_SR2_rap%d_pt%d.root",rapBin,ptBin),"RECREATE");

	RooHist* hpull_ctau = ctauFrame->pullHist("myHist","myCurve",kTRUE);
	hpull_ctau->SetMarkerSize(onia::markerSize_ML);
	for(int i=0;i<hpull_ctau->GetN();i++){
		hpull_ctau->SetPointEYlow(i,0.);
		hpull_ctau->SetPointEYhigh(i,0.);
		double x,y;
		hpull_ctau->GetPoint(i,x,y);
		pull->Fill(y);
	}
	pullFile->cd();
	pull->Write();
	pullFile->Close();

	double pullMean=pull->GetMean();
	double pullRMS=pull->GetRMS();
	double err_pullMean=pull->GetMeanError();
	double err_pullRMS=pull->GetRMSError();

	ctauFramePull->addPlotable(hpull_ctau,"P");


	cout<<"nBackgroundInSR2: "<<ws->var("var_ev")->getVal()*ws->var("var_fTotInSR2")->getVal()*ws->var("var_fracBackgroundInSR2")->getVal()<<endl;
	cout<<"nChic0InSR2: "<<ws->var("var_ev")->getVal()*ws->var("var_fTotInSR2")->getVal()*ws->var("var_fracChic0InSR2")->getVal()<<endl;
	cout<<"nChic1InSR2: "<<ws->var("var_ev")->getVal()*ws->var("var_fTotInSR2")->getVal()*ws->var("var_fracChic1InSR2")->getVal()<<endl;
	cout<<"nChic2InSR2: "<<ws->var("var_ev")->getVal()*ws->var("var_fTotInSR2")->getVal()*ws->var("var_fracChic2InSR2")->getVal()<<endl;

	double minCompFrac=0.03;

	bool plotBackground=false;
	bool plotChic0=false;
	bool plotChic1=false;
	bool plotChic2=false;
	if(ws->var("var_fracBackgroundInSR2")->getVal()>minCompFrac) plotBackground=true;
	if(ws->var("var_fracChic0InSR2")->getVal()>minCompFrac) plotChic0=true;
	if(ws->var("var_fracChic1InSR2")->getVal()>minCompFrac) plotChic1=true;
	if(ws->var("var_fracChic2InSR2")->getVal()>minCompFrac) plotChic2=true;

	if(plotBackground){
		cout<<"Plotting background"<<endl;
		fullPdf->plotOn(ctauFrame,
				Components("ML_background,ML_comb_background"),
				//Normalization(ws->var("var_ev")->getVal()*ws->var("var_fTotInSR2")->getVal(),2),
				LineStyle(onia::lineStyle_subComps_ML),
				LineColor(onia::colorBackground),
				LineWidth(onia::lineWidth_ML),
				ProjWData(*JpsictErr, *dataSR2Proj), NumCPU(1));
		cout<<"Plotting background finished"<<endl;
	}

	if(plotChic0){
		cout<<"Plotting Chic0"<<endl;
		fullPdf->plotOn(ctauFrame,
				Components("ML_chic0"),
				//Normalization(ws->var("var_ev")->getVal()*ws->var("var_fTotInSR2")->getVal(),2),
				LineStyle(onia::lineStyle_subComps_ML),
				LineColor(onia::colorChic0),
				LineWidth(onia::lineWidth_ML),
				ProjWData(*JpsictErr, *dataSR2Proj), NumCPU(1));
		cout<<"Plotting Chic0 finished"<<endl;
	}

	if(plotChic1){
		cout<<"Plotting Chic1"<<endl;
		fullPdf->plotOn(ctauFrame,
				Components("ML_chic1"),
				//Normalization(ws->var("var_ev")->getVal()*ws->var("var_fTotInSR2")->getVal(),2),
				LineStyle(onia::lineStyle_subComps_ML),
				LineColor(onia::colorChic1),
				LineWidth(onia::lineWidth_ML),
				ProjWData(*JpsictErr, *dataSR2Proj), NumCPU(1));
		cout<<"Plotting Chic1 finished"<<endl;
	}

	if(plotChic2){
		cout<<"Plotting Chic2"<<endl;
		fullPdf->plotOn(ctauFrame,
				Components("ML_chic2"),
				//Normalization(ws->var("var_ev")->getVal()*ws->var("var_fTotInSR2")->getVal(),2),
				LineStyle(onia::lineStyle_subComps_ML),
				LineColor(onia::colorChic2),
				LineWidth(onia::lineWidth_ML),
				ProjWData(*JpsictErr, *dataSR2Proj), NumCPU(1));
		cout<<"Plotting PR Chic2"<<endl;
		fullPdf->plotOn(ctauFrame,
				Components("ML_chic2_PR"),
				//Normalization(ws->var("var_ev")->getVal()*ws->var("var_fTotInSR2")->getVal(),2),
				LineStyle(onia::lineStyle_subCompPRsignal_ML),
				LineColor(onia::colorChic2),
				LineWidth(onia::lineWidth_PRsignal_ML),
				ProjWData(*JpsictErr, *dataSR2Proj), NumCPU(1));
		cout<<"Plotting NP Chic2"<<endl;
		fullPdf->plotOn(ctauFrame,
				Components("ML_chic2_NP"),
				//Normalization(ws->var("var_ev")->getVal()*ws->var("var_fTotInSR2")->getVal(),2),
				LineStyle(onia::lineStyle_subCompPRsignal_ML),
				LineColor(onia::colorChic2),
				LineWidth(onia::lineWidth_PRsignal_ML),
				ProjWData(*JpsictErr, *dataSR2Proj), NumCPU(1));
		cout<<"Plotting Chic2 finished"<<endl;
	}

	if(correctResolutionForPlotting){
		ws->var("ctResolution")->setVal(ws->var("ctResolution")->getVal()/resCorrFactor);
		ws->var("ctResolution2")->setVal(ws->var("ctResolution2")->getVal()/resCorrFactor);
	}


	double lineWidthRegions = 1.5;

	gStyle->SetLineStyleString(11,"12 50");

	int lineStyleRegions=11;

	TLine *linePRLow = new TLine(ws->var("var_PRMin")->getVal(), minY, ws->var("var_PRMin")->getVal(), maxY);
	TLine *linePRHigh = new TLine(ws->var("var_PRMax")->getVal(), minY, ws->var("var_PRMax")->getVal(), maxY);
	TLine *lineNPLow = new TLine(ws->var("var_NPMin")->getVal(), minY, ws->var("var_NPMin")->getVal(), maxY);
	TLine *lineNPHigh = new TLine(ws->var("var_NPMax")->getVal(), minY, ws->var("var_NPMax")->getVal(), maxY);

	linePRLow->SetLineWidth(lineWidthRegions);
	linePRLow->SetLineColor(onia::colorPR);
	linePRLow->SetLineStyle(lineStyleRegions);
	linePRHigh->SetLineWidth(lineWidthRegions);
	linePRHigh->SetLineColor(onia::colorPR);
	linePRHigh->SetLineStyle(lineStyleRegions);
	lineNPLow->SetLineWidth(lineWidthRegions);
	lineNPLow->SetLineColor(onia::colorNP);
	lineNPLow->SetLineStyle(lineStyleRegions);
	lineNPHigh->SetLineWidth(lineWidthRegions);
	lineNPHigh->SetLineColor(onia::colorNP);
	lineNPHigh->SetLineStyle(lineStyleRegions);


	if(ws->var("var_NPMin")->getVal()==ws->var("var_PRMax")->getVal())
		lineNPLow->SetLineWidth(lineWidthRegions*2);

	TH1* legend_Tot = data->createHistogram("legend_Tot",*Jpsict,Binning(50)) ; legend_Tot->SetLineColor(kBlue) ; legend_Tot->SetLineStyle(1) ; legend_Tot->SetLineWidth(2.) ;
	TH1* legend_Background = data->createHistogram("legend_Background",*Jpsict,Binning(50)) ; legend_Background->SetLineColor(onia::colorBackground) ; legend_Background->SetLineStyle(onia::lineStyle_subComps_ML) ; legend_Background->SetLineWidth(onia::lineWidth_ML) ;
	TH1* legend_Chic0 = data->createHistogram("legend_Chic0",*Jpsict,Binning(50)) ; legend_Chic0->SetLineColor(onia::colorChic0) ; legend_Chic0->SetLineStyle(onia::lineStyle_subComps_ML) ; legend_Chic0->SetLineWidth(onia::lineWidth_ML) ;
	TH1* legend_Chic1 = data->createHistogram("legend_Chic1",*Jpsict,Binning(50)) ; legend_Chic1->SetLineColor(onia::colorChic1) ; legend_Chic1->SetLineStyle(onia::lineStyle_subComps_ML) ; legend_Chic1->SetLineWidth(onia::lineWidth_ML) ;
	TH1* legend_Chic2 = data->createHistogram("legend_Chic2",*Jpsict,Binning(50)) ; legend_Chic2->SetLineColor(onia::colorChic2) ; legend_Chic2->SetLineStyle(onia::lineStyle_subComps_ML) ; legend_Chic2->SetLineWidth(onia::lineWidth_ML) ;

	TH1* legend_ChicPR = data->createHistogram("legend_ChicPR",*Jpsict,Binning(50)) ; legend_ChicPR->SetLineColor(onia::colorChic2) ; legend_ChicPR->SetLineStyle(onia::lineStyle_subCompPRsignal_ML) ; legend_ChicPR->SetLineWidth(onia::lineWidth_PRsignal_ML) ;

	TLegend* CtauLegend=new TLegend(0.13,0.71,0.23,0.91);
	CtauLegend->SetFillColor(kWhite);
	CtauLegend->SetFillStyle(0);
	CtauLegend->SetTextFont(42);
	CtauLegend->SetTextSize(0.035);
	CtauLegend->SetBorderSize(0.);
	CtauLegend->AddEntry(legend_Tot,"sum","l");
	if(plotBackground) CtauLegend->AddEntry(legend_Background,"BG","l");
	if(plotChic0) CtauLegend->AddEntry(legend_Chic0,"#chi_{c0}","l");
	if(plotChic1) CtauLegend->AddEntry(legend_Chic1,"#chi_{c1}","l");
	if(plotChic2){
		CtauLegend->AddEntry(legend_Chic2,"#chi_{c2}","l");
		CtauLegend->AddEntry(legend_ChicPR,"#chi^{PR}_{c2}, #chi^{NP}_{c2}","l");
	}


	double left=0.7, top=0.885, textSize=0.03;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	double step=textSize*1.6;


	gStyle->SetPadBottomMargin(0.08); //0.12
	gStyle->SetPadLeftMargin(0.09); //0.12
	gStyle->SetPadRightMargin(0.035); //0.05
	gStyle->SetPadTopMargin(0.05); //0.05

	TCanvas *c1;
	TPad *pad1;
	TPad *pad2;

	ctauFrame->SetMinimum(minY);
	ctauFrame->SetMaximum(maxY);

	for(int LinLog=0; LinLog<2; LinLog++){
		cout<<"LinLog "<<LinLog<<endl;

		c1=new TCanvas("c1","",1000,900);

		c1->cd();
		pad1 = new TPad("pad1","pad1",0.,0.,1.,0.3);
		pad1->SetGridy();
		pad1->SetBottomMargin(0.2);
		pad1->Draw();
		c1->cd();
		pad2 = new TPad("pad2","pad2",0.,0.3,1.,1.);
		pad2->Draw();


		//Sig
		pad2->cd(0);
		if(LinLog==0) pad2->SetLogy(false);
		if(LinLog==1) pad2->SetLogy(true);


		ctauFrame->Draw(); CtauLegend->Draw("same");
		lineNPLow->Draw("same");
		linePRLow->Draw("same");
		linePRHigh->Draw("same");


		left=0.54; top=0.885; textSize=0.030; latex->SetTextSize(textSize);
		if(rapBin==0) latex->DrawLatex(left,top,Form("|y%s| < %.1f",onia::KinParticleChar,onia::rapForPTRange[onia::kNbRapForPTBins]));
		else if(rapBin==1) latex->DrawLatex(left,top,Form("|y%s| < %.1f",onia::KinParticleChar,onia::rapForPTRange[rapBin]));
		else latex->DrawLatex(left,top,Form("%.1f < |y%s| < %.1f",onia::rapForPTRange[rapBin-1],onia::KinParticleChar,onia::rapForPTRange[rapBin]));
		top-=step;
		if(ptBin==0)
			latex->DrawLatex(left,top,Form("%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin],onia::KinParticleChar,onia::pTRange[rapBin][onia::kNbPTBins[rapBin]]));
		else
			latex->DrawLatex(left,top,Form("%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin-1],onia::KinParticleChar,onia::pTRange[rapBin][ptBin]));

		top-=step;
		latex->SetTextColor(kRed);
		latex->DrawLatex(left,top,Form("J/#psi SR, #chi SR2"));
		latex->SetTextColor(kBlack);
		top-=step;
		latex->DrawLatex(left,top,Form("#chi^{2}/ndf = %.1f / %d", chi2_Ctau, ndof_Ctau));
		top-=step;

	    RooRealVar var_chi2ndf_Lifetime_SR2("var_chi2ndf_Lifetime_SR2","var_chi2ndf_Lifetime_SR2",double(chi2_Ctau/ndof_Ctau)); if(!ws->var("var_chi2ndf_Lifetime_SR2")) ws->import(var_chi2ndf_Lifetime_SR2); else ws->var("var_chi2ndf_Lifetime_SR2")->setVal(double(chi2_Ctau/ndof_Ctau));


		textSize=0.03; latex->SetTextSize(textSize);
		left=0.805; top=0.885;
		double stepsizeTimes=1.9;
		latex->SetTextColor(onia::colorPR);
		latex->SetTextSize(textSize);
		latex->DrawLatex(left,top,Form("f^{PRSR2}_{PR#chi_{c2}}  =  %.3f",ws->var("var_fracPRChic2InPRSR2")->getVal()));
		top-=textSize*stepsizeTimes;
		latex->DrawLatex(left,top,Form("f^{PRSR2}_{NP#chi_{c2}}  =  %.3f",ws->var("var_fracNPChic2InPRSR2")->getVal()));
		top-=textSize*stepsizeTimes;
		latex->DrawLatex(left,top,Form("f^{PRSR2}_{BG}  =  %.3f",ws->var("var_fracBackgroundInPRSR2")->getVal()));
		top-=textSize*stepsizeTimes;
		if(plotChic1){
			latex->DrawLatex(left,top,Form("f^{PRSR2}_{PR#chi_{c1}}  =  %.3f",ws->var("var_fracPRChic1InPRSR2")->getVal()));
			top-=textSize*stepsizeTimes;
			latex->DrawLatex(left,top,Form("f^{PRSR2}_{NP#chi_{c1}}  =  %.3f",ws->var("var_fracNPChic1InPRSR2")->getVal()));
			top-=textSize*stepsizeTimes;
		}
		latex->SetTextColor(onia::colorNP);
		latex->SetTextSize(textSize);
		latex->DrawLatex(left,top,Form("f^{NPSR2}_{NP#chi_{c2}}  =  %.3f",ws->var("var_fracNPChic2InNPSR2")->getVal()));
		top-=textSize*stepsizeTimes;
		latex->DrawLatex(left,top,Form("f^{NPSR2}_{PR#chi_{c2}}  =  %.3f",ws->var("var_fracPRChic2InNPSR2")->getVal()));
		top-=textSize*stepsizeTimes;
		latex->DrawLatex(left,top,Form("f^{NPSR2}_{BG}  =  %.3f",ws->var("var_fracBackgroundInNPSR2")->getVal()));
		top-=textSize*stepsizeTimes;
		if(plotChic1){
			latex->DrawLatex(left,top,Form("f^{NPSR2}_{NP#chi_{c1}}  =  %.3f",ws->var("var_fracNPChic1InNPSR2")->getVal()));
			top-=textSize*stepsizeTimes;
			latex->DrawLatex(left,top,Form("f^{NPSR2}_{PR#chi_{c1}}  =  %.3f",ws->var("var_fracPRChic1InNPSR2")->getVal()));
			top-=textSize*stepsizeTimes;
		}
		latex->SetTextColor(1);


		pad1->cd(0); pad1->SetLogy(0);
		ctauFramePull->Draw();

		textSize=0.06; latex->SetTextSize(textSize);
		left=0.385; top=0.875;
		latex->SetTextSize(textSize);
		//latex->DrawLatex(left,top,Form("#mu_{pull}  =  %.3f #pm %.3f,  #sigma_{pull}  =  %.3f #pm %.3f",pullMean, err_pullMean, pullRMS, err_pullRMS));

		c1->cd();

		std::stringstream saveCtau;
		if(!zoom){
			if(LinLog==0) saveCtau << "Fit/chicFit/ctau_SR2_lin_rap" << rapBin << "_pt" << ptBin << ".pdf";
			if(LinLog==1) saveCtau << "Fit/chicFit/ctau_SR2_log_rap" << rapBin << "_pt" << ptBin << ".pdf";
		}
		else{
			if(LinLog==0) saveCtau << "Fit/chicFit/ctau_SR2_zoom_lin_rap" << rapBin << "_pt" << ptBin << ".pdf";
			if(LinLog==1) saveCtau << "Fit/chicFit/ctau_SR2_zoom_log_rap" << rapBin << "_pt" << ptBin << ".pdf";
		}
		c1->SaveAs(saveCtau.str().c_str());


	}

	delete c1;
	delete legend_Tot;
	delete legend_Background;
	delete legend_Chic0;
	delete legend_Chic1;
	delete legend_Chic2;

	double returnChi2Ndf=double(chi2_Ctau/ndof_Ctau);

	return returnChi2Ndf;

}



//==============================================
double plotLifetimeLSB(RooWorkspace *ws, int rapBin, int ptBin, int nState, bool SpeedPlotting, int nSpeedPlotting, bool zoom){
	int nbins=onia::LifetimePlotBins;
	TGaxis::SetMaxDigits(3);

	double PlotMin, PlotMax;
	if(zoom){
		PlotMin=onia::ctPlotMinZoom;
		PlotMax=onia::ctPlotMaxZoom;
	}
	else{
		PlotMin=onia::ctPlotMin;
		PlotMax=onia::ctPlotMax;
	}

	double binWidth=(PlotMax-PlotMin)/double(nbins)*1000;

	bool correctResolutionForPlotting=false;
	double resCorrFactor=1.0525;
	if(ptBin==0) resCorrFactor=1.0825;
	if(ptBin==1) resCorrFactor=1.0725;
	if(correctResolutionForPlotting){
		ws->var("ctResolution")->setVal(ws->var("ctResolution")->getVal()*resCorrFactor);
		ws->var("ctResolution2")->setVal(ws->var("ctResolution2")->getVal()*resCorrFactor);
	}


	RooRealVar *chicMass = ws->var("chicMass");
	RooRealVar *Jpsict = ws->var("Jpsict");
	RooRealVar *JpsictErr = ws->var("JpsictErr");
	assert( 0 != Jpsict );

	RooPlot *ctauFrame = ((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins),Range(PlotMin, PlotMax));
	assert ( 0 != ctauFrame );
	ctauFrame->SetName(Form("ctau_plot_rap%d_pt%d",rapBin,ptBin));
	ctauFrame->SetTitle("");
	ctauFrame->GetYaxis()->SetTitle(Form("Events / %1.0f micron",binWidth));
	ctauFrame->GetYaxis()->SetTitleOffset(1.3);
	ctauFrame->GetXaxis()->SetRangeUser(PlotMin, PlotMax);

	RooPlot *ctauFramePull = ((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins),Range(PlotMin, PlotMax));
	assert ( 0 != ctauFramePull );
	ctauFramePull->SetName(Form("pullctau_plot_rap%d_pt%d",rapBin,ptBin));
	ctauFramePull->SetTitle("");
	ctauFramePull->GetYaxis()->SetTitle("pull");
	ctauFramePull->GetXaxis()->SetTitleSize(0.08);
	ctauFramePull->GetYaxis()->SetTitleSize(0.08);
	ctauFramePull->GetXaxis()->SetLabelSize(0.08);
	ctauFramePull->GetYaxis()->SetLabelSize(0.08);
	ctauFramePull->GetYaxis()->SetTitleOffset(0.4);
	ctauFramePull->GetYaxis()->SetRangeUser(-5.99,5.99);
	ctauFramePull->GetXaxis()->SetRangeUser(PlotMin, PlotMax);


	RooAbsData *data= ws->data(Form("data_rap%d_pt%d_SR",rapBin,ptBin));
	assert ( 0 != data );
	RooFitResult* fitRlt = dynamic_cast<RooFitResult*>(ws->obj(Form("fitresult_rap%d_pt%d",rapBin,ptBin)));
	assert ( 0 != fitRlt);

	fitRlt->Print();


    //calculate ctau ranges
    double sig1MaxMass = ws->var("var_sig1MaxMass")->getVal();
    double sig1MinMass = ws->var("var_sig1MinMass")->getVal();
    double sig2MaxMass = ws->var("var_sig2MaxMass")->getVal();
    double sig2MinMass = ws->var("var_sig2MinMass")->getVal();
    double lsbMaxMass = ws->var("var_lsbMaxMass")->getVal();
    double lsbMinMass = ws->var("var_lsbMinMass")->getVal();
    double rsbMaxMass = ws->var("var_rsbMaxMass")->getVal();
    double rsbMinMass = ws->var("var_rsbMinMass")->getVal();


    // define data in different regions
    std::stringstream cutLSB;
    cutLSB << "chicMass > " << lsbMinMass << " && chicMass < " << lsbMaxMass;

    std::stringstream binNameLSB;
    binNameLSB  << "data_rap" << rapBin << "_pt" << ptBin << "_LSB";

    RooAbsData* dataLSB = data->reduce(Cut(cutLSB.str().c_str()));

    RooAbsData* dataLSBProj;
    if(SpeedPlotting){
    	//dataLSBProj = dataLSB->reduce(EventRange(0,nSpeedPlotting));
    	dataLSBProj = data->reduce(EventRange(0,nSpeedPlotting));
    }
    else{
    	//dataLSBProj = dataLSB->reduce(EventRange(0,1e10));
    	dataLSBProj = data->reduce(EventRange(0,1e10));
    }
    cout<<"number of events in dataLSB = "<<dataLSB->numEntries()<<endl;
    cout<<"number of events in dataLSBProj = "<<dataLSBProj->numEntries()<<endl;

	RooAbsPdf *fullPdf = ws->pdf("ML_fullModel_LSB");
	assert ( 0 != fullPdf );

	int nEntries = dataLSB->numEntries();

	double minY = 0.;
	double maxY = 0.;
	double enlargeYby=onia::enlargeYby_ML;

	dataLSB->plotOn(ctauFrame,MarkerSize(onia::markerSize_ML), Name("myHist"));
	maxY = ctauFrame->GetMaximum()*enlargeYby;
	minY = 5e-1;
	ctauFrame->GetYaxis()->SetRangeUser(minY,maxY);

	cout<<"Plotting sum"<<endl;
	fullPdf->plotOn(ctauFrame,
			//Normalization(ws->var("var_ev")->getVal()*ws->var("var_fTotInLSB")->getVal(),2),
			LineWidth(onia::lineWidth_ML),
			ProjWData(*JpsictErr, *dataLSBProj),
			NumCPU(1),
			Name("myCurve"));
	cout<<"Plotting sum finished"<<endl;

	//------get chi2------------
	int parsFit=(fitRlt->floatParsFinal()).getSize(); //this used the full p.d.f.
	int nBins_Ctau=ctauFrame->GetNbinsX();
	double chi2Pre_Ctau=ctauFrame->chiSquare(parsFit);  //reduced chi-squared = chi2/ndof
	int ndof_Ctau=nBins_Ctau-parsFit;  //num of degree of freedom
	double chi2_Ctau=chi2Pre_Ctau*ndof_Ctau;

	TH1F* pull = new TH1F("pull","pull distribution", 100,-10.,10.);
	gSystem->mkdir("Fit/root",kTRUE);gSystem->mkdir("Fit/chicFit",kTRUE);
	TFile *pullFile = new TFile(Form("Fit/root/pull_ctau_LSB_rap%d_pt%d.root",rapBin,ptBin),"RECREATE");

	RooHist* hpull_ctau = ctauFrame->pullHist("myHist","myCurve",kTRUE);
	hpull_ctau->SetMarkerSize(onia::markerSize_ML);
	for(int i=0;i<hpull_ctau->GetN();i++){
		hpull_ctau->SetPointEYlow(i,0.);
		hpull_ctau->SetPointEYhigh(i,0.);
		double x,y;
		hpull_ctau->GetPoint(i,x,y);
		pull->Fill(y);
	}
	pullFile->cd();
	pull->Write();
	pullFile->Close();

	double pullMean=pull->GetMean();
	double pullRMS=pull->GetRMS();
	double err_pullMean=pull->GetMeanError();
	double err_pullRMS=pull->GetRMSError();

	ctauFramePull->addPlotable(hpull_ctau,"P");


	cout<<"nBackgroundInLSB: "<<ws->var("var_ev")->getVal()*ws->var("var_fTotInLSB")->getVal()*ws->var("var_fracBackgroundInLSB")->getVal()<<endl;
	cout<<"nChic0InLSB: "<<ws->var("var_ev")->getVal()*ws->var("var_fTotInLSB")->getVal()*ws->var("var_fracChic0InLSB")->getVal()<<endl;
	cout<<"nChic1InLSB: "<<ws->var("var_ev")->getVal()*ws->var("var_fTotInLSB")->getVal()*ws->var("var_fracChic1InLSB")->getVal()<<endl;
	cout<<"nChic2InLSB: "<<ws->var("var_ev")->getVal()*ws->var("var_fTotInLSB")->getVal()*ws->var("var_fracChic2InLSB")->getVal()<<endl;


	double minCompFrac=0.03;

	bool plotBackground=false;
	bool plotChic0=false;
	bool plotChic1=false;
	bool plotChic2=false;
	if(ws->var("var_fracBackgroundInLSB")->getVal()>minCompFrac) plotBackground=true;
	if(ws->var("var_fracChic0InLSB")->getVal()>minCompFrac) plotChic0=true;
	if(ws->var("var_fracChic1InLSB")->getVal()>minCompFrac) plotChic1=true;
	if(ws->var("var_fracChic2InLSB")->getVal()>minCompFrac) plotChic2=true;



	if(plotBackground){
		cout<<"Plotting comb background"<<endl;
		fullPdf->plotOn(ctauFrame,
				Components("ML_comb_background"),
				//Normalization(ws->var("var_ev")->getVal()*ws->var("var_fTotInLSB")->getVal(),2),
				LineStyle(onia::lineStyle_subCompPRsignal_ML),
				LineColor(onia::colorBackground),
				LineWidth(onia::lineWidth_ML),
				ProjWData(*JpsictErr, *dataLSBProj), NumCPU(1));
		cout<<"Plotting comb background finished"<<endl;
		cout<<"Plotting background"<<endl;
		fullPdf->plotOn(ctauFrame,
				Components("ML_background"),
				//Normalization(ws->var("var_ev")->getVal()*ws->var("var_fTotInLSB")->getVal(),2),
				LineStyle(onia::lineStyle_subComps_ML),
				LineColor(onia::colorBackground),
				LineWidth(onia::lineWidth_ML),
				ProjWData(*JpsictErr, *dataLSBProj), NumCPU(1));
		cout<<"Plotting background finished"<<endl;
	}

	if(plotChic0){
		cout<<"Plotting Chic0"<<endl;
		fullPdf->plotOn(ctauFrame,
				Components("ML_chic0"),
				//Normalization(ws->var("var_ev")->getVal()*ws->var("var_fTotInLSB")->getVal(),2),
				LineStyle(onia::lineStyle_subComps_ML),
				LineColor(onia::colorChic0),
				LineWidth(onia::lineWidth_ML),
				ProjWData(*JpsictErr, *dataLSBProj), NumCPU(1));
		cout<<"Plotting Chic0 finished"<<endl;
	}

	if(plotChic1){
		cout<<"Plotting Chic1"<<endl;
		fullPdf->plotOn(ctauFrame,
				Components("ML_chic1"),
				//Normalization(ws->var("var_ev")->getVal()*ws->var("var_fTotInLSB")->getVal(),2),
				LineStyle(onia::lineStyle_subComps_ML),
				LineColor(onia::colorChic1),
				LineWidth(onia::lineWidth_ML),
				ProjWData(*JpsictErr, *dataLSBProj), NumCPU(1));
		cout<<"Plotting Chic1 finished"<<endl;
	}

	if(plotChic2){
		cout<<"Plotting Chic2"<<endl;
		fullPdf->plotOn(ctauFrame,
				Components("ML_chic2"),
				//Normalization(ws->var("var_ev")->getVal()*ws->var("var_fTotInLSB")->getVal(),2),
				LineStyle(onia::lineStyle_subComps_ML),
				LineColor(onia::colorChic2),
				LineWidth(onia::lineWidth_ML),
				ProjWData(*JpsictErr, *dataLSBProj), NumCPU(1));
		cout<<"Plotting Chic2 finished"<<endl;
	}

	if(correctResolutionForPlotting){
		ws->var("ctResolution")->setVal(ws->var("ctResolution")->getVal()/resCorrFactor);
		ws->var("ctResolution2")->setVal(ws->var("ctResolution2")->getVal()/resCorrFactor);
	}


	double lineWidthRegions = 1.5;

	gStyle->SetLineStyleString(11,"12 50");

	int lineStyleRegions=11;

	TLine *linePRLow = new TLine(ws->var("var_PRMin")->getVal(), minY, ws->var("var_PRMin")->getVal(), maxY);
	TLine *linePRHigh = new TLine(ws->var("var_PRMax")->getVal(), minY, ws->var("var_PRMax")->getVal(), maxY);
	TLine *lineNPLow = new TLine(ws->var("var_NPMin")->getVal(), minY, ws->var("var_NPMin")->getVal(), maxY);
	TLine *lineNPHigh = new TLine(ws->var("var_NPMax")->getVal(), minY, ws->var("var_NPMax")->getVal(), maxY);

	linePRLow->SetLineWidth(lineWidthRegions);
	linePRLow->SetLineColor(onia::colorPR);
	linePRLow->SetLineStyle(lineStyleRegions);
	linePRHigh->SetLineWidth(lineWidthRegions);
	linePRHigh->SetLineColor(onia::colorPR);
	linePRHigh->SetLineStyle(lineStyleRegions);
	lineNPLow->SetLineWidth(lineWidthRegions);
	lineNPLow->SetLineColor(onia::colorNP);
	lineNPLow->SetLineStyle(lineStyleRegions);
	lineNPHigh->SetLineWidth(lineWidthRegions);
	lineNPHigh->SetLineColor(onia::colorNP);
	lineNPHigh->SetLineStyle(lineStyleRegions);


	if(ws->var("var_NPMin")->getVal()==ws->var("var_PRMax")->getVal())
		lineNPLow->SetLineWidth(lineWidthRegions*2);

	TH1* legend_Tot = data->createHistogram("legend_Tot",*Jpsict,Binning(50)) ; legend_Tot->SetLineColor(kBlue) ; legend_Tot->SetLineStyle(1) ; legend_Tot->SetLineWidth(2.) ;
	TH1* legend_Background = data->createHistogram("legend_Background",*Jpsict,Binning(50)) ; legend_Background->SetLineColor(onia::colorBackground) ; legend_Background->SetLineStyle(onia::lineStyle_subComps_ML) ; legend_Background->SetLineWidth(onia::lineWidth_ML) ;
	TH1* legend_Comb_Background = data->createHistogram("legend_Comb_Background",*Jpsict,Binning(50)) ; legend_Comb_Background->SetLineColor(onia::colorBackground) ; legend_Comb_Background->SetLineStyle(onia::lineStyle_subCompPRsignal_ML) ; legend_Comb_Background->SetLineWidth(onia::lineWidth_ML) ;
	TH1* legend_Chic0 = data->createHistogram("legend_Chic0",*Jpsict,Binning(50)) ; legend_Chic0->SetLineColor(onia::colorChic0) ; legend_Chic0->SetLineStyle(onia::lineStyle_subComps_ML) ; legend_Chic0->SetLineWidth(onia::lineWidth_ML) ;
	TH1* legend_Chic1 = data->createHistogram("legend_Chic1",*Jpsict,Binning(50)) ; legend_Chic1->SetLineColor(onia::colorChic1) ; legend_Chic1->SetLineStyle(onia::lineStyle_subComps_ML) ; legend_Chic1->SetLineWidth(onia::lineWidth_ML) ;
	TH1* legend_Chic2 = data->createHistogram("legend_Chic2",*Jpsict,Binning(50)) ; legend_Chic2->SetLineColor(onia::colorChic2) ; legend_Chic2->SetLineStyle(onia::lineStyle_subComps_ML) ; legend_Chic2->SetLineWidth(onia::lineWidth_ML) ;

	TLegend* CtauLegend=new TLegend(0.13,0.71,0.23,0.91);
	CtauLegend->SetFillColor(kWhite);
	CtauLegend->SetFillStyle(0);
	CtauLegend->SetTextFont(42);
	CtauLegend->SetTextSize(0.035);
	CtauLegend->SetBorderSize(0.);
	CtauLegend->AddEntry(legend_Tot,"sum","l");
	if(plotBackground){
		CtauLegend->AddEntry(legend_Background,"J/#psi#gamma BG","l");
		CtauLegend->AddEntry(legend_Comb_Background,"#mu#mu#gamma BG","l");
	}
	if(plotChic0) CtauLegend->AddEntry(legend_Chic0,"#chi_{c0}","l");
	if(plotChic1) CtauLegend->AddEntry(legend_Chic1,"#chi_{c1}","l");
	if(plotChic2) CtauLegend->AddEntry(legend_Chic2,"#chi_{c2}","l");


	double left=0.7, top=0.885, textSize=0.03;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	double step=textSize*1.6;


	gStyle->SetPadBottomMargin(0.08); //0.12
	gStyle->SetPadLeftMargin(0.09); //0.12
	gStyle->SetPadRightMargin(0.035); //0.05
	gStyle->SetPadTopMargin(0.05); //0.05

	TCanvas *c1;
	TPad *pad1;
	TPad *pad2;

	ctauFrame->SetMinimum(minY);
	ctauFrame->SetMaximum(maxY);

	for(int LinLog=0; LinLog<2; LinLog++){
		cout<<"LinLog "<<LinLog<<endl;

		c1=new TCanvas("c1","",1000,900);

		c1->cd();
		pad1 = new TPad("pad1","pad1",0.,0.,1.,0.3);
		pad1->SetGridy();
		pad1->SetBottomMargin(0.2);
		pad1->Draw();
		c1->cd();
		pad2 = new TPad("pad2","pad2",0.,0.3,1.,1.);
		pad2->Draw();


		//Sig
		pad2->cd(0);
		if(LinLog==0) pad2->SetLogy(false);
		if(LinLog==1) pad2->SetLogy(true);


		ctauFrame->Draw(); CtauLegend->Draw("same");
		lineNPLow->Draw("same");
		linePRLow->Draw("same");
		linePRHigh->Draw("same");


		left=0.54; top=0.885; textSize=0.030; latex->SetTextSize(textSize);
		if(rapBin==0) latex->DrawLatex(left,top,Form("|y%s| < %.1f",onia::KinParticleChar,onia::rapForPTRange[onia::kNbRapForPTBins]));
		else if(rapBin==1) latex->DrawLatex(left,top,Form("|y%s| < %.1f",onia::KinParticleChar,onia::rapForPTRange[rapBin]));
		else latex->DrawLatex(left,top,Form("%.1f < |y%s| < %.1f",onia::rapForPTRange[rapBin-1],onia::KinParticleChar,onia::rapForPTRange[rapBin]));
		top-=step;
		if(ptBin==0)
			latex->DrawLatex(left,top,Form("%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin],onia::KinParticleChar,onia::pTRange[rapBin][onia::kNbPTBins[rapBin]]));
		else
			latex->DrawLatex(left,top,Form("%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin-1],onia::KinParticleChar,onia::pTRange[rapBin][ptBin]));

		top-=step;
		latex->SetTextColor(kRed);
		latex->DrawLatex(left,top,Form("J/#psi SR, #chi LSB"));
		latex->SetTextColor(kBlack);
		top-=step;
		latex->DrawLatex(left,top,Form("#chi^{2}/ndf = %.1f / %d", chi2_Ctau, ndof_Ctau));

	    RooRealVar var_chi2ndf_Lifetime_LSB("var_chi2ndf_Lifetime_LSB","var_chi2ndf_Lifetime_LSB",double(chi2_Ctau/ndof_Ctau)); if(!ws->var("var_chi2ndf_Lifetime_LSB")) ws->import(var_chi2ndf_Lifetime_LSB); else ws->var("var_chi2ndf_Lifetime_LSB")->setVal(double(chi2_Ctau/ndof_Ctau));


		textSize=0.03; latex->SetTextSize(textSize);
		left=0.775; top=0.885;
		double stepsizeTimes=1.9;
		latex->SetTextColor(onia::colorPR);
		latex->SetTextSize(textSize);
		latex->DrawLatex(left,top,Form("f^{PRLSB}_{BG}  =  %.3f",ws->var("var_fracBackgroundInPRLSB")->getVal()));
		top-=textSize*stepsizeTimes;
		if(plotChic0){
			latex->DrawLatex(left,top,Form("f^{PRLSB}_{PR#chi_{c0}}  =  %.3f",ws->var("var_fracPRChic0InPRLSB")->getVal()));
			top-=textSize*stepsizeTimes;
			latex->DrawLatex(left,top,Form("f^{PRLSB}_{NP#chi_{c0}}  =  %.3f",ws->var("var_fracNPChic0InPRLSB")->getVal()));
			top-=textSize*stepsizeTimes;
		}
		if(plotChic1){
			latex->DrawLatex(left,top,Form("f^{PRLSB}_{PR#chi_{c1}}  =  %.3f",ws->var("var_fracPRChic1InPRLSB")->getVal()));
			top-=textSize*stepsizeTimes;
			latex->DrawLatex(left,top,Form("f^{PRLSB}_{NP#chi_{c1}}  =  %.3f",ws->var("var_fracNPChic1InPRLSB")->getVal()));
			top-=textSize*stepsizeTimes;
		}
		if(plotChic2){
			latex->DrawLatex(left,top,Form("f^{PRLSB}_{PR#chi_{c2}}  =  %.3f",ws->var("var_fracPRChic2InPRLSB")->getVal()));
			top-=textSize*stepsizeTimes;
			latex->DrawLatex(left,top,Form("f^{PRLSB}_{NP#chi_{c2}}  =  %.3f",ws->var("var_fracNPChic2InPRLSB")->getVal()));
			top-=textSize*stepsizeTimes;
		}
		latex->SetTextColor(onia::colorNP);
		latex->SetTextSize(textSize);
		latex->DrawLatex(left,top,Form("f^{NPLSB}_{BG}  =  %.3f",ws->var("var_fracBackgroundInNPLSB")->getVal()));
		top-=textSize*stepsizeTimes;
		if(plotChic0){
			latex->DrawLatex(left,top,Form("f^{NPLSB}_{NP#chi_{c0}}  =  %.3f",ws->var("var_fracNPChic0InNPLSB")->getVal()));
			top-=textSize*stepsizeTimes;
			latex->DrawLatex(left,top,Form("f^{NPLSB}_{PR#chi_{c0}}  =  %.3f",ws->var("var_fracPRChic0InNPLSB")->getVal()));
			top-=textSize*stepsizeTimes;
		}
		if(plotChic1){
			latex->DrawLatex(left,top,Form("f^{NPLSB}_{NP#chi_{c1}}  =  %.3f",ws->var("var_fracNPChic1InNPLSB")->getVal()));
			top-=textSize*stepsizeTimes;
			latex->DrawLatex(left,top,Form("f^{NPLSB}_{PR#chi_{c1}}  =  %.3f",ws->var("var_fracPRChic1InNPLSB")->getVal()));
			top-=textSize*stepsizeTimes;
		}
		if(plotChic2){
			latex->DrawLatex(left,top,Form("f^{NPLSB}_{NP#chi_{c2}}  =  %.3f",ws->var("var_fracNPChic2InNPLSB")->getVal()));
			top-=textSize*stepsizeTimes;
			latex->DrawLatex(left,top,Form("f^{NPLSB}_{PR#chi_{c2}}  =  %.3f",ws->var("var_fracPRChic2InNPLSB")->getVal()));
			top-=textSize*stepsizeTimes;
		}
		latex->SetTextColor(1);


		pad1->cd(0); pad1->SetLogy(0);
		ctauFramePull->Draw();

		textSize=0.06; latex->SetTextSize(textSize);
		left=0.385; top=0.875;
		latex->SetTextSize(textSize);
		//latex->DrawLatex(left,top,Form("#mu_{pull}  =  %.3f #pm %.3f,  #sigma_{pull}  =  %.3f #pm %.3f",pullMean, err_pullMean, pullRMS, err_pullRMS));

		c1->cd();

		std::stringstream saveCtau;
		if(!zoom){
			if(LinLog==0) saveCtau << "Fit/chicFit/ctau_LSB_lin_rap" << rapBin << "_pt" << ptBin << ".pdf";
			if(LinLog==1) saveCtau << "Fit/chicFit/ctau_LSB_log_rap" << rapBin << "_pt" << ptBin << ".pdf";
		}
		else{
			if(LinLog==0) saveCtau << "Fit/chicFit/ctau_LSB_zoom_lin_rap" << rapBin << "_pt" << ptBin << ".pdf";
			if(LinLog==1) saveCtau << "Fit/chicFit/ctau_LSB_zoom_log_rap" << rapBin << "_pt" << ptBin << ".pdf";
		}
		c1->SaveAs(saveCtau.str().c_str());


	}

	delete c1;
	delete legend_Tot;
	delete legend_Background;
	delete legend_Chic0;
	delete legend_Chic1;
	delete legend_Chic2;

	double returnChi2Ndf=double(chi2_Ctau/ndof_Ctau);

	return returnChi2Ndf;

}




//==============================================
double plotLifetimeRSB(RooWorkspace *ws, int rapBin, int ptBin, int nState, bool SpeedPlotting, int nSpeedPlotting, bool zoom){
	int nbins=onia::LifetimePlotBins;
	TGaxis::SetMaxDigits(3);

	double PlotMin, PlotMax;
	if(zoom){
		PlotMin=onia::ctPlotMinZoom;
		PlotMax=onia::ctPlotMaxZoom;
	}
	else{
		PlotMin=onia::ctPlotMin;
		PlotMax=onia::ctPlotMax;
	}

	double binWidth=(PlotMax-PlotMin)/double(nbins)*1000;

	bool correctResolutionForPlotting=false;
	double resCorrFactor=1.0525;
	if(ptBin==0) resCorrFactor=1.0825;
	if(ptBin==1) resCorrFactor=1.0725;
	if(correctResolutionForPlotting){
		ws->var("ctResolution")->setVal(ws->var("ctResolution")->getVal()*resCorrFactor);
		ws->var("ctResolution2")->setVal(ws->var("ctResolution2")->getVal()*resCorrFactor);
	}


	RooRealVar *chicMass = ws->var("chicMass");
	RooRealVar *Jpsict = ws->var("Jpsict");
	RooRealVar *JpsictErr = ws->var("JpsictErr");
	assert( 0 != Jpsict );

	RooPlot *ctauFrame = ((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins),Range(PlotMin, PlotMax));
	assert ( 0 != ctauFrame );
	ctauFrame->SetName(Form("ctau_plot_rap%d_pt%d",rapBin,ptBin));
	ctauFrame->SetTitle("");
	ctauFrame->GetYaxis()->SetTitle(Form("Events / %1.0f micron",binWidth));
	ctauFrame->GetYaxis()->SetTitleOffset(1.3);
	ctauFrame->GetXaxis()->SetRangeUser(PlotMin, PlotMax);

	RooPlot *ctauFramePull = ((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins),Range(PlotMin, PlotMax));
	assert ( 0 != ctauFramePull );
	ctauFramePull->SetName(Form("pullctau_plot_rap%d_pt%d",rapBin,ptBin));
	ctauFramePull->SetTitle("");
	ctauFramePull->GetYaxis()->SetTitle("pull");
	ctauFramePull->GetXaxis()->SetTitleSize(0.08);
	ctauFramePull->GetYaxis()->SetTitleSize(0.08);
	ctauFramePull->GetXaxis()->SetLabelSize(0.08);
	ctauFramePull->GetYaxis()->SetLabelSize(0.08);
	ctauFramePull->GetYaxis()->SetTitleOffset(0.4);
	ctauFramePull->GetYaxis()->SetRangeUser(-5.99,5.99);
	ctauFramePull->GetXaxis()->SetRangeUser(PlotMin, PlotMax);


	RooAbsData *data= ws->data(Form("data_rap%d_pt%d_SR",rapBin,ptBin));
	assert ( 0 != data );
	RooFitResult* fitRlt = dynamic_cast<RooFitResult*>(ws->obj(Form("fitresult_rap%d_pt%d",rapBin,ptBin)));
	assert ( 0 != fitRlt);

	fitRlt->Print();


    //calculate ctau ranges
    double sig1MaxMass = ws->var("var_sig1MaxMass")->getVal();
    double sig1MinMass = ws->var("var_sig1MinMass")->getVal();
    double sig2MaxMass = ws->var("var_sig2MaxMass")->getVal();
    double sig2MinMass = ws->var("var_sig2MinMass")->getVal();
    double lsbMaxMass = ws->var("var_lsbMaxMass")->getVal();
    double lsbMinMass = ws->var("var_lsbMinMass")->getVal();
    double rsbMaxMass = ws->var("var_rsbMaxMass")->getVal();
    double rsbMinMass = ws->var("var_rsbMinMass")->getVal();


    // define data in different regions
    std::stringstream cutRSB;
    cutRSB << "chicMass > " << rsbMinMass << " && chicMass < " << rsbMaxMass;

    std::stringstream binNameRSB;
    binNameRSB  << "data_rap" << rapBin << "_pt" << ptBin << "_RSB";

    RooAbsData* dataRSB = data->reduce(Cut(cutRSB.str().c_str()));

    RooAbsData* dataRSBProj;
    if(SpeedPlotting){
    	//dataRSBProj = dataRSB->reduce(EventRange(0,nSpeedPlotting));
    	dataRSBProj = data->reduce(EventRange(0,nSpeedPlotting));
    }
    else{
    	//dataRSBProj = dataRSB->reduce(EventRange(0,1e10));
    	dataRSBProj = data->reduce(EventRange(0,1e10));
    }
    cout<<"number of events in dataRSB = "<<dataRSB->numEntries()<<endl;
    cout<<"number of events in dataRSBProj = "<<dataRSBProj->numEntries()<<endl;

	RooAbsPdf *fullPdf = ws->pdf("ML_fullModel_RSB");
	assert ( 0 != fullPdf );

	int nEntries = dataRSB->numEntries();

	double minY = 0.;
	double maxY = 0.;
	double enlargeYby=onia::enlargeYby_ML;

	dataRSB->plotOn(ctauFrame,MarkerSize(onia::markerSize_ML), Name("myHist"));
	maxY = ctauFrame->GetMaximum()*enlargeYby;
	minY = 5e-1;
	ctauFrame->GetYaxis()->SetRangeUser(minY,maxY);

	cout<<"Plotting sum"<<endl;
	fullPdf->plotOn(ctauFrame,
			//Normalization(ws->var("var_ev")->getVal()*ws->var("var_fTotInRSB")->getVal(),2),
			LineWidth(onia::lineWidth_ML),
			ProjWData(*JpsictErr, *dataRSBProj),
			NumCPU(1),
			Name("myCurve"));
	cout<<"Plotting sum finished"<<endl;

	//------get chi2------------
	int parsFit=(fitRlt->floatParsFinal()).getSize(); //this used the full p.d.f.
	int nBins_Ctau=ctauFrame->GetNbinsX();
	double chi2Pre_Ctau=ctauFrame->chiSquare(parsFit);  //reduced chi-squared = chi2/ndof
	int ndof_Ctau=nBins_Ctau-parsFit;  //num of degree of freedom
	double chi2_Ctau=chi2Pre_Ctau*ndof_Ctau;

	TH1F* pull = new TH1F("pull","pull distribution", 100,-10.,10.);
	gSystem->mkdir("Fit/root",kTRUE);gSystem->mkdir("Fit/chicFit",kTRUE);
	TFile *pullFile = new TFile(Form("Fit/root/pull_ctau_RSB_rap%d_pt%d.root",rapBin,ptBin),"RECREATE");

	RooHist* hpull_ctau = ctauFrame->pullHist("myHist","myCurve",kTRUE);
	hpull_ctau->SetMarkerSize(onia::markerSize_ML);
	for(int i=0;i<hpull_ctau->GetN();i++){
		hpull_ctau->SetPointEYlow(i,0.);
		hpull_ctau->SetPointEYhigh(i,0.);
		double x,y;
		hpull_ctau->GetPoint(i,x,y);
		pull->Fill(y);
	}
	pullFile->cd();
	pull->Write();
	pullFile->Close();

	double pullMean=pull->GetMean();
	double pullRMS=pull->GetRMS();
	double err_pullMean=pull->GetMeanError();
	double err_pullRMS=pull->GetRMSError();

	ctauFramePull->addPlotable(hpull_ctau,"P");


	cout<<"nBackgroundInRSB: "<<ws->var("var_ev")->getVal()*ws->var("var_fTotInRSB")->getVal()*ws->var("var_fracBackgroundInRSB")->getVal()<<endl;
	cout<<"nChic0InRSB: "<<ws->var("var_ev")->getVal()*ws->var("var_fTotInRSB")->getVal()*ws->var("var_fracChic0InRSB")->getVal()<<endl;
	cout<<"nChic1InRSB: "<<ws->var("var_ev")->getVal()*ws->var("var_fTotInRSB")->getVal()*ws->var("var_fracChic1InRSB")->getVal()<<endl;
	cout<<"nChic2InRSB: "<<ws->var("var_ev")->getVal()*ws->var("var_fTotInRSB")->getVal()*ws->var("var_fracChic2InRSB")->getVal()<<endl;

	double minCompFrac=0.03;

	bool plotBackground=false;
	bool plotChic0=false;
	bool plotChic1=false;
	bool plotChic2=false;
	if(ws->var("var_fracBackgroundInRSB")->getVal()>minCompFrac) plotBackground=true;
	if(ws->var("var_fracChic0InRSB")->getVal()>minCompFrac) plotChic0=true;
	if(ws->var("var_fracChic1InRSB")->getVal()>minCompFrac) plotChic1=true;
	if(ws->var("var_fracChic2InRSB")->getVal()>minCompFrac) plotChic2=true;

	cout<<"SumRSB = "<<ws->var("var_fracBackgroundInRSB")->getVal()+ws->var("var_fracChic0InRSB")->getVal()+ws->var("var_fracChic1InRSB")->getVal()+ws->var("var_fracChic2InRSB")->getVal()<<endl;
	cout<<"var_fracBackgroundInRSB = "<<ws->var("var_fracBackgroundInRSB")->getVal()<<endl;


	if(plotBackground){
		cout<<"Plotting comb background"<<endl;
		fullPdf->plotOn(ctauFrame,
				Components("ML_comb_background"),
				//Normalization(ws->var("var_ev")->getVal()*ws->var("var_fTotInRSB")->getVal(),2),
				LineStyle(onia::lineStyle_subCompPRsignal_ML),
				LineColor(onia::colorBackground),
				LineWidth(onia::lineWidth_ML),
				ProjWData(*JpsictErr, *dataRSBProj), NumCPU(1));
		cout<<"Plotting comb background finished"<<endl;
		cout<<"Plotting background"<<endl;
		fullPdf->plotOn(ctauFrame,
				Components("ML_background"),
				//Normalization(ws->var("var_ev")->getVal()*ws->var("var_fTotInRSB")->getVal(),2),
				LineStyle(onia::lineStyle_subComps_ML),
				LineColor(onia::colorBackground),
				LineWidth(onia::lineWidth_ML),
				ProjWData(*JpsictErr, *dataRSBProj), NumCPU(1));
		cout<<"Plotting background finished"<<endl;
	}

	if(plotChic0){
		cout<<"Plotting Chic0"<<endl;
		fullPdf->plotOn(ctauFrame,
				Components("ML_chic0"),
				//Normalization(ws->var("var_ev")->getVal()*ws->var("var_fTotInRSB")->getVal(),2),
				LineStyle(onia::lineStyle_subComps_ML),
				LineColor(onia::colorChic0),
				LineWidth(onia::lineWidth_ML),
				ProjWData(*JpsictErr, *dataRSBProj), NumCPU(1));
		cout<<"Plotting Chic0 finished"<<endl;
	}

	if(plotChic1){
		cout<<"Plotting Chic1"<<endl;
		fullPdf->plotOn(ctauFrame,
				Components("ML_chic1"),
				//Normalization(ws->var("var_ev")->getVal()*ws->var("var_fTotInRSB")->getVal(),2),
				LineStyle(onia::lineStyle_subComps_ML),
				LineColor(onia::colorChic1),
				LineWidth(onia::lineWidth_ML),
				ProjWData(*JpsictErr, *dataRSBProj), NumCPU(1));
		cout<<"Plotting Chic1 finished"<<endl;
	}

	if(plotChic2){
		cout<<"Plotting Chic2"<<endl;
		fullPdf->plotOn(ctauFrame,
				Components("ML_chic2"),
				//Normalization(ws->var("var_ev")->getVal()*ws->var("var_fTotInRSB")->getVal(),2),
				LineStyle(onia::lineStyle_subComps_ML),
				LineColor(onia::colorChic2),
				LineWidth(onia::lineWidth_ML),
				ProjWData(*JpsictErr, *dataRSBProj), NumCPU(1));
		cout<<"Plotting Chic2 finished"<<endl;
	}

	if(correctResolutionForPlotting){
		ws->var("ctResolution")->setVal(ws->var("ctResolution")->getVal()/resCorrFactor);
		ws->var("ctResolution2")->setVal(ws->var("ctResolution2")->getVal()/resCorrFactor);
	}


	double lineWidthRegions = 1.5;

	gStyle->SetLineStyleString(11,"12 50");

	int lineStyleRegions=11;

	TLine *linePRLow = new TLine(ws->var("var_PRMin")->getVal(), minY, ws->var("var_PRMin")->getVal(), maxY);
	TLine *linePRHigh = new TLine(ws->var("var_PRMax")->getVal(), minY, ws->var("var_PRMax")->getVal(), maxY);
	TLine *lineNPLow = new TLine(ws->var("var_NPMin")->getVal(), minY, ws->var("var_NPMin")->getVal(), maxY);
	TLine *lineNPHigh = new TLine(ws->var("var_NPMax")->getVal(), minY, ws->var("var_NPMax")->getVal(), maxY);

	linePRLow->SetLineWidth(lineWidthRegions);
	linePRLow->SetLineColor(onia::colorPR);
	linePRLow->SetLineStyle(lineStyleRegions);
	linePRHigh->SetLineWidth(lineWidthRegions);
	linePRHigh->SetLineColor(onia::colorPR);
	linePRHigh->SetLineStyle(lineStyleRegions);
	lineNPLow->SetLineWidth(lineWidthRegions);
	lineNPLow->SetLineColor(onia::colorNP);
	lineNPLow->SetLineStyle(lineStyleRegions);
	lineNPHigh->SetLineWidth(lineWidthRegions);
	lineNPHigh->SetLineColor(onia::colorNP);
	lineNPHigh->SetLineStyle(lineStyleRegions);


	if(ws->var("var_NPMin")->getVal()==ws->var("var_PRMax")->getVal())
		lineNPLow->SetLineWidth(lineWidthRegions*2);

	TH1* legend_Tot = data->createHistogram("legend_Tot",*Jpsict,Binning(50)) ; legend_Tot->SetLineColor(kBlue) ; legend_Tot->SetLineStyle(1) ; legend_Tot->SetLineWidth(2.) ;
	TH1* legend_Background = data->createHistogram("legend_Background",*Jpsict,Binning(50)) ; legend_Background->SetLineColor(onia::colorBackground) ; legend_Background->SetLineStyle(onia::lineStyle_subComps_ML) ; legend_Background->SetLineWidth(onia::lineWidth_ML) ;
	TH1* legend_Comb_Background = data->createHistogram("legend_Comb_Background",*Jpsict,Binning(50)) ; legend_Comb_Background->SetLineColor(onia::colorBackground) ; legend_Comb_Background->SetLineStyle(onia::lineStyle_subCompPRsignal_ML) ; legend_Comb_Background->SetLineWidth(onia::lineWidth_ML) ;
	TH1* legend_Chic0 = data->createHistogram("legend_Chic0",*Jpsict,Binning(50)) ; legend_Chic0->SetLineColor(onia::colorChic0) ; legend_Chic0->SetLineStyle(onia::lineStyle_subComps_ML) ; legend_Chic0->SetLineWidth(onia::lineWidth_ML) ;
	TH1* legend_Chic1 = data->createHistogram("legend_Chic1",*Jpsict,Binning(50)) ; legend_Chic1->SetLineColor(onia::colorChic1) ; legend_Chic1->SetLineStyle(onia::lineStyle_subComps_ML) ; legend_Chic1->SetLineWidth(onia::lineWidth_ML) ;
	TH1* legend_Chic2 = data->createHistogram("legend_Chic2",*Jpsict,Binning(50)) ; legend_Chic2->SetLineColor(onia::colorChic2) ; legend_Chic2->SetLineStyle(onia::lineStyle_subComps_ML) ; legend_Chic2->SetLineWidth(onia::lineWidth_ML) ;

	TLegend* CtauLegend=new TLegend(0.13,0.71,0.23,0.91);
	CtauLegend->SetFillColor(kWhite);
	CtauLegend->SetFillStyle(0);
	CtauLegend->SetTextFont(42);
	CtauLegend->SetTextSize(0.035);
	CtauLegend->SetBorderSize(0.);
	CtauLegend->AddEntry(legend_Tot,"sum","l");
	if(plotBackground){
		CtauLegend->AddEntry(legend_Background,"J/#psi#gamma BG","l");
		CtauLegend->AddEntry(legend_Comb_Background,"#mu#mu#gamma BG","l");
	}
	if(plotChic0) CtauLegend->AddEntry(legend_Chic0,"#chi_{c0}","l");
	if(plotChic1) CtauLegend->AddEntry(legend_Chic1,"#chi_{c1}","l");
	if(plotChic2) CtauLegend->AddEntry(legend_Chic2,"#chi_{c2}","l");


	double left=0.7, top=0.885, textSize=0.03;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	double step=textSize*1.6;


	gStyle->SetPadBottomMargin(0.08); //0.12
	gStyle->SetPadLeftMargin(0.09); //0.12
	gStyle->SetPadRightMargin(0.035); //0.05
	gStyle->SetPadTopMargin(0.05); //0.05

	TCanvas *c1;
	TPad *pad1;
	TPad *pad2;

	ctauFrame->SetMinimum(minY);
	ctauFrame->SetMaximum(maxY);

	for(int LinLog=0; LinLog<2; LinLog++){
		cout<<"LinLog "<<LinLog<<endl;

		c1=new TCanvas("c1","",1000,900);

		c1->cd();
		pad1 = new TPad("pad1","pad1",0.,0.,1.,0.3);
		pad1->SetGridy();
		pad1->SetBottomMargin(0.2);
		pad1->Draw();
		c1->cd();
		pad2 = new TPad("pad2","pad2",0.,0.3,1.,1.);
		pad2->Draw();


		//Sig
		pad2->cd(0);
		if(LinLog==0) pad2->SetLogy(false);
		if(LinLog==1) pad2->SetLogy(true);


		ctauFrame->Draw(); CtauLegend->Draw("same");
		lineNPLow->Draw("same");
		linePRLow->Draw("same");
		linePRHigh->Draw("same");


		left=0.54; top=0.885; textSize=0.030; latex->SetTextSize(textSize);
		if(rapBin==0) latex->DrawLatex(left,top,Form("|y%s| < %.1f",onia::KinParticleChar,onia::rapForPTRange[onia::kNbRapForPTBins]));
		else if(rapBin==1) latex->DrawLatex(left,top,Form("|y%s| < %.1f",onia::KinParticleChar,onia::rapForPTRange[rapBin]));
		else latex->DrawLatex(left,top,Form("%.1f < |y%s| < %.1f",onia::rapForPTRange[rapBin-1],onia::KinParticleChar,onia::rapForPTRange[rapBin]));
		top-=step;
		if(ptBin==0)
			latex->DrawLatex(left,top,Form("%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin],onia::KinParticleChar,onia::pTRange[rapBin][onia::kNbPTBins[rapBin]]));
		else
			latex->DrawLatex(left,top,Form("%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin-1],onia::KinParticleChar,onia::pTRange[rapBin][ptBin]));

		top-=step;
		latex->SetTextColor(kRed);
		latex->DrawLatex(left,top,Form("J/#psi SR, #chi RSB"));
		latex->SetTextColor(kBlack);
		top-=step;
		latex->DrawLatex(left,top,Form("#chi^{2}/ndf = %.1f / %d", chi2_Ctau, ndof_Ctau));
		top-=step;

	    RooRealVar var_chi2ndf_Lifetime_RSB("var_chi2ndf_Lifetime_RSB","var_chi2ndf_Lifetime_RSB",double(chi2_Ctau/ndof_Ctau)); if(!ws->var("var_chi2ndf_Lifetime_RSB")) ws->import(var_chi2ndf_Lifetime_RSB); else ws->var("var_chi2ndf_Lifetime_RSB")->setVal(double(chi2_Ctau/ndof_Ctau));


		textSize=0.03; latex->SetTextSize(textSize);
		left=0.775; top=0.885;
		double stepsizeTimes=1.9;
		latex->SetTextColor(onia::colorPR);
		latex->SetTextSize(textSize);
		latex->DrawLatex(left,top,Form("f^{PRRSB}_{BG}  =  %.3f",ws->var("var_fracBackgroundInPRRSB")->getVal()));
		top-=textSize*stepsizeTimes;
		if(plotChic2){
			latex->DrawLatex(left,top,Form("f^{PRRSB}_{PR#chi_{c2}}  =  %.3f",ws->var("var_fracPRChic2InPRRSB")->getVal()));
			top-=textSize*stepsizeTimes;
			latex->DrawLatex(left,top,Form("f^{PRRSB}_{NP#chi_{c2}}  =  %.3f",ws->var("var_fracNPChic2InPRRSB")->getVal()));
			top-=textSize*stepsizeTimes;
		}
		latex->SetTextColor(onia::colorNP);
		latex->SetTextSize(textSize);
		latex->DrawLatex(left,top,Form("f^{NPRSB}_{BG}  =  %.3f",ws->var("var_fracBackgroundInNPRSB")->getVal()));
		top-=textSize*stepsizeTimes;
		if(plotChic2){
			latex->DrawLatex(left,top,Form("f^{NPRSB}_{NP#chi_{c2}}  =  %.3f",ws->var("var_fracNPChic2InNPRSB")->getVal()));
			top-=textSize*stepsizeTimes;
			latex->DrawLatex(left,top,Form("f^{NPRSB}_{PR#chi_{c2}}  =  %.3f",ws->var("var_fracPRChic2InNPRSB")->getVal()));
			top-=textSize*stepsizeTimes;
		}
		latex->SetTextColor(1);


		pad1->cd(0); pad1->SetLogy(0);
		ctauFramePull->Draw();

		textSize=0.06; latex->SetTextSize(textSize);
		left=0.385; top=0.875;
		latex->SetTextSize(textSize);
		//latex->DrawLatex(left,top,Form("#mu_{pull}  =  %.3f #pm %.3f,  #sigma_{pull}  =  %.3f #pm %.3f",pullMean, err_pullMean, pullRMS, err_pullRMS));

		c1->cd();

		std::stringstream saveCtau;
		if(!zoom){
			if(LinLog==0) saveCtau << "Fit/chicFit/ctau_RSB_lin_rap" << rapBin << "_pt" << ptBin << ".pdf";
			if(LinLog==1) saveCtau << "Fit/chicFit/ctau_RSB_log_rap" << rapBin << "_pt" << ptBin << ".pdf";
		}
		else{
			if(LinLog==0) saveCtau << "Fit/chicFit/ctau_RSB_zoom_lin_rap" << rapBin << "_pt" << ptBin << ".pdf";
			if(LinLog==1) saveCtau << "Fit/chicFit/ctau_RSB_zoom_log_rap" << rapBin << "_pt" << ptBin << ".pdf";
		}
		c1->SaveAs(saveCtau.str().c_str());


	}

	delete c1;
	delete legend_Tot;
	delete legend_Background;
	delete legend_Chic0;
	delete legend_Chic1;
	delete legend_Chic2;

	double returnChi2Ndf=double(chi2_Ctau/ndof_Ctau);

	return returnChi2Ndf;

}




void latexFloatingLifetimePars(RooWorkspace *ws, TLatex* latex){


	double textSize, left, top;

	textSize=0.03; latex->SetTextSize(textSize);
	left=0.725; top=0.885;
	double stepsizeTimes=1.9;
	latex->SetTextSize(textSize);
	if(!ws->var("ctResolution")->isConstant()){
	latex->DrawLatex(left,top,Form("#sigma^{scale}_{l}  =  %.3f #pm %.3f",ws->var("ctResolution")->getVal(), ws->var("ctResolution")->getError()));
	top-=textSize*stepsizeTimes;
	}
	if(!ws->var("ctResolution2")->isConstant()){
	latex->DrawLatex(left,top,Form("#sigma^{scale2}_{l}  =  %.3f #pm %.3f",ws->var("ctResolution2")->getVal(), ws->var("ctResolution2")->getError()));
	top-=textSize*stepsizeTimes;
	}
	if(!ws->var("fracGauss2")->isConstant()){
	latex->DrawLatex(left,top,Form("f_{G_{2}}  =  %.3f #pm %.3f",ws->var("fracGauss2")->getVal(), ws->var("fracGauss2")->getError()));
	top-=textSize*stepsizeTimes;
	}
	if(!ws->var("NP_TauChic")->isConstant()){
	latex->DrawLatex(left,top,Form("#tau^{NP}_{#chi_{cJ}}  =  %.3f #pm %.3f mm",ws->var("NP_TauChic")->getVal(), ws->var("NP_TauChic")->getError()));
	top-=textSize*stepsizeTimes;
	}
	if(!ws->var("NP_TauBkg")->isConstant()){
	latex->DrawLatex(left,top,Form("#tau_{#psiBG}  =  %.3f #pm %.3f mm",ws->var("NP_TauBkg")->getVal(), ws->var("NP_TauBkg")->getError()));
	top-=textSize*stepsizeTimes;
	}
	//if(ws->var("FD_TauBkg")!=NULL)
	//if(!ws->var("FD_TauBkg")->isConstant()){
	//latex->DrawLatex(left,top,Form("#tau^{LS}_{BG}  =  %.3f #pm %.3f mm",ws->var("FD_TauBkg")->getVal(), ws->var("FD_TauBkg")->getError()));
	//top-=textSize*stepsizeTimes;
	//}
	//if(ws->var("DSD_TauBkg")!=NULL)
	//if(!ws->var("DSD_TauBkg")->isConstant()){
	//latex->DrawLatex(left,top,Form("#tau^{DS}_{BG}  =  %.3f #pm %.3f mm",ws->var("DSD_TauBkg")->getVal(), ws->var("DSD_TauBkg")->getError()));
	//top-=textSize*stepsizeTimes;
	//}
	if(!ws->var("fBkgNP")->isConstant()){
	latex->DrawLatex(left,top,Form("f^{#psiBG}_{NP}  =  %.3f #pm %.3f",ws->var("fBkgNP")->getVal(), ws->var("fBkgNP")->getError()));
	top-=textSize*stepsizeTimes;
	}
	if(ws->var("fBkgFD")!=NULL)
	if(!ws->var("fBkgFD")->isConstant()){
	latex->DrawLatex(left,top,Form("f^{BG}_{LS}  =  %.3f #pm %.3f",ws->var("fBkgFD")->getVal(), ws->var("fBkgFD")->getError()));
	top-=textSize*stepsizeTimes;
	}
	if(!ws->var("fracNP_chic0")->isConstant()){
	latex->DrawLatex(left,top,Form("f^{#chi_{c0}}_{NP}  =  %.3f #pm %.3f",ws->var("fracNP_chic0")->getVal(), ws->var("fracNP_chic0")->getError()));
	top-=textSize*stepsizeTimes;
	}
	if(!ws->var("fracNP_chic1")->isConstant()){
	latex->DrawLatex(left,top,Form("f^{#chi_{c1}}_{NP}  =  %.3f #pm %.3f",ws->var("fracNP_chic1")->getVal(), ws->var("fracNP_chic1")->getError()));
	top-=textSize*stepsizeTimes;
	}
	if(!ws->var("fracNP_chic2")->isConstant()){
	latex->DrawLatex(left,top,Form("f^{#chi_{c2}}_{NP}  =  %.3f #pm %.3f",ws->var("fracNP_chic2")->getVal(), ws->var("fracNP_chic2")->getError()));
	top-=textSize*stepsizeTimes;
	}


return;

}
