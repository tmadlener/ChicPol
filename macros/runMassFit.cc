#ifndef __CINT__
#endif

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include "TLatex.h"
using namespace std;

#include "massFit.cc"

Double_t paramMassRapParabola(Double_t *x, Double_t *par);

//========================================================
// code to read input arguments 
template<typename T>
void fromSplit(const std::string& key, const std::string &arg, T& out)
{
  const char delim = '=';
  // Skip if key or delimiter not there                                                                                                                                            
  if ((arg.find(key) == std::string::npos) ||
      (arg.find(delim) == std::string::npos))
    return;

  std::string skey, sval;
  std::stringstream sstr(arg);
  std::getline(sstr, skey, delim); // Dummy read to skip key                                                                                                                       
  std::getline(sstr, sval, delim); // Get value                                                                                                                                    
  T tout;
  if (!(std::istringstream(sval) >> std::boolalpha >> tout))
    return;
  out = tout;
  std::cout << std::boolalpha << skey << ": "  << out << std::endl;
}

// Special version for string without the conversion                                                                                                                                 
template<>
void fromSplit(const std::string& key, const std::string &arg, std::string &out)
{
  const char delim = '=';
  // Skip if key or delimiter not there                                                                                                                                            
  if ((arg.find(key) == std::string::npos) ||
      (arg.find(delim) == std::string::npos))
    return;
  std::string skey, sval;
  std::stringstream sstr(arg);
  std::getline(sstr, skey, delim); // Dummy read to skip key                                                                                                                       
  std::getline(sstr, sval, delim); // Get value                                                                                                                                    
  out = sval;
  std::cout << skey << ": "  << out << std::endl;
}


//===================================================
int main(int argc, char* argv[]){

  // Set defaults
    int rapMin = 999,
			 	rapMax = 999,
			 	ptMin = 999,
			 	ptMax = 999,
			 	nState = 999;

		bool fitMassPR = false,
				 fitMassNP = false,
				 MC = false;
	

    // Loop over argument list
    for (int i=1; i < argc; i++)
      {
	std::string arg = argv[i];
        fromSplit("rapMin", arg, rapMin);
        fromSplit("rapMax", arg, rapMax);
        fromSplit("ptMin", arg, ptMin);
        fromSplit("ptMax", arg, ptMax);
        fromSplit("nState", arg, nState);
        fromSplit("fitMassPR", arg, fitMassPR);
        fromSplit("fitMassNP", arg, fitMassNP);
        fromSplit("MC", arg, MC);
      }

    std::cout << "-----------------------\n"
	      << "Fitting mass for \n"
	      << "y bins " << rapMin << " - " << rapMax << "\n"
	      << "and pT bins "  << ptMin << " - " << ptMax << "\n"
	      << "-----------------------" << std::endl;

    const int nRap=rapMax-rapMin+1;
    const int nPt=ptMax-ptMin+1;

    double rapCenter_[nRap][nPt];
    double errrapCenter_[nRap][nPt];
    double MassRap_[nRap][nPt];
    double errMassRap_[nRap][nPt];
    double MeanMassRap_[nRap][nPt];
    double errMeanMassRap_[nRap][nPt];
    double AlphaMassRap_[nRap][nPt];
    double errAlphaMassRap_[nRap][nPt];

    int iRapCount=0;
    int iPtCount;
    for(int iRap = rapMin; iRap <= rapMax; iRap++){
  	  iPtCount=0;
      for(int iPT = ptMin; iPT <= ptMax; iPT++){
  		std::stringstream tempFrom;
  		tempFrom << "tmpFiles/backupWorkSpace/ws_createWorkspace_Chi_rap" << iRap << "_pt" << iPT << ".root";
  		const std::string infilenameFrom = tempFrom.str().c_str();

  		std::stringstream tempTo;
  		tempTo << "tmpFiles/backupWorkSpace/ws_MassFit_Jpsi_rap" << iRap << "_pt" << iPT << ".root";
  		const std::string infilenameTo = tempTo.str().c_str();

		gSystem->CopyFile(infilenameFrom.c_str(),infilenameTo.c_str(),kTRUE);

		massFit(infilenameTo.c_str(), iRap, iPT, nState, fitMassPR, fitMassNP, MC);


	    TFile *infile = new TFile(infilenameTo.c_str(), "UPDATE");
	    if(!infile){
	        std::cout << "Error: failed to open file with dataset" << std::endl;
	    }
	    RooWorkspace *ws=(RooWorkspace *)infile->Get("ws_masslifetime");
	    if(!ws){
	        std::cout << "Error: failed to open workspace " << std::endl;
	    }


  		rapCenter_[iRapCount][iPtCount]=(onia::rapForPTRange[iRap-1]+onia::rapForPTRange[iRap])/2.;
  		errrapCenter_[iRapCount][iPtCount]=(onia::rapForPTRange[iRap]-onia::rapForPTRange[iRap-1])/2.;
 		MassRap_[iRapCount][iPtCount]=ws->var("CBsigma_p0_jpsi")->getVal();
  		errMassRap_[iRapCount][iPtCount]=ws->var("CBsigma_p0_jpsi")->getError();
 		MeanMassRap_[iRapCount][iPtCount]=ws->var("CBmass_p0_jpsi")->getVal();
  		errMeanMassRap_[iRapCount][iPtCount]=ws->var("CBmass_p0_jpsi")->getError();
 		AlphaMassRap_[iRapCount][iPtCount]=ws->var("CBalpha_p0_jpsi")->getVal();
  		errAlphaMassRap_[iRapCount][iPtCount]=ws->var("CBalpha_p0_jpsi")->getError();

  		cout<<"iRap "<<iRap<<", iPT "<<iPT<<" -> sigma = "<<MassRap_[iRapCount][iPtCount]<<endl;
  		cout<<"iRap "<<iRap<<", iPT "<<iPT<<" -> mean = "<<MeanMassRap_[iRapCount][iPtCount]<<endl;
  		cout<<"iRap "<<iRap<<", iPT "<<iPT<<" -> alpha = "<<AlphaMassRap_[iRapCount][iPtCount]<<endl;

  		infile->Close();

	    iPtCount++;
      }
      iRapCount++;
    }

    bool fitRapDependence=false;

    double rapCenter[nRap];
    double errrapCenter[nRap];
    double MassRap[nRap];
    double errMassRap[nRap];
    double MeanMassRap[nRap];
    double errMeanMassRap[nRap];
    double AlphaMassRap[nRap];
    double errAlphaMassRap[nRap];


	iPtCount=0;
    for(int iPT = ptMin; iPT <= ptMax; iPT++){


       iRapCount=0;
       for(int iRap = rapMin; iRap <= rapMax; iRap++){
    	   rapCenter[iRapCount]=rapCenter_[iRapCount][iPtCount];
    	   errrapCenter[iRapCount]=errrapCenter_[iRapCount][iPtCount];
    	   MassRap[iRapCount]=MassRap_[iRapCount][iPtCount];
    	   errMassRap[iRapCount]=errMassRap_[iRapCount][iPtCount];
    	   MeanMassRap[iRapCount]=MeanMassRap_[iRapCount][iPtCount];
    	   errMeanMassRap[iRapCount]=errMeanMassRap_[iRapCount][iPtCount];
    	   AlphaMassRap[iRapCount]=AlphaMassRap_[iRapCount][iPtCount];
    	   errAlphaMassRap[iRapCount]=errAlphaMassRap_[iRapCount][iPtCount];

    	   iRapCount++;
       }


  //  int param_nRap=100;
    double param_rapMin=onia::rapForPTRange[0];
    double param_rapMax=onia::rapForPTRange[onia::kNbRapForPTBins+1];

  	TGraphErrors *graphMassRap = new TGraphErrors(nRap,rapCenter,MassRap,errrapCenter,errMassRap);
  	TGraphErrors *graphMeanMassRap = new TGraphErrors(nRap,rapCenter,MeanMassRap,errrapCenter,errMeanMassRap);
  	TGraphErrors *graphAlphaMassRap = new TGraphErrors(nRap,rapCenter,AlphaMassRap,errrapCenter,errAlphaMassRap);

  	graphMassRap->SetMarkerStyle(20);
  	graphMassRap->SetMarkerColor(kRed);
  	graphMeanMassRap->SetMarkerStyle(20);
  	graphMeanMassRap->SetMarkerColor(kGreen+2);
  	graphAlphaMassRap->SetMarkerStyle(20);
  	graphAlphaMassRap->SetMarkerColor(kMagenta);

  	char name[200];

  	sprintf(name, "fParabola");
  	TF1* fParabola = new TF1(name, paramMassRapParabola, 0., 1.5, 3);
  	double a=0.02;
  	double b=0.005;
  	double c=0.02;
  	fParabola->SetParameter(0,a);
  	fParabola->SetParameter(1,b);
  	fParabola->SetParameter(2,c);
  	if(fitRapDependence) graphMassRap->Fit(fParabola, "0", "", 0., 1.5);

  	fParabola->SetLineColor(kBlue);
  	fParabola->SetLineWidth(1.);

  	double a_res=fParabola->GetParameter(0);
  	double a_err=fParabola->GetParError(0);
  	double b_res=fParabola->GetParameter(1);
  	double b_err=fParabola->GetParError(1);
  	double c_res=fParabola->GetParameter(2);
  	double c_err=fParabola->GetParError(2);

  			TLegend* DMLegend=new TLegend(0.2,0.785,0.6,0.925);
  			DMLegend->SetFillColor(0);
  	//		DMLegend->SetTextFont(72);
  			DMLegend->SetTextSize(0.0245);
  			DMLegend->SetBorderSize(0);
  	//		DMLegend->SetMargin(0.135);
  			char DMLegendEntry[200];



  		  char DMyTitle[200];
  		  sprintf(DMyTitle,"dimuon mass resolution [GeV]");

  		  TCanvas* c2 = new TCanvas("c2","c2",1200,1100);
  		  gStyle->SetPalette(1);
  	 	  gPad->SetFillColor(kWhite);

  	 	  double MargLeft=0.125;//0.175;
  	 	  double MargRight=0.025;//0.1;
  	 	  double MargTop=0.05;//0.175;

  	 	  double MarkerSizeTom=1.75;

  	 	  gPad->SetLeftMargin(MargLeft);
  	      gPad->SetRightMargin(MargRight);
  	      gPad->SetTopMargin(MargTop);

  	      //TH1F* histAxis = new TH1F("","",100,param_rapMin, param_rapMax);

  		   TH1F *histAxis = new TH1F;

  			  histAxis->SetTitle(0);
  			  histAxis->SetStats(0);
  			  histAxis->SetMarkerStyle(25);
  			  histAxis->SetMarkerSize(MarkerSizeTom);
  			  histAxis->SetMarkerColor(kBlack);
  			  histAxis->SetLineColor(kBlack);
  			  histAxis->SetTitle(0);
  			  histAxis->SetStats(0);
  			  histAxis->SetMarkerColor(kRed);
  			  histAxis->SetMarkerStyle(20);
  			  histAxis->SetLineColor(kRed);
  			  histAxis->SetMarkerSize(MarkerSizeTom);

  		   histAxis = c2->DrawFrame(0,0.0,1.5,0.06);

  			  histAxis->SetXTitle("dimuon |y|");
  		 	  histAxis->SetYTitle(DMyTitle);
  			  histAxis->GetXaxis()->SetTitleOffset(1.2);
  			  histAxis->GetYaxis()->SetTitleOffset(1.5);


  		  histAxis->Draw("");

  			sprintf(DMLegendEntry,"Data");
  			DMLegend->AddEntry(graphMassRap,DMLegendEntry,"lp");
  			sprintf(DMLegendEntry,"Parametrization");
  			DMLegend->AddEntry(fParabola,DMLegendEntry,"l");


  			DMLegend->Draw("same");

  			graphMassRap->Draw("p,same");
  			fParabola->Draw("l,same");

  			double left=0.6, top=0.35, textSize=0.025;
  			TLatex *latex=new TLatex();
  			latex->SetTextFont(42);
  			latex->SetNDC(kTRUE);
  			latex->SetTextSize(textSize);
  			double step=textSize*1.3;

  			if(iPT==0)
  				latex->DrawLatex(left,top,Form("%.1f < p%s_{T} < %.1f GeV",onia::pTRange[0][iPT],onia::KinParticleChar,onia::pTRange[0][onia::kNbPTMaxBins]));
  			else
  				latex->DrawLatex(left,top,Form("%.1f < p%s_{T} < %.1f GeV",onia::pTRange[0][iPT-1],onia::KinParticleChar,onia::pTRange[0][iPT]));

  			top-=step;
  			top-=step;
  			latex->DrawLatex(left,top,Form("#sigma_{p0}  =  %.3f #pm %.3f GeV",a_res, a_err));
  			top-=step;
  			latex->DrawLatex(left,top,Form("#sigma_{p1}  =  %.3f #pm %.3f GeV",b_res, b_err));
  			top-=step;
  			latex->DrawLatex(left,top,Form("#sigma_{p2}  =  %.3f #pm %.3f GeV",c_res, c_err));

  			top-=step;
  			top-=step;
  			latex->DrawLatex(left,top,Form("#sigma(|y|)  =  #sigma_{p0} + #sigma_{p1}*|y| + #sigma_{p2}*|y|^{2}"));

  			char savename[200];
  			  c2->SetFrameBorderMode(0);
  		  	  sprintf(savename,"Figures/JpsiMassRap_pT%d.pdf",iPT);
  		  	  if(fitRapDependence) c2->SaveAs(savename);






	  		  	sprintf(name, "fParabolaMean");
	  		  	TF1* fParabolaMean = new TF1(name, paramMassRapParabola, 0., 1.5, 3);
  	  		  	a=3.094;
  	  		  	b=-0.002;
  	  		  	c=-0.002;
	  		  	fParabolaMean->FixParameter(0,a);
	  		  	fParabolaMean->SetParameter(1,b);
	  		  	fParabolaMean->SetParameter(2,c);
	  		  	if(fitRapDependence) graphMeanMassRap->Fit(fParabolaMean, "0", "", 0., 1.5);

	  		  	fParabolaMean->SetLineColor(kBlue);
	  		  	fParabolaMean->SetLineWidth(1.);

	  		  	a_res=fParabolaMean->GetParameter(0);
	  		  	a_err=fParabolaMean->GetParError(0);
	  		  	b_res=fParabolaMean->GetParameter(1);
	  		  	b_err=fParabolaMean->GetParError(1);
	  		  	c_res=fParabolaMean->GetParameter(2);
	  		  	c_err=fParabolaMean->GetParError(2);


  	  		  sprintf(DMyTitle,"pole mass #mu [GeV]");

  	  		  TCanvas* c3 = new TCanvas("c3","c3",1200,1100);
  	  		  gStyle->SetPalette(1);
  	  	 	  gPad->SetFillColor(kWhite);

  	  	 	  gPad->SetLeftMargin(MargLeft);
  	  	      gPad->SetRightMargin(MargRight);
  	  	      gPad->SetTopMargin(MargTop);


  	  		   TH1F *MeanhistAxis = new TH1F;

  	  			  MeanhistAxis->SetTitle(0);
  	  			  MeanhistAxis->SetStats(0);
  	  			  MeanhistAxis->SetMarkerStyle(25);
  	  			  MeanhistAxis->SetMarkerSize(MarkerSizeTom);
  	  			  MeanhistAxis->SetMarkerColor(kBlack);
  	  			  MeanhistAxis->SetLineColor(kBlack);
  	  			  MeanhistAxis->SetTitle(0);
  	  			  MeanhistAxis->SetStats(0);
  	  			  MeanhistAxis->SetMarkerColor(kRed);
  	  			  MeanhistAxis->SetMarkerStyle(20);
  	  			  MeanhistAxis->SetLineColor(kRed);
  	  			  MeanhistAxis->SetMarkerSize(MarkerSizeTom);

  	  		   MeanhistAxis = c2->DrawFrame(0,3.086,1.5,3.098);

  	  			  MeanhistAxis->SetXTitle("dimuon |y|");
  	  		 	  MeanhistAxis->SetYTitle(DMyTitle);
  	  			  MeanhistAxis->GetXaxis()->SetTitleOffset(1.2);
  	  			  MeanhistAxis->GetYaxis()->SetTitleOffset(1.685);


  	  		  MeanhistAxis->Draw("");


  	  			DMLegend->Draw("same");

  	  			graphMeanMassRap->Draw("p,same");
  	  			fParabolaMean->Draw("l,same");

  	  			left=0.6; top=0.35; textSize=0.025;
  	  			latex->SetTextFont(42);
  	  			latex->SetNDC(kTRUE);
  	  			latex->SetTextSize(textSize);
  	  			step=textSize*1.3;

  	  			if(iPT==0)
  	  				latex->DrawLatex(left,top,Form("%.1f < p%s_{T} < %.1f GeV",onia::pTRange[0][iPT],onia::KinParticleChar,onia::pTRange[0][onia::kNbPTMaxBins]));
  	  			else
  	  				latex->DrawLatex(left,top,Form("%.1f < p%s_{T} < %.1f GeV",onia::pTRange[0][iPT-1],onia::KinParticleChar,onia::pTRange[0][iPT]));

  	  			top-=step;
  	  			top-=step;
  	  			latex->DrawLatex(left,top,Form("#mu_{p0}  =  %.3f #pm %.3f GeV",a_res, a_err));
  	  			top-=step;
  	  			latex->DrawLatex(left,top,Form("#mu_{p1}  =  %.3f #pm %.3f GeV",b_res, b_err));
  	  			top-=step;
  	  			latex->DrawLatex(left,top,Form("#mu_{p2}  =  %.3f #pm %.3f GeV",c_res, c_err));

  	  			top-=step;
  	  			top-=step;
  	  			latex->DrawLatex(left,top,Form("#mu(|y|)  =  #mu_{p0} + #mu_{p1}*|y| + #mu_{p2}*|y|^{2}"));

  	  			  c3->SetFrameBorderMode(0);
  	  		  	  sprintf(savename,"Figures/JpsiMeanMassRap_pT%d.pdf",iPT);
  	  		  	  if(fitRapDependence) c3->SaveAs(savename);





  	  		  	sprintf(name, "fParabolaAlpha");
  	  		  	TF1* fParabolaAlpha = new TF1(name, paramMassRapParabola, 0., 1.5, 3);
	  		  	a=1.7;
	  		  	b=0.25;
	  		  	c=0.;
  	  		  	fParabolaAlpha->SetParameter(0,a);
  	  		  	fParabolaAlpha->SetParameter(1,b);
  	  		  	fParabolaAlpha->FixParameter(2,c);
  	  		  	if(fitRapDependence) graphAlphaMassRap->Fit(fParabolaAlpha, "0", "", 0., 1.5);

  	  		  	fParabolaAlpha->SetLineColor(kBlue);
  	  		  	fParabolaAlpha->SetLineWidth(1.);

  	  		  	a_res=fParabolaAlpha->GetParameter(0);
  	  		  	a_err=fParabolaAlpha->GetParError(0);
  	  		  	b_res=fParabolaAlpha->GetParameter(1);
  	  		  	b_err=fParabolaAlpha->GetParError(1);
  	  		  	c_res=fParabolaAlpha->GetParameter(2);
  	  		  	c_err=fParabolaAlpha->GetParError(2);

  	  	  		  sprintf(DMyTitle,"#alpha");

  	  	  		  TCanvas* c4 = new TCanvas("c4","c4",1200,1100);
  	  	  		  gStyle->SetPalette(1);
  	  	  	 	  gPad->SetFillColor(kWhite);

  	  	  	 	  gPad->SetLeftMargin(MargLeft);
  	  	  	      gPad->SetRightMargin(MargRight);
  	  	  	      gPad->SetTopMargin(MargTop);


  	  	  		   TH1F *AlphahistAxis = new TH1F;

  	  	  			  AlphahistAxis->SetTitle(0);
  	  	  			  AlphahistAxis->SetStats(0);
  	  	  			  AlphahistAxis->SetMarkerStyle(25);
  	  	  			  AlphahistAxis->SetMarkerSize(MarkerSizeTom);
  	  	  			  AlphahistAxis->SetMarkerColor(kBlack);
  	  	  			  AlphahistAxis->SetLineColor(kBlack);
  	  	  			  AlphahistAxis->SetTitle(0);
  	  	  			  AlphahistAxis->SetStats(0);
  	  	  			  AlphahistAxis->SetMarkerColor(kRed);
  	  	  			  AlphahistAxis->SetMarkerStyle(20);
  	  	  			  AlphahistAxis->SetLineColor(kRed);
  	  	  			  AlphahistAxis->SetMarkerSize(MarkerSizeTom);

  	  	  		   AlphahistAxis = c2->DrawFrame(0,1.2,1.5,2.6);

  	  	  			  AlphahistAxis->SetXTitle("dimuon |y|");
  	  	  		 	  AlphahistAxis->SetYTitle(DMyTitle);
  	  	  			  AlphahistAxis->GetXaxis()->SetTitleOffset(1.2);
  	  	  			  AlphahistAxis->GetYaxis()->SetTitleOffset(1.65);


  	  	  		  AlphahistAxis->Draw("");


  	  	  			DMLegend->Draw("same");

  	  	  			graphAlphaMassRap->Draw("p,same");
  	  	  			fParabolaAlpha->Draw("l,same");

  	  	  			left=0.6; top=0.35; textSize=0.025;
  	  	  			latex->SetTextFont(42);
  	  	  			latex->SetNDC(kTRUE);
  	  	  			latex->SetTextSize(textSize);
  	  	  			step=textSize*1.3;

  	  	  			if(iPT==0)
  	  	  				latex->DrawLatex(left,top,Form("%.1f < p%s_{T} < %.1f GeV",onia::pTRange[0][iPT],onia::KinParticleChar,onia::pTRange[0][onia::kNbPTMaxBins]));
  	  	  			else
  	  	  				latex->DrawLatex(left,top,Form("%.1f < p%s_{T} < %.1f GeV",onia::pTRange[0][iPT-1],onia::KinParticleChar,onia::pTRange[0][iPT]));

  	    			top-=step;
  	    			top-=step;
  	    			latex->DrawLatex(left,top,Form("#alpha_{p0}  =  %.3f #pm %.3f GeV",a_res, a_err));
  	    			top-=step;
  	    			latex->DrawLatex(left,top,Form("#alpha_{p1}  =  %.3f #pm %.3f GeV",b_res, b_err));

  	    			top-=step;
  	    			top-=step;
  	    			latex->DrawLatex(left,top,Form("#alpha(|y|)  =  #alpha_{p0} + #alpha_{p1}*|y|"));

  	  	  			  c4->SetFrameBorderMode(0);
  	  	  		  	  sprintf(savename,"Figures/JpsiAlphaMassRap_pT%d.pdf",iPT);
  	  	  		  	  if(fitRapDependence) c4->SaveAs(savename);

  		    iPtCount++;
  	    }


    return 0;
}

//==================================
Double_t paramMassRapParabola(Double_t *x, Double_t *par){
	//Double_t result=par[0]+x[0]*x[0]*par[1]+x[0]*x[0]*x[0]*par[2];
	Double_t result=par[0]+x[0]*par[1]+x[0]*x[0]*par[2];
	return result;
}
