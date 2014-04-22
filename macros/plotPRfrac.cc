#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TFrame.h"
#include "TCanvas.h"
#include "TGaxis.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TRandom.h"
#include "TVectorD.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TGraph.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"

void plotPRfrac(){

	gStyle->SetPadBottomMargin(0.11);
	gStyle->SetPadLeftMargin(0.08); //0.12
	gStyle->SetPadRightMargin(0.02); //0.05
	gStyle->SetPadTopMargin(0.05); //0.05

	const int RapBins = 2;
	double legendsize=0.04;

	TCanvas *c1=new TCanvas("c1","");
	c1->SetTickx();
	c1->SetTicky();
	gStyle->SetTitleFillColor(10);
	gStyle->SetTitleBorderSize(1);
	gStyle->SetOptFit(0);            
	gStyle->SetOptStat(1);
	gStyle->SetTitleFont(22);        
	gStyle->SetStatFont(22); 
	gStyle->SetStatColor(10);  
	gStyle->SetStatBorderSize(1);
	gStyle->SetLabelFont(22,"X");
	gStyle->SetLabelFont(22,"Y");
	gStyle->SetTitleXOffset(1.2);
	gStyle->SetTitleYOffset(1.2);    
	gStyle->SetHistLineWidth(2);
	gStyle->SetStatX(0.9);
	gStyle->SetStatY(0.9);
	gStyle->SetTitleX(0.15);
	gStyle->SetTitleY(0.96);

	char dataID[500];
	sprintf(dataID,"SetOfCuts11_ctauScen5_FracLSB-1_newMLfit_30Apr2013_correctfLSB");

	TFile *infile[2];
	TGraphAsymmErrors *graph_FracNP[2][RapBins]; 
	TGraphAsymmErrors *graph_FracPR[2][RapBins]; 
	TGraphAsymmErrors *graph_FracBG[2][RapBins];

	for(int ifile=0; ifile<2; ifile++){
		infile[ifile] = new TFile( Form("DataFiles/%s/Psi%dS/Fit/parameter/evaluateCtau/PR_3sigma/graphFrac_3.0sigma.root",dataID,ifile+1),"R");
		if(!infile[ifile]){cout<<"no input file"<<endl; return;}

		for(int rapBin = 1; rapBin < RapBins+1; rapBin++){
			graph_FracPR[ifile][rapBin-1] = (TGraphAsymmErrors*)infile[ifile]->Get(Form("graph_FracPR_%d",rapBin-1));
			graph_FracNP[ifile][rapBin-1] = (TGraphAsymmErrors*)infile[ifile]->Get(Form("graph_FracNP_%d",rapBin-1));
			graph_FracBG[ifile][rapBin-1] = (TGraphAsymmErrors*)infile[ifile]->Get(Form("graph_FracBG_%d",rapBin-1));

			graph_FracNP[ifile][rapBin-1] -> SetMarkerColor(kRed+1);
			graph_FracNP[ifile][rapBin-1] -> SetLineColor(kRed+1);

			graph_FracPR[ifile][rapBin-1] -> GetXaxis() -> SetTitle("#it{p}_{T} [GeV]");
			graph_FracNP[ifile][rapBin-1] -> GetXaxis() -> SetTitle("#it{p}_{T} [GeV]");
			graph_FracBG[ifile][rapBin-1] -> GetXaxis() -> SetTitle("#it{p}_{T} [GeV]");

			graph_FracPR[ifile][rapBin-1] -> GetYaxis() -> SetRangeUser(0.,1.2);
			graph_FracNP[ifile][rapBin-1] -> GetYaxis() -> SetRangeUser(0.,1.2);
			graph_FracBG[ifile][rapBin-1] -> GetYaxis() -> SetRangeUser(0.,1.2);
		}
	}

	graph_FracPR[0][0] -> SetMarkerStyle(20);
	graph_FracPR[1][0] -> SetMarkerStyle(24);

	graph_FracNP[0][0] -> SetMarkerStyle(21);
	graph_FracNP[1][0] -> SetMarkerStyle(25);

	graph_FracBG[0][0] -> SetMarkerStyle(29);
	graph_FracBG[1][0] -> SetMarkerStyle(30);
	graph_FracBG[0][0] -> SetMarkerSize(1.5);
	graph_FracBG[1][0] -> SetMarkerSize(1.5);

	double leftVal=0.57, topVal=0.88;
	double left=leftVal, top=topVal, textSize=0.045;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	double stepLatex=textSize*1.3;

	double blX = 0.75, blY = 0.45, trX = 0.93, trY = 0.6;
	TLegend* legend_rap1=new TLegend(blX,blY,trX,trY);
	legend_rap1->SetFillColor(kWhite);
	legend_rap1->SetTextFont(42);
	legend_rap1->SetTextSize(legendsize);
	legend_rap1->SetBorderSize(0.);
	legend_rap1->AddEntry(graph_FracPR[0][0],"Prompt","p");
	legend_rap1->AddEntry(graph_FracNP[0][0],"Nonprompt","p");
	legend_rap1->AddEntry(graph_FracBG[0][0],"Background","p");

	
	graph_FracPR[0][0]->GetXaxis()->SetTitleSize(0.043);
	graph_FracPR[0][0]->GetYaxis()->SetTitleSize(0.043);
	graph_FracPR[0][0]->GetXaxis()->SetTitleOffset(1.1);
	graph_FracPR[0][0]->GetYaxis()->SetTitleOffset(0.8);

	graph_FracPR[0][0]->Draw("AP");
	graph_FracNP[0][0]->Draw("P");
	graph_FracBG[0][0]->Draw("P");
	graph_FracPR[1][0]->Draw("P");
	graph_FracNP[1][0]->Draw("P");
	graph_FracBG[1][0]->Draw("P");

	legend_rap1->Draw();

	textSize=0.040;

	latex->DrawLatex(0.42,0.72,"J/#psi");
	latex->DrawLatex(0.42,0.54,"#psi(2S)");

	latex->DrawLatex(0.8,0.4,"|y| < 0.6");


	latex->SetTextFont(62);
	latex->SetTextSize(textSize);
	left=leftVal+0.15; top=topVal;
	stepLatex=textSize*1.3;
	//latex->DrawLatex(left,top, "CMS  Preliminary");
	latex->DrawLatex(left,top, "CMS");
	top -= stepLatex;
	latex->DrawLatex(left,top, "pp   #sqrt{s} = 7 TeV");
	top -= stepLatex;
	latex->DrawLatex(left,top, "L  =  4.9 fb^{-1}");


	c1->SaveAs("fraction_PsiNS_rap1.pdf");

	return ;
}
