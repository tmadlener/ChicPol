/*
 * CtauStudyPlots.C
 *
 *  Created on: May 23, 2014
 *      Author: valentinknuenz
 */



#include "Riostream.h"
#include "TSystem.h"
#include "TString.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TF2.h"
#include "TLatex.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraphAsymmErrors.h"
#include "TPad.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"


#include "../interface/commonVar.h"

using namespace onia;

void Plot2DHisto(TH2D* plothist, char xTitle[200], char yTitle[200], char savename[200], char tex[200]);
void Plot1DHisto(TH1D* plothist, char Title[200], char savename[200], char tex[200], double yMax);
void Plot1DHistoComp(TH1D* plothist_original, TH1D* plothist_cut, char Title[200], char savename[200], char tex1[200], char tex2[200], double yMax, bool Normalize);

void CtauStudyPlots(){



	  gROOT->Reset();
	  gROOT->SetBatch();

	  char saveDir[200];
	  char PlotID[200];
	  char savename[200];
	  sprintf(saveDir,"polFit/Figures/CtauStudy");
	  gSystem->mkdir(saveDir);
	  sprintf(PlotID,"2014May26_PlotLargerRegions");
	  sprintf(saveDir,"%s/%s",saveDir,PlotID);
	  gSystem->mkdir(saveDir);





	  	int iRap=1;
		int nBins=100;
		int nBins2D=60;
		const int nPT=onia::kNbPTMaxBins;

		double ctauBorder=0.05;
		double LxyBorder=0.25;
		double ctauerrBorder=0.05;
		double ctausigBorder=2.5;

		ctauerrBorder=0.15;//

		TH1D* h_ctau_2011[nPT+1];
		TH1D* h_ctau_2012[nPT+1];
		double h_ctau_min=-ctauBorder;
		double h_ctau_max=ctauBorder;

		TH1D* h_ctauerr_2011[nPT+1];
		TH1D* h_ctauerr_2012[nPT+1];
		double h_ctauerr_min=0.;
		double h_ctauerr_max=ctauerrBorder;

		TH1D* h_ctauerrpos_2011[nPT+1];
		TH1D* h_ctauerrpos_2012[nPT+1];
		double h_ctauerrpos_min=0.;
		double h_ctauerrpos_max=ctauerrBorder;

		TH1D* h_ctauerrneg_2011[nPT+1];
		TH1D* h_ctauerrneg_2012[nPT+1];
		double h_ctauerrneg_min=0.;
		double h_ctauerrneg_max=ctauerrBorder;

		TH1D* h_ctausig_2011[nPT+1];
		TH1D* h_ctausig_2012[nPT+1];
		double h_ctausig_min=-ctausigBorder;
		double h_ctausig_max=ctausigBorder;

		TH1D* h_Lxy_2011[nPT+1];
		TH1D* h_Lxy_2012[nPT+1];
		double h_Lxy_min=-LxyBorder;
		double h_Lxy_max=LxyBorder;

		//
		h_ctau_min=-0.3;
		h_ctau_max=1.;
		h_Lxy_min=-1.5;
		h_Lxy_max=5.;
		h_ctausig_min=-15;
		h_ctausig_max=50;
		h_ctau_min=-0.3;
		h_ctau_max=1.;



		TH2D* h_ctausig_vs_ctauerr_2011[nPT+1];
		TH2D* h_ctausig_vs_ctauerr_2012[nPT+1];
		TH2D* h_ctausig_vs_ctau_2011[nPT+1];
		TH2D* h_ctausig_vs_ctau_2012[nPT+1];
		TH2D* h_ctau_vs_ctauerr_2011[nPT+1];
		TH2D* h_ctau_vs_ctauerr_2012[nPT+1];

		for(int iPT = 0; iPT <= nPT+1; iPT++){

			h_ctau_2011[iPT] = new TH1D(Form("h_ctau_2011_pt%d",iPT), Form("h_ctau_2011_pt%d",iPT), nBins, h_ctau_min, h_ctau_max);
			h_ctau_2012[iPT] = new TH1D(Form("h_ctau_2012_pt%d",iPT), Form("h_ctau_2012_pt%d",iPT), nBins, h_ctau_min, h_ctau_max);
			h_ctauerr_2011[iPT] = new TH1D(Form("h_ctauerr_2011_pt%d",iPT), Form("h_ctauerr_2011_pt%d",iPT), nBins, h_ctauerr_min, h_ctauerr_max);
			h_ctauerr_2012[iPT] = new TH1D(Form("h_ctauerr_2012_pt%d",iPT), Form("h_ctauerr_2012_pt%d",iPT), nBins, h_ctauerr_min, h_ctauerr_max);
			h_ctauerrpos_2011[iPT] = new TH1D(Form("h_ctauerrpos_2011_pt%d",iPT), Form("h_ctauerrpos_2011_pt%d",iPT), nBins, h_ctauerrpos_min, h_ctauerrpos_max);
			h_ctauerrpos_2012[iPT] = new TH1D(Form("h_ctauerrpos_2012_pt%d",iPT), Form("h_ctauerrpos_2012_pt%d",iPT), nBins, h_ctauerrpos_min, h_ctauerrpos_max);
			h_ctauerrneg_2011[iPT] = new TH1D(Form("h_ctauerrneg_2011_pt%d",iPT), Form("h_ctauerrneg_2011_pt%d",iPT), nBins, h_ctauerrneg_min, h_ctauerrneg_max);
			h_ctauerrneg_2012[iPT] = new TH1D(Form("h_ctauerrneg_2012_pt%d",iPT), Form("h_ctauerrneg_2012_pt%d",iPT), nBins, h_ctauerrneg_min, h_ctauerrneg_max);
			h_ctausig_2011[iPT] = new TH1D(Form("h_ctausig_2011_pt%d",iPT), Form("h_ctausig_2011_pt%d",iPT), nBins, h_ctausig_min, h_ctausig_max);
			h_ctausig_2012[iPT] = new TH1D(Form("h_ctausig_2012_pt%d",iPT), Form("h_ctausig_2012_pt%d",iPT), nBins, h_ctausig_min, h_ctausig_max);

			h_Lxy_2011[iPT] = new TH1D(Form("h_Lxy_2011_pt%d",iPT), Form("h_Lxy_2011_pt%d",iPT), nBins, h_Lxy_min, h_Lxy_max);
			h_Lxy_2012[iPT] = new TH1D(Form("h_Lxy_2012_pt%d",iPT), Form("h_Lxy_2012_pt%d",iPT), nBins, h_Lxy_min, h_Lxy_max);

			h_ctausig_vs_ctauerr_2011[iPT] = new TH2D(Form("h_ctausig_vs_ctauerr_2011_pt%d",iPT), Form("h_ctausig_vs_ctauerr_2011_pt%d",iPT), nBins2D, h_ctausig_min, h_ctausig_max, nBins2D, h_ctauerr_min, h_ctauerr_max);
			h_ctausig_vs_ctauerr_2012[iPT] = new TH2D(Form("h_ctausig_vs_ctauerr_2012_pt%d",iPT), Form("h_ctausig_vs_ctauerr_2012_pt%d",iPT), nBins2D, h_ctausig_min, h_ctausig_max, nBins2D, h_ctauerr_min, h_ctauerr_max);
			h_ctausig_vs_ctau_2011[iPT] = new TH2D(Form("h_ctausig_vs_ctau_2011_pt%d",iPT), Form("h_ctausig_vs_ctau_2011_pt%d",iPT), nBins2D, h_ctausig_min, h_ctausig_max, nBins2D, h_ctau_min, h_ctau_max);
			h_ctausig_vs_ctau_2012[iPT] = new TH2D(Form("h_ctausig_vs_ctau_2012_pt%d",iPT), Form("h_ctausig_vs_ctau_2012_pt%d",iPT), nBins2D, h_ctausig_min, h_ctausig_max, nBins2D, h_ctau_min, h_ctau_max);
			h_ctau_vs_ctauerr_2011[iPT] = new TH2D(Form("h_ctau_vs_ctauerr_2011_pt%d",iPT), Form("h_ctau_vs_ctauerr_2011_pt%d",iPT), nBins2D, h_ctau_min, h_ctau_max, nBins2D, h_ctauerr_min, h_ctauerr_max);
			h_ctau_vs_ctauerr_2012[iPT] = new TH2D(Form("h_ctau_vs_ctauerr_2012_pt%d",iPT), Form("h_ctau_vs_ctauerr_2012_pt%d",iPT), nBins2D, h_ctau_min, h_ctau_max, nBins2D, h_ctauerr_min, h_ctauerr_max);

		}


		  TFile* inFileData2011 = new TFile("/scratch/knuenz/Polarization/RootInput/Psi/TTree_Onia2MuMu_v30_PromptRecoAB_10May2012_Jpsi.root");
		  TFile* inFileData2012 = new TFile("/scratch/knuenz/Polarization/RootInput/ChicPol/chic_rootuple_subFeb2014.root");

		  TTree* Tree2011;
		  TTree* Tree2012;

		  Tree2011=(TTree*)inFileData2011->Get("data");
		  Tree2012=(TTree*)inFileData2012->Get("rootuple/chicTree");


		  TLorentzVector  *Jpsi=0;
		  TLorentzVector  *lepP=0;
		  TLorentzVector  *lepN=0;
		  double ctau, ctauerr;

		  Tree2012->SetBranchAddress("dimuon_p4", &Jpsi);
		  Tree2012->SetBranchAddress("muonP_p4", &lepP);
		  Tree2012->SetBranchAddress("muonN_p4", &lepN);
		  Tree2012->SetBranchAddress("ctpv", &ctau);
		  Tree2012->SetBranchAddress("ctpv_error", &ctauerr);

		  Tree2011->SetBranchAddress("JpsiP", &Jpsi);
		  Tree2011->SetBranchAddress("muPosP", &lepP);
		  Tree2011->SetBranchAddress("muNegP", &lepN);
		  Tree2011->SetBranchAddress("Jpsict", &ctau);
		  Tree2011->SetBranchAddress("JpsictErr", &ctauerr);



			Long64_t nentries2011 = Tree2011->GetEntries();
			Long64_t nb2011 = 0;


			std::cout << "number of entries 2011 = " << nentries2011 << std::endl;

			Long64_t nentries2012 = Tree2012->GetEntries();
			Long64_t nb2012 = 0;

			std::cout << "number of entries 2012 = " << nentries2012 << std::endl;

			int nentries=200000;

			//nentries2011=nentries;
			//nentries2012=nentries;


			//loop over the events
			for (Long64_t jentry=0; jentry<nentries2011; jentry++) {

				if(jentry % 100000 == 0) std::cout << "event " << jentry << " of " << nentries2011 << std::endl;

				nb2011 = Tree2011->GetEntry(jentry);


				//if( !( HLT_Dimuon10_Jpsi_Barrel_v1 == 1 ||
				//       HLT_Dimuon10_Jpsi_Barrel_v2 == 1 ||
				//       HLT_Dimuon10_Jpsi_Barrel_v3 == 1 ||
				//       HLT_Dimuon10_Jpsi_Barrel_v5 == 1 ||
				//       HLT_Dimuon10_Jpsi_Barrel_v6 == 1 ||
				//       HLT_Dimuon10_Jpsi_Barrel_v9 == 1 )
				//       ) continue;
				//)

				for(int iPT = 0; iPT <= nPT+1; iPT++){


					Double_t ptMin;
					Double_t ptMax;
					if(iPT==0){
						ptMin = onia::pTRange[iRap][0];
						ptMax = onia::pTRange[iRap][onia::kNbPTBins[0]];
					} else{
						ptMin = onia::pTRange[iRap][iPT-1];
						ptMax = onia::pTRange[iRap][iPT];
					}

				    if(Jpsi->Pt()>ptMin && Jpsi->Pt()<ptMax){
						h_ctau_2011[iPT]->Fill(ctau);
						h_ctauerr_2011[iPT]->Fill(ctauerr);
						if(ctau>0&&ctau<ctauBorder) h_ctauerrpos_2011[iPT]->Fill(ctauerr);
						if(ctau<0&&ctau>-ctauBorder) h_ctauerrneg_2011[iPT]->Fill(ctauerr);
						h_ctausig_2011[iPT]->Fill(ctau/ctauerr);
						h_Lxy_2011[iPT]->Fill(ctau*Jpsi->Pt()/Jpsi->M());
						h_ctausig_vs_ctauerr_2011[iPT]->Fill(ctau/ctauerr,ctauerr);
						h_ctausig_vs_ctau_2011[iPT]->Fill(ctau/ctauerr,ctau);
						h_ctau_vs_ctauerr_2011[iPT]->Fill(ctau,ctauerr);
				    }

				}


			}
			cout<<"Finished Loop over "<<nentries2011<<" entries"<<endl;





			//loop over the events
			for (Long64_t jentry=0; jentry<nentries2012; jentry++) {

				if(jentry % 100000 == 0) std::cout << "event " << jentry << " of " << nentries2012 << std::endl;

				nb2012 = Tree2012->GetEntry(jentry);

				ctau*=10.;
				ctauerr*=10.;

				for(int iPT = 0; iPT <= nPT+1; iPT++){


					Double_t ptMin;
					Double_t ptMax;
					if(iPT==0){
						ptMin = onia::pTRange[iRap][0];
						ptMax = onia::pTRange[iRap][onia::kNbPTBins[0]];
					} else{
						ptMin = onia::pTRange[iRap][iPT-1];
						ptMax = onia::pTRange[iRap][iPT];
					}


				    if(Jpsi->Pt()>ptMin && Jpsi->Pt()<ptMax){
						h_ctau_2012[iPT]->Fill(ctau);
						h_ctauerr_2012[iPT]->Fill(ctauerr);
						if(ctau>0&&ctau<ctauBorder) h_ctauerrpos_2012[iPT]->Fill(ctauerr);
						if(ctau<0&&ctau>-ctauBorder) h_ctauerrneg_2012[iPT]->Fill(ctauerr);
						h_ctausig_2012[iPT]->Fill(ctau/ctauerr);
						h_Lxy_2012[iPT]->Fill(ctau*Jpsi->Pt()/Jpsi->M());
						h_ctausig_vs_ctauerr_2012[iPT]->Fill(ctau/ctauerr,ctauerr);
						h_ctausig_vs_ctau_2012[iPT]->Fill(ctau/ctauerr,ctau);
						h_ctau_vs_ctauerr_2012[iPT]->Fill(ctau,ctauerr);
				    }

				}


			}
			cout<<"Finished Loop over "<<nentries2012<<" entries"<<endl;







		char xTitle[200];
		char yTitle[200];

		char tex[200];
		char tex1[200];
		char tex2[200];
		double yMax=0.;//=CosThDist_Original->GetMaximum()*expandY;

		bool Normalize=true;


		for(int iPT = 0; iPT <= nPT+1; iPT++){

			sprintf(tex1,"2011 data");
			sprintf(tex2,"2012 data");

			sprintf(xTitle,"l");
			sprintf(savename,"%s/Comp_ctau_pT%d",saveDir, iPT);
			Normalize=true;
			Plot1DHistoComp(h_ctau_2011[iPT], h_ctau_2012[iPT], xTitle, savename, tex1, tex2, yMax, Normalize);

			sprintf(xTitle,"#sigma_{l}");
			sprintf(savename,"%s/Comp_ctauerr_pT%d",saveDir, iPT);
			Normalize=true;
			Plot1DHistoComp(h_ctauerr_2011[iPT], h_ctauerr_2012[iPT], xTitle, savename, tex1, tex2, yMax, Normalize);

			sprintf(xTitle,"l / #sigma_{l}");
			sprintf(savename,"%s/Comp_ctausig_pT%d",saveDir, iPT);
			Normalize=true;
			Plot1DHistoComp(h_ctausig_2011[iPT], h_ctausig_2012[iPT], xTitle, savename, tex1, tex2, yMax, Normalize);


			sprintf(xTitle,"L_{xy}");
			sprintf(savename,"%s/Comp_Lxy_pT%d",saveDir, iPT);
			Normalize=true;
			Plot1DHistoComp(h_Lxy_2011[iPT], h_Lxy_2012[iPT], xTitle, savename, tex1, tex2, yMax, Normalize);




			sprintf(tex1,Form("%1.2f < l < 0",-ctauBorder));
			sprintf(tex2,Form("0 < l < %1.2f",ctauBorder));

			sprintf(xTitle,"#sigma_{l}");
			sprintf(savename,"%s/Comp_ctauerr_posneg2011_pT%d",saveDir, iPT);
			Normalize=true;
			Plot1DHistoComp(h_ctauerrneg_2011[iPT], h_ctauerrpos_2011[iPT], xTitle, savename, tex1, tex2, yMax, Normalize);

			sprintf(xTitle,"#sigma_{l}");
			sprintf(savename,"%s/Comp_ctauerr_posneg2012_pT%d",saveDir, iPT);
			Normalize=true;
			Plot1DHistoComp(h_ctauerrneg_2012[iPT], h_ctauerrpos_2012[iPT], xTitle, savename, tex1, tex2, yMax, Normalize);





			char ptchar[200];
			char rapchar[200];
			char yearchar[200];

			if(iPT==0) sprintf(ptchar,"%.0f < p%s_{T} < %.0f GeV",onia::pTRange[1][iPT],onia::KinParticleChar,onia::pTRange[1][onia::kNbPTBins[1]]);
			else sprintf(ptchar,"%.0f < p%s_{T} < %.0f GeV",onia::pTRange[1][iPT-1],onia::KinParticleChar,onia::pTRange[1][iPT]);
			sprintf(rapchar,"|y%s| < %.1f",onia::KinParticleChar,onia::rapForPTRange[1]);


		sprintf(yearchar,"2011 data");

			sprintf(tex,"%s, %s, %s", yearchar, rapchar, ptchar);
			sprintf(xTitle,"l / #sigma_{l}");
			sprintf(yTitle,"#sigma_{l}");
			sprintf(savename,"%s/2D_ctausig_vs_ctauerr_2011_pT%d",saveDir, iPT);
			Plot2DHisto(h_ctausig_vs_ctauerr_2011[iPT], xTitle, yTitle, savename, tex);

			sprintf(tex,"%s, %s, %s", yearchar, rapchar, ptchar);
			sprintf(xTitle,"l / #sigma_{l}");
			sprintf(yTitle,"l");
			sprintf(savename,"%s/2D_ctausig_vs_ctau_2011_pT%d",saveDir, iPT);
			Plot2DHisto(h_ctausig_vs_ctau_2011[iPT], xTitle, yTitle, savename, tex);

			sprintf(tex,"%s, %s, %s", yearchar, rapchar, ptchar);
			sprintf(xTitle,"l");
			sprintf(yTitle,"#sigma_{l}");
			sprintf(savename,"%s/2D_ctau_vs_ctauerr_2011_pT%d",saveDir, iPT);
			Plot2DHisto(h_ctau_vs_ctauerr_2011[iPT], xTitle, yTitle, savename, tex);


		sprintf(yearchar,"2012 data");

			sprintf(tex,"%s, %s, %s", yearchar, rapchar, ptchar);
			sprintf(xTitle,"l / #sigma_{l}");
			sprintf(yTitle,"#sigma_{l}");
			sprintf(savename,"%s/2D_ctausig_vs_ctauerr_2012_pT%d",saveDir, iPT);
			Plot2DHisto(h_ctausig_vs_ctauerr_2012[iPT], xTitle, yTitle, savename, tex);

			sprintf(tex,"%s, %s, %s", yearchar, rapchar, ptchar);
			sprintf(xTitle,"l / #sigma_{l}");
			sprintf(yTitle,"l");
			sprintf(savename,"%s/2D_ctausig_vs_ctau_2012_pT%d",saveDir, iPT);
			Plot2DHisto(h_ctausig_vs_ctau_2012[iPT], xTitle, yTitle, savename, tex);

			sprintf(tex,"%s, %s, %s", yearchar, rapchar, ptchar);
			sprintf(xTitle,"l");
			sprintf(yTitle,"#sigma_{l}");
			sprintf(savename,"%s/2D_ctau_vs_ctauerr_2012_pT%d",saveDir, iPT);
			Plot2DHisto(h_ctau_vs_ctauerr_2012[iPT], xTitle, yTitle, savename, tex);


		}



		for(int iPT = 0; iPT <= nPT+1; iPT++){

			delete h_ctau_2011[iPT];
			delete h_ctau_2012[iPT];
			delete h_ctauerr_2011[iPT];
			delete h_ctauerr_2012[iPT];
			delete h_ctauerrpos_2011[iPT];
			delete h_ctauerrpos_2012[iPT];
			delete h_ctauerrneg_2011[iPT];
			delete h_ctauerrneg_2012[iPT];
			delete h_ctausig_2011[iPT];
			delete h_ctausig_2012[iPT];
			delete h_Lxy_2011[iPT];
			delete h_Lxy_2012[iPT];

		}

		inFileData2011->Close();
		inFileData2012->Close();

		delete inFileData2011;
		delete inFileData2012;

		return;


}


void Plot2DHisto(TH2D* plothist, char xTitle[200], char yTitle[200], char savename[200], char tex[200]){

	TGaxis::SetMaxDigits(3);

	TCanvas *c1 = new TCanvas("", "", 1200, 1000);
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gPad->SetFillColor(kWhite);
	gPad->SetLeftMargin(0.15);
	gPad->SetRightMargin(0.2);
	gPad->SetTopMargin(0.1);

	double yOffset=1.4;
	plothist->GetYaxis()->SetTitleOffset(yOffset);
	plothist->SetStats(0);
	plothist->SetTitle(0);
	plothist->GetYaxis()->SetTitle(yTitle);
	plothist->GetXaxis()->SetTitle(xTitle);
	plothist->SetMinimum(0.);

	plothist->Draw("colz");

	double left=0.25, top=0.935, textSize=0.04;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	latex->DrawLatex(left,top,tex);



	char linlog_savename[200];
	c1->SetLogz(false);
	sprintf(linlog_savename,"%s_lin.pdf",savename);
	c1->SaveAs(linlog_savename);
	c1->SetLogz(true);
	sprintf(linlog_savename,"%s_log.pdf",savename);
	c1->SaveAs(linlog_savename);

	TProfile* profileX = plothist->ProfileX("", 1, -1, "so");
	TProfile* profileY = plothist->ProfileY("", 1, -1, "so");

	profileX->Draw("same");

	sprintf(linlog_savename,"%s_log_AddProfile.pdf",savename);
	c1->SaveAs(linlog_savename);

	profileX->Draw();
	latex->DrawLatex(left,top,tex);
	sprintf(linlog_savename,"%s_log_ProfileX.pdf",savename);
	c1->SaveAs(linlog_savename);
	profileY->Draw();
	latex->DrawLatex(left,top,tex);
	sprintf(linlog_savename,"%s_log_ProfileY.pdf",savename);
	c1->SaveAs(linlog_savename);

	delete c1;
	delete latex;

}

void Plot1DHisto(TH1D* plothist, char Title[200], char savename[200], char tex[200], double yMax){

	TGaxis::SetMaxDigits(3);

	TCanvas *c1 = new TCanvas("", "", 1200, 1000);
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gPad->SetFillColor(kWhite);
	gPad->SetLeftMargin(0.15);
	gPad->SetRightMargin(0.1);
	gPad->SetTopMargin(0.1);

	double yOffset=1.4;
	plothist->GetYaxis()->SetTitleOffset(yOffset);
	plothist->SetStats(0);
	plothist->SetTitle(0);
	plothist->GetXaxis()->SetTitle(Title);
	plothist->GetYaxis()->SetTitle("Events");

	plothist->SetMinimum(0.);
	plothist->SetMaximum(yMax);

	plothist->SetMarkerStyle(20);
	plothist->SetMarkerColor(kGreen+2);

	plothist->Draw("e");

	double left=0.375, top=0.935, textSize=0.04;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	latex->DrawLatex(left,top,tex);

	c1->SaveAs(savename);

	delete c1;
	delete latex;

}

void Plot1DHistoComp(TH1D* plothist_original, TH1D* plothist_cut, char Title[200], char savename[200], char tex1[200], char tex2[200], double yMax, bool Normalize){

	TGaxis::SetMaxDigits(3);

	TCanvas *c1 = new TCanvas("", "", 1200, 1000);
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gPad->SetFillColor(kWhite);
	gPad->SetLeftMargin(0.15);
	gPad->SetRightMargin(0.1);
	gPad->SetTopMargin(0.1);

	double yOffset=1.6;
	plothist_original->GetYaxis()->SetTitleOffset(yOffset);
	plothist_original->SetStats(0);
	plothist_original->SetTitle(0);
	plothist_original->GetXaxis()->SetTitle(Title);
	plothist_original->GetYaxis()->SetTitle("Events");
	if(Normalize) plothist_original->GetYaxis()->SetTitle("arbitr. units");


	plothist_original->SetMinimum(0.);
	plothist_original->SetMaximum(yMax);

	if(Normalize){
		plothist_original->Scale(1./plothist_original->Integral());
		plothist_cut->Scale(1./plothist_cut->Integral());
		double buffmax_original=plothist_original->GetMaximum();
		double buffmax_cut=plothist_cut->GetMaximum();
		double buffmax=buffmax_original;
		if(buffmax_original<buffmax_cut) buffmax=buffmax_cut;
		plothist_original->SetMaximum(buffmax*1.3);
	}

	plothist_original->SetMarkerStyle(20);
	plothist_original->SetMarkerColor(kGreen+2);
	plothist_original->SetLineColor(kGreen+2);
	plothist_cut->SetMarkerStyle(20);
	plothist_cut->SetMarkerColor(kRed+2);
	plothist_cut->SetLineColor(kRed+2);

	plothist_original->Draw("");
	plothist_cut->Draw("same");

	double left=0.375, top=0.935, textSize=0.04;
	//double step=textSize*1.5;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	latex->SetTextColor(kGreen+2);
	latex->DrawLatex(left,top,tex1);
	left+=0.2;
	latex->SetTextColor(kRed+2);
	latex->DrawLatex(left,top,tex2);

	char linlog_savename[200];
	c1->SetLogy(false);
	sprintf(linlog_savename,"%s_lin.pdf",savename);
	c1->SaveAs(linlog_savename);
	c1->SetLogy(true);
	sprintf(linlog_savename,"%s_log.pdf",savename);
	c1->SaveAs(linlog_savename);

	delete c1;
	delete latex;

}

