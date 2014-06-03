/*
 * CtauErrModel.C
 *
 *  Created on: May 26, 2014
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

void CtauErrModel(){



	  gROOT->Reset();
	  gROOT->SetBatch();

	  char saveDir[200];
	  char PlotID[200];
	  char savename[200];
	  sprintf(saveDir,"polFit/Figures/CtauErrModel");
	  gSystem->mkdir(saveDir);
	  sprintf(PlotID,"2014May26_MoreLbins");
	  sprintf(saveDir,"%s/%s",saveDir,PlotID);
	  gSystem->mkdir(saveDir);





		int nBins=150;
		const int nPT=5;
		const int nRAP=2;
		const int nL=15;

		const double bordersPT[nPT+1] = {0., 12., 16., 20., 30., 100.};
		const double bordersRAP[nRAP+1] = {0., 0.6, 2.};
		const double bordersL[nL+1] = {onia::ctVarMin, -0.05, -0.03, -0.02, -0.015, -0.01, -0.005,  0., 0.005, 0.01, 0.015, 0.02, 0.03, 0.05, 0.1, onia::ctVarMax};


		double ctauerrBorder=0.15;
		TH1D* h_ctauerr_2011[nRAP+1][nPT+1][nL+1];
		TH1D* h_ctauerr_2012[nRAP+1][nPT+1][nL+1];
		double h_ctauerr_min=0.;
		double h_ctauerr_max=ctauerrBorder;


		for(int iRAP = 0; iRAP < nRAP+1; iRAP++){
			for(int iPT = 0; iPT < nPT+1; iPT++){
				for(int iL = 0; iL < nL+1; iL++){

					h_ctauerr_2011[iRAP][iPT][iL] = new TH1D(Form("h_ctauerr_2011_rap%d_pt%d_l%d",iRAP, iPT, iL), Form("h_ctauerr_2011_rap%d_pt%d_l%d",iRAP, iPT, iL), nBins, h_ctauerr_min, h_ctauerr_max);
					h_ctauerr_2012[iRAP][iPT][iL] = new TH1D(Form("h_ctauerr_2012_rap%d_pt%d_l%d",iRAP, iPT, iL), Form("h_ctauerr_2012_rap%d_pt%d_l%d",iRAP, iPT, iL), nBins, h_ctauerr_min, h_ctauerr_max);

				}
			}
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

		  Int_t HLT_Dimuon10_Jpsi_Barrel_v1;
		  Int_t HLT_Dimuon10_Jpsi_Barrel_v2;
		  Int_t HLT_Dimuon10_Jpsi_Barrel_v3;
		  Int_t HLT_Dimuon10_Jpsi_Barrel_v5;
		  Int_t HLT_Dimuon10_Jpsi_Barrel_v6;
		  Double_t JpsiVprob;

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
		  Tree2011->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v1", &HLT_Dimuon10_Jpsi_Barrel_v1);
		  Tree2011->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v2", &HLT_Dimuon10_Jpsi_Barrel_v2);
		  Tree2011->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v3", &HLT_Dimuon10_Jpsi_Barrel_v3);
		  Tree2011->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v5", &HLT_Dimuon10_Jpsi_Barrel_v5);
		  Tree2011->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v6", &HLT_Dimuon10_Jpsi_Barrel_v6);
		  Tree2011->SetBranchAddress("JpsiVprob", &JpsiVprob);





			Long64_t nentries2011 = Tree2011->GetEntries();
			Long64_t nb2011 = 0;


			std::cout << "number of entries 2011 = " << nentries2011 << std::endl;

			Long64_t nentries2012 = Tree2012->GetEntries();
			Long64_t nb2012 = 0;

			std::cout << "number of entries 2012 = " << nentries2012 << std::endl;

			int nentries=2000000;

			//nentries2011=nentries;
			//nentries2012=nentries;


			int nEventsPassingCuts=0;
			int nEventsFilled=0;
			//loop over the events
			for (Long64_t jentry=0; jentry<nentries2011; jentry++) {

				if(jentry % 100000 == 0) std::cout << "event " << jentry << " of " << nentries2011 << std::endl;

				nb2011 = Tree2011->GetEntry(jentry);

				if(JpsiVprob<0.01) continue;

				bool fiducial=true;

				if(TMath::Abs(lepP->Eta())<1.2 && lepP->Pt()<4.5) fiducial=false;
				if(TMath::Abs(lepP->Eta())>1.2 && TMath::Abs(lepP->Eta())<1.4 && lepP->Pt()<3.5) fiducial=false;
				if(TMath::Abs(lepP->Eta())>1.4 && TMath::Abs(lepP->Eta())<1.6 && lepP->Pt()<3.) fiducial=false;
				if(TMath::Abs(lepN->Eta())<1.2 && lepN->Pt()<4.5) fiducial=false;
				if(TMath::Abs(lepN->Eta())>1.2 && TMath::Abs(lepN->Eta())<1.4 && lepN->Pt()<3.5) fiducial=false;
				if(TMath::Abs(lepN->Eta())>1.4 && TMath::Abs(lepN->Eta())<1.6 && lepN->Pt()<3.) fiducial=false;

				if(!fiducial) continue;

				if( !( HLT_Dimuon10_Jpsi_Barrel_v1 == 1 ||
				       HLT_Dimuon10_Jpsi_Barrel_v2 == 1 ||
				       HLT_Dimuon10_Jpsi_Barrel_v3 == 1 ||
				       HLT_Dimuon10_Jpsi_Barrel_v5 == 1 ||
				       HLT_Dimuon10_Jpsi_Barrel_v6 == 1)
				       ) continue;

				bool inRegionOfInterest=false;
				if(Jpsi->Pt()>bordersPT[0] && Jpsi->Pt()<bordersPT[nPT] && TMath::Abs(Jpsi->Rapidity())>bordersRAP[0] && TMath::Abs(Jpsi->Rapidity())<bordersRAP[nRAP] && ctau>bordersL[0] && ctau<bordersL[nL])
					inRegionOfInterest=true;
				if(!inRegionOfInterest) continue;


				nEventsPassingCuts++;

				for(int iRAP = 0; iRAP < nRAP+1; iRAP++){
					for(int iPT = 0; iPT < nPT+1; iPT++){
						for(int iL = 0; iL < nL+1; iL++){



							Double_t ptMin;
							Double_t ptMax;
							if(iPT==0){
								ptMin = bordersPT[0];
								ptMax = bordersPT[nPT];
							} else{
								ptMin = bordersPT[iPT-1];
								ptMax = bordersPT[iPT];
							}
							Double_t rapMin;
							Double_t rapMax;
							if(iRAP==0){
								rapMin = bordersRAP[0];
								rapMax = bordersRAP[nRAP];
							} else{
								rapMin = bordersRAP[iRAP-1];
								rapMax = bordersRAP[iRAP];
							}
							Double_t lMin;
							Double_t lMax;
							if(iL==0){
								lMin = bordersL[0];
								lMax = bordersL[nL];
							} else{
								lMin = bordersL[iL-1];
								lMax = bordersL[iL];
							}


							if(Jpsi->Pt()>ptMin && Jpsi->Pt()<ptMax && TMath::Abs(Jpsi->Rapidity())>rapMin && TMath::Abs(Jpsi->Rapidity())<rapMax && ctau>lMin && ctau<lMax){
								h_ctauerr_2011[iRAP][iPT][iL]->Fill(ctauerr);
								nEventsFilled++;

						}

						}

					}
				}

			}
			cout<<"2011: Finished Loop over "<<nentries2011<<" entries"<<endl;

			cout<<"2011: "<<nEventsPassingCuts<<" entries passing cuts"<<endl;
			cout<<"2011: "<<nEventsFilled<<" entries filled"<<endl;




			nEventsPassingCuts=0;
			nEventsFilled=0;
			//loop over the events
			for (Long64_t jentry=0; jentry<nentries2012; jentry++) {

				if(jentry % 100000 == 0) std::cout << "event " << jentry << " of " << nentries2012 << std::endl;

				nb2012 = Tree2012->GetEntry(jentry);

				ctau*=10.;
				ctauerr*=10.;

				bool fiducial=true;

				if(TMath::Abs(lepP->Eta())<1.2 && lepP->Pt()<4.5) fiducial=false;
				if(TMath::Abs(lepP->Eta())>1.2 && TMath::Abs(lepP->Eta())<1.4 && lepP->Pt()<3.5) fiducial=false;
				if(TMath::Abs(lepP->Eta())>1.4 && TMath::Abs(lepP->Eta())<1.6 && lepP->Pt()<3.) fiducial=false;
				if(TMath::Abs(lepN->Eta())<1.2 && lepN->Pt()<4.5) fiducial=false;
				if(TMath::Abs(lepN->Eta())>1.2 && TMath::Abs(lepN->Eta())<1.4 && lepN->Pt()<3.5) fiducial=false;
				if(TMath::Abs(lepN->Eta())>1.4 && TMath::Abs(lepN->Eta())<1.6 && lepN->Pt()<3.) fiducial=false;

				if(!fiducial) continue;

				bool inRegionOfInterest=false;
				if(Jpsi->Pt()>bordersPT[0] && Jpsi->Pt()<bordersPT[nPT] && TMath::Abs(Jpsi->Rapidity())>bordersRAP[0] && TMath::Abs(Jpsi->Rapidity())<bordersRAP[nRAP] && ctau>bordersL[0] && ctau<bordersL[nL])
					inRegionOfInterest=true;
				if(!inRegionOfInterest) continue;


				nEventsPassingCuts++;

				for(int iRAP = 0; iRAP < nRAP+1; iRAP++){
					for(int iPT = 0; iPT < nPT+1; iPT++){
						for(int iL = 0; iL < nL+1; iL++){



							Double_t ptMin;
							Double_t ptMax;
							if(iPT==0){
								ptMin = bordersPT[0];
								ptMax = bordersPT[nPT];
							} else{
								ptMin = bordersPT[iPT-1];
								ptMax = bordersPT[iPT];
							}
							Double_t rapMin;
							Double_t rapMax;
							if(iRAP==0){
								rapMin = bordersRAP[0];
								rapMax = bordersRAP[nRAP];
							} else{
								rapMin = bordersRAP[iRAP-1];
								rapMax = bordersRAP[iRAP];
							}
							Double_t lMin;
							Double_t lMax;
							if(iL==0){
								lMin = bordersL[0];
								lMax = bordersL[nL];
							} else{
								lMin = bordersL[iL-1];
								lMax = bordersL[iL];
							}

							if(Jpsi->Pt()>ptMin && Jpsi->Pt()<ptMax && TMath::Abs(Jpsi->Rapidity())>rapMin && TMath::Abs(Jpsi->Rapidity())<rapMax && ctau>lMin && ctau<lMax){
								h_ctauerr_2012[iRAP][iPT][iL]->Fill(ctauerr);
								nEventsFilled++;
							}

						}

					}
				}


			}
			cout<<"2012: Finished Loop over "<<nentries2012<<" entries"<<endl;

			cout<<"2012: "<<nEventsPassingCuts<<" entries passing cuts"<<endl;
			cout<<"2012: "<<nEventsFilled<<" entries filled"<<endl;

			inFileData2011->Close();
			inFileData2012->Close();

			delete inFileData2011;
			delete inFileData2012;




		char xTitle[200];
		char yTitle[200];

		char tex[200];
		char tex1[200];
		char tex2[200];
		double yMax=0.;//=CosThDist_Original->GetMaximum()*expandY;

		bool Normalize=true;


		for(int iRAP = 0; iRAP < nRAP+1; iRAP++){
			for(int iPT = 0; iPT < nPT+1; iPT++){
				for(int iL = 0; iL < nL+1; iL++){

					//if(iRAP!=0 || iPT!=0) continue;

					cout<<h_ctauerr_2011[iRAP][iPT][iL]->GetName()<<" "<<h_ctauerr_2011[iRAP][iPT][iL]->GetEntries()<<endl;
					cout<<h_ctauerr_2012[iRAP][iPT][iL]->GetName()<<" "<<h_ctauerr_2012[iRAP][iPT][iL]->GetEntries()<<endl;


					sprintf(tex1,"2011 data");
					sprintf(tex2,"2012 data");

					sprintf(xTitle,"#sigma_{l}");
					sprintf(savename,"%s/Comp_ctauerr_rap%d_pT%d_l%d",saveDir, iRAP, iPT, iL);
					Plot1DHistoComp(h_ctauerr_2011[iRAP][iPT][iL], h_ctauerr_2012[iRAP][iPT][iL], xTitle, savename, tex1, tex2, yMax, Normalize);

					//cout<<"plotted "<<h_ctauerr_2011[iRAP][iPT][iL]->GetName()<<endl;

				}
			}
		}



		for(int iRAP = 0; iRAP < nRAP+1; iRAP++){
			for(int iPT = 0; iPT < nPT+1; iPT++){
				for(int iL = 0; iL < nL+1; iL++){

					if(iRAP!=0 || iPT!=0) continue;
					cout<<h_ctauerr_2011[iRAP][iPT][iL]->GetName()<<" "<<h_ctauerr_2011[iRAP][iPT][iL]->GetEntries()<<endl;
					cout<<h_ctauerr_2012[iRAP][iPT][iL]->GetName()<<" "<<h_ctauerr_2012[iRAP][iPT][iL]->GetEntries()<<endl;

				}
			}
		}
		for(int iRAP = 0; iRAP < nRAP+1; iRAP++){
			for(int iPT = 0; iPT < nPT+1; iPT++){
				for(int iL = 0; iL < nL+1; iL++){

					if(iL!=0 || iPT!=0) continue;
					cout<<h_ctauerr_2011[iRAP][iPT][iL]->GetName()<<" "<<h_ctauerr_2011[iRAP][iPT][iL]->GetEntries()<<endl;
					cout<<h_ctauerr_2012[iRAP][iPT][iL]->GetName()<<" "<<h_ctauerr_2012[iRAP][iPT][iL]->GetEntries()<<endl;

				}
			}
		}
		for(int iRAP = 0; iRAP < nRAP+1; iRAP++){
			for(int iPT = 0; iPT < nPT+1; iPT++){
				for(int iL = 0; iL < nL+1; iL++){

					if(iRAP!=0 || iL!=0) continue;
					cout<<h_ctauerr_2011[iRAP][iPT][iL]->GetName()<<" "<<h_ctauerr_2011[iRAP][iPT][iL]->GetEntries()<<endl;
					cout<<h_ctauerr_2012[iRAP][iPT][iL]->GetName()<<" "<<h_ctauerr_2012[iRAP][iPT][iL]->GetEntries()<<endl;

				}
			}
		}

		sprintf(savename,"%s/CtauErrModel_histograms.root",saveDir);
		TFile *outfile = new TFile(savename,"UPDATE");

		for(int iRAP = 0; iRAP < nRAP+1; iRAP++){
			for(int iPT = 0; iPT < nPT+1; iPT++){
				for(int iL = 0; iL < nL+1; iL++){

					outfile->cd();
					h_ctauerr_2011[iRAP][iPT][iL]->Write();
					h_ctauerr_2012[iRAP][iPT][iL]->Write();

				}
			}
		}


		outfile->Write();
		outfile->Close();
		delete outfile;
		outfile = NULL;

		for(int iRAP = 0; iRAP < nRAP+1; iRAP++){
			for(int iPT = 0; iPT < nPT+1; iPT++){
				for(int iL = 0; iL < nL+1; iL++){

					delete h_ctauerr_2011[iRAP][iPT][iL];
					delete h_ctauerr_2012[iRAP][iPT][iL];

				}
			}
		}




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

