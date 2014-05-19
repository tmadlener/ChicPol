/*
 * CutStudiesMCrec.C
 *
 *  Created on: May 13, 2014
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


#include "../../interface/commonVar.h"

using namespace onia;

void Plot2DHisto(TH2D* plothist, char xTitle[200], char yTitle[200], char savename[200], char tex[200]);
void Plot1DHisto(TH1D* plothist, char Title[200], char savename[200], char tex[200], double yMax);
void Plot1DHistoComp(TH1D* plothist_original, TH1D* plothist_cut, char Title[200], char savename[200], char tex1[200], char tex2[200], double yMax, bool Normalize);
void CalcCosthPhi(TLorentzVector* lepP, TLorentzVector* lepN, vector<double> &costh, vector<double> &phi);

void CutStudiesMCrec(){

	  bool useGenFile=false;
	  bool useRecFile=false;
	  bool useDataFile=true;


	  gROOT->Reset();
	  gROOT->SetBatch();

	  char saveDir[200];
	  char PlotID[200];
	  char savename[200];
	  sprintf(saveDir,"Figures/CutStudiesMCrec");
	  gSystem->mkdir(saveDir);
	  sprintf(PlotID,"2014May13");
	  sprintf(saveDir,"%s/%s",saveDir,PlotID);
	  gSystem->mkdir(saveDir);

	  const double pbeam_ = 7000.; // exact number irrelevant as long as pbeam >> Mprot
	  const double Mprot_ = 0.9382720;
	  const double Mlepton_ = 0.10566;  // (muon)
	  const double gPI_ = TMath::Pi();
	  const double Ebeam_ = sqrt( pbeam_*pbeam_ + Mprot_*Mprot_ );
	  TLorentzVector beam1_LAB_( 0., 0., pbeam_, Ebeam_ );
	  TLorentzVector beam2_LAB_( 0., 0., -pbeam_, Ebeam_ );

		const int nCuts=1;

		int nBinsCT=100;
		int nBinsPH=100;
		int nBinsCT2D=25;
		int nBinsPH2D=25;
		int nBinsPsiMass=100;
		int nBinsChiMass=100;

		TH2D* CosThPhiDist_Original;
		TH2D* CosThPhiDist_Cut;
		TH1D* CosThDist_Original;
		TH1D* CosThDist_Cut;
		TH1D* PhiDist_Original;
		TH1D* PhiDist_Cut;
		TH1D* PsiMassDist_Original;
		TH1D* PsiMassDist_Cut;
		TH1D* ChiMassDist_Original;
		TH1D* ChiMassDist_Cut;
		TH1D* ChiKVFMassDist_Original;

		CosThPhiDist_Original = new TH2D("CosThPhiDist_Original","CosThPhiDist_Original",nBinsCT2D, onia::cosTMin, onia::cosTMax, nBinsPH2D, onia::phiPolMin, onia::phiPolMax);
		CosThDist_Original = new TH1D("CosThDist_Original","CosThDist_Original",nBinsCT, onia::cosTMin, onia::cosTMax);
		PhiDist_Original = new TH1D("PhiDist_Original","PhiDist_Original",nBinsPH, onia::phiPolMin, onia::phiPolMax);
		PsiMassDist_Original = new TH1D("PsiMassDist_Original","PsiMassDist_Original",nBinsPsiMass, onia::massMin, onia::massMax);
		ChiMassDist_Original = new TH1D("ChiMassDist_Original","ChiMassDist_Original",nBinsChiMass, onia::chimassMin, onia::chimassMax);
		ChiKVFMassDist_Original = new TH1D("ChiKVFMassDist_Original","ChiKVFMassDist_Original",nBinsChiMass, onia::chimassMin, onia::chimassMax);

		CosThPhiDist_Cut = new TH2D(Form("CosThPhiDist_Cut"),Form("CosThPhiDist_Cut"),nBinsCT2D, onia::cosTMin, onia::cosTMax, nBinsPH2D, onia::phiPolMin, onia::phiPolMax);
		CosThDist_Cut = new TH1D(Form("CosThDist_Cut"), Form("CosThDist_Cut"),nBinsCT, onia::cosTMin, onia::cosTMax);
		PhiDist_Cut = new TH1D(Form("PhiDist_Cut"), Form("PhiDist_Cut"),nBinsPH, onia::phiPolMin, onia::phiPolMax);
		PsiMassDist_Cut = new TH1D("PsiMassDist_Cut","PsiMassDist_Cut",nBinsPsiMass, onia::massMin, onia::massMax);
		ChiMassDist_Cut = new TH1D("ChiMassDist_Cut","ChiMassDist_Cut",nBinsChiMass, onia::chimassMin, onia::chimassMax);



	  TFile* inFileGen = new TFile("/scratch/knuenz/Polarization/RootInput/ChicPol/chic_rootuple_MC_15M_gen.root");
	  TFile* inFileRec = new TFile("/scratch/knuenz/Polarization/RootInput/ChicPol/chic_rootuple_MC_15M_sel.root");
	  TFile* inFileData = new TFile("/scratch/knuenz/Polarization/RootInput/ChicPol/chic_rootuple_subFeb2014.root");


	  char savename_pre[200];
	  if(useGenFile) sprintf(savename_pre,"MCgen");
	  else if(useRecFile) sprintf(savename_pre,"MCrec");
	  else if(useDataFile) sprintf(savename_pre,"Data");

	  TTree* Tree;

	  if(useGenFile) Tree=(TTree*)inFileGen->Get("rootuple/GenParticlesTree");
	  else if(useRecFile) Tree=(TTree*)inFileRec->Get("rootuple/chicTree");
	  else if(useDataFile) Tree=(TTree*)inFileData->Get("rootuple/chicTree");

		Long64_t nentries = Tree->GetEntries();
		Long64_t nb = 0;

		std::cout << "number of entries = " << nentries << std::endl;

		TLorentzVector  *chi=0;
		TLorentzVector  *chi_rf=0;
		TLorentzVector  *Jpsi=0;
		TLorentzVector  *lepP=0;
		TLorentzVector  *lepN=0;
		TLorentzVector  *photon=0;
		double probFit1S;

		if(useGenFile){
			Tree->SetBranchAddress("chic_p4", &chi);
			Tree->SetBranchAddress("Jpsi_p4", &Jpsi);
			Tree->SetBranchAddress("photon_p4", &photon);
			Tree->SetBranchAddress("muP_p4", &lepP);
			Tree->SetBranchAddress("muM_p4", &lepN);
		}

		else if(useRecFile){
			Tree->SetBranchAddress("rf1S_chi_p4", &chi_rf);
			Tree->SetBranchAddress("chi_p4", &chi);
			Tree->SetBranchAddress("dimuon_p4", &Jpsi);
			Tree->SetBranchAddress("photon_p4", &photon);
			Tree->SetBranchAddress("muonP_p4", &lepP);
			Tree->SetBranchAddress("muonN_p4", &lepN);
			Tree->SetBranchAddress("probFit1S", &probFit1S);
		}
		else if(useDataFile){
			Tree->SetBranchAddress("rf1S_chi_p4", &chi_rf);
			Tree->SetBranchAddress("chi_p4", &chi);
			Tree->SetBranchAddress("dimuon_p4", &Jpsi);
			Tree->SetBranchAddress("photon_p4", &photon);
			Tree->SetBranchAddress("muonP_p4", &lepP);
			Tree->SetBranchAddress("muonN_p4", &lepN);
			Tree->SetBranchAddress("probFit1S", &probFit1S);
		}

		//nentries=1000;

		//loop over the events
		for (Long64_t jentry=0; jentry<nentries; jentry++) {

			//std::cout << "event " << jentry << " of " << nentries << std::endl;
			if(jentry % 100000 == 0) std::cout << "event " << jentry << " of " << nentries << std::endl;

			//Long64_t ientry = LoadTree(jentry);
			//if (ientry < 0) break;
			nb = Tree->GetEntry(jentry);

			vector<double> costh(3,0);
			vector<double> phi(3,0);
			CalcCosthPhi(lepP, lepN, costh, phi);



		    CosThPhiDist_Original->Fill(costh[2], phi[2]);
		    CosThDist_Original->Fill(costh[2]);
		    PhiDist_Original->Fill(phi[2]);
		    PsiMassDist_Original->Fill(Jpsi->M());
		    ChiMassDist_Original->Fill(chi->M()-Jpsi->M()+onia::MpsiPDG);
		    ChiKVFMassDist_Original->Fill(chi_rf->M());

		    if(probFit1S>0.01){
				CosThPhiDist_Cut->Fill(costh[2], phi[2]);
				CosThDist_Cut->Fill(costh[2]);
				PhiDist_Cut->Fill(phi[2]);
			    PsiMassDist_Cut->Fill(Jpsi->M());
			    ChiMassDist_Cut->Fill(chi->M()-Jpsi->M()+onia::MpsiPDG);
		    }

		}
		cout<<"Finished Loop over "<<nentries<<" entries"<<endl;
	  //Tree->Print();

		char xTitle[200];
		char yTitle[200];
		sprintf(xTitle,"cos#vartheta^{PX}");
		sprintf(yTitle,"#varphi^{PX} [deg]");

		char tex[200];
		char tex1[200];
		char tex2[200];
		double expandY=1.2;
		double yMax_Costh=CosThDist_Original->GetMaximum()*expandY;
		double yMax_Phi=PhiDist_Original->GetMaximum()*expandY;
		double yMax_PsiMass=PsiMassDist_Original->GetMaximum()*expandY;
		double yMax_ChiMass=ChiMassDist_Original->GetMaximum()*expandY;

		bool Normalize=false;

		sprintf(tex,"No cuts");
		sprintf(savename,"%s/%s_2D_CosThPhiDist_Original.pdf",saveDir, savename_pre);
		Plot2DHisto(CosThPhiDist_Original, xTitle, yTitle, savename, tex);
		sprintf(tex,"prob(KVF)>0.01");
		sprintf(savename,"%s/%s_2D_CosThPhiDist_Cut.pdf",saveDir, savename_pre);
		Plot2DHisto(CosThPhiDist_Cut, xTitle, yTitle, savename, tex);

		sprintf(tex1,"No cut");
		sprintf(tex2,"prob(KVF)>0.01");

		sprintf(savename,"%s/%s_1D_CosThDist_Comp.pdf",saveDir, savename_pre);
		Normalize=false;
		Plot1DHistoComp(CosThDist_Original, CosThDist_Cut, xTitle, savename, tex1, tex2, yMax_Costh, Normalize);
		sprintf(savename,"%s/%s_1D_PhiDist_Comp.pdf",saveDir, savename_pre);
		Normalize=false;
		Plot1DHistoComp(PhiDist_Original, PhiDist_Cut, yTitle, savename, tex1, tex2, yMax_Phi, Normalize);

		sprintf(savename,"%s/%s_1D_CosThDist_Comp_normalized.pdf",saveDir, savename_pre);
		Normalize=true;
		Plot1DHistoComp(CosThDist_Original, CosThDist_Cut, xTitle, savename, tex1, tex2, yMax_Costh, Normalize);
		sprintf(savename,"%s/%s_1D_PhiDist_Comp_normalized.pdf",saveDir, savename_pre);
		Normalize=true;
		Plot1DHistoComp(PhiDist_Original, PhiDist_Cut, yTitle, savename, tex1, tex2, yMax_Phi, Normalize);


		sprintf(xTitle,"M^{#psi}");
		sprintf(savename,"%s/%s_1D_PsiMassDist_Comp.pdf",saveDir, savename_pre);
		Normalize=false;
		Plot1DHistoComp(PsiMassDist_Original, PsiMassDist_Cut, xTitle, savename, tex1, tex2, yMax_PsiMass, Normalize);
		sprintf(xTitle,"M^{#chi}");
		sprintf(savename,"%s/%s_1D_ChiMassDist_Comp.pdf",saveDir, savename_pre);
		Normalize=false;
		Plot1DHistoComp(ChiMassDist_Original, ChiMassDist_Cut, xTitle, savename, tex1, tex2, yMax_ChiMass, Normalize);




		sprintf(tex1,"M^{Q}");
		sprintf(tex2,"M^{KVF}");
		sprintf(xTitle,"M^{#chi}");
		sprintf(savename,"%s/%s_1D_ChiMass_vs_KVFMass_Comp.pdf",saveDir, savename_pre);
		Normalize=false;
		Plot1DHistoComp(ChiMassDist_Original, ChiKVFMassDist_Original, xTitle, savename, tex1, tex2, yMax_ChiMass, Normalize);


		sprintf(tex1,"No cut");
		sprintf(tex2,"prob(KVF)>0.01");
		sprintf(xTitle,"M^{#psi}");
		sprintf(savename,"%s/%s_1D_PsiMassDist_Comp_normalized.pdf",saveDir, savename_pre);
		Normalize=true;
		Plot1DHistoComp(PsiMassDist_Original, PsiMassDist_Cut, xTitle, savename, tex1, tex2, yMax_PsiMass, Normalize);
		sprintf(xTitle,"M^{#chi}");
		sprintf(savename,"%s/%s_1D_ChiMassDist_Comp_normalized.pdf",saveDir, savename_pre);
		Normalize=true;
		Plot1DHistoComp(ChiMassDist_Original, ChiMassDist_Cut, xTitle, savename, tex1, tex2, yMax_ChiMass, Normalize);


		sprintf(tex1,"M^{Q}");
		sprintf(tex2,"M^{KVF}");
		sprintf(xTitle,"M^{#chi}");
		sprintf(savename,"%s/%s_1D_ChiMass_vs_KVFMass_Comp_normalized.pdf",saveDir, savename_pre);
		Normalize=true;
		Plot1DHistoComp(ChiMassDist_Original, ChiKVFMassDist_Original, xTitle, savename, tex1, tex2, yMax_ChiMass, Normalize);


		delete CosThPhiDist_Original;
		delete CosThDist_Original;
		delete PhiDist_Original;
		delete PsiMassDist_Original;
		delete ChiKVFMassDist_Original;

		delete CosThPhiDist_Cut;
		delete CosThDist_Cut;
		delete PhiDist_Cut;
		delete PsiMassDist_Cut;




//return;


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

	double left=0.15, top=0.935, textSize=0.04;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	latex->DrawLatex(left,top,tex);

	c1->SaveAs(savename);

}

void Plot1DHisto(TH1D* plothist, char Title[200], char savename[200], char tex[200], double yMax){

	TGaxis::SetMaxDigits(3);

	TCanvas *c1 = new TCanvas("", "", 1200, 1000);
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gPad->SetFillColor(kWhite);
	gPad->SetLeftMargin(0.15);
	gPad->SetRightMargin(0.05);
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

}

void Plot1DHistoComp(TH1D* plothist_original, TH1D* plothist_cut, char Title[200], char savename[200], char tex1[200], char tex2[200], double yMax, bool Normalize){

	TGaxis::SetMaxDigits(3);

	TCanvas *c1 = new TCanvas("", "", 1200, 1000);
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gPad->SetFillColor(kWhite);
	gPad->SetLeftMargin(0.15);
	gPad->SetRightMargin(0.05);
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
	double step=textSize*1.5;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	latex->SetTextColor(kGreen+2);
	latex->DrawLatex(left,top,tex1);
	left+=0.2;
	latex->SetTextColor(kRed+2);
	latex->DrawLatex(left,top,tex2);

	c1->SaveAs(savename);

}

void CalcCosthPhi(TLorentzVector* lepP, TLorentzVector* lepN, vector<double> &costh, vector<double> &phi){

	  const double pbeam_ = 7000.; // exact number irrelevant as long as pbeam >> Mprot
	  const double Mprot_ = 0.9382720;
	  const double Mlepton_ = 0.10566;  // (muon)
	  const double gPI_ = TMath::Pi();
	  const double Ebeam_ = sqrt( pbeam_*pbeam_ + Mprot_*Mprot_ );
	  TLorentzVector beam1_LAB_( 0., 0., pbeam_, Ebeam_ );
	  TLorentzVector beam2_LAB_( 0., 0., -pbeam_, Ebeam_ );


    double lepP_pT  = lepP->Pt();
    double lepN_pT  = lepN->Pt();

    double lepP_eta = lepP->PseudoRapidity();
    double lepN_eta = lepN->PseudoRapidity();

    // dilepton 4-vector:

    TLorentzVector dilepton = *lepP + *lepN;
    double pT   = dilepton.Pt();
    double rap  = dilepton.Rapidity();
    double mass = dilepton.M();

    // calculation of decay angles in three polarization frames

    // reference directions to calculate angles:

    TVector3 lab_to_dilep = -dilepton.BoostVector();

    TLorentzVector beam1_DILEP = beam1_LAB_;
    beam1_DILEP.Boost(lab_to_dilep);         // beam1 in the dilepton rest frame
    TLorentzVector beam2_DILEP = beam2_LAB_;
    beam2_DILEP.Boost(lab_to_dilep);         // beam2 in the dilepton rest frame

    TVector3 beam1_direction     = beam1_DILEP.Vect().Unit();
    TVector3 beam2_direction     = beam2_DILEP.Vect().Unit();
    TVector3 dilep_direction     = dilepton.Vect().Unit();
    TVector3 beam1_beam2_bisect  = ( beam1_direction - beam2_direction ).Unit();


    // all polarization frames have the same Y axis = the normal to the plane formed by
    // the directions of the colliding hadrons:

    TVector3 Yaxis = ( beam1_direction.Cross( beam2_direction ) ).Unit();

    // flip of y axis with rapidity:

    if ( rap < 0. ) Yaxis = - Yaxis;

    TVector3 perpendicular_to_beam = ( beam1_beam2_bisect.Cross( Yaxis ) ).Unit();


    // positive lepton in the dilepton rest frame:

    TLorentzVector lepton_DILEP = *lepP;
    lepton_DILEP.Boost(lab_to_dilep);

    // CS frame angles:

    TVector3 newZaxis = beam1_beam2_bisect;
    TVector3 newYaxis = Yaxis;
    TVector3 newXaxis = newYaxis.Cross( newZaxis );

    TRotation rotation;
    rotation.SetToIdentity();
    rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation.Invert();  // transforms coordinates from the "xyz" frame to the new frame
    TVector3 lepton_DILEP_rotated = lepton_DILEP.Vect();
    lepton_DILEP_rotated.Transform(rotation);

    double costh_CS = lepton_DILEP_rotated.CosTheta();
    double phi_CS   = lepton_DILEP_rotated.Phi() * 180. / gPI_;
    double phith_CS;
    if ( costh_CS < 0. ) phith_CS = phi_CS - 135.;
    if ( costh_CS > 0. ) phith_CS = phi_CS - 45.;
    if ( phith_CS < -180. ) phith_CS = 360. + phith_CS;


    // HELICITY frame angles:

    newZaxis = dilep_direction;
    newYaxis = Yaxis;
    newXaxis = newYaxis.Cross( newZaxis );

    rotation.SetToIdentity();
    rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation.Invert();
    lepton_DILEP_rotated = lepton_DILEP.Vect();
    lepton_DILEP_rotated.Transform(rotation);

    double costh_HX = lepton_DILEP_rotated.CosTheta();
    double phi_HX   = lepton_DILEP_rotated.Phi() * 180. / gPI_;
    double phith_HX;
    if ( costh_HX < 0. ) phith_HX = phi_HX - 135.;
    if ( costh_HX > 0. ) phith_HX = phi_HX - 45.;
    if ( phith_HX < -180. ) phith_HX = 360. + phith_HX;

		//PhiHX test
		//if(PhiHX_test) { if(phi_HX>80. || phi_HX<91.) continue; }

    // PERPENDICULAR HELICITY frame angles:

    newZaxis = perpendicular_to_beam;
    newYaxis = Yaxis;
    newXaxis = newYaxis.Cross( newZaxis );

    rotation.SetToIdentity();
    rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation.Invert();
    lepton_DILEP_rotated = lepton_DILEP.Vect();
    lepton_DILEP_rotated.Transform(rotation);

    double costh_PX = lepton_DILEP_rotated.CosTheta();
    double phi_PX   = lepton_DILEP_rotated.Phi() * 180. / gPI_;
    double phith_PX;
    if ( costh_PX < 0. ) phith_PX = phi_PX - 135.;
    if ( costh_PX > 0. ) phith_PX = phi_PX - 45.;
    if ( phith_PX < -180. ) phith_PX = 360. + phith_PX;

    costh[0]=costh_CS;
    costh[1]=costh_HX;
    costh[2]=costh_PX;

    phi[0]=phi_CS;
    phi[1]=phi_HX;
    phi[2]=phi_PX;

}
