/*
 * PhotonEffStudies.C
 *
 *  Created on: Apr 22, 2014
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

void Plot2DHisto(TH2D* plothist, char xTitle[200], char yTitle[200], char savename[200], char savenamelog[200], char tex[200], double yMax, bool addLog);
void Plot1DHisto(TH1D* plothist, char Title[200], char savename[200], char savenamelog[200], char tex[200], double yMax, bool addLog);
void CalcCosthPhi(TLorentzVector* lepP, TLorentzVector* lepN, vector<double> &costh, vector<double> &phi);

void PhotonEffStudies(){

	  gROOT->Reset();
	  gROOT->SetBatch();

	  char saveDir[200];
	  char PlotID[200];
	  char savename[200];
	  char savenamelog[200];
	  sprintf(saveDir,"Figures/PhotonEffStudies");
	  sprintf(PlotID,"2014April22");
	  sprintf(saveDir,"%s/%s",saveDir,PlotID);
	  gSystem->mkdir(saveDir);



		int nBinsCT=100;
		int nBinsPH=100;
		int nBinsCT2D=25;
		int nBinsPH2D=25;

		int nBinsPhotonPt=100;
		double min_PhotonPt=0;
		double max_PhotonPt=6.;

		int nBinsPhotonEta=100;
		double min_PhotonEta=-2.4;
		double max_PhotonEta=2.4;

		int nBinsPhotonPt2D=25;
		int nBinsPhotonEta2D=25;


		TH2D* CosThPhiDist_Gen;
		TH2D* PhotonPtEtaDist_Gen;

		TH1D* CosThDist_Gen;
		TH1D* PhiDist_Gen;
		TH1D* PhotonPtDist_Gen;
		TH1D* PhotonEtaDist_Gen;

		CosThPhiDist_Gen = new TH2D("CosThPhiDist_Gen","CosThPhiDist_Gen",nBinsCT2D, onia::cosTMin, onia::cosTMax, nBinsPH2D, onia::phiPolMin, onia::phiPolMax);
		PhotonPtEtaDist_Gen = new TH2D("PhotonPtEtaDist_Gen","PhotonPtEtaDist_Gen",nBinsPhotonPt2D, min_PhotonPt, max_PhotonPt, nBinsPhotonEta2D, min_PhotonEta, max_PhotonEta);

		CosThDist_Gen = new TH1D("CosThDist_Gen","CosThDist_Gen",nBinsCT, onia::cosTMin, onia::cosTMax);
		PhiDist_Gen = new TH1D("PhiDist_Gen","PhiDist_Gen",nBinsPH, onia::phiPolMin, onia::phiPolMax);
		PhotonPtDist_Gen = new TH1D("PhotonPtDist_Gen","PhotonPtDist_Gen",nBinsPhotonPt, min_PhotonPt, max_PhotonPt);
		PhotonEtaDist_Gen = new TH1D("PhotonEtaDist_Gen","PhotonEtaDist_Gen",nBinsPhotonEta, min_PhotonEta, max_PhotonEta);



		TH2D* CosThPhiDist_Rec;
		TH2D* PhotonPtEtaDist_Rec;

		TH1D* CosThDist_Rec;
		TH1D* PhiDist_Rec;
		TH1D* PhotonPtDist_Rec;
		TH1D* PhotonEtaDist_Rec;

		CosThPhiDist_Rec = new TH2D("CosThPhiDist_Rec","CosThPhiDist_Rec",nBinsCT2D, onia::cosTMin, onia::cosTMax, nBinsPH2D, onia::phiPolMin, onia::phiPolMax);
		PhotonPtEtaDist_Rec = new TH2D("PhotonPtEtaDist_Rec","PhotonPtEtaDist_Rec",nBinsPhotonPt2D, min_PhotonPt, max_PhotonPt, nBinsPhotonEta2D, min_PhotonEta, max_PhotonEta);

		CosThDist_Rec = new TH1D("CosThDist_Rec","CosThDist_Rec",nBinsCT, onia::cosTMin, onia::cosTMax);
		PhiDist_Rec = new TH1D("PhiDist_Rec","PhiDist_Rec",nBinsPH, onia::phiPolMin, onia::phiPolMax);
		PhotonPtDist_Rec = new TH1D("PhotonPtDist_Rec","PhotonPtDist_Rec",nBinsPhotonPt, min_PhotonPt, max_PhotonPt);
		PhotonEtaDist_Rec = new TH1D("PhotonEtaDist_Rec","PhotonEtaDist_Rec",nBinsPhotonEta, min_PhotonEta, max_PhotonEta);




	  TFile* inFileGen = new TFile("/scratch/knuenz/Polarization/RootInput/ChicPol/chic_rootuple_MC_15M_gen.root");
	  TFile* inFileRec = new TFile("/scratch/knuenz/Polarization/RootInput/ChicPol/chic_rootuple_MC_15M_sel.root");

	  TTree* GenTree=(TTree*)inFileGen->Get("rootuple/GenParticlesTree");
	  TTree* RecTree=(TTree*)inFileRec->Get("rootuple/chicTree");

		double fracOfStatToRun=1.0;




		TLorentzVector  *chi_Gen=0;
		TLorentzVector  *Jpsi_Gen=0;
		TLorentzVector  *lepP_Gen=0;
		TLorentzVector  *lepN_Gen=0;
		TLorentzVector  *photon_Gen=0;

		Long64_t nentries_Gen = GenTree->GetEntries();
		Long64_t nb_Gen = 0;

		std::cout << "number of Gen entries = " << nentries_Gen << std::endl;

		GenTree->SetBranchAddress("chic_p4", &chi_Gen);
		GenTree->SetBranchAddress("Jpsi_p4", &Jpsi_Gen);
		GenTree->SetBranchAddress("photon_p4", &photon_Gen);
		GenTree->SetBranchAddress("muP_p4", &lepP_Gen);
		GenTree->SetBranchAddress("muM_p4", &lepN_Gen);

		nentries_Gen*=fracOfStatToRun;

		//loop over the events
		for (Long64_t jentry=0; jentry<nentries_Gen; jentry++) {

			if(jentry % 100000 == 0) std::cout << "event " << jentry << " of " << nentries_Gen << std::endl;

			nb_Gen = GenTree->GetEntry(jentry);

			vector<double> costh(3,0);
			vector<double> phi(3,0);
			CalcCosthPhi(lepP_Gen, lepN_Gen, costh, phi);

		    CosThPhiDist_Gen->Fill(costh[2], phi[2]);
		    CosThDist_Gen->Fill(costh[2]);
		    PhiDist_Gen->Fill(phi[2]);

		    PhotonPtEtaDist_Gen->Fill(photon_Gen->Pt(),photon_Gen->Eta());
		    PhotonPtDist_Gen->Fill(photon_Gen->Pt());
		    PhotonEtaDist_Gen->Fill(photon_Gen->Eta());

		}
		cout<<"Finished Loop over "<<nentries_Gen<<" Gen entries"<<endl;


		TLorentzVector  *chi_Rec=0;
		TLorentzVector  *Jpsi_Rec=0;
		TLorentzVector  *lepP_Rec=0;
		TLorentzVector  *lepN_Rec=0;
		TLorentzVector  *photon_Rec=0;

		Long64_t nentries_Rec = RecTree->GetEntries();
		Long64_t nb_Rec = 0;

		std::cout << "number of Rec entries = " << nentries_Rec << std::endl;

		RecTree->SetBranchAddress("chi_p4", &chi_Rec);
		RecTree->SetBranchAddress("dimuon_p4", &Jpsi_Rec);
		RecTree->SetBranchAddress("photon_p4", &photon_Rec);
		RecTree->SetBranchAddress("muonP_p4", &lepP_Rec);
		RecTree->SetBranchAddress("muonN_p4", &lepN_Rec);

		nentries_Rec*=fracOfStatToRun;

		//loop over the events
		for (Long64_t jentry=0; jentry<nentries_Rec; jentry++) {

			if(jentry % 100000 == 0) std::cout << "event " << jentry << " of " << nentries_Rec << std::endl;

			nb_Rec = RecTree->GetEntry(jentry);

			vector<double> costh(3,0);
			vector<double> phi(3,0);
			CalcCosthPhi(lepP_Rec, lepN_Rec, costh, phi);

		    CosThPhiDist_Rec->Fill(costh[2], phi[2]);
		    CosThDist_Rec->Fill(costh[2]);
		    PhiDist_Rec->Fill(phi[2]);

		    PhotonPtEtaDist_Rec->Fill(photon_Rec->Pt(),photon_Rec->Eta());
		    PhotonPtDist_Rec->Fill(photon_Rec->Pt());
		    PhotonEtaDist_Rec->Fill(photon_Rec->Eta());

		}
		cout<<"Finished Loop over "<<nentries_Rec<<" Rec entries"<<endl;



		CosThPhiDist_Rec->Sumw2(true);
		CosThDist_Gen->Sumw2(true);
		PhiDist_Gen->Sumw2(true);
		PhotonPtDist_Gen->Sumw2(true);
		PhotonEtaDist_Gen->Sumw2(true);
		PhotonPtEtaDist_Gen->Sumw2(true);

		CosThPhiDist_Rec->Sumw2(true);
		CosThDist_Rec->Sumw2(true);
		PhiDist_Rec->Sumw2(true);
		PhotonPtDist_Rec->Sumw2(true);
		PhotonEtaDist_Rec->Sumw2(true);
		PhotonPtEtaDist_Rec->Sumw2(true);


		CosThPhiDist_Rec->Divide(CosThPhiDist_Gen);
		CosThDist_Rec->Divide(CosThDist_Gen);
		PhiDist_Rec->Divide(PhiDist_Gen);
		PhotonPtDist_Rec->Divide(PhotonPtDist_Gen);
		PhotonEtaDist_Rec->Divide(PhotonEtaDist_Gen);
		PhotonPtEtaDist_Rec->Divide(PhotonPtEtaDist_Gen);


		char xTitle[200];
		char yTitle[200];
		sprintf(xTitle,"cos#vartheta^{PX}");
		sprintf(yTitle,"#varphi^{PX} [deg]");

		char tex[200];
		double expandY=1.2;
		double yMax=1.5;

		sprintf(tex,"");

		sprintf(savename,"%s/2D_Eff_CosThPhi.pdf",saveDir);
		yMax=CosThPhiDist_Rec->GetMaximum()*expandY;
		Plot2DHisto(CosThPhiDist_Rec, xTitle, yTitle, savename, savenamelog, tex, yMax, false);
		sprintf(savename,"%s/1D_Eff_CosTh.pdf",saveDir);
		yMax=CosThDist_Rec->GetMaximum()*expandY;
		Plot1DHisto(CosThDist_Rec, xTitle, savename, savenamelog, tex, yMax, false);
		sprintf(savename,"%s/1D_Eff_Phi.pdf",saveDir);
		yMax=PhiDist_Rec->GetMaximum()*expandY;
		Plot1DHisto(PhiDist_Rec, yTitle, savename, savenamelog, tex, yMax, false);

		sprintf(xTitle,"photon p_{T}");
		sprintf(yTitle,"photon #eta");

		sprintf(savename,"%s/1D_Eff_PhotonPt.pdf",saveDir);
		sprintf(savenamelog,"%s/1D_Eff_PhotonPt_log.pdf",saveDir);
		yMax=PhotonPtDist_Rec->GetMaximum()*expandY;
		Plot1DHisto(PhotonPtDist_Rec, xTitle, savename, savenamelog, tex, yMax, true);
		sprintf(savename,"%s/1D_Eff_PhotonEta.pdf",saveDir);
		sprintf(savenamelog,"%s/1D_Eff_PhotonEta_log.pdf",saveDir);
		yMax=PhotonEtaDist_Rec->GetMaximum()*expandY;
		Plot1DHisto(PhotonEtaDist_Rec, yTitle, savename, savenamelog, tex, yMax, true);
		sprintf(savename,"%s/2D_Eff_PhotonPtEta.pdf",saveDir);
		sprintf(savenamelog,"%s/2D_Eff_PhotonPtEta_log.pdf",saveDir);
		yMax=PhotonPtEtaDist_Rec->GetMaximum()*expandY;
		Plot2DHisto(PhotonPtEtaDist_Rec, xTitle, yTitle, savename, savenamelog, tex, yMax, true);









//return;


}


void Plot2DHisto(TH2D* plothist, char xTitle[200], char yTitle[200], char savename[200], char savenamelog[200], char tex[200], double yMax, bool addLog){

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
	plothist->SetMaximum(yMax);

	plothist->Draw("colz");

	double left=0.15, top=0.935, textSize=0.04;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	latex->DrawLatex(left,top,tex);

	c1->SaveAs(savename);

	if(addLog){
		yMax=5e-2;
		plothist->SetMinimum(5e-6);
		plothist->SetMaximum(yMax);
		c1->SetLogz(true);
		c1->SaveAs(savenamelog);
	}


}

void Plot1DHisto(TH1D* plothist, char Title[200], char savename[200], char savenamelog[200], char tex[200], double yMax, bool addLog){

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
	plothist->GetYaxis()->SetTitle("#chi Efficiency");

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

	if(addLog){
		yMax=5e-2;
		plothist->SetMinimum(5e-6);
		plothist->SetMaximum(yMax);
		c1->SetLogy(true);
		c1->SaveAs(savenamelog);
	}

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
