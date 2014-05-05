/*
 * PhotonCutStudies.C
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

void Plot2DHisto(TH2D* plothist, char xTitle[200], char yTitle[200], char savename[200], char tex[200]);
void Plot1DHisto(TH1D* plothist, char Title[200], char savename[200], char tex[200], double yMax);
void CalcCosthPhi(TLorentzVector* lepP, TLorentzVector* lepN, vector<double> &costh, vector<double> &phi);

void PhotonCutStudies(){

	  gROOT->Reset();
	  gROOT->SetBatch();

	  char saveDir[200];
	  char PlotID[200];
	  char savename[200];
	  sprintf(saveDir,"Figures/PhotonCutStudies");
	  sprintf(PlotID,"tmp_2014April22");
	  sprintf(saveDir,"%s/%s",saveDir,PlotID);
	  gSystem->mkdir(saveDir);

	  const double pbeam_ = 7000.; // exact number irrelevant as long as pbeam >> Mprot
	  const double Mprot_ = 0.9382720;
	  const double Mlepton_ = 0.10566;  // (muon)
	  const double gPI_ = TMath::Pi();
	  const double Ebeam_ = sqrt( pbeam_*pbeam_ + Mprot_*Mprot_ );
	  TLorentzVector beam1_LAB_( 0., 0., pbeam_, Ebeam_ );
	  TLorentzVector beam2_LAB_( 0., 0., -pbeam_, Ebeam_ );

		const int nCuts=6;

		int nBinsCT=100;
		int nBinsPH=100;
		int nBinsCT2D=25;
		int nBinsPH2D=25;

		TH2D* CosThPhiDist_Original;
		TH2D* CosThPhiDist_Cut[nCuts];
		TH1D* CosThDist_Original;
		TH1D* CosThDist_Cut[nCuts];
		TH1D* PhiDist_Original;
		TH1D* PhiDist_Cut[nCuts];

		CosThPhiDist_Original = new TH2D("CosThPhiDist_Original","CosThPhiDist_Original",nBinsCT2D, onia::cosTMin, onia::cosTMax, nBinsPH2D, onia::phiPolMin, onia::phiPolMax);
		CosThDist_Original = new TH1D("CosThDist_Original","CosThDist_Original",nBinsCT, onia::cosTMin, onia::cosTMax);
		PhiDist_Original = new TH1D("PhiDist_Original","PhiDist_Original",nBinsPH, onia::phiPolMin, onia::phiPolMax);

		for (int iCut=0; iCut<nCuts; iCut++) {

			CosThPhiDist_Cut[iCut] = new TH2D(Form("CosThPhiDist_Cut%d",iCut+1),Form("CosThPhiDist_Cut%d",iCut+1),nBinsCT2D, onia::cosTMin, onia::cosTMax, nBinsPH2D, onia::phiPolMin, onia::phiPolMax);
			CosThDist_Cut[iCut] = new TH1D(Form("CosThDist_Cut%d",iCut+1), Form("CosThDist_Cut%d",iCut+1),nBinsCT, onia::cosTMin, onia::cosTMax);
			PhiDist_Cut[iCut] = new TH1D(Form("PhiDist_Cut%d",iCut+1), Form("PhiDist_Cut%d",iCut+1),nBinsPH, onia::phiPolMin, onia::phiPolMax);

		}


	  TFile* inFileGen = new TFile("/scratch/knuenz/Polarization/RootInput/ChicPol/chic_rootuple_MC_15M_gen.root");
	  TFile* inFileRec = new TFile("/scratch/knuenz/Polarization/RootInput/ChicPol/chic_rootuple_MC_15M_sel.root");

	  TTree* GenTree=(TTree*)inFileGen->Get("rootuple/GenParticlesTree");


		Long64_t nentries = GenTree->GetEntries();
		Long64_t nb = 0;

		std::cout << "number of entries = " << nentries << std::endl;

		TLorentzVector  *chi=0;
		TLorentzVector  *Jpsi=0;
		TLorentzVector  *lepP=0;
		TLorentzVector  *lepN=0;
		TLorentzVector  *photon=0;

		GenTree->SetBranchAddress("chic_p4", &chi);
		GenTree->SetBranchAddress("Jpsi_p4", &Jpsi);
		GenTree->SetBranchAddress("photon_p4", &photon);
		GenTree->SetBranchAddress("muP_p4", &lepP);
		GenTree->SetBranchAddress("muM_p4", &lepN);


		nentries=1000;

		//loop over the events
		for (Long64_t jentry=0; jentry<nentries; jentry++) {

			//std::cout << "event " << jentry << " of " << nentries << std::endl;
			if(jentry % 100000 == 0) std::cout << "event " << jentry << " of " << nentries << std::endl;

			//Long64_t ientry = LoadTree(jentry);
			//if (ientry < 0) break;
			nb = GenTree->GetEntry(jentry);

			vector<double> costh(3,0);
			vector<double> phi(3,0);
			CalcCosthPhi(lepP, lepN, costh, phi);



		    CosThPhiDist_Original->Fill(costh[2], phi[2]);
		    CosThDist_Original->Fill(costh[2]);
		    PhiDist_Original->Fill(phi[2]);

		    if(photon->Pt()>0.5){
				CosThPhiDist_Cut[0]->Fill(costh[2], phi[2]);
				CosThDist_Cut[0]->Fill(costh[2]);
				PhiDist_Cut[0]->Fill(phi[2]);
		    }
		    if(photon->Pt()>1.){
				CosThPhiDist_Cut[1]->Fill(costh[2], phi[2]);
				CosThDist_Cut[1]->Fill(costh[2]);
				PhiDist_Cut[1]->Fill(phi[2]);
		    }
		    if(photon->Pt()>2.){
				CosThPhiDist_Cut[2]->Fill(costh[2], phi[2]);
				CosThDist_Cut[2]->Fill(costh[2]);
				PhiDist_Cut[2]->Fill(phi[2]);
		    }
		    if(TMath::Abs(photon->Eta())<1.5){
				CosThPhiDist_Cut[3]->Fill(costh[2], phi[2]);
				CosThDist_Cut[3]->Fill(costh[2]);
				PhiDist_Cut[3]->Fill(phi[2]);
		    }
		    if(TMath::Abs(photon->Eta())<1.){
				CosThPhiDist_Cut[4]->Fill(costh[2], phi[2]);
				CosThDist_Cut[4]->Fill(costh[2]);
				PhiDist_Cut[4]->Fill(phi[2]);
		    }
		    if(TMath::Abs(photon->Eta())<0.5){
				CosThPhiDist_Cut[5]->Fill(costh[2], phi[2]);
				CosThDist_Cut[5]->Fill(costh[2]);
				PhiDist_Cut[5]->Fill(phi[2]);
		    }

		}
		cout<<"Finished Loop over "<<nentries<<" entries"<<endl;
	  //GenTree->Print();

		char xTitle[200];
		char yTitle[200];
		sprintf(xTitle,"cos#vartheta^{PX}");
		sprintf(yTitle,"#varphi^{PX} [deg]");

		char tex[200];
		double expandY=1.2;
		double yMax_Costh=CosThDist_Original->GetMaximum()*expandY;
		double yMax_Phi=PhiDist_Original->GetMaximum()*expandY;

		sprintf(tex,"No cuts");

		sprintf(savename,"%s/2D_CosThPhiDist_Original.pdf",saveDir);
		Plot2DHisto(CosThPhiDist_Original, xTitle, yTitle, savename, tex);
		sprintf(savename,"%s/1D_CosThDist_Original.pdf",saveDir);
		Plot1DHisto(CosThDist_Original, xTitle, savename, tex, yMax_Costh);
		sprintf(savename,"%s/1D_PhiDist_Original.pdf",saveDir);
		Plot1DHisto(PhiDist_Original, yTitle, savename, tex, yMax_Phi);

		for (int iCut=0; iCut<nCuts; iCut++) {

			if(iCut==0) sprintf(tex,"Photon p_{T} > 0.5 GeV");
			if(iCut==1) sprintf(tex,"Photon p_{T} > 1.0 GeV");
			if(iCut==2) sprintf(tex,"Photon p_{T} > 2.0 GeV");
			if(iCut==3) sprintf(tex,"Photon |#eta| < 1.5");
			if(iCut==4) sprintf(tex,"Photon |#eta| < 1.0");
			if(iCut==5) sprintf(tex,"Photon |#eta| < 0.5");

			yMax_Costh=CosThDist_Cut[iCut]->GetMaximum()*expandY;
			yMax_Phi=PhiDist_Cut[iCut]->GetMaximum()*expandY;

			sprintf(savename,"%s/2D_CosThPhiDist_photonCut%d.pdf",saveDir,iCut+1);
			Plot2DHisto(CosThPhiDist_Cut[iCut], xTitle, yTitle, savename, tex);
			sprintf(savename,"%s/1D_CosThDist_photonCut%d.pdf",saveDir,iCut+1);
			Plot1DHisto(CosThDist_Cut[iCut], xTitle, savename, tex, yMax_Costh);
			sprintf(savename,"%s/1D_PhiDist_photonCut%d.pdf",saveDir,iCut+1);
			Plot1DHisto(PhiDist_Cut[iCut], yTitle, savename, tex, yMax_Phi);

		}





		delete CosThPhiDist_Original;
		delete CosThDist_Original;
		delete PhiDist_Original;

		for (int iCut=0; iCut<nCuts; iCut++) {

			delete CosThPhiDist_Cut[iCut];
			delete CosThDist_Cut[iCut];
			delete PhiDist_Cut[iCut];

		}



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
