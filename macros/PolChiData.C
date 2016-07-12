#define PolChiData_cxx
#include "PolChiData.h"
#include "commonVar.h"
#include "effsAndCuts.h"

#include "TLorentzVector.h"
#include <TH2.h>
#include <TCanvas.h>

#include <string>
#include <iostream>
#include <sstream>

TH1F *Reco_StatEv;
TH1F *Reco_Onia_mass[onia::kNbPTMaxBins+1][onia::kNbRapForPTBins+1];
TH2F *Reco_Onia_rap_pT;
TH1F *Reco_Onia_pt[onia::kNbRapForPTBins+1];
TH1F *Reco_Onia_rap[onia::kNbPTMaxBins+1];

TTree *treeOut;
TLorentzVector *lepP, *photon, *jpsi, *chic, *lepN;
TLorentzVector *lepP_rf, *photon_rf, *jpsi_rf, *chic_rf, *lepN_rf;


void PolChiData::Loop(int nState, bool rejectCowboys, int FidCuts, bool MC, bool RequestTrigger, bool removeEta0p2_0p3,
                      bool cutDeltaREllDpt, bool correctCtau, bool useRefittedChic, bool cutDimuon10Gev, double muAccShift) {

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  Long64_t cutAtRecEvent = nentries;
  Long64_t count = 0;
  Long64_t nb = 0;

  std::cout << "number of entries = " << nentries << std::endl;

  // create branches for certain output variables
  double mQ = 0;
  double jpsict = 0;
  double jpsictErr = 0;
  double numPV=0;
  double probKVF=0;

  treeOut->Branch("Jpsict", &jpsict, "Jpsict/D");
  treeOut->Branch("JpsictErr", &jpsictErr, "JpsictErr/D");
  treeOut->Branch("mQ", &mQ, "mQ/D");
  treeOut->Branch("probKVF", &probFit1S, "probKVF/D");

  //fChain->Print();

  //nentries=1000;


  //TFile* infile = new TFile("/afs/hephy.at/scratch/k/knuenz/tmp/ChicPol/JpsictErr2011_15_20.root", "READ");
  //cout<<"opened file"<<endl;
  //TH1F* h_JpsictErr=(TH1F*)infile->Get("h");
  //h_JpsictErr->Print();
  //cout<<"opened hist"<<endl;

  //loop over the events
  for (Long64_t jentry=0; jentry<nentries; jentry++) {


    if(jentry % 100000 == 0) std::cout << "event " << jentry << " of " << nentries << std::endl;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);

    //fChain->GetEvent(1);


    // Define VARIABLES
    //TODO: change variable definitions once real file is available

    //4-momentum vectors
    lepP = muonP_p4;
    lepN = muonN_p4;
    jpsi = dimuon_p4;
    photon = photon_p4;
    chic = chi_p4;

    lepP_rf = rf1S_muonP_p4;
    lepN_rf = rf1S_muonN_p4;
    jpsi_rf = rf1S_dimuon_p4;
    photon_rf = rf1S_photon_p4;
    chic_rf = rf1S_chi_p4;

    // store the (redundant) non-refitted chic in the chic_rf Branch if desired in order to have it consistently used in later steps
    if (!useRefittedChic) {
      chic_rf = chi_p4;
    }

    //lepP -> 	SetXYZM(1.,1.,1.,1.);
    //lepN -> 	SetXYZM(1.,1.,1.,1.);
    //jpsi -> 	SetXYZM(1.,1.,1.,1.);
    //photon -> 	SetXYZM(1.,1.,1.,1.);
    //chic -> 	SetXYZM(1.,1.,1.,1.);

    //MuonVars
    double etaMuPos = lepP->PseudoRapidity();
    double etaMuNeg = lepN->PseudoRapidity();
    double pTMuPos = lepP->Pt();
    double pTMuNeg = lepN->Pt();
    double pMuPos = lepP->P();
    double pMuNeg = lepN->P();

    //DimuonVars
    double onia_mass = jpsi->M();
    double onia_pt = jpsi->Pt();
    double onia_P = jpsi->P();
    double onia_eta = jpsi->PseudoRapidity();
    double onia_rap = jpsi->Rapidity();
    double onia_phi = jpsi->Phi();
    double deltaPhi = lepN->Phi() - lepP->Phi();;
    if(deltaPhi > TMath::Pi()) deltaPhi -= 2.*TMath::Pi();
    else if(deltaPhi < -TMath::Pi()) deltaPhi += 2.*TMath::Pi();

    // select mass window, mass-resolution dependence on rapidity, values in MeV, valid for Ups1S 2012 data
    double massMin = 0, massMax = 0;
    double p0 = +62.62;
    double p1 = +56.30;
    double p2 = -20.77;
    //sigma = p0+p1*|y|^2+p2*|y|^3
    double sigmaRap = p0+p1*TMath::Abs(onia_rap)*TMath::Abs(onia_rap)+p2*TMath::Abs(onia_rap)*TMath::Abs(onia_rap)*TMath::Abs(onia_rap); //MeV
    sigmaRap/=1000.;//GeV
    //scale from ups1S down to Jpsi
    sigmaRap*=onia::MpsiPDG/onia::Mups1SPDG;
    // different mass region for different states
    massMin = onia::MpsiPDG-onia::nSigMass*sigmaRap;
    massMax = onia::MpsiPDG+onia::nSigMass*sigmaRap;

    // count all events
    Reco_StatEv->Fill(0.5);

    //Muon CUTS
    //apply fiducial cuts
    bool muonsInAcc = isMuonInAcceptance(FidCuts-1, pTMuPos, etaMuPos, muAccShift) && isMuonInAcceptance(FidCuts-1,pTMuNeg, etaMuNeg, muAccShift);

    if(!muonsInAcc) continue;
    // count events after fiducial cuts
    Reco_StatEv->Fill(1.5);


    //Dimuon CUTS

    if(onia_pt > 990.) continue;
    if (cutDimuon10Gev) {
      if (onia_pt < 10.) continue;
    }
    //if(jpsiVprob < onia::cut_vtxProb) continue;
    if(rejectCowboys)
      if(deltaPhi < 0.)  continue;

    //check the trigger flag: 0... no trigger, 1 ... triggered+matched, 3 ... triggered (HLT_DoubleMu0)
    //for a full list of accessible triggers, check https://espace.cern.ch/cms-quarkonia/onia-polarization/L1%20%20HLT/unprescaledTriggersVsRun.aspx
    int trigDecision = -99;

    //if(TMath::Abs(ctpv/ctpv_error)>2.5) continue;

    //TODO: add Jpsi vertex probability cut, add masscut (fit)

    // restrict to barrel for upsilon and jpsi
    //if(TMath::Abs(onia_rap) > onia::rap) continue;

    //if(onia_mass < massMin || onia_mass > massMax) continue;


    // count events after dimuon cuts
    Reco_StatEv->Fill(2.5);

    //Gamma CUTS

    if(photon->Pt()<onia::cut_gammapt) continue;
    if(TMath::Abs(photon->PseudoRapidity())>onia::cut_gammaeta) continue;
    if(conv_vertex<onia::cut_RconvMin || conv_vertex>onia::cut_RconvMax) continue;


    // count events after gamma cuts
    Reco_StatEv->Fill(3.5);

    //Chic CUTS
    if(dz>onia::cut_dz) continue;
    if(probFit1S<onia::cut_probFit) continue;

    // count events after chic cuts
    Reco_StatEv->Fill(4.5);

    numPV=numPrimaryVertices;
    mQ=chic->M()-jpsi->M()+onia::MpsiPDG;
    jpsict=ctpv*10; //convert cm in mm
    jpsictErr=ctpv_error*10; //convert cm in mm
    probKVF=probFit1S;


    double min=0.02;
    double s1=500;
    double s2=1.+s1*min*min;
    double ScaleFactorjpsictErr=(1+s1*(jpsictErr)*(jpsictErr))/s2;

    //if(jpsictErr>min)
    //jpsictErr*=ScaleFactorjpsictErr;

    //if(jentry % 2 == 0)
    //jpsictErr=h_JpsictErr->GetRandom();

    if(correctCtau){
      // double LifetimeCorrFactor = jpsi->Pt() / jpsi->M() * chic_rf->M() / chic->Pt();
      double LifetimeCorrFactor = jpsi->Pt() / jpsi->M() * onia::Mchi1PDG / chic->Pt(); // correction with PDG mass
      jpsict    = jpsict    * LifetimeCorrFactor ;
      jpsictErr = jpsictErr * LifetimeCorrFactor ;
    }

    treeOut->Fill();

    //remaining of the events will be used for the analysis
    count++;

  } // for loop over events

  std::cout << "number of reconstructed events: " << count << " of a total of " << nentries << " events" << std::endl;

} // void
