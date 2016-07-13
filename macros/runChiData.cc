#ifndef __CINT__
#endif

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
using namespace std;

#include "TChain.h"
#include "rootIncludes.inc"
#include "PolChiData.C"

#include "clarg_parsing.h"

void BookHistosReco(int nState);
void WriteHistosReco(const std::string &fNameOut);

/**
 ** NOTE:
 ** - This code does not fail gracefully, it simply skips input it doesn't understand
 ** - boolalpha is used: you must write true, false for booleans, not 0 or 1
 ** - inputTree1, inputTree2, inputTree3, inputTree4, inputTrees is now a vector instead
 ** - all input is on the form key=value (no spaces, no other order, only one argument allowed)
 ** - inputTree= will give a segfault (empty value for string)
 ** - inputTree1=xyz will still work, but the internal ordering in the vector doesnt see the 1
 **/

//===================================================
int main(int argc, char* argv[]){

  // Set defaults
  bool
    rejectCowboys = true,
    MC=false,
    RequestTrigger=false,
    removeEta0p2_0p3=false,
    cutDeltaREllDpt=false;
  bool rejectSeagulls = false;
  bool correctCtau  = false;
  bool useRefittedChic = true;
  bool cutDimuon10Gev = false;
  std::vector<std::string> inputTrees;
  int
    FidCuts = 999,
    nState = 999;

  double muAccShift = 0.; // overall shift that is applied to the single muon acceptance pT cuts
  // Loop over argument list
  for (int i=1; i < argc; i++)
    {
      std::string arg = argv[i];
      fromSplit("rejectCowboys", arg, rejectCowboys);
      fromSplit("MC", arg, MC);
      fromSplit("RequestTrigger", arg, RequestTrigger);
      fromSplit("removeEta0p2_0p3", arg, removeEta0p2_0p3);
      fromSplit("cutDeltaREllDpt", arg, cutDeltaREllDpt);
      fromSplit("FidCuts", arg, FidCuts);
      fromSplit("nState", arg, nState);
      fromSplit("useRefittedChic", arg, useRefittedChic);
      std::string str;
      fromSplit("inputTree", arg, str);
      fromSplit("correctCtau", arg, correctCtau);
      if (!str.empty()) inputTrees.push_back(str);
      fromSplit("cutDimuon10Gev", arg, cutDimuon10Gev);
      fromSplit("muAccShift", arg, muAccShift);
      fromSplit("rejectSeagullss", arg, rejectSeagulls);
    }

  std::cout << "-----------------------------------\n" <<
    "nState = 1,2,3...Upsilon\n" <<
    "nState = 4...JPsi\n" <<
    "nState = 5...PsiPrime\n" <<
    "nState = 6...Chic1\n" <<
    "nState = 7...Chic2\n" <<
    "-----------------------------------" << std::endl;

  // get input files
  TChain *chain = new TChain("rootuple/chicTree");
  for (int i = 0; i < inputTrees.size(); i++){
    chain->Add(inputTrees[i].c_str());
  }
  TTree *tree = chain;


  // define output file
  const std::string fNameOut = "tmpFiles/selEvents_data.root";
  TFile *fOut = fOut = new TFile(fNameOut.c_str(), "RECREATE");

  PolChiData treeReco(tree);
  BookHistosReco(nState);
  printf("after booking of histo\n");
  //tree->GetEvent(1);
  treeReco.Loop(nState, rejectCowboys, FidCuts, MC, RequestTrigger, removeEta0p2_0p3, cutDeltaREllDpt, correctCtau, useRefittedChic,
                cutDimuon10Gev, muAccShift, rejectSeagulls);
  printf("writing out the histograms\n");
  WriteHistosReco(fNameOut.c_str());

  fOut->Close();

  return 0;
}
//==========================================
void BookHistosReco(int nState){

  //mass
  int nBinsMass = 320;
  double  massMin = onia::massMin;
  double massMax = onia::massMax;

  //pt
  int nBinsPt = 1000;
  double pTMin = 0., pTMax = 100.;
  //rap
  int nBinsRap = 100;
  double rapMin = -2.5, rapMax = 2.5;

  //statistics
  Reco_StatEv = new TH1F("Reco_StatEv", "", 12, 0., 12.);

  // pt vs y
  std::string name1 = "Reco_Onia_rap_pt";
  std::string title1 = ";y(#mu#mu);p_{T}^{#mu#mu} [GeV]";
  //Reco_Onia_rap_pT = new TH2F(name1.c_str(), title1.c_str(), nBinsRap, rapMin, rapMax, nBinsPt, pTMin, pTMax);
  //Reco_Onia_rap_pT->Sumw2();

  // book histograms for different pt and y bins
  for(int iRap = 0; iRap < onia::kNbRapForPTBins+1; iRap++){

    //pT
    std::stringstream namePt;
    namePt << "Reco_Onia_pt_rap" << iRap;
    std::string titlePt = ";p_{T}^{#mu#mu} [GeV]";
    //Reco_Onia_pt[iRap]  = new TH1F(namePt.str().c_str(), titlePt.c_str(), nBinsPt, pTMin, pTMax);
    //Reco_Onia_pt[iRap]->Sumw2();

    for(int iPT = 0; iPT < onia::kNbPTMaxBins+1; iPT++){

      //Mass:
      std::stringstream name;
      name << "Reco_Onia_mass_rap" << iRap << "_pT" << iPT;
      std::string title = ";M [GeV]";
      //Reco_Onia_mass[iPT][iRap] = new TH1F(name.str().c_str(), title.c_str(), nBinsMass, massMin, massMax);
      //Reco_Onia_mass[iPT][iRap]->Sumw2();

    } // iRap
  } // iPT


  for(int iPT = 0; iPT < onia::kNbPTMaxBins+1; iPT++){
    //rap
    std::stringstream nameRap;
    nameRap << "Reco_Onia_rap_pT" << iPT;
    std::string titleRap = ";y(#mu#mu)";
    //Reco_Onia_rap[iPT]  = new TH1F(nameRap.str().c_str(), titleRap.c_str(), nBinsRap, rapMin,rapMax);
    //Reco_Onia_rap[iPT]->Sumw2();
  } // iPT

  //prepare the branches for the output tree
  treeOut = new TTree ("selectedData", "selected events");

  lepP = new TLorentzVector();
  lepN = new TLorentzVector();
  jpsi = new TLorentzVector();
  photon = new TLorentzVector();
  chic = new TLorentzVector();

  lepP_rf = new TLorentzVector();
  lepN_rf = new TLorentzVector();
  jpsi_rf = new TLorentzVector();
  photon_rf = new TLorentzVector();
  chic_rf = new TLorentzVector();

  treeOut->Branch("lepP", "TLorentzVector", &lepP);
  treeOut->Branch("lepN", "TLorentzVector", &lepN);
  treeOut->Branch("jpsi", "TLorentzVector", &jpsi);
  treeOut->Branch("photon", "TLorentzVector", &photon);
  treeOut->Branch("chic", "TLorentzVector", &chic);

  treeOut->Branch("lepP_rf", "TLorentzVector", &lepP_rf);
  treeOut->Branch("lepN_rf", "TLorentzVector", &lepN_rf);
  treeOut->Branch("jpsi_rf", "TLorentzVector", &jpsi_rf);
  treeOut->Branch("photon_rf", "TLorentzVector", &photon_rf);
  treeOut->Branch("chic_rf", "TLorentzVector", &chic_rf);

}

//==========================================
void WriteHistosReco(const std::string &fNameOut){

  treeOut->Write();

  Reco_StatEv->Write();
  //Reco_Onia_rap_pT->Write();
  //
  //for(int iRap = 0; iRap < onia::kNbRapForPTBins+1; iRap++){
  //	Reco_Onia_pt[iRap]->Write();
  //
  //  for(int iPT = 0; iPT < onia::kNbPTMaxBins+1; iPT++){
  //    Reco_Onia_mass[iPT][iRap]->Write();
  //  }
  //}
  //
  //for(int iPTBin = 0; iPTBin < onia::kNbPTMaxBins+1; iPTBin++)
  //  Reco_Onia_rap[iPTBin]->Write();
}
