#ifndef __CINT__

#include "rootIncludes.inc"
#include "clarg_parsing.h"
#include "commonVar.h"

#include "runMC_helperStructs.h"

#include "TChain.h"

//stl
#include <iostream>
#include <string>
#include <vector>
#include <sstream>

// some forward declarations, definitions follow at the bottom of this file
void BookHistosReco(int nState, OutHistograms& outHistos);
void WriteHistosReco(const std::string& fNameOut, OutTree& outTree, OutHistograms& outHistos);

/**
 * This is the main file for creating the selEvents_data.root file as e.g. runData or runChiData does for data.
 * It is built similarly to those and as them does not fail gracefully but assumes that the appropriate inputs
 * are provided.
 */
int main(int argc, char* argv[])
{
  // default values
  int nState = 999;
  int FidCuts = 999;
  std::vector<std::string> inputTrees;
  bool requestTrigger = true;
  bool rejectCowboys = false;
  bool removeEta0p2_0p3 = false;
  bool cutDeltaREllDpt = false;

  for (int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    fromSplit("nState", arg, nState);
    fromSplit("FidCuts", arg, FidCuts);
    fromSplit("RequestTrigger", arg, requestTrigger);
    fromSplit("rejectCowboys", arg, rejectCowboys);
    fromSplit("removeEta0p2_0p3", arg, removeEta0p2_0p3);
    fromSplit("cutDeltaREllDpt", arg, cutDeltaREllDpt);

    std::string str;
    fromSplit("inputTree", arg, str);
    if (!str.empty()) inputTrees.push_back(str);
  }

  TChain* chain = new TChain("data"); // "rootuple/chicTree" for chic
  for (size_t i = 0; i < inputTrees.size(); ++i) {
    chain->Add(inputTrees[i].c_str());
  }
  TTree* tree = chain;

  // define output file
  const std::string fNameOut = "tmpFiles/selEvents_data.root";
  TFile* fOut = new TFile(fNameOut.c_str(), "RECREATE");
  OutTree outTree("selectedData", "selected events");
  OutHistograms outHistos; // WARNING: raw pointers pointing to nothing!

  BookHistosReco(nState, outHistos);
  PolDataMC polData(tree);
  polData.Loop(outTree, outHistos, nState, FidCuts, requestTrigger, rejectCowboys, removeEta0p2_0p3, cutDeltaREllDpt);
  WriteHistosReco(fNameOut, outTree, outHistos);

  fOut->Close();
  return 0;
}

void BookHistosReco(int nState, OutHistograms& outHistos)
{
  // mass
  int nBinsMass = 320;
  double massMin = onia::massMin;
  double massMax = onia::massMax;

  // pt
  int nBinsPt = 1000;
  double pTMin = 0., pTMax = 100.;

  // rap
  int nBinsRap = 100;
  double rapMin = -2.5, rapMax = 2.5;

  // statistics
  outHistos.Reco_StatEv = new TH1F("Reco_StateEv", "", 12, 0., 12.);

  // pt vs y
  const std::string nameRapPt = "Reco_Onia_rap_pt";
  const std::string titleRapPt = ";y(#mu#mu);p_{T}^{#mu#mu} [GeV]";
  outHistos.Reco_Onia_rap_pT = new TH2F(nameRapPt.c_str(), titleRapPt.c_str(), nBinsRap, rapMin, rapMax, nBinsPt, pTMin, pTMax);
  outHistos.Reco_Onia_rap_pT->Sumw2();

  // book histograms for different pt and y bins
  for (int iRap = 0; iRap < onia::kNbRapForPTBins+1; ++iRap) {
    // pT
    std::stringstream namePt;
    namePt << "Reco_Onia_pt_rap" << iRap;
    const std::string titlePt = ";p_{T}^{#mu#mu} [GeV]";
    outHistos.Reco_Onia_pt[iRap] = new TH1F(namePt.str().c_str(), titlePt.c_str(), nBinsPt, pTMin, pTMax);
    outHistos.Reco_Onia_pt[iRap]->Sumw2();

    for (int iPt = 0; iPt < onia::kNbPTMaxBins+1; ++iPt) {
      std::stringstream nameRap;
      nameRap << "Reco_Onia_mass_rap" << iRap << "_pT" << iPt;
      const std::string titleRap = ";M [GeV]";
      outHistos.Reco_Onia_mass[iPt][iRap] = new TH1F(nameRap.str().c_str(), titleRap.c_str(), nBinsMass, massMin, massMax);
      outHistos.Reco_Onia_mass[iPt][iRap]->Sumw2();
    }
  }

  for (int iPt = 0; iPt < onia::kNbPTMaxBins+1; ++iPt) {
    std::stringstream nameRap;
    nameRap << "Reco_Onia_rap_pT" << iPt;
    const std::string titleRap = ";y(#mu#mu)";
    outHistos.Reco_Onia_rap[iPt] = new TH1F(nameRap.str().c_str(), titleRap.c_str(), nBinsRap, rapMin, rapMax);
    outHistos.Reco_Onia_rap[iPt]->Sumw2();
  }
  // prepare the branches for the output tree -> already done in constructor of OutTree
}

void WriteHistosReco(const std::string& fNameOut, OutTree& outTree, OutHistograms& outHistos)
{
  outTree.tree->Write();

  outHistos.Reco_StatEv->Write();
  outHistos.Reco_Onia_rap_pT->Write();

  for (int iRap = 0; iRap < onia::kNbRapForPTBins+1; ++iRap) {
    outHistos.Reco_Onia_pt[iRap]->Write();
    for (int iPt = 0; iPt < onia::kNbPTMaxBins+1; ++iPt) {
      outHistos.Reco_Onia_mass[iPt][iRap]->Write();
    }
  }

  for (int iPt = 0; iPt < onia::kNbPTMaxBins+1; ++iPt) {
    outHistos.Reco_Onia_rap[iPt]->Write();
  }
}

#endif
