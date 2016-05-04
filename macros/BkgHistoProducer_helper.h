#ifndef BKGHISTOPRODUCER_HELPER_H__
#define BKGHISTOPRODUCER_HELPER_H__

#include "rootIncludes.inc"

#include <vector>
#include <iostream> // this can probably removed after develpment

/** struct containing the variables that are stored in the data root file. */
struct BkgHistoRootVars {
  /** default constructor. Initialize all pointers to nullptr (as soon as c++11 is used), and double values to -1.*/
  BkgHistoRootVars() : lepP(NULL), lepN(NULL), jpsi(NULL), chic(NULL), chic_rf(NULL), jpsict(-1.), mQ(-1) {;}
  /** destructor. Delete all 4 vectors. */
  ~BkgHistoRootVars();
  TLorentzVector* lepP; /**< 4 momentum of positive muon. */
  TLorentzVector* lepN; /**< 4 momentum of negative muon. */
  TLorentzVector* jpsi; /**< 4 momentum of Jpsi. */
  TLorentzVector* chic; /**< 4 momentum of chic. */
  TLorentzVector* chic_rf; /**< 4 momentum of refitted chic. */
  double jpsict; /**< jpsi ct variable. */
  double mQ; /**< mQ variable. */
};

struct BkgHistoCosThetaHists {
  /** constructor. Creates the necessary number of histos (one for each frame) and initializes them to nullptr.*/
  BkgHistoCosThetaHists(const int nHists);
  /** destructor. Deletes all present histograms.*/
  ~BkgHistoCosThetaHists();
  /** store hists to root file. (all that are not NULL!) */
  void storeToFile(TFile* file);

  std::vector<TH2D*> hBG_L; /**< LSB histogram. */
  std::vector<TH2D*> hBG_R; /**< RSB histogram. */
  std::vector<TH2D*> hBG; /**< BG histogram. */
  std::vector<TH2D*> hNPBG; /**< NonPrompt BG histogram. */
  std::vector<TH2D*> hNPS; /**< NonPrompt signal histogram. */
  std::vector<TH2D*> hBGinNP_L; /**< BG in NP histogram for LSB. */
  std::vector<TH2D*> hBGinNP_R; /**< BG in NP histogram for RSB. */
  std::vector<TH2D*> hBGinNP; /**< BG in NP histogram. */
  std::vector<TH2D*> hTBG; /**< (?) histogram. */
  std::vector<TH2D*> hSR_L; /**< signal region LSB (?). */
  std::vector<TH2D*> hSR_R; /**< signal region RSB (?).*/
  std::vector<TH2D*> hSR; /**< signal region (?). */
};

// ================================================================================
//                               IMPLEMENTATION
// ================================================================================
BkgHistoRootVars::~BkgHistoRootVars()
{
  std::cout << "DESTRUCTOR OF BkgHistoRootVars" << std::endl;
  delete lepP;
  delete lepN;
  delete jpsi;
  delete chic;
  delete chic_rf;
}

BkgHistoCosThetaHists::BkgHistoCosThetaHists(const int nHists)
{
  std::cout << "CONSTRUCTOR OF BkgHistoCosThetaHists" << std::endl;
  // TODO: replace NULL with nullptr after migragtion to c++11
  for (int iH = 0; iH < nHists; ++iH) {
    hBG_L.push_back(NULL);
    hBG_R.push_back(NULL);
    hBG.push_back(NULL);
    hNPBG.push_back(NULL);
    hNPS.push_back(NULL);
    hBGinNP_L.push_back(NULL);
    hBGinNP_R.push_back(NULL);
    hBGinNP.push_back(NULL);
    hTBG.push_back(NULL);
    hSR_L.push_back(NULL);
    hSR_R.push_back(NULL);
    hSR.push_back(NULL);
  }
}

BkgHistoCosThetaHists::~BkgHistoCosThetaHists()
{
  std::cout << "DESTRUCTOR OF BkgHistoCosThetaHists" << std::endl;
  for (size_t iH = 0; iH < hBG_L.size(); ++iH) { // every vector should contain the same number of TH2Ds
    delete hBG_L[iH];
    delete hBG_R[iH];
    delete hBG[iH];
    delete hNPBG[iH];
    delete hNPS[iH];
    delete hBGinNP_L[iH];
    delete hBGinNP_R[iH];
    delete hBGinNP[iH];
    delete hTBG[iH];
    delete hSR_L[iH];
    delete hSR_R[iH];
    delete hSR[iH];
  }
}

void BkgHistoCosThetaHists::storeToFile(TFile* file)
{
  std::cout << "STORE TO FILE OF BkgHistoCosThetaHists" << std::endl;
  file->cd();
  for (size_t iH = 0; iH < hBG_L.size(); ++iH) { // every vector should contain the same number of TH2Ds
    if (hBG_L[iH]) hBG_L[iH]->Write();
    if (hBG_R[iH]) hBG_R[iH]->Write();
    if (hBG[iH]) hBG[iH]->Write();
    if (hNPBG[iH]) hNPBG[iH]->Write();
    if (hNPS[iH]) hNPS[iH]->Write();
    if (hBGinNP_L[iH]) hBGinNP_L[iH]->Write();
    if (hBGinNP_R[iH]) hBGinNP_R[iH]->Write();
    if (hBGinNP[iH]) hBGinNP[iH]->Write();
    if (hTBG[iH]) hTBG[iH]->Write();
    if (hSR_L[iH]) hSR_L[iH]->Write();
    if (hSR_R[iH]) hSR_R[iH]->Write();
    if (hSR[iH]) hSR[iH]->Write();
  }
}

#endif
