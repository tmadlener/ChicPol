#ifndef PolChiData_h
#define PolChiData_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TLorentzVector.h"

class PolChiData {
public :
  TTree *fChain;   //!pointer to the analyzed TTree or TChain
  int fCurrent; //!current Tree number in a TChain


  TLorentzVector  *chi_p4, *dimuon_p4, *muonP_p4, *muonN_p4, *photon_p4;
  TLorentzVector  *rf1S_chi_p4, *rf1S_dimuon_p4, *rf1S_muonP_p4, *rf1S_muonN_p4, *rf1S_photon_p4;

  Double_t ele_lowerPt_pt;
  Double_t ele_higherPt_pt;
  Double_t ctpv;
  Double_t ctpv_error;
  Double_t pi0_abs_mass;
  Double_t psi1S_nsigma;
  Double_t conv_vertex;
  Double_t dz;
  Double_t numPrimaryVertices;
  Int_t trigger;
  Double_t probFit1S;

  PolChiData(TTree *tree=0);
  virtual ~PolChiData();
  virtual int Cut(Long64_t entry);
  virtual int GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void Init(TTree *tree);
  virtual void Loop(int selDimuType, bool rejectCowboys, int FidCuts, bool MC, bool RequestTrigger, bool removeEta0p2_0p3,
                    bool cutDeltaREllDpt, bool correctCtau, bool useRefittedChic, bool cutDimuon10Gev, double muAccShift = 0.,
                    bool rejectSeagulls = false);
  virtual bool Notify();
  virtual void Show(Long64_t entry = -1);
};

#endif

#ifdef PolChiData_cxx
PolChiData::PolChiData(TTree *tree)
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("TTree_Onia2MuMu_V8_PromptReco_v4.root");

    if (!f) {
      f = new TFile("TTree_Onia2MuMu_V8_PromptReco_v4.root");
    }

    tree = (TTree*)gDirectory->Get("data");
  }

  Init(tree);
}

PolChiData::~PolChiData()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

int PolChiData::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t PolChiData::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void PolChiData::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;

  chi_p4 = 0;
  dimuon_p4 = 0;
  muonP_p4 = 0;
  muonN_p4 = 0;
  photon_p4 = 0;
  rf1S_chi_p4 = 0;
  rf1S_dimuon_p4 = 0;
  rf1S_muonP_p4 = 0;
  rf1S_muonN_p4 = 0;
  rf1S_photon_p4 = 0;

  fChain->SetBranchAddress("ele_lowerPt_pt",&ele_lowerPt_pt);
  fChain->SetBranchAddress("ele_higherPt_pt",&ele_higherPt_pt);
  fChain->SetBranchAddress("ctpv",&ctpv);
  fChain->SetBranchAddress("ctpv_error",&ctpv_error);
  fChain->SetBranchAddress("pi0_abs_mass",&pi0_abs_mass);
  fChain->SetBranchAddress("psi1S_nsigma",&psi1S_nsigma);
  fChain->SetBranchAddress("conv_vertex",&conv_vertex);
  fChain->SetBranchAddress("dz",&dz);
  fChain->SetBranchAddress("numPrimaryVertices",&numPrimaryVertices);
  fChain->SetBranchAddress("trigger",&trigger);
  fChain->SetBranchAddress("probFit1S",&probFit1S);

  fChain->SetBranchAddress("chi_p4", &chi_p4);
  fChain->SetBranchAddress("dimuon_p4", &dimuon_p4);
  fChain->SetBranchAddress("muonP_p4", &muonP_p4);
  fChain->SetBranchAddress("muonN_p4", &muonN_p4);
  fChain->SetBranchAddress("photon_p4", &photon_p4);

  fChain->SetBranchAddress("rf1S_chi_p4", &rf1S_chi_p4);
  fChain->SetBranchAddress("rf1S_dimuon_p4", &rf1S_dimuon_p4);
  fChain->SetBranchAddress("rf1S_muonP_p4", &rf1S_muonP_p4);
  fChain->SetBranchAddress("rf1S_muonN_p4", &rf1S_muonN_p4);
  fChain->SetBranchAddress("rf1S_photon_p4", &rf1S_photon_p4);

  Notify();
}

Bool_t PolChiData::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void PolChiData::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
int PolChiData::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef PolChiData_cxx
