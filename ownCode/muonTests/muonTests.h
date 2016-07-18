//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul 14 17:03:56 2016 by ROOT version 5.32/00
// from TTree data/CMSSW Quarkonia J/psi Polarization+Trigger Tree
// found on file: /data/users/ferraioc/Polarization/2012ppOniaData/r2012A_DoubMu_jpsi_v8.root
//////////////////////////////////////////////////////////

// TODO: changes necessary for chic input files (different/no triggers, different branch names, etc...)

#ifndef muonTests_h
#define muonTests_h

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TLorentzVector.h"

class muonTests {
public :
  TTree *fChain;   //!pointer to the analyzed TTree or TChain
  int fCurrent; //!current Tree number in a TChain

  // Declaration of leaf types
  TLorentzVector  *onia, *muPos, *muNeg;
  int eventNb;
  int runNb;
  int lumiBlock;
  int nPriVtx;

  double Jpsict;
  double JpsictErr;
  double JpsiMassErr;
  double JpsiVprob;
  double JpsiDistM1;
  double JpsiDphiM1;
  double JpsiDrM1;
  double JpsiDistM2;
  double JpsiDphiM2;
  double JpsiDrM2;

  // trigger
  // Jpsi
  int HLT_Dimuon8_Jpsi_v3;
  int HLT_Dimuon8_Jpsi_v4;
  int HLT_Dimuon8_Jpsi_v5;
  int HLT_Dimuon8_Jpsi_v6;
  int HLT_Dimuon8_Jpsi_v7;
  int HLT_Dimuon10_Jpsi_v3;
  int HLT_Dimuon10_Jpsi_v4;
  int HLT_Dimuon10_Jpsi_v5;
  int HLT_Dimuon10_Jpsi_v6;


  muonTests(TTree *tree=0);
  virtual ~muonTests();
  virtual int Cut(Long64_t entry);
  virtual int GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef muonTests_cxx
muonTests::muonTests(TTree *tree) : fChain(0)
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/afs/hephy.at/data/ikraetschmer01/ikraetschmer/TnP2012/InputFiles/validation/onia2MuMu_tree_validation.root");
    // TFile* f = (TFile*)gROOT->GetListOfFiles()->FindObject("/scratch/knuenz/Polarization/RootInput/ChicPol/chic_rootuple_subFeb2014.root"); // have to do some restructuring for this!
    if (!f || !f->IsOpen()) {
      f = new TFile("/afs/hephy.at/data/ikraetschmer01/ikraetschmer/TnP2012/InputFiles/validation/onia2MuMu_tree_validation.root");
      // f = new TFile("/scratch/knuenz/Polarization/RootInput/ChicPol/chic_rootuple_subFeb2014.root"); // same here
    }
    f->GetObject("data",tree);

  }

  Init(tree);
}

muonTests::~muonTests()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t muonTests::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t muonTests::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void muonTests::Init(TTree *tree)
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

  onia = 0;
  muNeg = 0;
  muPos = 0;

  fChain->SetBranchAddress("JpsiP", &onia);
  fChain->SetBranchAddress("muNegP", &muNeg);
  fChain->SetBranchAddress("muPosP", &muPos);

  fChain->SetBranchAddress("eventNb", &eventNb);
  fChain->SetBranchAddress("runNb", &runNb);
  fChain->SetBranchAddress("lumiBlock", &lumiBlock);
  fChain->SetBranchAddress("nPriVtx", &nPriVtx);

  fChain->SetBranchAddress("Jpsict", &Jpsict);
  fChain->SetBranchAddress("JpsictErr", &JpsictErr);
  fChain->SetBranchAddress("JpsiMassErr", &JpsiMassErr);
  fChain->SetBranchAddress("JpsiVprob", &JpsiVprob);
  fChain->SetBranchAddress("JpsiDistM1", &JpsiDistM1);
  fChain->SetBranchAddress("JpsiDphiM1", &JpsiDphiM1);
  fChain->SetBranchAddress("JpsiDrM1", &JpsiDrM1);
  fChain->SetBranchAddress("JpsiDistM2", &JpsiDistM2);
  fChain->SetBranchAddress("JpsiDphiM2", &JpsiDphiM2);
  fChain->SetBranchAddress("JpsiDrM2", &JpsiDrM2);


  fChain->SetBranchAddress("HLT_Dimuon8_Jpsi_v3", &HLT_Dimuon8_Jpsi_v3);
  fChain->SetBranchAddress("HLT_Dimuon8_Jpsi_v4", &HLT_Dimuon8_Jpsi_v4);
  fChain->SetBranchAddress("HLT_Dimuon8_Jpsi_v5", &HLT_Dimuon8_Jpsi_v5);
  fChain->SetBranchAddress("HLT_Dimuon8_Jpsi_v6", &HLT_Dimuon8_Jpsi_v6);
  fChain->SetBranchAddress("HLT_Dimuon8_Jpsi_v7", &HLT_Dimuon8_Jpsi_v7);
  fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_v3", &HLT_Dimuon10_Jpsi_v3);
  fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_v4", &HLT_Dimuon10_Jpsi_v4);
  fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_v5", &HLT_Dimuon10_Jpsi_v5);
  fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_v6", &HLT_Dimuon10_Jpsi_v6);


  Notify();
}

Bool_t muonTests::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void muonTests::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t muonTests::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef muonTests_cxx
