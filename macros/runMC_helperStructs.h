#ifndef RUNMC_HELPER_STRUCTS_H__
#define RUNMC_HELPER_STRUCTS_H__

#include "rootIncludes.inc"
#include "commonVar.h"
#include "effsAndCuts.h" // NOTE: make sure to include the right version

#include <vector>

/** helper struct to contain the output branches. */
struct OutTree {
  TTree* tree;
  TLorentzVector* lepP;
  TLorentzVector* lepN;
  TLorentzVector* jpsi;

  double jpsict;
  double jpsictErr;
  double jpsiMassErr;
  double jpsiVprob;

  OutTree(const std::string& branchname, const std::string& desc) :
    tree(new TTree(branchname.c_str(), desc.c_str())), lepP(NULL), lepN(NULL), jpsi(NULL),
    jpsict(0.), jpsictErr(0.), jpsiMassErr(0.), jpsiVprob(0.)
  {
    tree->Branch("lepP", "TLorentzVector", &lepP);
    tree->Branch("lepN", "TLorentzVector", &lepN);
    tree->Branch("jpsi", "TLorentzVector", &jpsi);

    tree->Branch("Jpsict", &jpsict, "Jpsict/D");
    tree->Branch("JpsictErr", &jpsictErr, "JpsictErr/D");
    tree->Branch("JpsiMassErr", &jpsiMassErr, "JpsiMassErr/D");
    tree->Branch("JpsiVprob", &jpsiVprob, "JpsiVprob/D");
  }
};

/** helper struct to contain the output histograms instead of having them lying around as global variables.*/
struct OutHistograms {
  TH1F* Reco_StatEv;
  TH1F* Reco_Onia_mass[onia::kNbPTMaxBins+1][onia::kNbRapForPTBins+1];
  TH2F* Reco_Onia_rap_pT;
  TH1F* Reco_Onia_pt[onia::kNbRapForPTBins+1];
  TH1F* Reco_Onia_rap[onia::kNbPTMaxBins+1];
};

/** the same as PolData and PolChiData except for MC. */
struct PolDataMC {
  TTree* fChain;
  int fCurrent;

  TLorentzVector* onia;
  TLorentzVector* muPos;
  TLorentzVector* muNeg;

  int eventNb;
  int runNb;
  int lumiBlock;
  int nPriVtx;

  double Jpsict;
  double JpsictErr;
  double JpsiMassErr;
  double JpsiVprob;

  // triggers
  std::vector<int> HLT_Dimuon8_Jpsi; /**< all different versions in one vector. */
  std::vector<int> HLT_Dimuon10_Jpsi; /**< all different versions in one vector. */

  std::vector<char> Dimuon8_vs; /**< used versions of the Trigger. */
  std::vector<char> Dimuon10_vs; /**< used versions of the Trigger. */

  PolDataMC(TTree* tree = NULL);
  virtual ~PolDataMC();
  virtual int Cut(Long64_t entry);
  virtual int GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void Init(TTree* tree);
  virtual void Loop(OutTree& outTree, OutHistograms& outHistos, int nState, int FidCuts, bool RequestTrigger,
                    bool rejectCowboys, bool removeEta0p2_0p3, bool cutDeltaREllDpt);
  virtual bool Notify();
  virtual void Show(Long64_t entry = -1);
  bool jpsiTriggerDecision(); /**< returns true if ANY of the triggers fired. */
  void fillHistograms(OutHistograms& outHistos, double rap, double pt, double mass); /** fill the global histograms. */
  void fillTree(OutTree& tree); /**< fill the outtree. */
};

PolDataMC::PolDataMC(TTree* tree)
{
  // NOTE: tmadlener: have to do this in this way because initializer lists come around in c++11 only
  const char v8[] = {'3', '4', '5', '6', '7'};
  Dimuon8_vs = std::vector<char>(v8, v8 + sizeof(v8) / sizeof(v8[0]));
  const char v10[] = {'3', '4', '5', '6'};
  Dimuon10_vs = std::vector<char>(v10, v10 + sizeof(v10) / sizeof(v10[0]));

  // tmadlener: skipped the check for a nullptr that is present in PolData and PolChiData
  Init(tree);
}

PolDataMC::~PolDataMC()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

int PolDataMC::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t PolDataMC::LoadTree(Long64_t entry)
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

void PolDataMC::Init(TTree* tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  if (!tree) return;
  fChain = tree;
  fCurrent = -1;

  onia = muPos = muNeg = NULL; // NOTE: potentially dangerous

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

  for (size_t i = 0; i < Dimuon8_vs.size(); ++i) {
    const std::string branchname = std::string("HLT_Dimuon8_Jpsi_v") + Dimuon8_vs[i];
    HLT_Dimuon8_Jpsi.push_back(0);
    fChain->SetBranchAddress(branchname.c_str(), &HLT_Dimuon8_Jpsi[i]);
  }

  for (size_t i = 0; i < Dimuon8_vs.size(); ++i) {
    const std::string branchname = std::string("HLT_Dimuon10_Jpsi_v") + Dimuon10_vs[i];
    HLT_Dimuon10_Jpsi.push_back(0);
    fChain->SetBranchAddress(branchname.c_str(), &HLT_Dimuon10_Jpsi[i]);
  }

}

Bool_t PolDataMC::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PolDataMC::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

int PolDataMC::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

void PolDataMC::Loop(OutTree& outTree, OutHistograms& outHistos, int nState, int FidCuts, bool RequestTrigger,
                     bool rejectCowboys, bool removeEta0p2_0p3, bool cutDeltaREllDpt)
{
  if (!fChain) return;

  Long64_t nentries = fChain->GetEntries();
  // Long64_t cutAtRecEvent = nentries;
  Long64_t count = 0;
  Long64_t nb = 0;

  std::cout << "number of entries = " << nentries << std::endl;

  // nentries = 250000; // for development
  for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
    if (jentry % 100000 == 0) std::cout << "event " << jentry << " of " << nentries << std::endl;

    if (LoadTree(jentry) < 0) break;
    nb = fChain->GetEntry(jentry);

    // process only reconstructed events
    if (onia->Pt() > 990.) continue;
    if (JpsiVprob < 0.01) continue;

    // count all events
    outHistos.Reco_StatEv->Fill(0.5);

    // different trigger for different particles
    // Jpsi trigger paths
    if (nState == 4) {
      if (RequestTrigger && !jpsiTriggerDecision()) continue; // use lazy evaluation to save some time
    }

    // count events after trigger
    outHistos.Reco_StatEv->Fill(1.5);

    double onia_rap = onia->Rapidity();
    if (TMath::Abs(onia_rap) > onia::rap) continue;
    outHistos.Reco_StatEv->Fill(2.5); // count events after rapidity cut

    double deltaPhi = muNeg->Phi() - muPos->Phi();
    if (deltaPhi > TMath::Pi()) deltaPhi -= 2. * TMath::Pi();
    else if (deltaPhi < -TMath::Pi()) deltaPhi += 2.*TMath::Pi();

    if (rejectCowboys) {
      if (deltaPhi < 0.) continue;
      outHistos.Reco_StatEv->Fill(3.5); // count events after cowboy rejection
    }

    double massMin = onia::massMin;
    double massMax = onia::massMax;

    double onia_mass = onia->M();
    if (onia_mass < massMin || onia_mass > massMax) continue;
    outHistos.Reco_StatEv->Fill(4.5); // events after mass cut

    double etaMuPos = muPos->Eta();
    double etaMuNeg = muNeg->Eta();
    double pTMuPos = muPos->Pt();
    double pTMuNeg = muNeg->Pt();

    // NOTE: make sure to provide the appropriate version of this function!
    if (!isMuonInAcceptance(FidCuts-1, pTMuPos, etaMuPos) || !isMuonInAcceptance(FidCuts-1, pTMuNeg, etaMuNeg))
      continue;
    outHistos.Reco_StatEv->Fill(5.5); // count events after fiducial cuts


    if (removeEta0p2_0p3) {
      if ( (TMath::Abs(etaMuPos) > 0.2 && TMath::Abs(etaMuPos) < 0.3) ||
           (TMath::Abs(etaMuNeg) > 0.2 && TMath::Abs(etaMuNeg) < 0.3) ) {
        continue;
      }
      outHistos.Reco_StatEv->Fill(6.5);
    }
    double onia_pt = onia->Pt();

    if (cutDeltaREllDpt) {
      double deltaEta = etaMuNeg - etaMuPos;
      double deltaPt = pTMuNeg - pTMuPos;
      double deltaREll = TMath::Sqrt(deltaEta*deltaEta + TMath::Power(1.2*deltaPhi, 2));
      double deltaREllDpt = TMath::Sqrt(deltaREll*deltaREll + TMath::Power(0.00157*deltaPt, 2));

      if (onia_pt > 35. && onia_pt < 40. && deltaREllDpt < 0.18) continue;
      if (onia_pt > 40. && onia_pt < 50. && deltaREllDpt < 0.16) continue;
      if (onia_pt > 50. && onia_pt < 70. && deltaREllDpt < 0.14) continue;

      outHistos.Reco_StatEv->Fill(7.5);
    }

    fillHistograms(outHistos, onia_rap, onia_pt, onia_mass);
    fillTree(outTree);

    count++;
  } // end event loop

  std::cout << "number of reconstructed events: " << count << " of a total of " << nentries << " events" << std::endl;
}

bool PolDataMC::jpsiTriggerDecision()
{
  //check the trigger flag: 0... no trigger, 1 ... triggered+matched, 3 ... triggered (HLT_DoubleMu0)
  //for a full list of accessible triggers, check https://espace.cern.ch/cms-quarkonia/onia-polarization/L1%20%20HLT/unprescaledTriggersVsRun.aspx

  typedef std::vector<int>::const_iterator vecIt;
  for (vecIt it = HLT_Dimuon8_Jpsi.begin(); it != HLT_Dimuon8_Jpsi.end(); ++it) {
    if (*it == 1) return true;
  }
  for (vecIt it = HLT_Dimuon10_Jpsi.begin(); it != HLT_Dimuon10_Jpsi.end(); ++it) {
    if (*it == 1) return true;
  }
  return false;
}

void PolDataMC::fillHistograms(OutHistograms& outHistos, double rap, double pt, double mass)
{
  // filling mass, pt and y histograms
  // indices [0] contain all events while [1], etc. show events according to the pt and y bin

  outHistos.Reco_Onia_rap_pT->Fill(rap, pt);
  outHistos.Reco_Onia_mass[0][0]->Fill(mass);
  outHistos.Reco_Onia_rap[0]->Fill(rap);
  outHistos.Reco_Onia_pt[0]->Fill(pt);

  for (int iRap = 0; iRap < onia::kNbRapForPTBins; ++iRap) {
    // events integrated in pt in different rapidity bins
    if (TMath::Abs(rap) > onia::rapForPTRange[iRap] && TMath::Abs(rap) < onia::rapForPTRange[iRap + 1]) {
      outHistos.Reco_Onia_pt[iRap+1]->Fill(pt);
      outHistos.Reco_Onia_mass[0][iRap+1]->Fill(mass);

      for (int iPt = 0; iPt < onia::kNbPTMaxBins; ++iPt) {
        // events for different pt and y bins
        if (pt > onia::pTRange[0][iPt] && pt < onia::pTRange[0][iPt+1]) {
          outHistos.Reco_Onia_mass[iPt+1][iRap+1]->Fill(mass);
        }
      }
    }
  }

  for (int iPt = 0; iPt < onia::kNbPTMaxBins; ++iPt) {
    // events integrated in rapidity for different pT bins
    if (pt > onia::pTRange[0][iPt] && pt < onia::pTRange[0][iPt+1]) {
      outHistos.Reco_Onia_rap[iPt+1]->Fill(rap);
      outHistos.Reco_Onia_mass[iPt+1][0]->Fill(mass);
    }
  }
}

void PolDataMC::fillTree(OutTree& outTree)
{
  outTree.lepP = muPos;
  outTree.lepN = muNeg;
  outTree.jpsi = onia;
  outTree.jpsict = Jpsict;
  outTree.jpsictErr = JpsictErr;
  outTree.jpsiMassErr = JpsiMassErr;
  outTree.jpsiVprob = JpsiVprob;

  outTree.tree->Fill();
}
#endif
