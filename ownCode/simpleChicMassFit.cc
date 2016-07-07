#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooPolynomial.h"

#include <iostream>
#include <string>

#include "commonVar.h"

/** check if the event is in the kinematic regions for the analysis */
bool kinAccept(const double pt, const double rap)
{
  if (pt >= onia::pTRange[0][0] && pt < onia::pTRange[0][onia::kNbPTBins[0]] &&
      TMath::Abs(rap) >= onia::rapForPTRange[0] && TMath::Abs(rap) < onia::rapForPTRange[onia::kNbRapForPTBins]) {
    return true;
  }

  return false;
}

/** check whether the event would be accepted by the cuts applied at the createWorkspace stage. */
bool wsAccept(const TLorentzVector* jpsi, const TLorentzVector* chic, const double ct, const double ctErr)
{
  if (jpsi->M() > onia::massMin && jpsi->M() < onia::massMax &&
      jpsi->Pt() < 1000. && jpsi->Pt() > 0. &&
      TMath::Abs(jpsi->Rapidity()) < 2. &&
      chic->M() > onia::chimassMin && chic->M() < onia::chimassMax &&
      chic->Pt() < 100. && chic->Pt() > 0. &&
      TMath::Abs(chic->Rapidity()) < onia::chirap &&
      ct > onia::ctVarMin && ct < onia::ctVarMax &&
      ctErr > 0.0001 && ctErr < 1.) {

    return kinAccept(chic->Pt(), chic->Rapidity());
  }

  return false;
}

/** creates all necessary objects in the workspace for the mass only fit. */
void setupMassShapes(RooWorkspace* ws)
{
  // assuming MC == false in chiMassLifetimeFit (different ranges for the widths of the CB functions)
  ws->factory("RooCBShape::M_chic1(chicMass, CBmass1[3.5,3.45,3.54], CBsigma1[0.008,0.003,0.02], CBalpha1[0.6,0.2,1.1], CBn[2.5,1.8,3.2])");
  RooFormulaVar PES("PES", Form("(@0-%f)/%f", onia::MpsiPDG, onia::Mchi1PDG-onia::MpsiPDG), RooArgList(*ws->var("CBmass1")));
  ws->import(PES);

  RooFormulaVar CBmass0("CBmass0", Form("@0*%f+%f", onia::Mchi0PDG-onia::MpsiPDG, onia::MpsiPDG), RooArgList(*ws->function("PES")));
  RooFormulaVar CBsigma0("CBsigma0", Form("@0*%f", (onia::Mchi0PDG-onia::MpsiPDG)/(onia::Mchi1PDG-onia::MpsiPDG)), RooArgList(*ws->var("CBsigma1")));
  ws->import(CBmass0);
  ws->import(CBsigma0);

  RooFormulaVar CBmass2("CBmass2",Form("@0*%f+%f",onia::Mchi2PDG-onia::MpsiPDG, onia::MpsiPDG),RooArgList(*ws->function("PES")));
  RooFormulaVar CBsigma2("CBsigma2",Form("@0*%f",(onia::Mchi2PDG-onia::MpsiPDG)/(onia::Mchi1PDG-onia::MpsiPDG)),RooArgList(*ws->var("CBsigma1")));
  ws->import(CBmass2);
  ws->import(CBsigma2);

  RooFormulaVar CBn2("CBn2", Form("@0"), RooArgList(*ws->var("CBn")));
  ws->import(CBn2);
  ws->factory("RooCBShape::M_chic2(chicMass, CBmass2, CBsigma2, CBalpha2[0.6,0.2,1.1], CBn2)");

  RooFormulaVar CBalpha0("CBalpha0",Form("(@0+@1)/2."),RooArgList(*ws->var("CBalpha1"),*ws->var("CBalpha2")));
  ws->import(CBalpha0);

  ws->factory("RooVoigtian::M_chic0(chicMass,CBmass0,CBsigma0, CBwidth0[0.0104])");
}

/** setup the background functions (mass only) */
void setupBackground(RooWorkspace* ws)
{
  ws->factory("q01S[3.1,3.0,3.2]");
  ws->factory("alpha1[0.6,0.,5.0]");
  ws->factory("beta1[-2.5,-10.,0.]");
  RooFormulaVar a1("a1", "TMath::Abs(@0-@1)", RooArgList(*ws->var("chicMass"), *ws->var("q01S")));
  RooFormulaVar b1("b1","@0*(@1-@2)",RooArgList(*ws->var("beta1"),*ws->var("chicMass"), *ws->var("q01S")));
  RooFormulaVar signum1("signum1","(TMath::Sign(-1.,@0-@1)+1)/2.",RooArgList(*ws->var("chicMass"),*ws->var("q01S")));

  ws->factory("BK_p1[0.,-1.,1.]");
  ws->factory("BK_p2[0.,-1.,1.]");

  RooAbsPdf* M_background = new RooPolynomial("M_background", "M_background", *ws->var("chicMass"), RooArgList(*ws->var("BK_p1"), *ws->var("BK_p2")));
  ws->import(*M_background);
}

/** setup different factors that are needed. */
void setupMiscFactors(RooWorkspace* ws)
{
  ws->factory("fracSignal_chic1[0.7,0.6,0.8]");
  ws->factory("fracSignal_chic0[0.03,0.,0.01]");
  ws->factory("fracBackground[0.6,0.,1.]");

  // TODO: make this non hardcoded
  ws->factory("jpsi_fBkg[0.0182726]");
}

/**
 * TODO: for MC closure have at the moment ten events more by this than by createWorkspace!
 * -> LOOK INTO THIS!
 */
int main(int argc, char* argv[])
{
  std::cout << "Opening file: " << argv[1] << std::endl;
  TFile* inFile = TFile::Open(argv[1], "READ");
  TTree* inTree = (TTree*) inFile->Get("selectedData");

  TLorentzVector* chic_rf = NULL;
  inTree->SetBranchAddress("chic_rf", &chic_rf);
  TLorentzVector* jpsi = NULL;
  inTree->SetBranchAddress("jpsi", &jpsi);
  double ct;
  inTree->SetBranchAddress("Jpsict", &ct);
  double ctErr;
  inTree->SetBranchAddress("JpsictErr", &ctErr);

  RooRealVar mass("chicMass", "M^{#chi} [GeV]", onia::chimassMin, onia::chimassMax);
  RooArgList massList(mass);
  RooDataSet* massSet = new RooDataSet("chicMassSet", "", massList);

  for (int i = 0; i < inTree->GetEntries(); ++i) {
    inTree->GetEntry(i);

    if (!wsAccept(jpsi, chic_rf, ct, ctErr)) continue;

    mass.setVal(chic_rf->M());
    massSet->add(massList);
  }

  // massSet->Print();

  RooWorkspace* ws = new RooWorkspace("ws_masslifetime");
  ws->import(*massSet);

  setupMassShapes(ws);
  setupBackground(ws);
  setupMassShapes(ws);

  ws->factory("SUM::M_signal(fracSignal_chic1*M_chic1, fracSignal_chic0*M_chic0, M_chic2)");
  ws->factory("SUM:M_fullModel(fracBackground*M_background, jpsi_fBkg*M_background, M_signal)");

  ws->writeToFile("FitWorkspace.root");
  delete chic_rf; delete jpsi; delete massSet; delete ws;
  return 0;
}
