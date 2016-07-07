#include "../macros/bkgHistos_helper.h"
// #include "../interface/commonVar.h"g

#include "TFile.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"

#include <iostream>

int main(int argc, char* argv[])
{
  TFile* file = TFile::Open(argv[1]);
  RooWorkspace* ws = static_cast<RooWorkspace*>(file->Get("ws_masslifetime"));

  std::cout << argv[1] << std::endl;

  // ws->Print();

  // double CBmass1 = getVarVal(ws, "CBmass1");
  double CBmass0 = getVarVal(ws, "CBmass0");
  double CBmass2 = getVarVal(ws, "CBmass2");

  // std::cout << "=============================" << std::endl;
  // std::cout << "CBmass0 = " << CBmass0 << std::endl;
  // std::cout << "CBmass1 = " << CBmass1 << std::endl;
  // std::cout << "CBmass2 = " << CBmass2 << std::endl;
  // std::cout << "------------------------------" << std::endl;
  // std::cout << "CBmass0 / CBmass1 = " << CBmass0 / CBmass1 << std::endl;
  // std::cout << "CBmass0 / CBmass2 = " << CBmass0 / CBmass2 << std::endl;
  // std::cout << "CBmass1 / CBmass2 = " << CBmass1 / CBmass2 << std::endl;
  // std::cout << "------------------------------" << std::endl;
  // std::cout << "CBmass0 - CBmass1 = " << CBmass0 - CBmass1 << std::endl;
  // std::cout << "CBmass0 - CBmass2 = " << CBmass0 - CBmass2 << std::endl;
  // std::cout << "CBmass1 - CBmass2 = " << CBmass1 - CBmass2 << std::endl;

  double CBsigma0 = getVarVal(ws, "CBsigma0");
  double CBsigma1 = getVarVal(ws, "CBsigma1");
  double CBsigma2 = getVarVal(ws, "CBsigma2");

  // std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
  // std::cout << "CBsigma0 = " << CBsigma0 << std::endl;
  // std::cout << "CBsigma1 = " << CBsigma1 << std::endl;
  // std::cout << "CBsigma2 = " << CBsigma2 << std::endl;
  // std::cout << "------------------------------" << std::endl;
  // std::cout << "CBsigma0 / CBsigma1 = " << CBsigma0 / CBsigma1 << std::endl;
  // std::cout << "CBsigma0 / CBsigma2 = " << CBsigma0 / CBsigma2 << std::endl;
  // std::cout << "CBsigma2 / CBsigma1 = " << CBsigma2 / CBsigma1 << std::endl;

  // double CBn = getVarVal(ws, "CBn");
  // double CBn2 = getVarVal(ws, "CBn2");

  // std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl;
  // std::cout << "CBn = " << CBn << std::endl;
  // std::cout << "CBn2 = " << CBn2 << std::endl;

  double Qchic0 = onia::Mchi0PDG - onia::MpsiPDG;
  double Qchic1 = onia::Mchi1PDG - onia::MpsiPDG;
  double Qchic2 = onia::Mchi2PDG - onia::MpsiPDG;
  double pes = getVarVal(ws, "PES");

  // std::cout << "------------------------------" << std::endl;
  // std::cout << "(M0 - Mpsi)/Mpsi = " << (CBmass0 - onia::MpsiPDG) / pes << ", Qchic0 = " << Qchic0 << std::endl;
  // std::cout << "(M2 - Mpsi)/Mpsi = " << (CBmass2 - onia::MpsiPDG) / pes << ", Qchic2 = " << Qchic2 << std::endl;
  // std::cout << "##############################" << std::endl;

  double mQchic0 = (CBmass0 - onia::MpsiPDG) / pes;
  double mQchic2 = (CBmass2 - onia::MpsiPDG) / pes;

  double Q2to1 = Qchic2 / Qchic1;
  double Q0to1 = Qchic0 / Qchic1;

  double sQ2to1 = CBsigma2 / CBsigma1;
  double sQ0to1 = CBsigma0 / CBsigma1;

  std::cout << "------------------------------" << std::endl;
  std::cout << "Convergence: " << getVarVal(ws, "var_covQualHesse") << " (Hesse), " << getVarVal(ws, "var_covQualMigrad") << " (Migrad)" << std::endl;
  std::cout << "Qchic0 (fitted) = " << mQchic0 << ", (constrained) = " << Qchic0 << ", rel. err. [%] = " << (mQchic0 - Qchic0) / Qchic0 * 100 << std::endl;
  std::cout << "Qchic2 (fitted) = " << mQchic2 << ", (constrained) = " << Qchic2 << ", rel. err. [%] = " << (mQchic2 - Qchic2) / Qchic2 * 100 << std::endl;

  std::cout << "Q0 / Q1 (fitted) = " << sQ0to1 << ", (constrained) = " << Q0to1 << ", rel. err. [%] = " << (sQ0to1 - Q0to1) / Q0to1 * 100 << std::endl;
  std::cout << "Q2 / Q1 (fitted) = " << sQ2to1 << ", (constrained) = " << Q2to1 << ", rel. err. [%] = " << (sQ2to1 - Q2to1) / Q2to1 * 100 << std::endl;
  return 0;
}
