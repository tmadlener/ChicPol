#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "RooWorkspace.h"

#include "../macros/BkgHistoProducer_helper.h"
#include "commonVar.h"
// #include "string_stuff.h"
#include "vector_stuff.h"

#include <string>
#include <vector>
#include <iostream>
#include <sstream>

/** helper struct. */
struct EventContent {
  EventContent(const std::string& dataFile); /** constructor. */
  // ~EventContent(); /**< destructor. not defined at the moment. */
  void getEvent(int i) { m_inTree->GetEntry(i); }
  void getWorkspace(const std::string& name);
  void getRegionValues();
  bool inPtRapRange(int iPt, int iRap);
  bool inRangeM(const BkgHistoRange& range) { return range.accept(m_chic->M()); }
  bool inRangeL(const BkgHistoRange& range) { return range.accept(m_jpsict); }

  // members data
  TLorentzVector* m_jpsi;
  TLorentzVector* m_chic;
  double m_jpsict;

  TFile* m_fitFile;
  RooWorkspace* m_ws;

  TFile* m_dataFile;
  TTree* m_inTree;

  NamedVarStore<double> m_fitVars;
};

EventContent::EventContent(const std::string& dataFile) : m_jpsi(NULL), m_chic(NULL), m_fitFile(NULL), m_ws(NULL)
{
  m_dataFile = TFile::Open(dataFile.c_str());
  if (!m_dataFile) std::cerr << "File Openining failure" << std::endl;
  m_inTree = static_cast<TTree*>(m_dataFile->Get("selectedData"));
  if (!m_inTree) std::cerr << "TTree retrieval error" << std::endl;

  m_inTree->SetBranchAddress("jpsi", &m_jpsi);
  m_inTree->SetBranchAddress("chic_rf", &m_chic);
  m_inTree->SetBranchAddress("Jpsict", &m_jpsict);
}

void EventContent::getWorkspace(const std::string& name)
{
  m_fitFile = TFile::Open(name.c_str());
  if (!m_fitFile) std::cerr << "Fit File Opening failure" << std::endl;
  m_ws = static_cast<RooWorkspace*>(m_fitFile->Get("ws_masslifetime"));
}

void EventContent::getRegionValues()
{
  m_fitVars.setFromWS(m_ws, "CBsigma_p0_jpsi", "p0");
  m_fitVars.setFromWS(m_ws, "CBsigma_p1_jpsi", "p1");
  m_fitVars.setFromWS(m_ws, "CBsigma_p2_jpsi", "p2");
  m_fitVars.setFromWS(m_ws, "CBmass_p0_jpsi", "CBmass_p0");

  m_fitVars.setFromWS(m_ws, "var_sig1MinMass", "massMinSR1");
  m_fitVars.setFromWS(m_ws, "var_sig2MinMass", "massMinSR2");
  m_fitVars.setFromWS(m_ws, "var_sig1MaxMass", "massMaxSR1");
  m_fitVars.setFromWS(m_ws, "var_sig2MaxMass", "massMaxSR2");

  m_fitVars.setFromWS(m_ws, "var_PRMin", "PRmin");
  m_fitVars.setFromWS(m_ws, "var_PRMax", "PRmax");
}

bool EventContent::inPtRapRange(int iPt, int iRap)
{
  double pt = m_chic->Pt();
  double absY = TMath::Abs(m_chic->Rapidity());

  return ( pt >= onia::pTRange[iRap][iPt-1] && pt < onia::pTRange[iRap][iPt] &&
           absY >= onia::rapForPTRange[iRap-1] && absY < onia::rapForPTRange[iRap] );
}

/** hacky, but works and this has to win no beauty contest. */
enum hists {
  all = 0, /**< all jpsis in (rap, pt) bin. */
  jpsi = 1, /**< all jpsis after applying jpsi mass cut. */
  chicM1 = 2, /**< all jpsis after applying jpsi mass cut + chic mass cut for SR1. */
  chicM2 = 3, /**< all jpsis after applying jpsi mass cut + chic mass cut for SR2. */
  chicL = 4, /**< all jpsis after applying jpsi mass cut + chic lifetime cut. */
  chicM1L = 5, /**< all jpsis after applying all cuts of the analysis for SR2. */
  chicM2L = 6 /**< all jpsis after applying all cuts of the analysis for SR2. */
};



struct OutputFile {
  OutputFile(int iPt, int iRap);
  void writeAndCloseFile();
  void setAxis();

  TFile* m_filePtr;
  std::vector<TH1D*> m_hists;
};

OutputFile::OutputFile(int iPt, int iRap)
{
  std::stringstream filename;
  filename << "jpsiMassHistos_rap" << iRap << "_pt" << iPt << ".root";
  m_filePtr = new TFile(filename.str().c_str(), "RECREATE");
  m_hists.push_back(new TH1D("all", /*"all J/#psi"*/"", 50, 2.95, 3.25));
  m_hists.push_back(new TH1D("jpsi", /*"J/#psi after jpsi mass cut"*/"", 50, 2.95, 3.25));
  m_hists.push_back(new TH1D("chicM1", /*"J/#psi after chic mass cut for SR1"*/"", 50, 2.95, 3.25));
  m_hists.push_back(new TH1D("chicM2", /*"J/#psi after chic mass cut for SR2"*/"", 50, 2.95, 3.25));
  m_hists.push_back(new TH1D("chicL", /*"J/#psi after chic lifetime cut"*/"", 50., 2.95, 3.25));
  m_hists.push_back(new TH1D("chicML1", /*"J/#psi after chic mass (SR1) & lifetime cut"*/"", 50, 2.95, 3.25));
  m_hists.push_back(new TH1D("chicML2", /*"J/#psi after chic mass (SR2) & lifetime cut"*/"", 50, 2.95, 3.25));

  setAxis();
}

void OutputFile::writeAndCloseFile()
{
  m_filePtr->cd();
  for (size_t i = 0; i < m_hists.size(); ++i) {
    m_hists[i]->Write();
  }
  m_filePtr->Write();
  m_filePtr->Close();
}

void OutputFile::setAxis()
{
  for (size_t i = 0; i < m_hists.size(); ++i) {
    m_hists[i]->GetXaxis()->SetTitle("M_{J/#psi}");
  }
}


/////////////////////////
// PROGRAM STARTS HERE //
/////////////////////////
/** extract the pt and rap bin numbers from the name of the workspace. NOTE: A LOT OF ASSUMPTIONS GO IN HERE!!*/
std::pair<int, int> getPtRapFromName(const std::string& name)
{
  std::vector<std::string> wsNameParts = splitString(name, '_');
  int iRap = atoi(wsNameParts[wsNameParts.size()-2].substr(3).c_str());

  wsNameParts = splitString(wsNameParts.back(), '.');
  int iPt = atoi(wsNameParts[0].substr(2).c_str());

  return std::make_pair(iPt, iRap);
}

inline bool jpsiMassAccepted(const TLorentzVector* jpsi, double mass, double p0, double p1, double p2, double y)
{
  return TMath::Abs(jpsi->M() - mass) < onia::nSigMass * rapSigma(p0, p1, p2, y);
}

/** setup the Legend. NOTE: you take ownership of the legend! */
TLegend* setupLegend(double x1, double y1, double x2, double y2)
{
  TLegend* leg = new TLegend(x1, y1, x2, y2);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(72);
  leg->SetTextSize(0.025);
  leg->SetBorderSize(1);

  return leg;
}

TH1D* divide(TH1D* h1, TH1D* h2)
{
  TH1D* h = static_cast<TH1D*>(h1->Clone());
  h->Divide(h1,h2);
  return h;
}

void makeRatioPlot(TH1D* h1, TH1D* h2, const std::string& filename)
{
  TCanvas c1("pCanvas", "pCanvas", 500, 500);
  TH1D* h = divide(h1,h2);
  c1.cd();
  h->Draw();
  c1.SaveAs(filename.c_str());
}

void makePdfs(const OutputFile& file)
{
  gStyle->SetOptStat(0);

  std::string filenameBase = file.m_filePtr->GetName();
  filenameBase = filenameBase.substr(0, filenameBase.length() - 5);

  makeRatioPlot(file.m_hists[chicM1L], file.m_hists[jpsi], filenameBase + "_PRSR1_Jpsi.pdf");
  makeRatioPlot(file.m_hists[chicM1], file.m_hists[jpsi], filenameBase + "_SR1_Jpsi.pdf");
  makeRatioPlot(file.m_hists[chicL], file.m_hists[jpsi], filenameBase + "_PR_Jpsi.pdf");
  makeRatioPlot(file.m_hists[chicM2], file.m_hists[jpsi], filenameBase + "_SR2_Jpsi.pdf");
  makeRatioPlot(file.m_hists[chicM2L], file.m_hists[jpsi], filenameBase + "_PRSR2_Jpsi.pdf");

  TCanvas c1("plotCanvas", "Jpsi mass distributions", 500, 500);
  c1.cd();

  file.m_hists[chicM1L]->DrawNormalized();
  file.m_hists[chicM1]->SetLineColor(kRed);
  file.m_hists[chicM1]->DrawNormalized("same");
  file.m_hists[chicL]->SetLineColor(kGreen);
  file.m_hists[chicL]->DrawNormalized("same");
  file.m_hists[jpsi]->SetLineColor(kBlack);
  file.m_hists[jpsi]->DrawNormalized("same");

  TLegend* leg = setupLegend(0.1, 0.7, 0.4, 0.9);
  leg->AddEntry(file.m_hists[chicM1L], "Analysis cuts", "l");
  leg->AddEntry(file.m_hists[chicM1], "chic mass cut SR1", "l");
  leg->AddEntry(file.m_hists[chicL], "chic lifetime cut", "l");
  leg->AddEntry(file.m_hists[jpsi], "jpsi (signal)", "l");

  leg->Draw();
  std::string savename = filenameBase + "_SR1.pdf";
  c1.SaveAs(savename.c_str());

  c1.Clear();
  c1.cd();

  file.m_hists[chicM2L]->DrawNormalized();
  file.m_hists[chicM2]->SetLineColor(kRed);
  file.m_hists[chicM2]->DrawNormalized("same");
  file.m_hists[chicL]->SetLineColor(kGreen);
  file.m_hists[chicL]->DrawNormalized("same");
  file.m_hists[jpsi]->SetLineColor(kBlack);
  file.m_hists[jpsi]->DrawNormalized("same");

  leg = setupLegend(0.1, 0.7, 0.4, 0.9);
  leg->AddEntry(file.m_hists[chicM2L], "Analysis cuts", "l");
  leg->AddEntry(file.m_hists[chicM2], "chic mass cut SR2", "l");
  leg->AddEntry(file.m_hists[chicL], "chic lifetime cut", "l");
  leg->AddEntry(file.m_hists[jpsi], "jpsi (signal)", "l");

  leg->Draw();
  savename = filenameBase + "_SR2.pdf";
  c1.SaveAs(savename.c_str());
}


/** "main function. */
void plotJpsiDistributions(const std::string& dataFile, const std::vector<std::string>& inputWorkspaces)
{
  EventContent events(dataFile);

  for (size_t iInput = 0; iInput < inputWorkspaces.size(); ++iInput) {
    events.getWorkspace(inputWorkspaces[iInput]);
    std::pair<int, int> ptRap = getPtRapFromName(inputWorkspaces[iInput]);
    OutputFile outFile(ptRap.first, ptRap.second);

    events.getRegionValues();
    double p0 = events.m_fitVars["p0"];
    double p1 = events.m_fitVars["p1"];
    double p2 = events.m_fitVars["p2"];
    double CBmass_p0 = events.m_fitVars["CBmass_p0"];
    const BkgHistoRange sr1Range(events.m_fitVars["massMinSR1"], events.m_fitVars["massMaxSR1"]);
    const BkgHistoRange sr2Range(events.m_fitVars["massMinSR2"], events.m_fitVars["massMaxSR2"]);
    const BkgHistoRange prRange(events.m_fitVars["PRmin"], events.m_fitVars["PRmax"]);

    unsigned allAcc = 0;
    unsigned jpsiAcc = 0;
    unsigned m1Acc = 0;
    unsigned m2Acc = 0;
    unsigned lAcc = 0;
    unsigned m1lAcc = 0;
    unsigned m2lAcc = 0;
    for (int i = 0; i < events.m_inTree->GetEntries(); ++i) {
      events.getEvent(i);
      if (!events.inPtRapRange(ptRap.first, ptRap.second)) continue;
      double m = events.m_jpsi->M();
      outFile.m_hists[all]->Fill(m);
      allAcc++;

      if (jpsiMassAccepted(events.m_jpsi, CBmass_p0, p0, p1, p2, events.m_jpsi->Rapidity())) {
        outFile.m_hists[jpsi]->Fill(m);
        jpsiAcc++;
        if (events.inRangeM(sr1Range)) {
          outFile.m_hists[chicM1]->Fill(m);
          m1Acc++;
        }
        if (events.inRangeM(sr2Range)) {
          outFile.m_hists[chicM2]->Fill(m);
          m2Acc++;
        }
        if (events.inRangeL(prRange)) {
          outFile.m_hists[chicL]->Fill(m);
          lAcc++;
        }
        if (events.inRangeM(sr1Range) && events.inRangeL(prRange)) {
          outFile.m_hists[chicM1L]->Fill(m);
          m1lAcc++;
        }
        if (events.inRangeM(sr2Range) && events.inRangeL(prRange)) {
          outFile.m_hists[chicM2L]->Fill(m);
          m2lAcc++;
        }
      }
    }

    std::cout << "=============== " << ptRap.first << ", " << ptRap.second << " ===============" << std::endl;
    std::cout << "Jpsi / all " << jpsiAcc << " / " << allAcc << std::endl;
    std::cout << "SR1 (mass) " << m1Acc << ", SR2 (mass) " << m2Acc << std::endl;
    std::cout << "lifetime " << lAcc << std::endl;
    std::cout << "final PRSR1 " << m1lAcc << ", PRSR2 " << m2lAcc << std::endl;
    std::cout << "====================================" << std::endl;

    makePdfs(outFile);
    outFile.writeAndCloseFile();
  }
}


#ifndef __CINT__
int main (int argc, char* argv[])
{
  std::vector<std::string> inputWorkspaces;
  for (int i = 2; i < argc; ++i) inputWorkspaces.push_back(argv[i]);

  plotJpsiDistributions(argv[1], inputWorkspaces);
}
#endif
