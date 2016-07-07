#include "TFile.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"

#include <string>
#include <vector>
#include <iostream>
#include <sstream>

/** get all efficiencies with the name bas effname and ending between etaMin and etaMax from file. */
template<typename EffType>
std::vector<EffType*> getEffsFromFile(TFile* file, const std::string& effname, int etaMin, int etaMax)
{
  std::vector<EffType*> vec;
  for (int i = etaMin; i < etaMax; ++i) {
    std::stringstream name;
    name << effname << i;
    vec.push_back(static_cast<EffType*>(file->Get(name.str().c_str())));
  }
  return vec;
}

/** set marker and line color of passed function to color. */
template<typename FuncType>
void setColor(FuncType* f, int color)
{
  f->SetMarkerColor(color);
  f->SetLineColor(color);
}

/** setup the Legend. NOTE: you take ownership of the legend! */
TLegend* setupLegend(double x1, double y1, double x2, double y2)
{
  TLegend* leg = new TLegend(x1, y1, x2, y2);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(72);
  leg->SetTextSize(0.035);
  leg->SetBorderSize(1);

  return leg;
}

template<typename FuncType>
void plotAndAddToLegend(FuncType* f, TLegend* leg, const std::string& legEntry, TCanvas* c)
{
  f->Draw();
  leg->AddEntry(f, legEntry.c_str(), "ple");
  c->Update();
}

/** plot efficiencies. takes two input filenames. fn1 = hybrid + MCTRUTH, fn2 = sigmoid. */
void plotEff(const std::string& fn1, const std::string& fn2)
{
  TFile* f1 = TFile::Open(fn1.c_str(), "READ");
  TFile* f2 = TFile::Open(fn2.c_str(), "READ");

  std::vector<TGraphAsymmErrors*> hybridEffs = getEffsFromFile<TGraphAsymmErrors>(f1, "gEffHybrid_DATA_PT_AETA", 0, 9);
  std::vector<TGraphAsymmErrors*> mcTruthEffs = getEffsFromFile<TGraphAsymmErrors>(f1, "gEff_MCTRUTH_PT_AETA", 0, 9);
  std::vector<TF1*> sigEffs = getEffsFromFile<TF1>(f2, "fitTotEff_DATA_pt_etaBin", 0, 9);

  for (size_t i = 0; i < 9; ++i) { // WARNING: hardcoded loop variables!
    TCanvas* c = new TCanvas("plotting Canvas", "", 500, 500);
    c->cd();

    std::stringstream title;
    title << "eff_comparison_eta_" << i;
    TH1F* plotHisto = c->DrawFrame(0., 0., 70., 1.1);
    plotHisto->SetXTitle("p_{T} [GeV/c]");
    plotHisto->SetYTitle("single muon eff");
    plotHisto->SetTitle(title.str().c_str());

    setColor(hybridEffs[i], kBlue);
    setColor(mcTruthEffs[i], kGreen);
    setColor(sigEffs[i], kRed);

    TLegend* plotLeg = setupLegend(0.72, 0.1, 0.95, 0.35);
    plotAndAddToLegend(sigEffs[i], plotLeg, "sigmoid", c);
    plotAndAddToLegend(hybridEffs[i], plotLeg, "hybrid", c);
    plotAndAddToLegend(mcTruthEffs[i], plotLeg, "mcTruth", c);


    plotLeg->Draw();
    c->Update();
    title << ".pdf";
    c->SaveAs(title.str().c_str());

    c->Close();
    delete c;
    delete plotLeg;
  }
}

#ifndef __CINT__
int main (int argc, char* argv[])
{
  plotEff(argv[1], argv[2]);
}
#endif
