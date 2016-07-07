#include "../interface/commonVar.h"

#include "compPol_helper.h"

#include "TFile.h"
// #include "TTree.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"

#include <string>
#include <vector>
#include <iostream>
#include <map>

const std::string rapBin = "rap1";
const std::vector<int> markerStyle = {24, 25, 26};
const std::vector<int> markerColor = {8, 9, 2, 6, 1};
const int markerSize = 8;

// UGLY: global variables for access across different functions
// gets set in main by a call to setupGlobals;
std::string yAxisTitle;
std::vector<std::string> graphNames;
std::vector<std::string> legendEntries;
std::vector<std::string> frames;

/** create and setup the canvas for plotting. NOTE: you take ownership of the returned pointer! */
TCanvas* setupCanvas(const std::string& identifier, const std::string& title)
{
  TCanvas* canv = new TCanvas(identifier.c_str(), title.c_str(), 1000, 800);
  canv->SetFillColor(kWhite);
  canv->SetGrid();
  canv->GetFrame()->SetFillColor(kWhite);
  canv->GetFrame()->SetBorderSize(0);
  canv->SetRightMargin(0.05);

  return canv;
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

void setStyle()
{
  gStyle->SetPalette(1,0);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.15);

  gStyle->SetTickLength(-0.02, "xyz");
  gStyle->SetLabelOffset(0.02, "x");
  gStyle->SetLabelOffset(0.02, "y");
  gStyle->SetTitleOffset(1.3, "x");
  gStyle->SetTitleOffset(1.4, "y");
  gStyle->SetTitleFillColor(kWhite);
}

/** main function plotting the same histogram of the two different files. */
void compBoundaries(const std::vector<std::string>& inputfns, const std::string& ofn)
{
  std::vector<TFile*> inputFiles;
  // for (const auto& fn : inputfns) { // no range-based for loops pre C++11
  for (size_t i=0; i < inputfns.size(); ++i) {
    const std::string& fn = inputfns[i];
    inputFiles.push_back(TFile::Open(fn.c_str()));
  }

  std::vector<std::vector<TGraphAsymmErrors*> > graphs;
  for (size_t iFile = 0; iFile < inputFiles.size(); ++iFile) {
    std::vector<TGraphAsymmErrors*> fileGraphs;
    for (size_t iGraph = 0; iGraph < graphNames.size(); ++iGraph) {
      const std::string& graphID = graphNames[iGraph];
      fileGraphs.push_back( (TGraphAsymmErrors*) inputFiles[iFile]->Get(graphID.c_str()) );
    }
    if(!fileGraphs.empty()) graphs.push_back(fileGraphs);
  }

  TCanvas* compCanvas = setupCanvas("compCanvas", "compCanvas");
  compCanvas->cd();

  TLegend* plotLegend = setupLegend(0.72, 0.1, 0.95, 0.35); // 0.825,0.75,0.95,0.9

  TH1F* plotHisto = compCanvas->DrawFrame(onia::pTRange[1][0], -1.1, onia::pTRange[1][5], 1.1);
  plotHisto->SetXTitle("p_{T} [GeV/c]");
  plotHisto->SetYTitle(yAxisTitle.c_str());
  plotHisto->GetYaxis()->SetTitleOffset(1.5);

  for (size_t iFile = 0; iFile < graphs.size(); ++iFile) {
    for (size_t iGraph = 0; iGraph < graphs[iFile].size(); ++iGraph) {
      TGraphAsymmErrors* graph = graphs[iFile][iGraph];
      graph->SetMarkerColor(markerColor[iFile]);
      graph->SetMarkerStyle(markerStyle[iGraph]);
      graph->SetMarkerSize(1.5);
      graph->SetLineColor(markerColor[iFile]);

      graph->Draw("P");
      std::string legEntry = frames[iGraph] + " " + legendEntries[iFile];
      plotLegend->AddEntry(graph, legEntry.c_str(), "ple");
      compCanvas->Update();
    }
  }

  TLatex* text = new TLatex(onia::pTRange[1][5]*0.22, -0.88, " |y^{#chi}| < 1.2");
  text->SetTextSize(0.035);
  text->Draw("same");

  plotLegend->Draw();
  compCanvas->Update();
  compCanvas->SaveAs(ofn.c_str());

  compCanvas->Close();
  delete compCanvas;

  // WARNING: I don't know for which of these objects ROOT has actually taken ownership! Deleting might cause crashes
  // delete plotHisto;
  // delete plotLegend;
  // delete text;
}

#ifndef __CINT__
int main (int argc, char* const argv[])
{
  std::vector<std::string> args;
  for (int i = 1; i < argc; ++i) { args.push_back( std::string(argv[i])); }
  cl_args arguments(args, rapBin);

  yAxisTitle = arguments.yAxisTitle;
  graphNames = arguments.graphs;
  legendEntries = arguments.legendEntries;
  frames = arguments.frames;

  compBoundaries(arguments.ifn, arguments.ofn);
  return 0;
}
#endif
