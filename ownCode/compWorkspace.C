#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"

#include <string>
#include <vector>
#include <iostream>

void compWorkspace(const std::vector<std::string>& ifns, const std::string& ofn, const std::string& wsname,
                   const std::string& variable, const std::string& dataset)
{
  std::vector<TFile*> inputFiles;
  for (size_t i = 0; i < ifns.size(); ++i) { inputFiles.push_back(TFile::Open(ifns[i].c_str())); }
  std::cout << "=============== FILES END ===============" << std::endl;

  std::vector<RooWorkspace*> workSpaces;
  for (size_t i = 0; i < inputFiles.size(); ++i) {
    if (inputFiles[i]) {
      workSpaces.push_back( (RooWorkspace*) inputFiles[i]->Get(wsname.c_str()) );
    } else {
      std::cerr << "invalid TFile pointer: " << i << std::endl;
    }
  }
  std::cout << "=============== WORKSPACES END ===============" << std::endl;

  std::vector<RooRealVar*> variables;
  for (size_t i = 0; i < workSpaces.size(); ++i) {
    if (workSpaces[i]) {
      variables.push_back( workSpaces[i]->var(variable.c_str()) );
    } else {
      std::cerr << "invalid RooWorkspace pointer: " << i << std::endl;
    }
  }
  std::cout << "=============== VARIABLES END ===============" << std::endl;

  std::vector<RooDataSet*> datasets;
  for (size_t i = 0; i < workSpaces.size(); ++i) {
    if (workSpaces[i]) {
      datasets.push_back( (RooDataSet*) workSpaces[i]->data(dataset.c_str()) );
    } else {
      std::cerr << "invalid RooWorkspace pointer: " << i << std::endl;
    }
  }
  std::cout << "=============== DATASETS END ===============" << std::endl;

  RooPlot* plotFrame = variables[0]->frame();
  std::vector<RooPlot*> plots;
  for (size_t i = 0; i < variables.size(); ++i) {
    // plots.push_back( variables[i]->frame() );
    // datasets[i]->plotOn(plots.back());
    datasets[i]->plotOn(plotFrame, RooFit::MarkerColor(i+1), RooFit::LineColor(i+1));
  }

  TCanvas can("canvas", "bla", 1000,800);
  // can.Divide(1,2);
  // can.cd(1); plots[0]->Draw();
  // can.cd(2); plots[1]->Draw();

  // TLegend* leg = new TLegend(0.4,0.6,0.89,0.89);
  // leg->AddEntry(plots[0], "first", "l");
  // leg->AddEntry(plots[1], "second", "l");

  can.cd();
  plotFrame->Draw();
  // leg->Draw();
  can.SaveAs(ofn.c_str());
}

#ifndef __CINT__
int main(int argc, char* argv[])
{
  // TODO: proper clarg parsing
  compWorkspace({std::string(argv[1]), std::string(argv[2])}, std::string(argv[3]), "ws_masslifetime",
                std::string(argv[4]), "jpsi_data_rap1_pt4");

  return 0;
}
#endif
