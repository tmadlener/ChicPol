/*
 * runPlotJpsiMassLifetime.cc
 *
 *  Created on: Apr 23, 2014
 *      Author: valentinknuenz
 */




#ifndef __CINT__
#endif

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
using namespace std;

#include "clarg_parsing.h"

#include "PlotJpsiMassLifetime.cc"

//===================================================
int main(int argc, char* argv[]){
  // Set defaults
  int
    rapMin = 999,
    rapMax = 999,
    ptMin = 999,
    ptMax = 999,
    nState = 999,
    PlottingJpsi = 999;

  // Loop over argument list
  for (int i=1; i < argc; i++)
    {
      std::string arg = argv[i];
      fromSplit("rapMin", arg, rapMin);
      fromSplit("rapMax", arg, rapMax);
      fromSplit("ptMin", arg, ptMin);
      fromSplit("ptMax", arg, ptMax);
      fromSplit("nState", arg, nState);
      fromSplit("PlottingJpsi", arg, PlottingJpsi);
    }

  std::cout << "-----------------------\n"
            << "Plot \n"
            << "y bins " << rapMin << " - " << rapMax << "\n"
            << "and pT bins "  << ptMin << " - " << ptMax << "\n"
            << "-----------------------" << std::endl;

  for(int iRap = rapMin; iRap <= rapMax; iRap++){
    for(int iPT = ptMin; iPT <= ptMax; iPT++){

      std::stringstream temp;
      std::stringstream temp2;
      temp << "tmpFiles/backupWorkSpace/ws_MassLifetimeFit_Jpsi_rap" << iRap << "_pt" << iPT << ".root";

      bool MassLifetimeFile=true;
      TFile *infileTest = TFile::Open(temp.str().c_str());
      if(!infileTest){
        std::cout << "Error: failed to open file with dataset" << std::endl;
        temp2 << "tmpFiles/backupWorkSpace/ws_MassFit_Jpsi_rap" << iRap << "_pt" << iPT << ".root";
        MassLifetimeFile=false;
      }
      else infileTest->Close();

      const std::string infilename = temp.str().c_str();
      const std::string infilename2 = temp2.str().c_str();

      if(MassLifetimeFile) PlotMassLifetime(infilename.c_str(), iRap, iPT, nState, PlottingJpsi);
      else PlotMassLifetime(infilename2.c_str(), iRap, iPT, nState, PlottingJpsi);

    }
  }
  return 0;
}
