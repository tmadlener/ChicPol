/*
 * runDefineRegionsAndFractions.cc
 *
 *  Created on: Jan 30, 2014
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

#include "DefineRegionsAndFractions.cc"

//===================================================
int main(int argc, char* argv[]){

  // Set defaults
  int rapMin = 999,
    rapMax = 999,
    ptMin = 999,
    ptMax = 999,
    nState = 999,
    rapFixTo = 999,
    ptFixTo = 999;

  bool runChiMassFitOnly = false;
  bool FixRegionsToInclusiveFit = false;
  bool doFractionUncer=false;
  double chic1MassLow = -1;
  double chic2MassLow = -1;
  double chic1MassHigh = -1;
  double chic2MassHigh = -1;
  double nSigPR = -1;
  double nSigNP = -1;
  
  // Loop over argument list
  for (int i=1; i < argc; i++)
    {
      std::string arg = argv[i];
      fromSplit("runChiMassFitOnly", arg, runChiMassFitOnly);
      fromSplit("doFractionUncer", arg, doFractionUncer);
      fromSplit("FixRegionsToInclusiveFit", arg, FixRegionsToInclusiveFit);
      fromSplit("rapMin", arg, rapMin);
      fromSplit("rapMax", arg, rapMax);
      fromSplit("ptMin", arg, ptMin);
      fromSplit("ptMax", arg, ptMax);
      fromSplit("nState", arg, nState);
      fromSplit("rapFixTo", arg, rapFixTo);
      fromSplit("ptFixTo", arg, ptFixTo);
      fromSplit("chic1MassLow", arg, chic1MassLow);
      fromSplit("chic2MassLow", arg, chic2MassLow);
      fromSplit("chic1MassHigh", arg, chic1MassHigh);
      fromSplit("chic2MassHigh", arg, chic2MassHigh);
      fromSplit("nSigPR", arg, nSigPR);
      fromSplit("nSigNP", arg, nSigNP);
    }

  std::cout << "-----------------------\n"
            << "Defining regions and calculating fractions for \n"
            << "y bins " << rapMin << " - " << rapMax << "\n"
            << "and pT bins "  << ptMin << " - " << ptMax << "\n"
            << "-----------------------" << std::endl;

  for(int iRap = rapMin; iRap <= rapMax; iRap++){
    for(int iPT = ptMin; iPT <= ptMax; iPT++){

      std::stringstream tempFrom;
      tempFrom << "tmpFiles/backupWorkSpace/ws_MassLifetimeFit_Chi_rap" << iRap << "_pt" << iPT << ".root";
      const std::string infilenameFrom = tempFrom.str().c_str();

      std::stringstream tempTo;
      tempTo << "tmpFiles/backupWorkSpace/ws_DefineRegionsAndFractions_Chi_rap" << iRap << "_pt" << iPT << ".root";
      const std::string infilenameTo = tempTo.str().c_str();

      gSystem->CopyFile(infilenameFrom.c_str(),infilenameTo.c_str(),kTRUE);

      if(FixRegionsToInclusiveFit&&iRap==rapFixTo&&iPT==ptFixTo) FixRegionsToInclusiveFit=false;
      DefineRegionsAndFractions(infilenameTo.c_str(), iRap, iPT, nState, runChiMassFitOnly, FixRegionsToInclusiveFit, rapFixTo, ptFixTo, doFractionUncer, chic1MassLow, chic1MassHigh, chic2MassLow, chic2MassHigh, nSigPR, nSigNP);

    }
  }

  return 0;
}
