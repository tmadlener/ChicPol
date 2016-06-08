/*
 * runChiMassLifetimeFit.cc
 *
 *  Created on: Jan 23, 2014
 *      Author: valentinknuenz
 */

#ifndef __CINT__
#endif

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
using namespace std;

#include "chiMassLifetimeFit.cc"
#include "TROOT.h"

#include "clarg_parsing.h"

//===================================================
int main(int argc, char* argv[]){

  // Set defaults
  int rapMin = 999,
    rapMax = 999,
    ptMin = 999,
    ptMax = 999,
    nState = 999;

  bool runChiMassFitOnly = false;
  bool MC = false;
  bool MCclosure = false;
  
  // Loop over argument list
  for (int i=1; i < argc; i++)
    {
      std::string arg = argv[i];
      fromSplit("runChiMassFitOnly", arg, runChiMassFitOnly);
      fromSplit("rapMin", arg, rapMin);
      fromSplit("rapMax", arg, rapMax);
      fromSplit("ptMin", arg, ptMin);
      fromSplit("ptMax", arg, ptMax);
      fromSplit("nState", arg, nState);
      fromSplit("MC", arg, MC);
      fromSplit("MCclosure", arg, MCclosure);
    }

  std::cout << "-----------------------\n"
            << "Fitting mass/lifetime for \n"
            << "y bins " << rapMin << " - " << rapMax << "\n"
            << "and pT bins "  << ptMin << " - " << ptMax << "\n"
            << "-----------------------" << std::endl;

  for(int iRap = rapMin; iRap <= rapMax; iRap++){
    for(int iPT = ptMin; iPT <= ptMax; iPT++){

      std::stringstream tempFrom;
      // tempFrom << "tmpFiles/backupWorkSpace/ws_MassLifetimeFit_Jpsi_rap" << iRap << "_pt" << iPT << ".root";
      tempFrom << "tmpFiles/backupWorkSpace/ws_MassFit_Jpsi_rap" << iRap << "_pt" << iPT << ".root";
      const std::string infilenameFrom = tempFrom.str().c_str();


      std::stringstream tempTo;
      tempTo << "tmpFiles/backupWorkSpace/ws_MassLifetimeFit_Chi_rap" << iRap << "_pt" << iPT << ".root";
      const std::string infilenameTo = tempTo.str().c_str();

      cout<<"copy file "<<infilenameFrom.c_str()<<"to "<<infilenameTo.c_str()<<endl;
      std::cout << "CopyFile: " << gSystem->CopyFile(infilenameFrom.c_str(),infilenameTo.c_str(),kTRUE) << std::endl;
      cout<<"copy file finished"<<endl;

      chiMassLifetimeFit(infilenameTo.c_str(), iRap, iPT, nState, runChiMassFitOnly, MC, MCclosure);

    }
  }

  return 0;
}
