#ifndef __CINT__
#endif

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
using namespace std;

#include "bkgHistos_helper.h"
#include "calculatePar.cc"
// #include "calcPol.C" // already in bkgHistos_helper.h
#include "bkgHistos.C"
#include "bkgHistos_chi.C"

#include "clarg_parsing.h"
//========================================================
// code to read input arguments

//===================================================
int main(int argc, char* argv[]){

  // Set defaults
  int
    rapMin = 999,
    rapMax = 999,
    ptMin = 999,
    ptMax = 999,
    nState = 999,
    ctauScen = 999,
    FracLSB = 999;
  bool
    MC = false,
    doCtauUncer = false,
    PolLSB = false,
    PolRSB = false,
    PolNP = false,
    forceBinning = false,
    folding = false,
    normApproach = false,
    scaleFracBg = false,
    subtractNP = false;
  std::string polDataPath;
  bool useRefittedChic = true;

  // Loop over argument list
  for (int i=1; i < argc; i++) {
    std::string arg = argv[i];
    fromSplit("rapMin", arg, rapMin);
    fromSplit("rapMax", arg, rapMax);
    fromSplit("ptMin", arg, ptMin);
    fromSplit("ptMax", arg, ptMax);
    fromSplit("nState", arg, nState);
    fromSplit("ctauScen", arg, ctauScen);
    fromSplit("FracLSB", arg, FracLSB);
    fromSplit("MC", arg, MC);
    fromSplit("doCtauUncer", arg, doCtauUncer);
    fromSplit("PolLSB", arg, PolLSB);
    fromSplit("PolRSB", arg, PolRSB);
    fromSplit("PolNP", arg, PolNP);
    fromSplit("forceBinning", arg, forceBinning);
    fromSplit("folding", arg, folding);
    fromSplit("normApproach", arg, normApproach);
    fromSplit("scaleFracBg", arg, scaleFracBg);
    fromSplit("polDataPath", arg, polDataPath);
    fromSplit("subtractNP", arg, subtractNP);
    fromSplit("useRefittedChic", arg, useRefittedChic);
  }

  std::cout << "-----------------------\n"
            << "Creating background model for \n"
            << "y bins " << rapMin << " - " << rapMax << "\n"
            << "and pT bins "  << ptMin << " - " << ptMax << "\n"
            << "-----------------------" << std::endl;

  for(int iRap = rapMin; iRap <= rapMax; iRap++){
    for(int iPT = ptMin; iPT <= ptMax; iPT++){

      if (nState == 4 || nState == 5) {
        std::stringstream temp;
        temp << "tmpFiles/fit_Psi" << nState-3 << "S_rap" << iRap << "_pt" << iPT << ".root";
        std::string infilename = temp.str().c_str();

        bkgHistos(infilename.c_str(), iRap, iPT, nState, folding, MC, doCtauUncer, PolLSB, PolRSB, PolNP,
                  ctauScen, FracLSB, forceBinning, normApproach, scaleFracBg, polDataPath.c_str());
      }
      else if (nState == 6) {
        std::stringstream temp;
        temp << "tmpFiles/backupWorkSpace/ws_DefineRegionsAndFractions_Chi_rap" << iRap << "_pt" << iPT << ".root";
        std::string infilename = temp.str().c_str();

        bkgHistos_chi(infilename.c_str(), iRap, iPT, folding, MC, PolLSB, PolRSB, PolNP, FracLSB, normApproach, subtractNP, useRefittedChic);
      }
    }
  }

  return 0;
}
