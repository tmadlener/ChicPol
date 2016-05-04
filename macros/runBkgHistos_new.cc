// temp file for development and compiling

#include "BkgHistoProducerFactory.h"
#include "BkgHistoProducer.h"

#include "clarg_parsing.h"

#include <string>
#include <sstream>

int main(int argc, char* argv[])
{
  int nState = -1;
  int rapMin = -1;
  int rapMax = -1;
  int ptMin = -1;
  int ptMax = -1;

  for (int iArg = 0; iArg < argc; ++iArg) {
    const std::string arg = std::string(argv[iArg]);
    fromSplit("nState", arg, nState);
    fromSplit("rapMin", arg, rapMin);
    fromSplit("rapMax", arg, rapMax);
    fromSplit("ptMin", arg, ptMin);
    fromSplit("ptMax", arg, ptMax);
  }

  BkgHistoProducerFactory factory;
  IBkgHistoProducer* bkgHistoProducer = factory.createBkgHistoProducer(nState);
  if (!bkgHistoProducer) {
    return -1; // exit early with error state -1
  }

  std::string infileBase = "tmpFiles/backupWorkSpace/ws_";
  if (nState == 4 || nState == 5) {
    infileBase += "MassFit_Jpsi";
  } else if (nState == 6) {
    // infileBase += "DefineRegionsAndFractions_Chi"; // for after dev
    infileBase += "MassFit_Jpsi"; // for development
  }

  for (int iRap = rapMin; iRap <= rapMax; ++iRap) {
    for (int iPt = ptMin; iPt <= ptMax; ++iPt) {
      std::stringstream infilename;
      infilename << infileBase << "_rap" << iRap << "_pt" << iPt << ".root";
      bkgHistoProducer->initialize(infilename.str(), iRap, iPt);
      bkgHistoProducer->fillHistos();
      bkgHistoProducer->storeHistos();
    }
  }

  delete bkgHistoProducer;

  return 0;
}
