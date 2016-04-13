#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include "TROOT.h"

#include "clarg_parsing.h"

#include "createWorkspace.C"

using namespace onia;


//=====================================================================
int main(int argc, char* argv[]){

  // set default values
  int nState = 999;
  bool correctCtau  = false;
  bool drawRapPt2D  = false;
  bool useRefittedMass = true;

  // Loop over argument list
  for (int i=1; i < argc; i++) {
    std::string arg = argv[i];
    fromSplit("nState", arg, nState);
    fromSplit("correctCtau", arg, correctCtau);
    fromSplit("drawRapPt2D", arg, drawRapPt2D);
    fromSplit("useRefittedChic", arg, useRefittedMass);
  }

  const std::string infilename = "tmpFiles/selEvents_data.root";
  createWorkspace(infilename, nState, correctCtau, drawRapPt2D, useRefittedMass);

  return 0;
}
