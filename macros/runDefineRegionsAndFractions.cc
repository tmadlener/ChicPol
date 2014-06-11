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

#include "DefineRegionsAndFractions.cc"

//========================================================
// code to read input arguments
template<typename T>
void fromSplit(const std::string& key, const std::string &arg, T& out)
{
  const char delim = '=';
  // Skip if key or delimiter not there
  if ((arg.find(key) == std::string::npos) ||
      (arg.find(delim) == std::string::npos))
    return;

  std::string skey, sval;
  std::stringstream sstr(arg);
  std::getline(sstr, skey, delim); // Dummy read to skip key
  std::getline(sstr, sval, delim); // Get value
  T tout;
  if (!(std::istringstream(sval) >> std::boolalpha >> tout))
    return;
  out = tout;
  std::cout << std::boolalpha << skey << ": "  << out << std::endl;
}

// Special version for string without the conversion
template<>
void fromSplit(const std::string& key, const std::string &arg, std::string &out)
{
  const char delim = '=';
  // Skip if key or delimiter not there
  if ((arg.find(key) == std::string::npos) ||
      (arg.find(delim) == std::string::npos))
    return;
  std::string skey, sval;
  std::stringstream sstr(arg);
  std::getline(sstr, skey, delim); // Dummy read to skip key
  std::getline(sstr, sval, delim); // Get value
  out = sval;
  std::cout << skey << ": "  << out << std::endl;
}


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
			DefineRegionsAndFractions(infilenameTo.c_str(), iRap, iPT, nState, runChiMassFitOnly, FixRegionsToInclusiveFit, rapFixTo, ptFixTo, doFractionUncer);

      }
    }

    return 0;
}






