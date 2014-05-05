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

#include "PlotJpsiMassLifetime.cc"

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
    	    TFile *infileTest = new TFile(temp.str().c_str(), "UPDATE");
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
