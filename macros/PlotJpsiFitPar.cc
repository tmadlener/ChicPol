/*
 * PlotJpsiFitPar.cc
 *
 *  Created on: Apr 29, 2014
 *      Author: valentinknuenz
 */

#include <iostream>
#include <string>
#include <sstream>
#include "calculatePar.cc"
#include "TBox.h"
#include "TMarker.h"
using namespace std;
using namespace RooFit;
using namespace onia;

void PlotJpsiFitPar(int  nState=4, int rapMin=0, int rapMax=0, int ptMin=0, int ptMax=0, int rapFixTo=0, int ptFixTo=0, bool AddInclusiveResult=false, double Ymin=0., double Ymax=0., char ParName[200]="default", char ParTitle[200]="default", char Folder[200]="default", char SaveName[200]="default", bool logY=false, int WhatKindOfParIndex=5);
void PlotJpsiRegionFracs(int  nState=4, int rapBin=0, int ptMin=0, int ptMax=0, double Ymin=0., double Ymax=0., char Folder[200]="default", char SaveName[200]="default", bool logY=false, int RegionCode=999);

double legendsize=0.035;

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

//=====================================================================
int main(int argc, char* argv[]){
  // set default values
  int nState = 999;
  bool doCtauUncer = false;
  bool AddInclusiveResult = false;
  int rapMin = 999,
    rapMax = 999,
    ptMin = 999,
    ptMax = 999,
    rapFixTo = 999,
    ptFixTo = 999;

  // Loop over argument list
  for (int i=1; i < argc; i++){
    std::string arg = argv[i];
    fromSplit("nState", arg, nState);
    fromSplit("doCtauUncer", arg, doCtauUncer);
    fromSplit("AddInclusiveResult", arg, AddInclusiveResult);
    fromSplit("rapMin", arg, rapMin);
    fromSplit("rapMax", arg, rapMax);
    fromSplit("ptMin", arg, ptMin);
    fromSplit("ptMax", arg, ptMax);
    fromSplit("rapFixTo", arg, rapFixTo);
    fromSplit("ptFixTo", arg, ptFixTo);
  }

  TFile *outfileCreate;

  char predirname[2000];
  char dirname[2000];
  sprintf(predirname,"Fit/parameter");
  gSystem->mkdir(predirname);
  sprintf(dirname,"%s/jpsi_mass",predirname);
  gSystem->mkdir(dirname);
  outfileCreate = new TFile(Form("%s/Par_chi.root",dirname,"RECREATE"));
  outfileCreate->Close();
  sprintf(dirname,"%s/jpsi_lifetime",predirname);
  gSystem->mkdir(dirname);
  outfileCreate = new TFile(Form("%s/Par_chi.root",dirname,"RECREATE"));
  outfileCreate->Close();
  sprintf(dirname,"%s/jpsi_fitqual",predirname);
  gSystem->mkdir(dirname);
  outfileCreate = new TFile(Form("%s/Par_chi.root",dirname,"RECREATE"));
  outfileCreate->Close();
  sprintf(dirname,"%s/jpsi_fractions",predirname);
  gSystem->mkdir(dirname);
  outfileCreate = new TFile(Form("%s/Par_chi.root",dirname,"RECREATE"));
  outfileCreate->Close();


  double Ymin, Ymax;
  char ParName[200];
  char ParTitle[200];
  char SaveName[200];
  char Folder[200];
  bool logY=false;


  sprintf(Folder, "jpsi_fractions");  Ymin = 0.; Ymax = 1.; logY=false; int rapBin=1;
  sprintf(SaveName, "default");

  for(int rapBin = rapMin; rapBin < rapMax+1; rapBin++){
    for(int RegionCode=1;RegionCode<4;RegionCode++)
      PlotJpsiRegionFracs(nState, rapBin, ptMin, ptMax, Ymin, Ymax, Folder, SaveName, logY, RegionCode);
  }


  int WhatKindOfParIndex=999;

  WhatKindOfParIndex=0; sprintf(SaveName, "CBmass_p0_jpsi"); sprintf(ParName, "CBmass_p0_jpsi"); sprintf(ParTitle,"#mu_{#psi} [GeV]"); sprintf(Folder, "jpsi_mass");  Ymin = 3.093; Ymax = 3.0945; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
  WhatKindOfParIndex=0; sprintf(SaveName, "CBmass_p1_jpsi"); sprintf(ParName, "CBmass_p1_jpsi"); sprintf(ParTitle,"#mu_{p1} [GeV]"); sprintf(Folder, "jpsi_mass");  Ymin = -0.01; Ymax = 0.01; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
  WhatKindOfParIndex=0; sprintf(SaveName, "CBmass_p2_jpsi"); sprintf(ParName, "CBmass_p2_jpsi"); sprintf(ParTitle,"#mu_{p2} [GeV]"); sprintf(Folder, "jpsi_mass");  Ymin = -0.01; Ymax = 0.01; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);

  WhatKindOfParIndex=0; sprintf(SaveName, "CBsigma_p0_jpsi"); sprintf(ParName, "CBsigma_p0_jpsi"); sprintf(ParTitle,"#sigma^{p0}_{#psi} [GeV]"); sprintf(Folder, "jpsi_mass");  Ymin = 0.017; Ymax = 0.024; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
  WhatKindOfParIndex=0; sprintf(SaveName, "CBsigma_p1_jpsi"); sprintf(ParName, "CBsigma_p1_jpsi"); sprintf(ParTitle,"#sigma^{p1}_{#psi} [GeV]"); sprintf(Folder, "jpsi_mass");  Ymin = -0.015; Ymax = 0.015; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
  WhatKindOfParIndex=0; sprintf(SaveName, "CBsigma_p2_jpsi"); sprintf(ParName, "CBsigma_p2_jpsi"); sprintf(ParTitle,"#sigma^{p2}_{#psi} [GeV]"); sprintf(Folder, "jpsi_mass");  Ymin = 0.005; Ymax = 0.02; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);

  WhatKindOfParIndex=0; sprintf(SaveName, "CBalpha_p0_jpsi"); sprintf(ParName, "CBalpha_p0_jpsi"); sprintf(ParTitle,"#alpha^{CB}_{#psi}"); sprintf(Folder, "jpsi_mass");  Ymin = 1.4; Ymax = 2.0; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
  WhatKindOfParIndex=0; sprintf(SaveName, "CBalpha_p1_jpsi"); sprintf(ParName, "CBalpha_p1_jpsi"); sprintf(ParTitle,"#alpha_{p1}"); sprintf(Folder, "jpsi_mass");  Ymin = 0.; Ymax = 0.8; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);



  WhatKindOfParIndex=123321; sprintf(SaveName, "bkgLambda_jpsi"); sprintf(ParName, "bkgLambda_jpsi"); sprintf(ParTitle,"#lambda_{BG} [GeV^{-1}]"); sprintf(Folder, "jpsi_mass");  Ymin = 0.; Ymax = 5.; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
  WhatKindOfParIndex=0; sprintf(SaveName, "fracBkg_jpsi"); sprintf(ParName, "fracBkg_jpsi"); sprintf(ParTitle,"f^{#psi}_{#mu#muBG}"); sprintf(Folder, "jpsi_mass");  Ymin = 0.05; Ymax = 0.12; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
  WhatKindOfParIndex=0; sprintf(SaveName, "CBn_jpsi"); sprintf(ParName, "CBn_jpsi"); sprintf(ParTitle,"n^{CB}_{#psi}"); sprintf(Folder, "jpsi_mass");  Ymin = 1.4; Ymax = 3.6; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
  WhatKindOfParIndex=3; sprintf(SaveName, "var_massres"); sprintf(ParName, "var_massres"); sprintf(ParTitle,"#sigma^{eff.}_{#mu#mu} [GeV]"); sprintf(Folder, "jpsi_mass");  Ymin = 0.022; Ymax = 0.032; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);


  WhatKindOfParIndex=999; sprintf(SaveName, "Chi2ndf_Lifetime_RSB"); sprintf(ParName, "var_chi2ndf_Lifetime_RSB"); sprintf(ParTitle,"#chi^{2} / ndf (lifetime RSB)"); sprintf(Folder, "jpsi_fitqual");  Ymin = 0.; Ymax = 5.; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
  WhatKindOfParIndex=999; sprintf(SaveName, "Chi2ndf_Lifetime_LSB"); sprintf(ParName, "var_chi2ndf_Lifetime_LSB"); sprintf(ParTitle,"#chi^{2} / ndf (lifetime LSB)"); sprintf(Folder, "jpsi_fitqual");  Ymin = 0.; Ymax = 5.; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
  WhatKindOfParIndex=999; sprintf(SaveName, "Chi2ndf_Lifetime_SR"); sprintf(ParName, "var_chi2ndf_Lifetime_SR"); sprintf(ParTitle,"#chi^{2} / ndf (lifetime SR)"); sprintf(Folder, "jpsi_fitqual");  Ymin = 0.; Ymax = 5.; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
  WhatKindOfParIndex=999; sprintf(SaveName, "Chi2ndf_Mass"); sprintf(ParName, "var_jpsi_chi2ndf_mass"); sprintf(ParTitle,"#chi^{2} / ndf (mass)"); sprintf(Folder, "jpsi_fitqual");  Ymin = 0.; Ymax = 5.; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);


  WhatKindOfParIndex=0; sprintf(SaveName, "Sig_ctResolution"); sprintf(ParName, "jpsi_ctResolution"); sprintf(ParTitle,"#sigma^{scale}_{l}"); sprintf(Folder, "jpsi_lifetime");  Ymin = 0.7; Ymax = 1.1; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
  WhatKindOfParIndex=0; sprintf(SaveName, "Sig_ctResolution2"); sprintf(ParName, "jpsi_ctResolution2"); sprintf(ParTitle,"#sigma^{scale2}_{l}"); sprintf(Folder, "jpsi_lifetime");  Ymin = 0.8; Ymax = 2.5; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
  WhatKindOfParIndex=0; sprintf(SaveName, "Sig_fracGauss2"); sprintf(ParName, "jpsi_fracGauss2"); sprintf(ParTitle,"f_{G_{2}}"); sprintf(Folder, "jpsi_lifetime");  Ymin = 0.; Ymax = 0.45; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
  WhatKindOfParIndex=0; sprintf(SaveName, "Sig_fPrompt"); sprintf(ParName, "jpsi_fPrompt"); sprintf(ParTitle,"f^{#psiSR}_{PR}"); sprintf(Folder, "jpsi_lifetime");  Ymin = 0.; Ymax = 1.; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
  WhatKindOfParIndex=3; sprintf(SaveName, "Sig_fNonPrompt"); sprintf(ParName, "jpsi_fNonPrompt"); sprintf(ParTitle,"f^{#psiSR}_{NP}"); sprintf(Folder, "jpsi_lifetime");  Ymin = 0.; Ymax = 1.; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
  WhatKindOfParIndex=0; sprintf(SaveName, "Sig_nonPromptTau"); sprintf(ParName, "jpsi_nonPromptTau"); sprintf(ParTitle,"#tau_{NP} [mm]"); sprintf(Folder, "jpsi_lifetime");  Ymin = 0.25; Ymax = 0.5; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);

  WhatKindOfParIndex=2; sprintf(SaveName, "BK_fracBkg_jpsi"); sprintf(ParName, "jpsi_fBkg"); sprintf(ParTitle,"f^{#psiSR}_{BG}"); sprintf(Folder, "jpsi_lifetime");  Ymin = 0.; Ymax = 0.07; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
  WhatKindOfParIndex=0; sprintf(SaveName, "BK_DSD_TauBkg"); sprintf(ParName, "jpsi_bkgTauDSD"); sprintf(ParTitle,"#tau^{DS}_{BG} [mm]"); sprintf(Folder, "jpsi_lifetime");  Ymin = 0.; Ymax = 0.06; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
  WhatKindOfParIndex=0; sprintf(SaveName, "BK_FD_TauBkg"); sprintf(ParName, "jpsi_bkgTauFD"); sprintf(ParTitle,"#tau^{LS}_{BG} [mm]"); sprintf(Folder, "jpsi_lifetime");  Ymin = 0.; Ymax = 0.4; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
  WhatKindOfParIndex=0; sprintf(SaveName, "BK_NP_TauBkg"); sprintf(ParName, "jpsi_bkgTauSSD"); sprintf(ParTitle,"#tau^{RS}_{BG} [mm]"); sprintf(Folder, "jpsi_lifetime");  Ymin = 0.3; Ymax = 0.45; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
  WhatKindOfParIndex=3; sprintf(SaveName, "BK_fBkgDS"); sprintf(ParName, "jpsi_fBkgDSD"); sprintf(ParTitle,"f_{BG}^{DS}"); sprintf(Folder, "jpsi_lifetime");  Ymin = 0.; Ymax = 0.4; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
  WhatKindOfParIndex=0; sprintf(SaveName, "BK_fBkgRS"); sprintf(ParName, "jpsi_fBkgSSDR"); sprintf(ParTitle,"f_{BG}^{RS}"); sprintf(Folder, "jpsi_lifetime");  Ymin = 0.5; Ymax = 0.9; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
  WhatKindOfParIndex=0; sprintf(SaveName, "BK_fBkgLS"); sprintf(ParName, "jpsi_fBkgSSDL"); sprintf(ParTitle,"f_{BG}^{LS}"); sprintf(Folder, "jpsi_lifetime");  Ymin = 0.; Ymax = 0.1; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);

  WhatKindOfParIndex=3; sprintf(SaveName, "fLSBpsi"); sprintf(ParName, "var_fLSBpsi"); sprintf(ParTitle,"f_{LSB}^{#psi}"); sprintf(Folder, "jpsi_fractions");  Ymin = 0.; Ymax = 1.; logY=false;
  PlotJpsiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);

  //jpsi_fNonPrompt
  //jpsi_fNonPromptLSB
  //jpsi_fNonPromptRSB


  return 0;
}


//=============================================
void PlotJpsiFitPar(int  nState, int rapMin, int rapMax, int ptMin, int ptMax, int rapFixTo, int ptFixTo, bool AddInclusiveResult, double Ymin, double Ymax, char ParName[200], char ParTitle[200], char Folder[200], char SaveName[200], bool logY, int WhatKindOfParIndex){
  int RapBins = rapMax-rapMin+1,
    PtBins  = ptMax-ptMin+1;

  std::stringstream savePath;
  savePath << "Fit/parameter";
  gSystem->mkdir(savePath.str().c_str(),kTRUE);
  savePath << "/" << Folder;
  gSystem->mkdir(savePath.str().c_str(),kTRUE);

  cout<<"Directory: "<<savePath.str().c_str()<<endl;

  double Xmin = 0.,  Xmax = 0.;
  Xmin = 5.001;    Xmax = 54.999;
  //Xmin = 5.001;    Xmax = 84.999;







  bool changeSign=false;
  if(WhatKindOfParIndex==123321){
    WhatKindOfParIndex=0;
    changeSign=true;
    cout<<"changeSign of this parameter: "<<ParName<<endl;
  }

  double pTmean[RapBins][PtBins];
  double pTmean_Err[RapBins][PtBins];
  double pTmean_Incl=-999.;
  double pTmean_Incl_Err=0.;

  double Par[RapBins][PtBins];
  double Par_Err[RapBins][PtBins];
  double Par_Incl=-999.;
  double Par_Incl_Err=0.;

  TFile *inFile;
  RooWorkspace *ws;
  RooDataSet *data;
  char inName[200];
  bool functionVar=false;


  if(AddInclusiveResult){
    cout<<"rapIncl "<<rapFixTo<<"   ptIncl "<<ptFixTo<<endl;

    std::stringstream infileNameStream;
    infileNameStream << "tmpFiles/backupWorkSpace/ws_MassLifetimeFit_Jpsi_rap" << rapFixTo << "_pt" << ptFixTo << ".root";
    const std::string infileName = infileNameStream.str().c_str();

    TFile *infile = new TFile(infileName.c_str(), "READ");
    if(!infile){
      std::cout << "Error: failed to open file with dataset" << std::endl;
    }
    RooWorkspace *ws=(RooWorkspace *)infile->Get("ws_masslifetime");
    if(!ws){
      std::cout << "Error: failed to open workspace " << std::endl;
    }

    //cout<<"infileName: "<<infileName.c_str()<<endl;

    pTmean_Incl = ws->var("var_chicMeanPt")->getVal();
    pTmean_Incl_Err = ws->var("var_chicMeanPt")->getError();
    if(KinParticleChi == false){
      pTmean_Incl = ws->var("var_jpsiMeanPt")->getVal();
      pTmean_Incl_Err = ws->var("var_jpsiMeanPt")->getError();
    }
                
    if(ws->var(ParName)==NULL){
      if(ws->function(ParName)==NULL) return;
      else functionVar=true;
    }
    if(!functionVar){
      Par_Incl = ws->var(ParName)->getVal();
      Par_Incl_Err = ws->var(ParName)->getError();
    }
    else{
      Par_Incl = ws->function(ParName)->getVal();
      Par_Incl_Err = 0.;
    }

    if(changeSign) Par_Incl*=-1.;

    cout<<"Inclusive pTmean: "<<pTmean_Incl<<" +-"<<pTmean_Incl_Err<<endl;
    cout<<"Inclusive "<<ParName<<": "<<Par_Incl<<" +-"<<Par_Incl_Err<<endl;

  }

  int nFixedBins=0;

  for(int rapBin = rapMin; rapBin < rapMax+1; rapBin++){
    int rapArrayIndex=rapBin-rapMin;
    for(int ptBin = ptMin; ptBin < ptMax+1; ptBin++){
      int ptArrayIndex=ptBin-ptMin;

      cout<<"rap "<<rapBin<<"   pt "<<ptBin<<endl;

      std::stringstream infileNameStream;
      infileNameStream << "tmpFiles/backupWorkSpace/ws_MassLifetimeFit_Jpsi_rap" << rapBin << "_pt" << ptBin << ".root";
      const std::string infileName = infileNameStream.str().c_str();

      TFile *infile = new TFile(infileName.c_str(), "READ");
      if(!infile){
        std::cout << "Error: failed to open file with dataset" << std::endl;
      }
      RooWorkspace *ws=(RooWorkspace *)infile->Get("ws_masslifetime");
      if(!ws){
        std::cout << "Error: failed to open workspace " << std::endl;
      }

      //cout<<"infileName: "<<infileName.c_str()<<endl;

      pTmean[rapArrayIndex][ptArrayIndex] = ws->var("var_chicMeanPt")->getVal();
      pTmean_Err[rapArrayIndex][ptArrayIndex] = ws->var("var_chicMeanPt")->getError();
      if(KinParticleChi == false){
        pTmean[rapArrayIndex][ptArrayIndex] = ws->var("var_jpsiMeanPt")->getVal();
        pTmean_Err[rapArrayIndex][ptArrayIndex] = ws->var("var_jpsiMeanPt")->getError();
      }

      if(ws->var(ParName)==NULL){
        if(ws->function(ParName)==NULL) return;
        else{
          cout<<"ws->var==NULL but ws->function!=0 -> function=true"<<endl;
          functionVar=true;
        }
      }
      if(!functionVar){
        Par[rapArrayIndex][ptArrayIndex] = ws->var(ParName)->getVal();
        Par_Err[rapArrayIndex][ptArrayIndex] = ws->var(ParName)->getError();
      }
      else{
        Par[rapArrayIndex][ptArrayIndex] = ws->function(ParName)->getVal();
        Par_Err[rapArrayIndex][ptArrayIndex] = 0.;
      }

      if(WhatKindOfParIndex==0 || WhatKindOfParIndex==1){
        if(!functionVar){
          if(ws->var(ParName)->isConstant()){
            WhatKindOfParIndex=1;
            nFixedBins++;
          }
        }
        if(functionVar){
          WhatKindOfParIndex=1;
          nFixedBins++;
        }
      }

      if(changeSign) Par[rapArrayIndex][ptArrayIndex]*=-1.;

      cout<<"pTmean: "<<pTmean[rapArrayIndex][ptArrayIndex]<<" +-"<<pTmean_Err[rapArrayIndex][ptArrayIndex]<<endl;
      cout<<ParName<<": "<<Par[rapArrayIndex][ptArrayIndex]<<" +-"<<Par_Err[rapArrayIndex][ptArrayIndex]<<endl;

    }
  }

  cout<<"read all parameters"<<endl;

  cout<<"nFixedBins "<<nFixedBins<<endl;
  cout<<"ptMax-ptMin+1 "<<ptMax-ptMin+1<<endl;

  if(WhatKindOfParIndex==1 && nFixedBins<ptMax-ptMin+1) WhatKindOfParIndex=4;

  int box_Incl_FillStyle=3001;
  int box_Incl_FillColor=416-7;
  int marker_Incl_Style=20;
  int marker_Incl_Color=416+2;
  TGraphErrors *graph_Incl_Phantom = new TGraphErrors();
  graph_Incl_Phantom->SetFillStyle(box_Incl_FillStyle);
  graph_Incl_Phantom->SetFillColor(box_Incl_FillColor);
  graph_Incl_Phantom->SetLineWidth(0.);
  graph_Incl_Phantom->SetLineColor(box_Incl_FillColor);
  graph_Incl_Phantom->SetMarkerStyle(marker_Incl_Style);
  graph_Incl_Phantom->SetMarkerColor(marker_Incl_Color);

  TGraphErrors *graph_Par[RapBins], graph_Par_Incl;

  double box_yMin, box_yMax;
  box_yMin=Par_Incl-Par_Incl_Err;
  box_yMax=Par_Incl+Par_Incl_Err;
  if(box_yMin<Ymin) box_yMin=Ymin;
  if(box_yMax>Ymax) box_yMax=Ymax;

  TBox* box_Par_Incl = new TBox( onia::pTRange[0][ptMin-1], box_yMin, onia::pTRange[0][ptMax], box_yMax);
  box_Par_Incl->SetFillStyle(box_Incl_FillStyle);
  box_Par_Incl->SetFillColor(box_Incl_FillColor);
  box_Par_Incl->SetLineWidth(0.);

  TMarker* marker_Par_Incl = new TMarker(pTmean_Incl, Par_Incl, marker_Incl_Style);
  marker_Par_Incl->SetMarkerColor(marker_Incl_Color);

  for(int rapBin = rapMin; rapBin < rapMax+1; rapBin++){
    int rapArrayIndex=rapBin-rapMin;
    graph_Par[rapArrayIndex] = new TGraphErrors(PtBins, pTmean[rapArrayIndex], Par[rapArrayIndex], pTmean_Err[rapArrayIndex], Par_Err[rapArrayIndex]);
  }

  double bottomMarg=0.11;
  double leftMarg=0.15;
  double rightMarg=0.02;
  double topMarg=0.02;

  gStyle->SetPadBottomMargin(bottomMarg); //0.12
  gStyle->SetPadLeftMargin(leftMarg); //0.12
  gStyle->SetPadRightMargin(rightMarg); //0.05
  gStyle->SetPadTopMargin(topMarg); //0.05

  double blX = 0.7, trX = 1.-rightMarg-0.05;
  double blY = 0.775, trY = 1.-topMarg-0.05;
  TLegend* legend=new TLegend(blX,blY,trX,trY);
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(legendsize);
  legend->SetBorderSize(0.);
  if(AddInclusiveResult){
    legend->AddEntry(graph_Incl_Phantom,Form("Incl. result"),"fp");
  }
  for(int rapBin = rapMin; rapBin < rapMax+1; rapBin++){
    int rapArrayIndex=rapBin-rapMin;
    if(rapBin==0) legend->AddEntry(graph_Par[rapArrayIndex],Form("|y%s| < %1.1f", onia::KinParticleChar, onia::rapForPTRange[onia::kNbRapForPTBins]),"lp");
    else if(rapBin==1) legend->AddEntry(graph_Par[rapArrayIndex],Form("|y%s| < %1.1f", onia::KinParticleChar, onia::rapForPTRange[1]),"lp");
    else if(rapBin>1) legend->AddEntry(graph_Par[rapArrayIndex],Form("%1.1f < |y%s| < %1.1f", onia::rapForPTRange[rapBin-1], onia::KinParticleChar, onia::rapForPTRange[rapBin]),"lp");
    //if(rapBin==0) legend->AddEntry(graph_Par[rapArrayIndex],Form("no Punzi terms"),"lp");
    //else if(rapBin==1) legend->AddEntry(graph_Par[rapArrayIndex],Form("with Punzi terms"),"lp");
  }


  TCanvas *c1=new TCanvas("c1","");
  c1->SetTickx();
  c1->SetTicky();
  //c1->SetGridx();
  //c1->SetGridy();
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);
  gStyle->SetTitleFont(22);
  gStyle->SetStatFont(22);
  gStyle->SetStatColor(10);
  gStyle->SetStatBorderSize(1);
  gStyle->SetLabelFont(22,"X");
  gStyle->SetLabelFont(22,"Y");
  gStyle->SetTitleXOffset(1.2);
  gStyle->SetTitleYOffset(0.825);
  gStyle->SetTitleYSize(0.08);
  gStyle->SetTitleXSize(0.04);
  gStyle->SetHistLineWidth(2);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  gStyle->SetTitleX(0.15);
  gStyle->SetTitleY(0.96);


  for(int rapBin = rapMin; rapBin < rapMax+1; rapBin++){
    int rapArrayIndex=rapBin-rapMin;
    graph_Par[rapArrayIndex]->SetTitle("");
    graph_Par[rapArrayIndex]->GetXaxis()->SetTitle(Form("p%s_{T} (GeV)",onia::KinParticleChar));
    graph_Par[rapArrayIndex]->GetYaxis()->SetTitle(ParTitle);
    graph_Par[rapArrayIndex]->GetXaxis()->SetLimits(Xmin, Xmax);
    graph_Par[rapArrayIndex]->GetYaxis()->SetRangeUser(Ymin, Ymax);
    graph_Par[rapArrayIndex]->SetMarkerStyle(onia::marker_rapForPTBins[rapArrayIndex]);
    graph_Par[rapArrayIndex]->SetMarkerColor(onia::colour_rapForPTBins[rapArrayIndex+2]);
    graph_Par[rapArrayIndex]->SetLineColor(onia::colour_rapForPTBins[rapArrayIndex+2]);
  }

  //double left=0.43, top=0.8+0.05, textSize=0.055;
  //if(nState==5) left=0.41;
  double left=0.45, top=0.87, textSize=0.055;
  if(nState==5) left=0.43;
  TLatex *latex=new TLatex();
  latex->SetTextFont(42);
  latex->SetNDC(kTRUE);
  latex->SetTextSize(textSize);
  double step=textSize*1.3;
  ///////////////////

  graph_Par[0]->Draw("AP");

  if(AddInclusiveResult){
    box_Par_Incl->Draw( "same" );
    marker_Par_Incl->Draw("same");
  }

  for(int rapBin = rapMin; rapBin < rapMax+1; rapBin++){
    int rapArrayIndex=rapBin-rapMin;
    graph_Par[rapArrayIndex]->Draw("Psame");
  }



  char WhatKindOfPar[200];
  if(WhatKindOfParIndex==0) sprintf(WhatKindOfPar,"Floating par.");
  if(WhatKindOfParIndex==1) sprintf(WhatKindOfPar,"Fixed par.");
  if(WhatKindOfParIndex==2) sprintf(WhatKindOfPar,"Constrained par.");
  if(WhatKindOfParIndex==3) sprintf(WhatKindOfPar,"Derived par.");
  if(WhatKindOfParIndex==4) sprintf(WhatKindOfPar,"Part. fixed par.");

  left=0.65; top=2*bottomMarg;
  latex->SetTextColor(kRed);
  if(WhatKindOfParIndex<5) latex->DrawLatex(left,top,Form(WhatKindOfPar));




  legend->Draw();

  if(logY) c1->SetLogy(true);
  else c1->SetLogy(false);

  c1->SaveAs(Form("%s/Par_%s.pdf",savePath.str().c_str(),SaveName));


  ///
  TFile *outfile  = new TFile(Form("%s/Par_chi.root",savePath.str().c_str()),"UPDATE");
  for(int rapBin = rapMin; rapBin < rapMax+1; rapBin++){
    int rapArrayIndex=rapBin-rapMin;
    graph_Par[rapArrayIndex]->SetName(Form("graph_%s_rap%d",ParName,rapBin)); graph_Par[rapArrayIndex]->Write();
    cout<<"added TGraph "<<Form("graph_%s_rap%d",ParName,rapBin)<<" to outputfile "<<Form("%s/Par_chi.root",savePath.str().c_str())<<endl;
  }
  outfile->Close();

  return;
}


void PlotJpsiRegionFracs(int  nState, int rapBin, int ptMin, int ptMax, double Ymin, double Ymax, char Folder[200], char SaveName[200], bool logY, int RegionCode){

  const int nContributions = 4;
  const int PtBins = ptMax-ptMin+1;

  std::stringstream savePath;
  savePath << "Fit/parameter";
  gSystem->mkdir(savePath.str().c_str(),kTRUE);
  savePath << "/" << Folder;
  gSystem->mkdir(savePath.str().c_str(),kTRUE);

  cout<<"Directory: "<<savePath.str().c_str()<<endl;

  double Xmin = 0.,  Xmax = 0.;
  Xmin = 5.001;    Xmax = 54.999;
  //Xmin = 5.001;    Xmax = 84.999;

  char RegName[200];
  if(RegionCode==1) sprintf(RegName,"LSB");
  if(RegionCode==2) sprintf(RegName,"SR");
  if(RegionCode==3) sprintf(RegName,"RSB");

  char ParName[200];






  double pTmean[PtBins];
  double pTmean_Err[PtBins];

  double Frac_BG[PtBins];
  double Frac_BG_Err[PtBins];
  double Frac_Jpsi[PtBins];
  double Frac_Jpsi_Err[PtBins];
  double Frac_PRJpsi[PtBins];
  double Frac_PRJpsi_Err[PtBins];
  double Frac_NPJpsi[PtBins];
  double Frac_NPJpsi_Err[PtBins];


  TFile *inFile;
  RooWorkspace *ws;
  RooDataSet *data;
  char inName[200];



  for(int ptBin = ptMin; ptBin < ptMax+1; ptBin++){
    int ptArrayIndex=ptBin-ptMin;

    cout<<"   pt "<<ptBin<<endl;

    std::stringstream infileNameStream;
    infileNameStream << "tmpFiles/backupWorkSpace/ws_MassLifetimeFit_Jpsi_rap" << rapBin << "_pt" << ptBin << ".root";
    const std::string infileName = infileNameStream.str().c_str();

    TFile *infile = new TFile(infileName.c_str(), "READ");
    if(!infile){
      std::cout << "Error: failed to open file with dataset" << std::endl;
    }
    RooWorkspace *ws=(RooWorkspace *)infile->Get("ws_masslifetime");
    if(!ws){
      std::cout << "Error: failed to open workspace " << std::endl;
    }

    //cout<<"infileName: "<<infileName.c_str()<<endl;

    pTmean[ptArrayIndex] = ws->var("var_chicMeanPt")->getVal();
    pTmean_Err[ptArrayIndex] = ws->var("var_chicMeanPt")->getError();
    if(KinParticleChi == false){
      pTmean[ptArrayIndex] = ws->var("var_jpsiMeanPt")->getVal();
      pTmean_Err[ptArrayIndex] = ws->var("var_jpsiMeanPt")->getError();
    }

    sprintf(ParName,"var_frac_jpsi_BGIn%s",RegName);
    Frac_BG[ptArrayIndex] = ws->var(ParName)->getVal();
    Frac_BG_Err[ptArrayIndex] = ws->var(ParName)->getError();
    sprintf(ParName,"var_frac_jpsi_SigIn%s",RegName);
    Frac_Jpsi[ptArrayIndex] = ws->var(ParName)->getVal();
    Frac_Jpsi_Err[ptArrayIndex] = ws->var(ParName)->getError();

    if(RegionCode==2) sprintf(ParName,"jpsi_fPrompt");
    else sprintf(ParName,"jpsi_fPrompt%s",RegName);
    Frac_PRJpsi[ptArrayIndex] = ws->function(ParName)->getVal();
    Frac_PRJpsi_Err[ptArrayIndex] = 0.;
    if(RegionCode==2) sprintf(ParName,"jpsi_fNonPrompt");
    else sprintf(ParName,"jpsi_fNonPrompt%s",RegName);
    Frac_NPJpsi[ptArrayIndex] = ws->function(ParName)->getVal();
    Frac_NPJpsi_Err[ptArrayIndex] = 0.;

  }


  cout<<"read all parameters"<<endl;



  TGraphErrors *graph_Frac_BG, *graph_Frac_Jpsi, *graph_Frac_PRJpsi, *graph_Frac_NPJpsi;

  graph_Frac_BG = new TGraphErrors(PtBins, pTmean, Frac_BG, pTmean_Err, Frac_BG_Err);
  graph_Frac_Jpsi = new TGraphErrors(PtBins, pTmean, Frac_Jpsi, pTmean_Err, Frac_Jpsi_Err);
  graph_Frac_PRJpsi = new TGraphErrors(PtBins, pTmean, Frac_PRJpsi, pTmean_Err, Frac_PRJpsi_Err);
  graph_Frac_NPJpsi = new TGraphErrors(PtBins, pTmean, Frac_NPJpsi, pTmean_Err, Frac_NPJpsi_Err);

  double markersize[nContributions]={1., 1., 1., 1.};
  int markerstyle[nContributions]={25, 20, 24, 24};
  int markercolor[nContributions]={onia::ColorMuMuBG, onia::ColorSumJpsiSignal, onia::ColorPRJpsi, onia::ColorNPJpsi};


  //if(RegionCode==1 || RegionCode==3 ) markersize[0]*=1.5;
  //if(RegionCode==2) markersize[1]*=1.5;


  graph_Frac_BG->SetMarkerStyle(markerstyle[0]);
  graph_Frac_BG->SetMarkerSize(markersize[0]);
  graph_Frac_BG->SetMarkerColor(markercolor[0]);
  graph_Frac_BG->SetLineColor(markercolor[0]);

  graph_Frac_Jpsi->SetMarkerStyle(markerstyle[1]);
  graph_Frac_Jpsi->SetMarkerSize(markersize[1]);
  graph_Frac_Jpsi->SetMarkerColor(markercolor[1]);
  graph_Frac_Jpsi->SetLineColor(markercolor[1]);

  graph_Frac_PRJpsi->SetMarkerStyle(markerstyle[2]);
  graph_Frac_PRJpsi->SetMarkerSize(markersize[2]);
  graph_Frac_PRJpsi->SetMarkerColor(markercolor[2]);
  graph_Frac_PRJpsi->SetLineColor(markercolor[2]);

  graph_Frac_NPJpsi->SetMarkerStyle(markerstyle[3]);
  graph_Frac_NPJpsi->SetMarkerSize(markersize[3]);
  graph_Frac_NPJpsi->SetMarkerColor(markercolor[3]);
  graph_Frac_NPJpsi->SetLineColor(markercolor[3]);


  double bottomMarg=0.11;
  double leftMarg=0.15;
  double rightMarg=0.02;
  double topMarg=0.02;

  gStyle->SetPadBottomMargin(bottomMarg); //0.12
  gStyle->SetPadLeftMargin(leftMarg); //0.12
  gStyle->SetPadRightMargin(rightMarg); //0.05
  gStyle->SetPadTopMargin(topMarg); //0.05

  double blX = 0.775, trX = 1.-rightMarg-0.05;
  double blY = 0.5, trY = 1.-topMarg-0.05;
  TLegend* legend=new TLegend(blX,blY,trX,trY);
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(legendsize);
  legend->SetBorderSize(0.);

  legend->AddEntry(graph_Frac_BG,"Bkg","lp");
  legend->AddEntry(graph_Frac_Jpsi,"PR+NP J/#psi","lp");
  legend->AddEntry(graph_Frac_PRJpsi,"PR J/#psi","lp");
  legend->AddEntry(graph_Frac_NPJpsi,"NP J/#psi","lp");



  TCanvas *c1=new TCanvas("c1","");
  c1->SetTickx();
  c1->SetTicky();
  //c1->SetGridx();
  //c1->SetGridy();
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);
  gStyle->SetTitleFont(22);
  gStyle->SetStatFont(22);
  gStyle->SetStatColor(10);
  gStyle->SetStatBorderSize(1);
  gStyle->SetLabelFont(22,"X");
  gStyle->SetLabelFont(22,"Y");
  gStyle->SetTitleXOffset(1.2);
  gStyle->SetTitleYOffset(1.8);
  gStyle->SetTitleYSize(0.04);
  gStyle->SetTitleXSize(0.04);
  gStyle->SetHistLineWidth(2);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  gStyle->SetTitleX(0.15);
  gStyle->SetTitleY(0.96);


  char ParTitle[200];
  sprintf(ParTitle,"Fractions in %s",RegName);

  graph_Frac_BG->SetTitle("");
  graph_Frac_BG->GetXaxis()->SetTitle(Form("p%s_{T} (GeV)",onia::KinParticleChar));
  graph_Frac_BG->GetYaxis()->SetTitle(ParTitle);
  graph_Frac_BG->GetXaxis()->SetLimits(Xmin, Xmax);
  graph_Frac_BG->GetYaxis()->SetRangeUser(Ymin, Ymax);

  double left=0.45, top=0.87, textSize=0.055;
  if(nState==5) left=0.43;
  TLatex *latex=new TLatex();
  latex->SetTextFont(42);
  latex->SetNDC(kTRUE);
  latex->SetTextSize(textSize);
  double step=textSize*1.3;
  ///////////////////

  graph_Frac_BG->Draw("AP");
  graph_Frac_Jpsi->Draw("Psame");
  graph_Frac_PRJpsi->Draw("Psame");
  graph_Frac_NPJpsi->Draw("Psame");


  legend->Draw();

  if(logY) c1->SetLogy(true);
  else c1->SetLogy(false);

  c1->SaveAs(Form("%s/Fractions_In_%s_rap%d.pdf",savePath.str().c_str(),RegName,rapBin));


  ///
  TFile *outfile  = new TFile(Form("%s/Par_chi.root",savePath.str().c_str()),"UPDATE");

  graph_Frac_BG->SetName(Form("graph_Frac_BG_in_%s",RegName)); graph_Frac_BG->Write();
  graph_Frac_Jpsi->SetName(Form("graph_Frac_Jpsi_in_%s",RegName)); graph_Frac_Jpsi->Write();

  outfile->Close();

  return;

}
