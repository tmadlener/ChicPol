#ifndef COMMON_VAR_H__
#define COMMON_VAR_H__

#include "TLorentzVector.h"
#include "TMath.h"

namespace onia{

  // beam energy in GeV
  // const double pbeam = 3500.; // 2011
  const double pbeam = 4000.; // 2012
  // masses
  const double Mprot = 0.9382720;
  const double muMass = 0.105658;
  const double MpsiPDG = 3.096916;
  const double Mups1SPDG = 9.46030;
  const double Mchi0PDG = 3.41475;
  const double Mchi1PDG = 3.51066;
  const double Mchi2PDG = 3.55620;
  const double Ebeam = sqrt( pbeam*pbeam + Mprot*Mprot );
  const TLorentzVector beam1_LAB( 0., 0., pbeam, Ebeam );
  const TLorentzVector beam2_LAB( 0., 0., -pbeam, Ebeam );

  // dimuon mass ranges
  const double massMin = 2.85;
  const double massMax = 3.3;
  const double nSigMass = 3.0;
  const double nSigBkgLow = 4.0;
  const double nSigBkgHigh = 3.5;

  //Measure polarization as function of which particle-kinemarics
  bool KinParticleChi = true;
  const char *KinParticleChar = "^{#chi}";
  bool KinParticleChiButJpsiRap = false;
  //bool KinParticleChi = false;
  //const char *KinParticleChar = "^{#psi}";

  //chi mass ranges
  const double chimassMin = 3.325;
  const double chimassMax = 3.725;
  //When running on 2011 Jpsi data:
  //const double chimassMin = 2.;
  //const double chimassMax = 4.;

  const double massChiSBMin = chimassMin;
  const double massChiSBMax = chimassMax;


  //const double nSigChi1MassLow = 8.0;
  //const double nSigChi1MassHigh = 3.0;
  //const double nSigChi2MassLow = 6.0;
  //const double nSigChi2MassHigh = 3.0;
  //const double nSigChiBkgLow = 10.0;
  //const double nSigChiBkgHigh = 3.5;

  const double fChi1MassLow = 0.1;
  const double fChi1MassHigh = 0.01;
  const double fChi2MassLow = 0.1;
  const double fChi2MassHigh = 0.01;
  const double fChiBkgLow = 0.05;
  const double fChiBkgHigh = 0.001;

  //lifetime ranges
  const double ctVarMin = -0.29375;
  const double ctVarMax = 1.00625;
  const double ctPlotMin = ctVarMin;
  const double ctPlotMax = ctVarMax;
  const double ctPlotMinZoom = -0.1;
  const double ctPlotMaxZoom = 0.25;
  const double nSigLifetimePR = 3.0;
  const double nSigLifetimeNP = 3.0;


  // Binning
  //const int kNbRapForPTBins = 2;
  //double rapForPTRange[kNbRapForPTBins+1] = {0., 1.2, 1.5};
  //double rapRange[2*kNbRapForPTBins+1] = {-1.5, -1.2, 0., 1.2, 1.5};
  const int kNbRapForPTBins = 1;
  double rapForPTRange[kNbRapForPTBins+1] = {0., 1.2};
  double rapRange[2*kNbRapForPTBins+1] = {-1.2, 0., 1.2};

#ifndef USE_CHIC_BINNING
#warning "USE_CHIC_BINNING is not defined. Setting it to 1"
#define USE_CHIC_BINNING 1
#endif
  //chic
#if USE_CHIC_BINNING == 1
  const int kNbPTMaxBins = 5;
#endif
#if USE_CHIC_BINNING == 2
  const int kNbPTMaxBins = 4; // chic2 binning
#endif
  const int kNbPTBins[kNbRapForPTBins+1] = {kNbPTMaxBins, kNbPTMaxBins};//all y, y1
  double pTRange[kNbRapForPTBins+1][kNbPTMaxBins+1] = {
#if USE_CHIC_BINNING == 1
    // chic1 binning (standard)
    {10., 15., 20., 25., 30., 50.},//all rapidities
    {10., 15., 20., 25., 30., 50.}};//forward rapidity
#endif
#if USE_CHIC_2_BINNING == 2
    // chic2 binning (last two pt bins merged)
    {10., 15., 20., 25., 50.}, // all rapidities
    {10., 15., 20., 25., 50.}}; // forward rapidity
#endif
  //chic1
  //{10., 14., 18., 22., 30., 50.},//all rapidities
  //{10., 14., 18., 22., 30., 50.}};
  //chic2
  //{10., 14., 20., 25., 50.},//all rapidities
  //{10., 14., 20., 25., 50.}};

  //// Binning
  //const int kNbRapForPTBins = 15;
  //Double_t rapForPTRange[kNbRapForPTBins+1] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5}; //2 September 2011
  ////study the negative and positive rapidity sides separately
  //Double_t rapRange[2*kNbRapForPTBins+1] = {-1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5};
  //
  ////Jpsi
  //const int kNbPTMaxBins = 6;
  //const int kNbPTBins[kNbRapForPTBins+1] = {kNbPTMaxBins, kNbPTMaxBins, kNbPTMaxBins, kNbPTMaxBins, kNbPTMaxBins, kNbPTMaxBins, kNbPTMaxBins, kNbPTMaxBins, kNbPTMaxBins, kNbPTMaxBins, kNbPTMaxBins, kNbPTMaxBins, kNbPTMaxBins, kNbPTMaxBins, kNbPTMaxBins};//all y, y1
  //double pTRange[kNbRapForPTBins+1][kNbPTMaxBins+1] = {
  //  {10., 15., 20., 25., 30., 50., 100.},//all rapidities
  //  {10., 15., 20., 25., 30., 50., 100.},//mid-rapidity
  //  {10., 15., 20., 25., 30., 50., 100.},//mid-rapidity
  //  {10., 15., 20., 25., 30., 50., 100.},//mid-rapidity
  //  {10., 15., 20., 25., 30., 50., 100.},//mid-rapidity
  //  {10., 15., 20., 25., 30., 50., 100.},//mid-rapidity
  //  {10., 15., 20., 25., 30., 50., 100.},//mid-rapidity
  //  {10., 15., 20., 25., 30., 50., 100.},//mid-rapidity
  //  {10., 15., 20., 25., 30., 50., 100.},//mid-rapidity
  //  {10., 15., 20., 25., 30., 50., 100.},//mid-rapidity
  //  {10., 15., 20., 25., 30., 50., 100.},//mid-rapidity
  //  {10., 15., 20., 25., 30., 50., 100.},//mid-rapidity
  //  {10., 15., 20., 25., 30., 50., 100.},//mid-rapidity
  //  {10., 15., 20., 25., 30., 50., 100.},//mid-rapidity
  //  {10., 15., 20., 25., 30., 50., 100.}};//forward rapidity

  //number of reference frames
  const int kNbFrames = 6;
  const char *frameLabel[kNbFrames] = {"CS", "HX", "PHX", "sGJ", "GJ1", "GJ2"};
  enum {CS, HX, PHX, sGJ, GJ1, GJ2};


  //polarization variables
  const int kNbPolVar = 2; //cosTheta, phi
  enum {cosThPol,phiPol};
  //cosTheta
  const int kNbBinsCosT = 16;
  double cosTMin = -1., cosTMax = 1.;
  //phi for pol.
  const int kNbBinsPhiPol = 16;
  double phiPolMin = -180., phiPolMax = 180.;



  //some make up to use the same colour and marker for each pT and rapidity bin
  //in every plotting macro:
  int colour_pT[] = {1, 2, 3, 4, 6, 7, 8, 49, 38, 46, 12, 40, 50};
  int marker_pT[] = {20, 21, 25, 22, 23, 26, 27, 28, 29, 30, 20, 20, 50};
  int colour_rapForPTBins[] = {1, 30, 4, 2, 3, kMagenta+1};
  int marker_rapForPTBins[] = {20, 25, 21, 20, 22, 29};

  //Chic plots:
  int colorBackground=1;
  int colorChic0=901;
  int colorChic1=417;
  int colorChic2=632;
  int colorPR=616;
  int colorNP=600;
  //Jpsi plots:
  int colorJpsi=633;
  int ColorMuMuBG=903;
  int ColorSumJpsiSignal=419;
  int ColorPRJpsi=600;
  int ColorNPJpsi=632;


  double lineWidth_ML=2.;
  double markerSize_ML=0.8;
  double enlargeYby_ML=1.2;
  int lineStyle_subComps_ML=7;
  int lineStyle_subCompPRsignal_ML=3;
  double lineWidth_PRsignal_ML=1.;

  int LifetimePlotBins=104;
  int ChicMassPlotBins=80;

  const int nRegions=4;
  int colorRegions[nRegions+1]={600,920+1,colorChic1,colorChic2,920+3};


  // Dimuon cuts

  double cut_vtxProb=0.01;
  double rap = 20.; //dimuon rap

  // Gamma cuts

  double cut_gammapt = 0.;
  double cut_gammaeta=15.;
  double cut_RconvMin = 1.5;
  double cut_RconvMax = 200;

  // Chi cuts

  double cut_dz=1.35;//1.
  double cut_probFit = 0.; //0.01;
  double chirap = rapForPTRange[kNbRapForPTBins]; //chi rap

}

#endif
