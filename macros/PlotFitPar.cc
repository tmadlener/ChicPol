#include <iostream>
#include <string>
#include <sstream>
#include "calculatePar.cc"
#include "TBox.h"
#include "TMarker.h"
using namespace std;
using namespace RooFit;
using namespace onia;

void PlotChiFitPar(int  nState=4, int rapMin=0, int rapMax=0, int ptMin=0, int ptMax=0, int rapFixTo=0, int ptFixTo=0, bool AddInclusiveResult=false, double Ymin=0., double Ymax=0., char ParName[200]="default", char ParTitle[200]="default", char Folder[200]="default", char SaveName[200]="default", bool logY=false, int WhatKindOfParIndex=5);
void PlotChiRegionFracs(int  nState=4, int rapBin=0, int ptMin=0, int ptMax=0, double Ymin=0., double Ymax=0., char Folder[200]="default", char SaveName[200]="default", bool logY=false, int RegionCode=999);

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
    sprintf(dirname,"%s/chic_mass",predirname);
    gSystem->mkdir(dirname);
	outfileCreate = new TFile(Form("%s/Par_chi.root",dirname,"RECREATE"));
	outfileCreate->Close();
    sprintf(dirname,"%s/chic_lifetime",predirname);
    gSystem->mkdir(dirname);
	outfileCreate = new TFile(Form("%s/Par_chi.root",dirname,"RECREATE"));
	outfileCreate->Close();
    sprintf(dirname,"%s/chic_fitqual",predirname);
    gSystem->mkdir(dirname);
	outfileCreate = new TFile(Form("%s/Par_chi.root",dirname,"RECREATE"));
	outfileCreate->Close();
    sprintf(dirname,"%s/chic_regions",predirname);
    gSystem->mkdir(dirname);
	outfileCreate = new TFile(Form("%s/Par_chi.root",dirname,"RECREATE"));
	outfileCreate->Close();
    sprintf(dirname,"%s/chic_fractions",predirname);
    gSystem->mkdir(dirname);
	outfileCreate = new TFile(Form("%s/Par_chi.root",dirname,"RECREATE"));
	outfileCreate->Close();
    sprintf(dirname,"%s/chic_additionals",predirname);
    gSystem->mkdir(dirname);
	outfileCreate = new TFile(Form("%s/Par_chi.root",dirname,"RECREATE"));
	outfileCreate->Close();


	double Ymin, Ymax;
	char ParName[200];
	char ParTitle[200];
	char SaveName[200];
	char Folder[200];
	bool logY=false;

	int WhatKindOfParIndex=999;

	sprintf(Folder, "chic_fractions");  Ymin = 0.; Ymax = 1.; logY=false; int rapBin=1;
	sprintf(SaveName, "default");

	for(int rapBin = rapMin; rapBin < rapMax+1; rapBin++){
		for(int RegionCode=1;RegionCode<9;RegionCode++)
			PlotChiRegionFracs(nState, rapBin, ptMin, ptMax, Ymin, Ymax, Folder, SaveName, logY, RegionCode);
	}

	WhatKindOfParIndex=500; sprintf(SaveName, "bfrac_Jpsi"); sprintf(ParName, "fBkgNP"); sprintf(ParTitle,"J/#psi B-fraction"); sprintf(Folder, "chic_additionals");  Ymin = 0.; Ymax = 1.; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=500; sprintf(SaveName, "bfrac_Jpsi_incl"); sprintf(ParName, "jpsi_bfraction"); sprintf(ParTitle,"J/#psi B-fraction"); sprintf(Folder, "chic_additionals");  Ymin = 0.; Ymax = 1.; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=501; sprintf(SaveName, "bfrac_Chic1"); sprintf(ParName, "fracNP_chic1"); sprintf(ParTitle,"#chi_{c1} B-fraction"); sprintf(Folder, "chic_additionals");  Ymin = 0.; Ymax = 1.; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=502; sprintf(SaveName, "bfrac_Chic2"); sprintf(ParName, "fracNP_chic2"); sprintf(ParTitle,"#chi_{c2} B-fraction"); sprintf(Folder, "chic_additionals");  Ymin = 0.; Ymax = 1.; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);

	WhatKindOfParIndex=3; sprintf(SaveName, "Sig_CBmass0"); sprintf(ParName, "CBmass0"); sprintf(ParTitle,"m_{#chi_{c0}} [GeV]"); sprintf(Folder, "chic_mass");  Ymin = 3.375; Ymax = 3.425; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=0; sprintf(SaveName, "Sig_CBmass1"); sprintf(ParName, "CBmass1"); sprintf(ParTitle,"m_{#chi_{c1}} [GeV]"); sprintf(Folder, "chic_mass");  Ymin = 3.5025; Ymax = 3.5075; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=3; sprintf(SaveName, "Sig_CBmass2"); sprintf(ParName, "CBmass2"); sprintf(ParTitle,"m_{#chi_{c2}} [GeV]"); sprintf(Folder, "chic_mass");  Ymin = 3.5425; Ymax = 3.5575; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=3; sprintf(SaveName, "Sig_CBsigma0"); sprintf(ParName, "CBsigma0"); sprintf(ParTitle,"#sigma_{#chi_{c0}} [GeV]"); sprintf(Folder, "chic_mass");  Ymin = 0.004; Ymax = 0.012; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=0; sprintf(SaveName, "Sig_CBsigma1"); sprintf(ParName, "CBsigma1"); sprintf(ParTitle,"#sigma_{#chi_{c1}} [GeV]"); sprintf(Folder, "chic_mass");  Ymin = 0.004; Ymax = 0.012; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=3; sprintf(SaveName, "Sig_CBsigma2"); sprintf(ParName, "CBsigma2"); sprintf(ParTitle,"#sigma_{#chi_{c2}} [GeV]"); sprintf(Folder, "chic_mass");  Ymin = 0.004; Ymax = 0.012; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=0; sprintf(SaveName, "BK_p1"); sprintf(ParName, "BK_p1"); sprintf(ParTitle,"p_{1}^{BG}"); sprintf(Folder, "chic_mass");  Ymin = -0.3; Ymax = -0.2; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=0; sprintf(SaveName, "BK_p2"); sprintf(ParName, "BK_p2"); sprintf(ParTitle,"p_{2}^{BG}"); sprintf(Folder, "chic_mass");  Ymin = -1e-2; Ymax = 1e-1;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=0; sprintf(SaveName, "Sig_CBalpha1"); sprintf(ParName, "CBalpha1"); sprintf(ParTitle,"#alpha^{CB}_{#chi_{c1}}"); sprintf(Folder, "chic_mass");  Ymin = 0.3; Ymax = 1.3; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=0; sprintf(SaveName, "Sig_CBalpha2"); sprintf(ParName, "CBalpha2"); sprintf(ParTitle,"#alpha^{CB}_{#chi_{c2}}"); sprintf(Folder, "chic_mass");  Ymin = 0.3; Ymax = 1.3; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=0; sprintf(SaveName, "Sig_CBn"); sprintf(ParName, "CBn"); sprintf(ParTitle,"n^{CB}_{#chi_{c1}}"); sprintf(Folder, "chic_mass");  Ymin = 1.; Ymax = 5.; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=0; sprintf(SaveName, "Sig_CBn2"); sprintf(ParName, "CBn2"); sprintf(ParTitle,"n^{CB}_{#chi_{c2}}"); sprintf(Folder, "chic_mass");  Ymin = 1.; Ymax = 5.; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=0; sprintf(SaveName, "fracBackground"); sprintf(ParName, "fracBackground"); sprintf(ParTitle,"f^{#chi}_{#psiBG}"); sprintf(Folder, "chic_mass");  Ymin = 0.6; Ymax = 0.8; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=0; sprintf(SaveName, "fracSignal_chic0"); sprintf(ParName, "fracSignal_chic0"); sprintf(ParTitle,"f^{Sig}_{#chi_{c0}}"); sprintf(Folder, "chic_mass");  Ymin = 0.; Ymax = 0.06; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=0; sprintf(SaveName, "fracSignal_chic1"); sprintf(ParName, "fracSignal_chic1"); sprintf(ParTitle,"f^{Sig}_{#chi_{c1}}"); sprintf(Folder, "chic_mass");  Ymin = 0.6; Ymax = 0.85; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=3; sprintf(SaveName, "fracSignal_chic2"); sprintf(ParName, "fracSignal_chic2"); sprintf(ParTitle,"f^{Sig}_{#chi_{c2}}"); sprintf(Folder, "chic_mass");  Ymin = 0.25; Ymax = 0.35; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	//WhatKindOfParIndex=0; sprintf(SaveName, "NumEvE"); sprintf(ParName, "NumEvE"); sprintf(ParTitle,"n^{Tot}_{fit}"); sprintf(Folder, "chic_mass");  Ymin = 1.; Ymax = 3e5; logY=true;
	//PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=3; sprintf(SaveName, "Sig_PES"); sprintf(ParName, "PES"); sprintf(ParTitle,"photon energy scale"); sprintf(Folder, "chic_mass");  Ymin = 0.975; Ymax = 1.; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);



	WhatKindOfParIndex=0; sprintf(SaveName, "BK_NP_TauBkg"); sprintf(ParName, "NP_TauBkg"); sprintf(ParTitle,"#tau_{#psiBG} [mm]"); sprintf(Folder, "chic_lifetime");  Ymin = 0.25; Ymax = 0.5; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=0; sprintf(SaveName, "NP_TauChic"); sprintf(ParName, "NP_TauChic"); sprintf(ParTitle,"#tau^{NP}_{#chi_{cJ}} [mm]"); sprintf(Folder, "chic_lifetime");  Ymin = 0.25; Ymax = 0.5; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=0; sprintf(SaveName, "Sig_ctResolution"); sprintf(ParName, "ctResolution"); sprintf(ParTitle,"#sigma^{scale}_{l}"); sprintf(Folder, "chic_lifetime");  Ymin = 0.7; Ymax = 1.1; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=0; sprintf(SaveName, "Sig_ctResolution2"); sprintf(ParName, "ctResolution2"); sprintf(ParTitle,"#sigma^{scale2}_{l}"); sprintf(Folder, "chic_lifetime");  Ymin = 0.8; Ymax = 2.5; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=0; sprintf(SaveName, "Sig_fracGauss2"); sprintf(ParName, "fracGauss2"); sprintf(ParTitle,"f_{G_{2}}"); sprintf(Folder, "chic_lifetime");  Ymin = 0.; Ymax = 1.; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=0; sprintf(SaveName, "BK_fBkgNP"); sprintf(ParName, "fBkgNP"); sprintf(ParTitle,"f^{#psiBG}_{NP}"); sprintf(Folder, "chic_lifetime");  Ymin = 0.4; Ymax = 0.7; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=0; sprintf(SaveName, "fracNP_chic0"); sprintf(ParName, "fracNP_chic0"); sprintf(ParTitle,"f^{#chi_{c0}}_{NP}"); sprintf(Folder, "chic_lifetime");  Ymin = 0.; Ymax = 1.; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=0; sprintf(SaveName, "fracNP_chic1"); sprintf(ParName, "fracNP_chic1"); sprintf(ParTitle,"f^{#chi_{c1}}_{NP}"); sprintf(Folder, "chic_lifetime");  Ymin = 0.15; Ymax = 0.35; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=0; sprintf(SaveName, "fracNP_chic2"); sprintf(ParName, "fracNP_chic2"); sprintf(ParTitle,"f^{#chi_{c2}}_{NP}"); sprintf(Folder, "chic_lifetime");  Ymin = 0.; Ymax = 0.2; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);

	WhatKindOfParIndex=4; sprintf(SaveName, "covQualMigrad"); sprintf(ParName, "var_covQualMigrad"); sprintf(ParTitle,"Migrad covQual"); sprintf(Folder, "chic_fitqual");  Ymin = 0.; Ymax = 40.; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=4; sprintf(SaveName, "covQualHesse"); sprintf(ParName, "var_covQualHesse"); sprintf(ParTitle,"Hesse covQual"); sprintf(Folder, "chic_fitqual");  Ymin = 0.; Ymax = 40.; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);


	WhatKindOfParIndex=4; sprintf(SaveName, "Chi2ndf_Lifetime_RSB"); sprintf(ParName, "var_chi2ndf_Lifetime_RSB"); sprintf(ParTitle,"#chi^{2} / ndf (lifetime RSB)"); sprintf(Folder, "chic_fitqual");  Ymin = 0.; Ymax = 3.; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=4; sprintf(SaveName, "Chi2ndf_Lifetime_LSB"); sprintf(ParName, "var_chi2ndf_Lifetime_LSB"); sprintf(ParTitle,"#chi^{2} / ndf (lifetime LSB)"); sprintf(Folder, "chic_fitqual");  Ymin = 0.; Ymax = 3.; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=4; sprintf(SaveName, "Chi2ndf_Lifetime_SR2"); sprintf(ParName, "var_chi2ndf_Lifetime_SR2"); sprintf(ParTitle,"#chi^{2} / ndf (lifetime SR2)"); sprintf(Folder, "chic_fitqual");  Ymin = 0.; Ymax = 3.; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=4; sprintf(SaveName, "Chi2ndf_Lifetime_SR1"); sprintf(ParName, "var_chi2ndf_Lifetime_SR1"); sprintf(ParTitle,"#chi^{2} / ndf (lifetime SR1)"); sprintf(Folder, "chic_fitqual");  Ymin = 0.; Ymax = 3.; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=4; sprintf(SaveName, "Chi2ndf_Lifetime"); sprintf(ParName, "var_chi2ndf_Lifetime"); sprintf(ParTitle,"#chi^{2} / ndf (lifetime)"); sprintf(Folder, "chic_fitqual");  Ymin = 0.; Ymax = 3.; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=4; sprintf(SaveName, "Chi2ndf_Mass"); sprintf(ParName, "var_chi2ndf_Mass"); sprintf(ParTitle,"#chi^{2} / ndf (mass)"); sprintf(Folder, "chic_fitqual");  Ymin = 0.; Ymax = 3.; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);


	WhatKindOfParIndex=3; sprintf(SaveName, "Def_Mass_SR1Max"); sprintf(ParName, "var_sig1MaxMass"); sprintf(ParTitle,"M^{SR1}_{high} [GeV]"); sprintf(Folder, "chic_regions");  Ymin = 3.41; Ymax = 3.58; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=3; sprintf(SaveName, "Def_Mass_SR1Min"); sprintf(ParName, "var_sig1MinMass"); sprintf(ParTitle,"M^{SR1}_{low} [GeV]"); sprintf(Folder, "chic_regions");  Ymin = 3.41; Ymax = 3.58; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=3; sprintf(SaveName, "Def_Mass_SR2Max"); sprintf(ParName, "var_sig2MaxMass"); sprintf(ParTitle,"M^{SR2}_{high} [GeV]"); sprintf(Folder, "chic_regions");  Ymin = 3.41; Ymax = 3.58; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=3; sprintf(SaveName, "Def_Mass_SR2Min"); sprintf(ParName, "var_sig2MinMass"); sprintf(ParTitle,"M^{SR2}_{low} [GeV]"); sprintf(Folder, "chic_regions");  Ymin = 3.41; Ymax = 3.58; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=3; sprintf(SaveName, "Def_Mass_LSBMax"); sprintf(ParName, "var_lsbMaxMass"); sprintf(ParTitle,"M^{LSB}_{high} [GeV]"); sprintf(Folder, "chic_regions");  Ymin = 3.41; Ymax = 3.58; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=3; sprintf(SaveName, "Def_Mass_RSBMin"); sprintf(ParName, "var_rsbMinMass"); sprintf(ParTitle,"M^{RSB}_{low} [GeV]"); sprintf(Folder, "chic_regions");  Ymin = 3.41; Ymax = 3.58; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);

	WhatKindOfParIndex=3; sprintf(SaveName, "Def_Ctau_PRMin"); sprintf(ParName, "var_PRMin"); sprintf(ParTitle,"l^{PR}_{low} [mm]"); sprintf(Folder, "chic_regions");  Ymin = -0.1; Ymax = 0.1; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=3; sprintf(SaveName, "Def_Ctau_PRMax"); sprintf(ParName, "var_PRMax"); sprintf(ParTitle,"l^{PR}_{high} [mm]"); sprintf(Folder, "chic_regions");  Ymin = -0.1; Ymax = 0.1; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=3; sprintf(SaveName, "Def_Ctau_NPMin"); sprintf(ParName, "var_NPMin"); sprintf(ParTitle,"l^{NP}_{low} [mm]"); sprintf(Folder, "chic_regions");  Ymin = -0.1; Ymax = 0.1; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);

	WhatKindOfParIndex=3; sprintf(SaveName, "Def_Ctau_Resolution"); sprintf(ParName, "var_ctres"); sprintf(ParTitle,"#sigma_{l} [mm]"); sprintf(Folder, "chic_regions");  Ymin = 0.; Ymax = 0.03; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);

	WhatKindOfParIndex=3; sprintf(SaveName, "N_PRChic1_in_PRSR1"); sprintf(ParName, "var_nPRChic1InPRSR1"); sprintf(ParTitle,"N^{PRSR1}_{PR#chi_{c1}}"); sprintf(Folder, "chic_regions");  Ymin = 700; Ymax = 5e4; logY=true;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=3; sprintf(SaveName, "N_PRChic2_in_PRSR2"); sprintf(ParName, "var_nPRChic2InPRSR2"); sprintf(ParTitle,"N^{PRSR2}_{PR#chi_{c2}}"); sprintf(Folder, "chic_regions");  Ymin = 700; Ymax = 5e4; logY=true;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);

	WhatKindOfParIndex=3; sprintf(SaveName, "ChicMeanAbsRap"); sprintf(ParName, "var_chicMeanAbsRap"); sprintf(ParTitle,"mean |y^{#chi}|"); sprintf(Folder, "chic_regions");  Ymin = 0.4; Ymax = 0.7; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);

	WhatKindOfParIndex=3; sprintf(SaveName, "ratio_PR_chic2_over_chic1"); sprintf(ParName, "var_ratio_PR_chic2_over_chic1"); sprintf(ParTitle,"N(#chi^{PR}_{c2}) / N(#chi^{PR}_{c1})"); sprintf(Folder, "chic_additionals");  Ymin = 0.; Ymax = 1.; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);
	WhatKindOfParIndex=3; sprintf(SaveName, "ratio_NP_chic2_over_chic1"); sprintf(ParName, "var_ratio_NP_chic2_over_chic1"); sprintf(ParTitle,"N(#chi^{NP}_{c2}) / N(#chi^{NP}_{c1})"); sprintf(Folder, "chic_additionals");  Ymin = 0.; Ymax = 1.; logY=false;
	PlotChiFitPar(nState, rapMin, rapMax, ptMin, ptMax, rapFixTo, ptFixTo, AddInclusiveResult, Ymin, Ymax, ParName, ParTitle, Folder, SaveName, logY, WhatKindOfParIndex);





	return 0;
}


//=============================================
void PlotChiFitPar(int  nState, int rapMin, int rapMax, int ptMin, int ptMax, int rapFixTo, int ptFixTo, bool AddInclusiveResult, double Ymin, double Ymax, char ParName[200], char ParTitle[200], char Folder[200], char SaveName[200], bool logY, int WhatKindOfParIndex){
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





	bool bfracJpsi=false;
	bool bfracChic1=false;
	bool bfracChic2=false;
	bool bfrac=false;

	if(WhatKindOfParIndex==500) bfracJpsi=true;
	if(WhatKindOfParIndex==501) bfracChic1=true;
	if(WhatKindOfParIndex==502) bfracChic2=true;

	if(bfracJpsi || bfracChic1 || bfracChic2) bfrac=true;




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

	char meanPtchar[200];
	sprintf(meanPtchar, "var_chicMeanPt_JpsiSR");
	if(bfracJpsi) sprintf(meanPtchar, "var_jpsiMeanPt_JpsiSR");

	if(AddInclusiveResult){
		cout<<"rapIncl "<<rapFixTo<<"   ptIncl "<<ptFixTo<<endl;

		std::stringstream infileNameStream;
		infileNameStream << "tmpFiles/backupWorkSpace/ws_DefineRegionsAndFractions_Chi_rap" << rapFixTo << "_pt" << ptFixTo << ".root";
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

		pTmean_Incl = ws->var(meanPtchar)->getVal();
		pTmean_Incl_Err = ws->var(meanPtchar)->getError();


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
		cout<<"Inclusive pTmean: "<<pTmean_Incl<<" +-"<<pTmean_Incl_Err<<endl;
		cout<<"Inclusive "<<ParName<<": "<<Par_Incl<<" +-"<<Par_Incl_Err<<endl;

	}

	for(int rapBin = rapMin; rapBin < rapMax+1; rapBin++){
		int rapArrayIndex=rapBin-rapMin;
		for(int ptBin = ptMin; ptBin < ptMax+1; ptBin++){
			int ptArrayIndex=ptBin-ptMin;

			cout<<"rap "<<rapBin<<"   pt "<<ptBin<<endl;

			std::stringstream infileNameStream;
			infileNameStream << "tmpFiles/backupWorkSpace/ws_DefineRegionsAndFractions_Chi_rap" << rapBin << "_pt" << ptBin << ".root";
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

			pTmean[rapArrayIndex][ptArrayIndex] = ws->var(meanPtchar)->getVal();
			pTmean_Err[rapArrayIndex][ptArrayIndex] = ws->var(meanPtchar)->getError();


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

			if(WhatKindOfParIndex==0){
				if(!functionVar)
					if(ws->var(ParName)->isConstant()) WhatKindOfParIndex=1;
				if(functionVar)
					WhatKindOfParIndex=1;
			}

			cout<<"pTmean: "<<pTmean[rapArrayIndex][ptArrayIndex]<<" +-"<<pTmean_Err[rapArrayIndex][ptArrayIndex]<<endl;
			cout<<ParName<<": "<<Par[rapArrayIndex][ptArrayIndex]<<" +-"<<Par_Err[rapArrayIndex][ptArrayIndex]<<endl;

		}
	}

	cout<<"read all parameters"<<endl;


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

	if(bfrac){
		blX-=0.2;
		trX-=0.2;
	}

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
		char legendentry[200];
		if(rapBin==0) sprintf(legendentry,"|y%s| < %1.1f", onia::KinParticleChar, onia::rapForPTRange[onia::kNbRapForPTBins]);
		else if(rapBin==1) sprintf(legendentry,"|y%s| < %1.1f", onia::KinParticleChar, onia::rapForPTRange[1]);
		else if(rapBin>1) sprintf(legendentry,"%1.1f < |y%s| < %1.1f", onia::rapForPTRange[rapBin-1], onia::KinParticleChar, onia::rapForPTRange[rapBin]);
		if(bfrac){
			sprintf(legendentry,"%s, BPH-13-001, #sqrt{s} = 8TeV", legendentry, onia::rapForPTRange[rapBin-1], onia::KinParticleChar, onia::rapForPTRange[rapBin]);
		}

		if(rapBin==0) legend->AddEntry(graph_Par[rapArrayIndex],legendentry,"lp");
		else if(rapBin==1) legend->AddEntry(graph_Par[rapArrayIndex],legendentry,"lp");
		else if(rapBin>1) legend->AddEntry(graph_Par[rapArrayIndex],legendentry,"lp");
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
		if(bfracJpsi) graph_Par[rapArrayIndex]->GetXaxis()->SetTitle(Form("p^{#psi}_{T} (GeV)"));
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








	if(bfracJpsi){

		cout<<"bfracJpsi"<<endl;
		// Jpsi bfraction results //
		const int PTBins1 = 10, PTBins2 = 6;
		double pTBinRap1[PTBins1+1] = {8., 9., 10., 11., 12., 13.5, 15., 18., 30., 45., 70.}; //rap 0.0 - 0.9
		double pTBinRap2[PTBins2+1] = {8., 9., 10., 12., 15., 30., 45.}; //rap 0.9 - 1.2

		double pTBinRap1_weight[PTBins1] = { 8.569, 9.515, 10.494, 11.485, 12.711, 14.206, 16.329, 21.77, 35.03, 52.7}; //rap 0.0 - 0.9
		double pTBinRap2_weight[PTBins2] = {8.524, 9.491, 10.921, 13.312, 18.95, 34.7}; ////rap 0.9 - 1.2

		double BFracRap1[PTBins1] = { 0.271, 0.285, 0.320, 0.345, 0.373, 0.417, 0.454, 0.535, 0.633, 0.646 };
		double BFracErrRap1[PTBins1] = { 0.006, 0.005, 0.005, 0.005, 0.005, 0.006, 0.006, 0.006, 0.015, 0.038 };
		double BFracRap2[PTBins1] = { 0.265, 0.281, 0.324, 0.380, 0.481, 0.616 };
		double BFracErrRap2[PTBins1] = { 0.007, 0.007, 0.006, 0.006, 0.007, 0.029 };


		double pTRap1[PTBins1];
		double pTRap1ErrLow[PTBins1];
		double pTRap1ErrHigh[PTBins1];
		double pTRap2[PTBins2];
		double pTRap2ErrLow[PTBins2];
		double pTRap2ErrHigh[PTBins2];
		for(int pt=0; pt<PTBins1; pt++){
			pTRap1[pt] = pTBinRap1_weight[pt];
			pTRap1ErrLow[pt] = 0.;
			pTRap1ErrHigh[pt] = 0.;
		}
		for(int pt=0; pt<PTBins2; pt++){
			pTRap2[pt] = pTBinRap2_weight[pt];
			pTRap2ErrLow[pt] = 0.;
			pTRap2ErrHigh[pt] = 0.;
		}

		TGraphAsymmErrors *graph_BFracRap1 = new TGraphAsymmErrors(PTBins1,
				pTRap1, BFracRap1, pTRap1ErrLow, pTRap1ErrHigh, BFracErrRap1, BFracErrRap1);
		TGraphAsymmErrors *graph_BFracRap2 = new TGraphAsymmErrors(PTBins2,
				pTRap2, BFracRap2, pTRap2ErrLow, pTRap2ErrHigh, BFracErrRap2, BFracErrRap2);

		graph_BFracRap1->SetMarkerStyle(24);
		graph_BFracRap1->SetMarkerColor(onia::colour_rapForPTBins[5]);
		graph_BFracRap1->SetLineColor(onia::colour_rapForPTBins[5]);

		graph_BFracRap2->SetMarkerStyle(24);
		graph_BFracRap2->SetMarkerColor(onia::colour_rapForPTBins[3]);
		graph_BFracRap2->SetLineColor(onia::colour_rapForPTBins[3]);

		legend->AddEntry(graph_BFracRap1,"|y^{#psi}| < 0.9, BPH-10-014, #sqrt{s} = 7TeV","lp");
		legend->AddEntry(graph_BFracRap2,"0.9 < |y^{#psi}| < 1.2, BPH-10-014, #sqrt{s} = 7TeV","lp");

		graph_BFracRap1->Draw("Psame");
		graph_BFracRap2->Draw("Psame");

	}



	if(bfracChic1 || bfracChic2){

		if(bfracChic1) cout<<"bfracChic1"<<endl;
		if(bfracChic2) cout<<"bfracChic2"<<endl;

		//ChicJ B-fractions from ATLAS
		//ATLAS-CONF-2013-095
		//pTBins (no mean pT given), for jpsi rapidity |y|<0.75

		const int PTBinsChi = 5;
		double pTBinChi[PTBinsChi+1] = {12., 14., 16., 18., 22., 30.};
		double pTBinChiWeight[PTBinsChi] = {13., 15., 17., 20., 26.};

		double BFracChi1[PTBinsChi] = {0.22, 0.24, 0.24, 0.27, 0.34};
		double BFracChi2[PTBinsChi] = {0.10, 0.09, 0.09, 0.11, 0.12};

		double BFracStatErrChi1[PTBinsChi] = {0.02, 0.02, 0.03, 0.03, 0.03};
		double BFracSystErrChi1[PTBinsChi] = {0.01, 0.01, 0.03, 0.01, 0.01};

		double BFracStatErrChi2[PTBinsChi] = {0.04, 0.04, 0.04, 0.05, 0.06};
		double BFracSystErrChi2[PTBinsChi] = {0.01, 0.02, 0.02, 0.02, 0.03};

		double BFracErrChi1[PTBinsChi];
		double BFracErrChi2[PTBinsChi];

		double pTChi1[PTBinsChi];
		double pTChi1ErrLow[PTBinsChi];
		double pTChi1ErrHigh[PTBinsChi];
		double pTChi2[PTBinsChi];
		double pTChi2ErrLow[PTBinsChi];
		double pTChi2ErrHigh[PTBinsChi];
		for(int pt=0; pt<PTBinsChi; pt++){
			pTChi1[pt] = pTBinChiWeight[pt];
			pTChi1ErrLow[pt] = 0.;
			pTChi1ErrHigh[pt] = 0.;
			BFracErrChi1[pt] = TMath::Sqrt(BFracStatErrChi1[pt]*BFracStatErrChi1[pt]+BFracSystErrChi1[pt]*BFracSystErrChi1[pt]);
		}
		for(int pt=0; pt<PTBinsChi; pt++){
			pTChi2[pt] = pTBinChiWeight[pt];
			pTChi2ErrLow[pt] = 0.;
			pTChi2ErrHigh[pt] = 0.;
			BFracErrChi2[pt] = TMath::Sqrt(BFracStatErrChi2[pt]*BFracStatErrChi2[pt]+BFracSystErrChi2[pt]*BFracSystErrChi2[pt]);
		}

		TGraphAsymmErrors *graph_BFracChi1 = new TGraphAsymmErrors(PTBinsChi,
				pTChi1, BFracChi1, pTChi1ErrLow, pTChi1ErrHigh, BFracErrChi1, BFracErrChi1);
		TGraphAsymmErrors *graph_BFracChi2 = new TGraphAsymmErrors(PTBinsChi,
				pTChi2, BFracChi2, pTChi2ErrLow, pTChi2ErrHigh, BFracErrChi2, BFracErrChi2);


		graph_BFracChi1->SetMarkerStyle(24);
		graph_BFracChi1->SetMarkerColor(onia::colour_rapForPTBins[5]);
		graph_BFracChi1->SetLineColor(onia::colour_rapForPTBins[5]);

		graph_BFracChi2->SetMarkerStyle(24);
		graph_BFracChi2->SetMarkerColor(onia::colour_rapForPTBins[3]);
		graph_BFracChi2->SetLineColor(onia::colour_rapForPTBins[3]);

		if(bfracChic1){
			legend->AddEntry(graph_BFracChi1,"|y^{#psi}| < 0.75, ATLAS, #sqrt{s} = 7TeV","lp");
			graph_BFracChi1->Draw("Psame");
		}

		if(bfracChic2){
			legend->AddEntry(graph_BFracChi2,"|y^{#psi}| < 0.75, ATLAS, #sqrt{s} = 7TeV","lp");
			graph_BFracChi2->Draw("Psame");
		}

	}










	char WhatKindOfPar[200];
	if(WhatKindOfParIndex==0) sprintf(WhatKindOfPar,"Floating par.");
	if(WhatKindOfParIndex==1) sprintf(WhatKindOfPar,"Fixed par.");
	if(WhatKindOfParIndex==2) sprintf(WhatKindOfPar,"Constrained par.");
	if(WhatKindOfParIndex==3) sprintf(WhatKindOfPar,"Derived par.");

	left=0.65; top=2*bottomMarg;
    latex->SetTextColor(kRed);
    if(WhatKindOfParIndex<4) latex->DrawLatex(left,top,Form(WhatKindOfPar));




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


void PlotChiRegionFracs(int  nState, int rapBin, int ptMin, int ptMax, double Ymin, double Ymax, char Folder[200], char SaveName[200], bool logY, int RegionCode){

	const int nContributions = 7;
	const int PtBins = ptMax-ptMin+1;

	std::stringstream savePath;
	savePath << "Fit/parameter";
	gSystem->mkdir(savePath.str().c_str(),kTRUE);
	savePath << "/" << Folder;
	gSystem->mkdir(savePath.str().c_str(),kTRUE);

	cout<<"Directory: "<<savePath.str().c_str()<<endl;

	double Xmin = 0.,  Xmax = 0.;
	Xmin = 5.001;    Xmax = 54.999;

	char RegName[200];
	if(RegionCode==1) sprintf(RegName,"PRLSB");
	if(RegionCode==2) sprintf(RegName,"PRSR1");
	if(RegionCode==3) sprintf(RegName,"PRSR2");
	if(RegionCode==4) sprintf(RegName,"PRRSB");
	if(RegionCode==5) sprintf(RegName,"NPLSB");
	if(RegionCode==6) sprintf(RegName,"NPSR1");
	if(RegionCode==7) sprintf(RegName,"NPSR2");
	if(RegionCode==8) sprintf(RegName,"NPRSB");

	char ParName[200];






	double pTmean[PtBins];
	double pTmean_Err[PtBins];

	double Frac_BG[PtBins];
	double Frac_BG_Err[PtBins];
	double Frac_PR_Chic0[PtBins];
	double Frac_PR_Chic0_Err[PtBins];
	double Frac_PR_Chic1[PtBins];
	double Frac_PR_Chic1_Err[PtBins];
	double Frac_PR_Chic2[PtBins];
	double Frac_PR_Chic2_Err[PtBins];
	double Frac_NP_Chic0[PtBins];
	double Frac_NP_Chic0_Err[PtBins];
	double Frac_NP_Chic1[PtBins];
	double Frac_NP_Chic1_Err[PtBins];
	double Frac_NP_Chic2[PtBins];
	double Frac_NP_Chic2_Err[PtBins];


	TFile *inFile;
	RooWorkspace *ws;
	RooDataSet *data;
	char inName[200];



		for(int ptBin = ptMin; ptBin < ptMax+1; ptBin++){
			int ptArrayIndex=ptBin-ptMin;

			cout<<"   pt "<<ptBin<<endl;

			std::stringstream infileNameStream;
			infileNameStream << "tmpFiles/backupWorkSpace/ws_DefineRegionsAndFractions_Chi_rap" << rapBin << "_pt" << ptBin << ".root";
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

			char meanPtchar[200];
			sprintf(meanPtchar, "var_chicMeanPt_Chic%s",RegName);

			pTmean[ptArrayIndex] = ws->var(meanPtchar)->getVal();
			pTmean_Err[ptArrayIndex] = ws->var(meanPtchar)->getError();
			sprintf(ParName,"var_fracBackgroundIn%s",RegName);
			Frac_BG[ptArrayIndex] = ws->var(ParName)->getVal();
			Frac_BG_Err[ptArrayIndex] = ws->var(ParName)->getError();
			sprintf(ParName,"var_fracPRChic0In%s",RegName);
			Frac_PR_Chic0[ptArrayIndex] = ws->var(ParName)->getVal();
			Frac_PR_Chic0_Err[ptArrayIndex] = ws->var(ParName)->getError();
			sprintf(ParName,"var_fracNPChic0In%s",RegName);
			Frac_NP_Chic0[ptArrayIndex] = ws->var(ParName)->getVal();
			Frac_NP_Chic0_Err[ptArrayIndex] = ws->var(ParName)->getError();
			sprintf(ParName,"var_fracPRChic1In%s",RegName);
			Frac_PR_Chic1[ptArrayIndex] = ws->var(ParName)->getVal();
			Frac_PR_Chic1_Err[ptArrayIndex] = ws->var(ParName)->getError();
			sprintf(ParName,"var_fracNPChic1In%s",RegName);
			Frac_NP_Chic1[ptArrayIndex] = ws->var(ParName)->getVal();
			Frac_NP_Chic1_Err[ptArrayIndex] = ws->var(ParName)->getError();
			sprintf(ParName,"var_fracPRChic2In%s",RegName);
			Frac_PR_Chic2[ptArrayIndex] = ws->var(ParName)->getVal();
			Frac_PR_Chic2_Err[ptArrayIndex] = ws->var(ParName)->getError();
			sprintf(ParName,"var_fracNPChic2In%s",RegName);
			Frac_NP_Chic2[ptArrayIndex] = ws->var(ParName)->getVal();
			Frac_NP_Chic2_Err[ptArrayIndex] = ws->var(ParName)->getError();

		}


	cout<<"read all parameters"<<endl;



	TGraphErrors *graph_Frac_BG, *graph_Frac_PR_Chic0, *graph_Frac_NP_Chic0, *graph_Frac_PR_Chic1, *graph_Frac_NP_Chic1, *graph_Frac_PR_Chic2, *graph_Frac_NP_Chic2;


	graph_Frac_BG = new TGraphErrors(PtBins, pTmean, Frac_BG, pTmean_Err, Frac_BG_Err);
	graph_Frac_PR_Chic0 = new TGraphErrors(PtBins, pTmean, Frac_PR_Chic0, pTmean_Err, Frac_PR_Chic0_Err);
	graph_Frac_NP_Chic0 = new TGraphErrors(PtBins, pTmean, Frac_NP_Chic0, pTmean_Err, Frac_NP_Chic0_Err);
	graph_Frac_PR_Chic1 = new TGraphErrors(PtBins, pTmean, Frac_PR_Chic1, pTmean_Err, Frac_PR_Chic1_Err);
	graph_Frac_NP_Chic1 = new TGraphErrors(PtBins, pTmean, Frac_NP_Chic1, pTmean_Err, Frac_NP_Chic1_Err);
	graph_Frac_PR_Chic2 = new TGraphErrors(PtBins, pTmean, Frac_PR_Chic2, pTmean_Err, Frac_PR_Chic2_Err);
	graph_Frac_NP_Chic2 = new TGraphErrors(PtBins, pTmean, Frac_NP_Chic2, pTmean_Err, Frac_NP_Chic2_Err);

	double markersize[nContributions]={1., 1., 1., 1., 1., 1., 1.};
	int markerstyle[nContributions]={25, 20, 24, 20, 24, 20, 24};
	int markercolor[nContributions]={1, onia::colorChic0, onia::colorChic0, onia::colorChic1, onia::colorChic1, onia::colorChic2, onia::colorChic2};


	if(RegionCode==1 || RegionCode==4 || RegionCode==5 || RegionCode==8) markersize[0]*=1.5;
	if(RegionCode==2) markersize[3]*=1.5;
	if(RegionCode==3) markersize[5]*=1.5;
	if(RegionCode==6) markersize[4]*=1.5;
	if(RegionCode==7) markersize[6]*=1.5;


	graph_Frac_BG->SetMarkerStyle(markerstyle[0]);
	graph_Frac_BG->SetMarkerSize(markersize[0]);
	graph_Frac_BG->SetMarkerColor(markercolor[0]);
	graph_Frac_BG->SetLineColor(markercolor[0]);

	graph_Frac_PR_Chic0->SetMarkerStyle(markerstyle[1]);
	graph_Frac_PR_Chic0->SetMarkerSize(markersize[1]);
	graph_Frac_PR_Chic0->SetMarkerColor(markercolor[1]);
	graph_Frac_PR_Chic0->SetLineColor(markercolor[1]);
	graph_Frac_NP_Chic0->SetMarkerStyle(markerstyle[2]);
	graph_Frac_NP_Chic0->SetMarkerSize(markersize[2]);
	graph_Frac_NP_Chic0->SetMarkerColor(markercolor[2]);
	graph_Frac_NP_Chic0->SetLineColor(markercolor[2]);

	graph_Frac_PR_Chic1->SetMarkerStyle(markerstyle[3]);
	graph_Frac_PR_Chic1->SetMarkerSize(markersize[3]);
	graph_Frac_PR_Chic1->SetMarkerColor(markercolor[3]);
	graph_Frac_PR_Chic1->SetLineColor(markercolor[3]);
	graph_Frac_NP_Chic1->SetMarkerStyle(markerstyle[4]);
	graph_Frac_NP_Chic1->SetMarkerSize(markersize[4]);
	graph_Frac_NP_Chic1->SetMarkerColor(markercolor[4]);
	graph_Frac_NP_Chic1->SetLineColor(markercolor[4]);

	graph_Frac_PR_Chic2->SetMarkerStyle(markerstyle[5]);
	graph_Frac_PR_Chic2->SetMarkerSize(markersize[5]);
	graph_Frac_PR_Chic2->SetMarkerColor(markercolor[5]);
	graph_Frac_PR_Chic2->SetLineColor(markercolor[5]);
	graph_Frac_NP_Chic2->SetMarkerStyle(markerstyle[6]);
	graph_Frac_NP_Chic2->SetMarkerSize(markersize[6]);
	graph_Frac_NP_Chic2->SetMarkerColor(markercolor[6]);
	graph_Frac_NP_Chic2->SetLineColor(markercolor[6]);



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
	legend->AddEntry(graph_Frac_PR_Chic0,"PR #chi_{c0}","lp");
	legend->AddEntry(graph_Frac_NP_Chic0,"NP #chi_{c0}","lp");
	legend->AddEntry(graph_Frac_PR_Chic1,"PR #chi_{c1}","lp");
	legend->AddEntry(graph_Frac_NP_Chic1,"NP #chi_{c1}","lp");
	legend->AddEntry(graph_Frac_PR_Chic2,"PR #chi_{c2}","lp");
	legend->AddEntry(graph_Frac_NP_Chic2,"NP #chi_{c2}","lp");



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
	graph_Frac_PR_Chic0->Draw("Psame");
	graph_Frac_NP_Chic0->Draw("Psame");
	graph_Frac_PR_Chic1->Draw("Psame");
	graph_Frac_NP_Chic1->Draw("Psame");
	graph_Frac_PR_Chic2->Draw("Psame");
	graph_Frac_NP_Chic2->Draw("Psame");


	legend->Draw();

	if(logY) c1->SetLogy(true);
	else c1->SetLogy(false);

	c1->SaveAs(Form("%s/Fractions_In_%s_rap%d.pdf",savePath.str().c_str(),RegName,rapBin));


	///
	TFile *outfile  = new TFile(Form("%s/Par_chi.root",savePath.str().c_str()),"UPDATE");

	graph_Frac_BG->SetName(Form("graph_Frac_BG_in_%s",RegName)); graph_Frac_BG->Write();
	graph_Frac_PR_Chic0->SetName(Form("graph_Frac_PR_Chic0_in_%s",RegName)); graph_Frac_PR_Chic0->Write();
	graph_Frac_NP_Chic0->SetName(Form("graph_Frac_NP_Chic0_in_%s",RegName)); graph_Frac_NP_Chic0->Write();
	graph_Frac_PR_Chic1->SetName(Form("graph_Frac_PR_Chic1_in_%s",RegName)); graph_Frac_PR_Chic1->Write();
	graph_Frac_NP_Chic1->SetName(Form("graph_Frac_NP_Chic1_in_%s",RegName)); graph_Frac_NP_Chic1->Write();
	graph_Frac_PR_Chic2->SetName(Form("graph_Frac_PR_Chic2_in_%s",RegName)); graph_Frac_PR_Chic2->Write();
	graph_Frac_NP_Chic2->SetName(Form("graph_Frac_NP_Chic2_in_%s",RegName)); graph_Frac_NP_Chic2->Write();

	outfile->Close();

	return;

}
