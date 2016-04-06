#include "rootIncludes.inc"
#include "commonVar.h"
#include "ToyMC.h"
#include "effsAndCuts.h"

#include "TSystem.h"
#include "TROOT.h"
#include "polGen.C"
#include "polRec.C"
#include "polFit.C"
#include "polPlot.C"
#include "TGraphAsymmErrors.h"
#include "TFile.h"

#include <time.h>

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

//====================================
int main(int argc, char* argv[]) {

	int nGenerations=999;
	int polScenSig=9999;
	int frameSig=999;
	int polScenBkg=999;
	int frameBkg=999;
	int rapBinMin=999;
	int rapBinMax=999;
	int ptBinMin=999;
	int ptBinMax=999;
	int nEff=999;
	int nRhoFactor=999;
	int nDileptonEff=999;
	int nRecEff=999;
	int nRecDileptonEff=999;
	int nRecRhoFactor=999;
	int FidCuts=999;
	int nSample=999;
	int ConstEvents=999;
	int nSkipGen=999;
	int ThisGen=999;
	int MPValgo=0;
	int nState=0;
	int nAmap=999;
	int nDenominatorAmap=999;
        int batch = 0,
            statBGF = 0,
            statBG = 0,
            statRho = 0,
            cutDelta = 0,
            nSigma = 0;

	bool ConstEvents_ = false,
            gen = false,
            rec = false,
            fit = false,
            plot = false,
            RealData = false,
            UseDifferingEff = false,
            MCeff = false,
            MCReceff = false,
            MCDileptoneff = false,
            MCDileptonReceff = false,
            scalePlots = false,
            NewAccCalc = false,
            deletePseudoData = false,
            useAmapApproach = false,
            useBatch = false,
            StatVarTotBGfraction = false,
            StatVarTotBGmodel = false,
            StatVarRho = false,
            cutDeltaREllDpt = false,
            PlotForPaper = false;

	double n_sigmas_signal = 3.;

        std::string storagedir, //Storage Directory
            basedir, //Code Directory
            JobID,
            TreeID;

        char *realdatadir;

        //Loop over argument list
        for (int i=1; i < argc; i++)
        {
            std::string arg = argv[i];
            fromSplit("nState", arg, nState);
            fromSplit("ptBinMin", arg, ptBinMin);
            fromSplit("ptBinMax", arg, ptBinMax);
            fromSplit("rapBinMin", arg, rapBinMin);
            fromSplit("rapBinMax", arg, rapBinMax);
            fromSplit("nGenerations", arg, nGenerations);
            fromSplit("frameSig", arg, frameSig);
            fromSplit("polScenSig", arg, polScenSig);
            fromSplit("frameBkg", arg, frameBkg);
            fromSplit("polScenBkg", arg, polScenBkg);
            fromSplit("nEff", arg, nEff);
            fromSplit("nDiEff", arg, nDileptonEff);
            fromSplit("nRhoFactor", arg, nRhoFactor);
            fromSplit("FidCuts", arg, FidCuts);
            fromSplit("nSample", arg, nSample);
            fromSplit("nSigma", arg, nSigma);
            fromSplit("ConstEvents", arg, ConstEvents);
            fromSplit("nSkipGen", arg, nSkipGen);
            fromSplit("ThisGen", arg, ThisGen);
            fromSplit("nAmap", arg, nAmap);
            fromSplit("nDenominatorAmap", arg, nDenominatorAmap);
            fromSplit("MPValgo", arg, MPValgo);
            fromSplit("UseBatch", arg, batch);
            if(batch == 1) useBatch = true;
            fromSplit("StatVarTotBGfraction", arg, statBGF);
            if(statBGF == 1) StatVarTotBGfraction = true;
            fromSplit("StatVarTotBGmodel", arg, statBG);
            if(statBG == 1) StatVarTotBGmodel = true;
            fromSplit("StatVarRho", arg, statRho);
            if(statRho == 1) StatVarRho = true;
            fromSplit("cutDeltaREllDpt", arg, cutDelta);
            if(cutDelta == 1) cutDeltaREllDpt = true;

            fromSplit("NewAccCalc", arg, NewAccCalc);
            fromSplit("UseConstEv", arg, ConstEvents_);
            //if(ConstEvents_ == true) cout<<"use constant number of reconstructed events"<<endl;
            fromSplit("deletePseudoData", arg, deletePseudoData);
            fromSplit("gen", arg, gen);
            fromSplit("rec", arg, rec);
            fromSplit("fit", arg, fit);
            fromSplit("plot", arg, plot);
            fromSplit("UseDifferingEff", arg, UseDifferingEff);
            fromSplit("nRecEff", arg, nRecEff);
            fromSplit("nRecDiEff", arg, nRecDileptonEff);
            fromSplit("nRecRhoFactor", arg, nRecRhoFactor);
            fromSplit("scalePlots", arg, scalePlots);
            fromSplit("useAmapApproach", arg, useAmapApproach);
            fromSplit("UseMCeff", arg, MCeff);
            fromSplit("UseMCReceff", arg, MCReceff);
            fromSplit("UseMCDileptoneff", arg, MCDileptoneff);
            fromSplit("UseMCDileptonReceff", arg, MCDileptonReceff);
            fromSplit("PlotForPaper", arg, PlotForPaper);

            fromSplit("JobID", arg, JobID);
            fromSplit("basedir", arg, basedir);
            fromSplit("storagedir", arg, storagedir);
            fromSplit("realdatadir", arg, realdatadir);
            fromSplit("TreeID", arg, TreeID);
        }

	double mass_signal_peak = 0, mass_signal_sigma = 0;
	if(nState < 4 && rapBinMin==1) mass_signal_sigma=0.075;
	else if(nState < 4 && rapBinMin==2) mass_signal_sigma=0.1;
	else if(nState == 4 && rapBinMin==1) mass_signal_sigma=0.025;
	else if(nState == 4 && rapBinMin==2) mass_signal_sigma=0.035;
	else if(nState == 5 && rapBinMin==1) mass_signal_sigma=0.026;
	else if(nState == 5 && rapBinMin==2) mass_signal_sigma=0.038;
	else if(nState == 5 && rapBinMin==3) mass_signal_sigma=0.048;
	else if(nState == 6) mass_signal_sigma=0.006;
	else if(nState == 7) mass_signal_sigma=0.006;

	if(nState==1) mass_signal_peak=9.5;
	else if(nState==2) mass_signal_peak=10.;
	else if(nState==3) mass_signal_peak=10.4;
	else if(nState==4) mass_signal_peak=3.097;
	else if(nState==5) mass_signal_peak=3.686;
	else if(nState==6) mass_signal_peak=3.51;
        else if(nState==7) mass_signal_peak=3.556;

        // get polarization values 
	double lambda_theta_sig_ = 0, lambda_phi_sig_ = 0, lambda_thetaphi_sig_ = 0;
	double lambda_theta_bkg_ = ToyMC::ScenarioBkg[0][polScenBkg-1];
	double lambda_phi_bkg_ = ToyMC::ScenarioBkg[1][polScenBkg-1];
	double lambda_thetaphi_bkg_ = ToyMC::ScenarioBkg[2][polScenBkg-1];

        // inject polarization from data
	bool injectRealDataPolarization = false;
	if(polScenSig==999) injectRealDataPolarization=true;
	if(injectRealDataPolarization){

            cout<<"injectRealDataPolarization"<<endl;            
            std::stringstream filename;
            filename << "/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/Systematics/TotalSyst/May20Centrals_CentralsFromAlteredPPDMay20_1SigmaStatError_FracHighCorrected/TGraphResults_Psi%" << nState << "S.root";
            cout << filename.str().c_str() << endl;
            TFile *infile1 = new TFile(filename.str().c_str(),"READ");

            int nRapBins = 2; 
            if(nState == 5) nRapBins=3;
            char GraphName[100];
            
            for(int rapBin = 1; rapBin < nRapBins+1; rapBin++){

                TGraphAsymmErrors *graph_lth, *graph_lph, *graph_ltp, *graph_lthstar, *graph_lphstar, *graph_ltilde;

                if(frameSig==1)  {
                    cout<<GraphName<<endl;
                    sprintf(GraphName,"lth_CS_rap%d",rapBin);
                    graph_lth = (TGraphAsymmErrors*) infile1->Get(GraphName);
                    sprintf(GraphName,"lph_CS_rap%d",rapBin);
                    graph_lph = (TGraphAsymmErrors*) infile1->Get(GraphName);
                    sprintf(GraphName,"ltp_CS_rap%d",rapBin);
                    graph_ltp = (TGraphAsymmErrors*) infile1->Get(GraphName);
                    sprintf(GraphName,"lthstar_CS_rap%d",rapBin);
                    graph_lthstar = (TGraphAsymmErrors*) infile1->Get(GraphName);
                    sprintf(GraphName,"lphstar_CS_rap%d",rapBin);
                    graph_lphstar = (TGraphAsymmErrors*) infile1->Get(GraphName);
                    sprintf(GraphName,"ltilde_CS_rap%d",rapBin);
                    graph_ltilde = (TGraphAsymmErrors*) infile1->Get(GraphName);
                }
                else if(frameSig==2)  {
                    sprintf(GraphName,"lth_HX_rap%d",rapBin);
                    graph_lth = (TGraphAsymmErrors*) infile1->Get(GraphName);
                    sprintf(GraphName,"lph_HX_rap%d",rapBin);
                    graph_lph = (TGraphAsymmErrors*) infile1->Get(GraphName);
                    sprintf(GraphName,"ltp_HX_rap%d",rapBin);
                    graph_ltp = (TGraphAsymmErrors*) infile1->Get(GraphName);
                    sprintf(GraphName,"lthstar_HX_rap%d",rapBin);
                    graph_lthstar = (TGraphAsymmErrors*) infile1->Get(GraphName);
                    sprintf(GraphName,"lphstar_HX_rap%d",rapBin);
                    graph_lphstar = (TGraphAsymmErrors*) infile1->Get(GraphName);
                    sprintf(GraphName,"ltilde_HX_rap%d",rapBin);
                    graph_ltilde = (TGraphAsymmErrors*) infile1->Get(GraphName);
                }
                else if(frameSig==3)  {
                    sprintf(GraphName,"lth_PX_rap%d",rapBin);
                    graph_lth = (TGraphAsymmErrors*) infile1->Get(GraphName);
                    sprintf(GraphName,"lph_PX_rap%d",rapBin);
                    graph_lph = (TGraphAsymmErrors*) infile1->Get(GraphName);
                    sprintf(GraphName,"ltp_PX_rap%d",rapBin);
                    graph_ltp = (TGraphAsymmErrors*) infile1->Get(GraphName);
                    sprintf(GraphName,"lthstar_PX_rap%d",rapBin);
                    graph_lthstar = (TGraphAsymmErrors*) infile1->Get(GraphName);
                    sprintf(GraphName,"lphstar_PX_rap%d",rapBin);
                    graph_lphstar = (TGraphAsymmErrors*) infile1->Get(GraphName);
                    sprintf(GraphName,"ltilde_PX_rap%d",rapBin);
                    graph_ltilde = (TGraphAsymmErrors*) infile1->Get(GraphName);
                }
                
                cout<<"TGraphs of all parameters loaded for frame "<<frameSig<<endl;

                double ptCentre, lth_lmean, lph_lmean, ltp_lmean;
                graph_lth->GetPoint(ptBinMin-1,ptCentre,lth_lmean);
                graph_lph->GetPoint(ptBinMin-1,ptCentre,lph_lmean);
                graph_ltp->GetPoint(ptBinMin-1,ptCentre,ltp_lmean);
                cout<<"Values of all parameters loaded for frame "<<frameSig<<endl;

                if(rapBin==rapBinMin){
                    lambda_theta_sig_ = lth_lmean;
                    lambda_phi_sig_ = lph_lmean;
                    lambda_thetaphi_sig_ = ltp_lmean;
                }
                
            }

            cout<<"Using Real Data Results as Input for toyMC-samples, rap"<<rapBinMin<<"_pT"<<ptBinMin<<", injected in frame "<<frameSig<<endl;
            cout<<"lth = "<<lambda_theta_sig_<<endl;
            cout<<"lph = "<<lambda_phi_sig_<<endl;
            cout<<"ltp = "<<lambda_thetaphi_sig_<<endl;
            
	}
        // inject polarization values from file
        else{
            lambda_theta_sig_ = ToyMC::ScenarioSig[0][polScenSig-1];
            lambda_phi_sig_ = ToyMC::ScenarioSig[1][polScenSig-1];
            lambda_thetaphi_sig_ = ToyMC::ScenarioSig[2][polScenSig-1];
	}

	if(!UseDifferingEff) {nRecEff=nEff; nRecDileptonEff=nDileptonEff; MCReceff=MCeff; MCDileptonReceff=MCDileptoneff; nRecRhoFactor=nRhoFactor;}

	Char_t *OutputDirectory;
        Char_t *TreeBinID_dataFile = "ToyMC";
	Char_t *TreeBinID;
	double f_BG;
	int n_events;
	char basestruct[1000],substruct[1000], dirstruct[1000], rapptstruct[1000], filenameFrom[1000], filenameTo[1000] , tmpfilename[1000], TreeBinID_[1000], effDir[1000];

	sprintf(basestruct,"%s/%s",storagedir.c_str(),JobID.c_str());
        gSystem->mkdir(basestruct);
	sprintf(substruct,"%s/Sig_frame%dscen%d_Bkg_frame%dscen%d",basestruct,frameSig,polScenSig,frameBkg,polScenBkg); 
        if(!RealData) gSystem->mkdir(substruct);
	sprintf(effDir,"%s/macros/polFit/EffFiles",basedir.c_str());

	cout<<"storagedir: "<<storagedir.c_str()<<endl;
	cout<<"JobID: "<<JobID.c_str()<<endl;
	cout<<"basestruct: "<<basestruct<<endl;
	cout<<"substruct: "<<substruct<<endl;

	time_t seconds; seconds = time (NULL); double time_0=seconds; double time_1;

	int iRap = rapBinMin;
	int iPt = ptBinMin;
	sprintf(TreeBinID_,"%s_rap%d_pT%d",TreeID.c_str(),iRap,iPt);
	char TreeBinID_dataFileChar[200];
	sprintf(TreeBinID_dataFileChar,"%s",TreeBinID_);
	TreeBinID_dataFile=TreeBinID_dataFileChar;
	if(useBatch) sprintf(TreeBinID_,"Fit%d_%s_rap%d_pT%d",ThisGen,TreeID.c_str(),iRap,iPt);
	TreeBinID=TreeBinID_;

	cout<<TreeBinID_dataFile<<endl;
	cout<<TreeBinID<<endl;

	double ptlow=onia::pTRange[iRap][iPt-1];
	double pthigh=onia::pTRange[iRap][iPt];
	double raplow=onia::rapForPTRange[iRap-1];
	double raphigh=onia::rapForPTRange[iRap];

	sprintf(rapptstruct,"%s/rap%d_pT%d",substruct,iRap,iPt); if(!RealData) gSystem->mkdir(rapptstruct);

	if(gen) {sprintf(filenameFrom,"%s/polGen.C",substruct);					sprintf(filenameTo,"%s/polGen.C",rapptstruct); 							gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);}
	if(rec) {sprintf(filenameFrom,"%s/polRec.C",substruct);					sprintf(filenameTo,"%s/polRec.C",rapptstruct); 							gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);}
	if(fit) {sprintf(filenameFrom,"%s/polFit.C",substruct);					sprintf(filenameTo,"%s/polFit.C",rapptstruct); 							gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);}
	if(plot) {sprintf(filenameFrom,"%s/polPlot.C",substruct);				sprintf(filenameTo,"%s/polPlot.C",rapptstruct); 						gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);}

	sprintf(filenameFrom,"%s/effsAndCuts.h",substruct);					sprintf(filenameTo,"%s/effsAndCuts.h",rapptstruct); 							gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);

	gSystem->cd(rapptstruct);

	//if(gen) gROOT->ProcessLine(".L polGen.C+");
	//if(rec) gROOT->ProcessLine(".L polRec.C+");
	//if(fit) gROOT->ProcessLine(".L polFit.C+");
	//if(plot) gROOT->ProcessLine(".L polPlot.C+");

	/// Extract number of signal and background events to be generated, as well as f_BG to be generated to result in desired effective f_BG:
	int nTargetEvents;
	nTargetEvents = ToyMC::numEvents[iRap-1][iPt-1];

	if(gen){

            int numEvCheck = 500000;
            f_BG = ToyMC::fracBackground[iRap-1][iPt-1];
            sprintf(tmpfilename,"%s/data.root",rapptstruct);
            TFile* dataFile = new TFile(tmpfilename, "READ");
            cout<<"f_BG: "<<f_BG<<endl;
            
            if(dataFile->Get("isBGdistribution")==NULL || dataFile == NULL){
                OutputDirectory=rapptstruct;
                polGen(raplow,raphigh,ptlow,pthigh,mass_signal_peak,mass_signal_sigma,n_sigmas_signal,numEvCheck,f_BG,lambda_theta_sig_,lambda_phi_sig_,lambda_thetaphi_sig_,lambda_theta_bkg_,lambda_phi_bkg_,lambda_thetaphi_bkg_,frameSig,frameBkg,-999,nState,OutputDirectory);
                std::cout << "gen finished" << std::endl;
                if(rec)polRec(raplow,raphigh,ptlow,pthigh,mass_signal_peak,mass_signal_sigma,n_sigmas_signal,nRecEff,nRecDileptonEff,nRecRhoFactor,FidCuts,OutputDirectory, true, effDir, MCReceff, MCDileptonReceff, iRap, iPt, useAmapApproach, nAmap, nDenominatorAmap);
                sprintf(tmpfilename,"%s/genData.root",rapptstruct);			gSystem->Unlink(tmpfilename);
                sprintf(tmpfilename,"%s/GenResults.root",rapptstruct);		gSystem->Unlink(tmpfilename);
            }

            sprintf(tmpfilename,"%s/data.root",rapptstruct);
            dataFile = new TFile(tmpfilename, "READ");
            TH1D* isBG_distribution = (TH1D*)dataFile->Get("isBGdistribution");
            
            double sigFact = isBG_distribution->GetBinContent(1)/(numEvCheck*(1-f_BG));
            double bkgFact = isBG_distribution->GetBinContent(2)/(numEvCheck*f_BG);
            
            dataFile->Close();
            
            
            if(ConstEvents_) nTargetEvents = ConstEvents;
            
            n_events = nTargetEvents/sigFact+(nTargetEvents/(1-f_BG)-nTargetEvents)/bkgFact;
            
            f_BG = (n_events-nTargetEvents/sigFact)/n_events;
            cout<<"n_events: "<<n_events<<endl;
            
	}
	/// Start actual Generation and Fits:
	cout<<"Start actual Generation and Fits....."<<endl;

	int iGen = ThisGen;
	int nTotalFits = nSkipGen+nGenerations;

	seconds = time (NULL); time_1=seconds;
	//sprintf(dirstruct,"%s/Generation%d",rapptstruct,iGen+nSkipGen); if(!RealData) gSystem->mkdir(dirstruct);
	sprintf(dirstruct,"%s/Generation%d",rapptstruct,iGen); if(!RealData) gSystem->mkdir(dirstruct);
	OutputDirectory=dirstruct;
	if(RealData) OutputDirectory=basestruct;

	cout<<"nState: "<<nState<<endl;
	cout<<"OutputDirectory: "<<OutputDirectory<<endl;
	cout<<"basestruct: "<<basestruct<<endl;

	if(gen) polGen(raplow,raphigh,ptlow,pthigh,mass_signal_peak,mass_signal_sigma,n_sigmas_signal,n_events,f_BG,lambda_theta_sig_,lambda_phi_sig_,lambda_thetaphi_sig_,lambda_theta_bkg_,lambda_phi_bkg_,lambda_thetaphi_bkg_,frameSig,frameBkg,iGen,nState,OutputDirectory);
	if(rec) polRec(raplow,raphigh,ptlow,pthigh,mass_signal_peak,mass_signal_sigma,n_sigmas_signal,nRecEff,nRecDileptonEff,nRecRhoFactor,FidCuts,OutputDirectory, false, effDir, MCReceff, MCDileptonReceff, iRap, iPt, useAmapApproach, nAmap, nDenominatorAmap, StatVarTotBGfraction, StatVarRho);
	if(fit) polFit(nSample,FidCuts, nEff, nDileptonEff, nRhoFactor, OutputDirectory, realdatadir, TreeBinID, TreeBinID_dataFile, RealData, effDir, MCeff, MCDileptoneff, iRap, iPt, NewAccCalc, MPValgo, useAmapApproach, nAmap, nDenominatorAmap, StatVarTotBGfraction, StatVarTotBGmodel, StatVarRho, cutDeltaREllDpt);
	if(plot) polPlot(OutputDirectory, TreeBinID, RealData, MPValgo, scalePlots, nTotalFits, nState, ptlow, pthigh, raplow, raphigh, PlotForPaper);

	//sprintf(dirstruct,"%s/Generation%d",rapptstruct,iGen+nSkipGen);
	sprintf(dirstruct,"%s/Generation%d",rapptstruct,iGen);
	if(deletePseudoData){
		sprintf(tmpfilename,"%s/genData.root",dirstruct);			gSystem->Unlink(tmpfilename);
		sprintf(tmpfilename,"%s/data.root",dirstruct);				gSystem->Unlink(tmpfilename);
		sprintf(tmpfilename,"%s/efficiency.root",dirstruct);		gSystem->Unlink(tmpfilename);
	}

	seconds = time (NULL);

	if(fit) cout<<"Proccessing time for this generation: "<<seconds-time_1<<" s"<<" (Corresponding to "<<(seconds-time_1)/60<<" min)"<<endl;
	if(fit) cout<<"Per signal event in final sample:     "<<(seconds-time_1)/nTargetEvents*1000<<" ms"<<endl;

	if(gen)  {sprintf(tmpfilename,"%s/polGen.C",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polGen_C.d",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polGen_C.so",rapptstruct);			gSystem->Unlink(tmpfilename);}
	if(rec)  {sprintf(tmpfilename,"%s/polRec.C",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polRec_C.d",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polRec_C.so",rapptstruct);			gSystem->Unlink(tmpfilename);}
	if(fit)  {sprintf(tmpfilename,"%s/polFit.C",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polFit_C.d",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polFit_C.so",rapptstruct);			gSystem->Unlink(tmpfilename);}
	if(plot)  {sprintf(tmpfilename,"%s/polPlot.C",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polPlot_C.d",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polPlot_C.so",rapptstruct);			gSystem->Unlink(tmpfilename);}
	if(fit && RealData)  {sprintf(tmpfilename,"%s/polFit_C.d",basestruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polFit_C.so",basestruct);			gSystem->Unlink(tmpfilename);}
	if(plot && RealData)  {sprintf(tmpfilename,"%s/polPlot_C.d",basestruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polPlot_C.so",basestruct);			gSystem->Unlink(tmpfilename);}

	sprintf(tmpfilename,"%s/effsAndCuts.h",rapptstruct);			gSystem->Unlink(tmpfilename);

	return 0;
}


