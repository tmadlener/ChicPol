#include "../../interface/rootIncludes.inc"

#include <string>
#include <iostream>
#include <sstream>
using namespace std;

void plotSingleMuPtEta() {
	TGaxis::SetMaxDigits(3);

	char directory[200];
	char Fig_directory[200];
	sprintf(Fig_directory,"FigBuffer/KinematicDist");
	gSystem->mkdir(Fig_directory,true);

	char DataID[200];
	//sprintf(DataID,"SetOfCuts11_ctauScen5_FracLSB-1_newMLfit_30Apr2013_correctfLSB");
	sprintf(DataID,"SetOfCuts11_ctauScen5_MC");
	//sprintf(Fig_directory,"%s/%s/Psi%dS",Fig_directory,DataID,nState-3);
	//gSystem->mkdir(Fig_directory,true);

	char filenamePre[200];
	sprintf(filenamePre, "/afs/ihep.ac.cn/users/z/zhangll/fs/work/polarization/PsiPol2011/macros/DataFiles");

	int pTmin=2;
	int pTmax=5;
	int rapmin=1;
	int rapmax=1;
	int nStateMin=4;
	int nStateMax=4;

	for(int nState=nStateMin; nState<=nStateMax; nState++){
		if(nState==4) { pTmin=3; pTmax=12; rapmin=1; rapmax=1;}
		if(nState==5) { pTmin=2; pTmax=5; rapmin=1; rapmax=3;}

	  sprintf(Fig_directory,"%s/%s/Psi%dS",Fig_directory,DataID,nState-3);
	  gSystem->mkdir(Fig_directory,true);

		cout << "nState: " << nState << "\n"
			<< "pTmin: " << pTmin << "\n"
			<< "pTmax: " << pTmax << "\n"
			<< "rapmin: " << rapmin << "\n"
			<< "rapmax: " << rapmax << endl;

		for(int irap=rapmin; irap<=rapmax; irap++){
			for(int ipt=pTmin; ipt<=pTmax; ipt++){

				char filename[500];
				sprintf(filename,"%s/%s/Psi%dS/tmpFiles/data_Psi%dS_rap%d_pT%d.root",filenamePre, DataID, 
						nState-3, nState-3, irap, ipt);
				cout << "filename: " << filename << endl;
				TFile *infile = new TFile(filename,"R");
				if(!infile) { cout<<"failed to open file.."<<endl; return; }

				TTree *tree = (TTree*) infile -> Get("selectedData") ;
				TLorentzVector* lepP  = new TLorentzVector();  tree -> SetBranchAddress( "lepP",   &lepP ); 
				TLorentzVector* lepN  = new TLorentzVector();  tree -> SetBranchAddress( "lepN",   &lepN );
				TLorentzVector* JpsiP = new TLorentzVector();  tree -> SetBranchAddress( "JpsiP",  &JpsiP );


				int nbinsRap = 80, nbinsPt = 100 ;
				double rapMin = -1.8, rapMax = 1.8 ;
				double ptMin = 0., ptMax = 50. ;
				if(ipt==12) ptMax = 70.;
				if(irap==1) { rapMin = -1.; rapMax = 1.; }

				TH2D* singleMu_EtaPt ; 
				TH1D* singleMu_Eta; TH1D* singleMu_Pt;
			  singleMu_EtaPt	= new TH2D( "singleMu_EtaPt", "singleMu_EtaPt", nbinsRap, rapMin, rapMax, nbinsPt, ptMin, ptMax);
				singleMu_Eta    = new TH1D( "singleMu_Eta", "singleMu_Eta", nbinsRap, rapMin, rapMax);
				singleMu_Pt     = new TH1D( "singleMu_Pt", "singleMu_Pt", nbinsPt, ptMin, ptMax);

				int n_events = tree -> GetEntries();

				int n_step = n_events/5; 
				int n_step_=1;

				for ( int i_event = 1; i_event <= n_events; i_event++ ) {
					  if (i_event%n_step == 0) {cout << n_step_*20 <<" % "<<endl; n_step_++;}

						tree->GetEvent( i_event-1 );

						double lepP_pT  = lepP->Pt();
						double lepN_pT  = lepN->Pt();

						double lepP_eta = lepP->PseudoRapidity();
						double lepN_eta = lepN->PseudoRapidity();

						// dilepton 4-vector:
						TLorentzVector dilepton = *lepP + *lepN;
						double pT   = dilepton.Pt();
						double rap  = dilepton.Rapidity();
						double mass = dilepton.M();

						singleMu_Eta -> Fill( lepP_eta ); 
						singleMu_Eta -> Fill( lepN_eta ); 

						singleMu_Pt -> Fill ( lepP_pT );
						singleMu_Pt -> Fill ( lepN_pT );

						singleMu_EtaPt -> Fill ( lepP_eta, lepP_pT );
						singleMu_EtaPt -> Fill ( lepN_eta, lepN_pT );

						//if(lepP_eta < 0.) cout << "lepP_eta: " << lepP_eta << endl;
						//if(lepP_pT > 50. ) cout << "lepP_pT: " << lepP_pT << endl;
						//if(lepN_pT > 50. ) cout << "lepN_pT: " << lepN_pT << endl;

				}//i_event


				gStyle->SetPadBottomMargin(0.12);
				gStyle->SetPadLeftMargin(0.12);
				gStyle->SetPadRightMargin(0.02);
				gStyle->SetPadTopMargin(0.05); 

				TCanvas* plotCanvas = new TCanvas("plotCanvas","plotCanvas",600,500);
				gStyle->SetPalette(1);
				gPad->SetFillColor(kWhite);
				gPad->SetRightMargin(0.17);

				singleMu_EtaPt->SetYTitle("p_{T}^{#mu} [GeV]");
				singleMu_EtaPt->SetXTitle("|#eta^{#mu}|");
				singleMu_EtaPt->SetTitle(0);
				singleMu_EtaPt->SetStats(0);
				singleMu_EtaPt->GetXaxis()->SetTitleOffset(1.2);
				singleMu_EtaPt->GetYaxis()->SetTitleOffset(1.1);

				singleMu_Eta->SetXTitle("|#eta^{#mu}|");
				singleMu_Eta->SetTitle(0);
				singleMu_Eta->SetStats(0);
				singleMu_Eta->GetXaxis()->SetTitleOffset(1.2);
				singleMu_Eta->GetYaxis()->SetTitleOffset(1.1);

				singleMu_Pt->SetXTitle("p_{T}^{#mu} [GeV]");
				singleMu_Pt->SetTitle(0);
				singleMu_Pt->SetStats(0);
				singleMu_Pt->GetXaxis()->SetTitleOffset(1.2);
				singleMu_Pt->GetYaxis()->SetTitleOffset(1.1);

				singleMu_EtaPt->Draw("colz");
				plotCanvas->Print(Form("%s/singleMu_EtaPt_Psi%dS_rap%d_pT%d.pdf",Fig_directory, nState-3, irap, ipt));

				singleMu_Eta->Draw("pe");
				plotCanvas->Print(Form("%s/singleMu_Eta_Psi%dS_rap%d_pT%d.pdf",Fig_directory, nState-3, irap, ipt));

				singleMu_Pt->Draw("pe");
				plotCanvas->Print(Form("%s/singleMu_Pt_Psi%dS_rap%d_pT%d.pdf",Fig_directory, nState-3, irap, ipt));

			}//ipt
		}//irap

	}//nState

}
