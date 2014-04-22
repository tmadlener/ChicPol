#include "../interface/rootIncludes.inc"

#include <string>
#include <iostream>
#include <sstream>
using namespace std;

void muRapPt(int nState=5) {
	TGaxis::SetMaxDigits(3);

	char filename[500];
	sprintf(filename, "DataFiles/SetOfCuts0_FrameworkTest_5Dec2012/Psi%dS/tmpFiles/selEvents_data.root", nState-3);
	TFile *f = TFile::Open(filename);
	TTree *tree = (TTree*)f->Get("selectedData");

	TH2D* rapPt_muP;
	TH2D* rapPt_muN;
	TH2D* rapPt_mu;
	//rapPt_muP = new TH2D( "rapPt_muP", "rapPt_muP",50,0.,2.5,100,0.,10.);
	//rapPt_muN = new TH2D( "rapPt_muN", "rapPt_muN",50,0.,2.5,100,0.,10.);
	rapPt_muP = new TH2D( "rapPt_muP", "rapPt_muP",40,0.,2.5,100,0.,10.);
	rapPt_muN = new TH2D( "rapPt_muN", "rapPt_muN",40,0.,2.5,100,0.,10.);
	if(nState==4){
		rapPt_muP = new TH2D( "rapPt_muP", "rapPt_muP",40,0.,2.0,100,0.,10.);
		rapPt_muN = new TH2D( "rapPt_muN", "rapPt_muN",40,0.,2.0,100,0.,10.);
	}

	double MassMin;
	double MassMax;

	if(nState==4){
		MassMin=3.011;//massPsi1S-onia::nSigMass*sigma1S;
		MassMax=3.174;//massPsi1S+onia::nSigMass*sigma1S;
		// sigma  27.2 MeV
		// mean 3.093 GeV
	}
	if(nState==5){
		////pT > 7
		//MassMin=3.576;//massPsi2S-onia::nSigMass*sigma2S;
		//MassMax=3.786;//massPsi2S+onia::nSigMass*sigma2S;
		//pT > 10
		MassMin=3.578;//massPsi2S-onia::nSigMass*sigma2S;
		MassMax=3.784;//massPsi2S+onia::nSigMass*sigma2S;
		// sigma 34.9 MeV // pT > 7
		// sigma 34.3 MeV // pT > 10
		// mean 3.681 GeV
	}

	std::stringstream cutMassRapPtP, cutMassRapPtN;

	cutMassRapPtP<<"(JpsiP.M() > " << MassMin << " && JpsiP.M() < "<< MassMax << ")"
		<< "&& JpsiP.Pt() > 10. && TMath::Abs(JpsiP.Rapidity()) < 1.5"
		<< "&& lepP.Pt() < 10. && TMath::Abs(lepP.Eta()) < 2.5 ";
	cutMassRapPtN<<"(JpsiP.M() > " << MassMin << " && JpsiP.M() < "<< MassMax << ")"
		<< "&& JpsiP.Pt() > 10. && TMath::Abs(JpsiP.Rapidity()) < 1.5"
		<< "&& lepN.Pt() < 10. && TMath::Abs(lepN.Eta()) < 2.5 ";
	if(nState==4){
		cutMassRapPtP<<"(JpsiP.M() > " << MassMin << " && JpsiP.M() < "<< MassMax << ")"
			<< "&& JpsiP.Pt() > 10. && TMath::Abs(JpsiP.Rapidity()) < 1.2"
			<< "&& lepP.Pt() < 10. && TMath::Abs(lepP.Eta()) < 2.0 ";
		cutMassRapPtN<<"(JpsiP.M() > " << MassMin << " && JpsiP.M() < "<< MassMax << ")"
			<< "&& JpsiP.Pt() > 10. && TMath::Abs(JpsiP.Rapidity()) < 1.2"
			<< "&& lepN.Pt() < 10. && TMath::Abs(lepN.Eta()) < 2.0 ";
	}

	cout<<"produce rapPt_muP.."<<endl;
	tree->Draw("lepP.Pt():TMath::Abs(lepP.Eta())>>rapPt_muP",cutMassRapPtP.str().c_str(),"colz");
	cout<<"produce rapPt_muN.."<<endl;
	tree->Draw("lepN.Pt():TMath::Abs(lepN.Eta())>>rapPt_muN",cutMassRapPtN.str().c_str(),"colz");
	cout<<"done!"<<endl;

	rapPt_mu = (TH2D*)rapPt_muP->Clone();
	rapPt_mu->Add(rapPt_muN);

	gStyle->SetPadBottomMargin(0.12);
	//gStyle->SetPadLeftMargin(0.08); //0.12
	gStyle->SetPadRightMargin(0.02); //0.05
	//gStyle->SetPadTopMargin(0.02); //0.05

	gStyle->SetPadLeftMargin(0.12); //0.12
	gStyle->SetPadTopMargin(0.05); //0.05

	//TCanvas* c2 = new TCanvas("c2","c2",700,500);
	TCanvas* c2 = new TCanvas("c2","c2",600,500);

	rapPt_mu->SetYTitle("p_{T}^{#mu} [GeV]");
	rapPt_mu->SetXTitle("|#eta^{#mu}|");
	gStyle->SetPalette(1);
	gPad->SetFillColor(kWhite);
	rapPt_mu->SetTitle(0);
	rapPt_mu->SetStats(0);
	gPad->SetRightMargin(0.17);
	rapPt_mu->GetXaxis()->SetTitleOffset(1.2);
	rapPt_mu->GetYaxis()->SetTitleOffset(1.1);

	rapPt_mu->Draw("colz");
	if(nState==4)
	 	c2->Print(Form("singleMu2DDist_psi%dS.pdf",nState-3));
	else
	 	c2->Print(Form("singleMu2DDist_psi%dS_rap1p5.pdf",nState-3));
}
