#include "rootIncludes.inc"
#include "commonVar.h"

#include <string>
#include <iostream>
#include <sstream>

using namespace RooFit;

void createWorkspace(const std::string &infilename, int nState, bool correctCtau, bool drawRapPt2D){
	gROOT->SetStyle("Plain");
	gStyle->SetTitleBorderSize(0);

	// Set some strings
	const std::string workspacename = "ws_masslifetime",
				treename = "selectedData";

	// Get the tree from the data file
	TFile *f = TFile::Open(infilename.c_str());
	TTree *tree = (TTree*)f->Get(treename.c_str());

	// Set branch addresses in tree to be able to import tree to roofit
	TLorentzVector* chic = new TLorentzVector;
	tree->SetBranchAddress("chic",&chic);
	TLorentzVector* chic_rf = new TLorentzVector;
	tree->SetBranchAddress("chic_rf",&chic_rf);
	TLorentzVector* jpsi = new TLorentzVector;
	tree->SetBranchAddress("jpsi",&jpsi);
	double lifetime = 0;
	tree->SetBranchAddress("Jpsict",&lifetime);
	double lifetimeErr = 0;
	tree->SetBranchAddress("JpsictErr",&lifetimeErr);

	char lifetimeTitle[200];

	sprintf(lifetimeTitle,"l^{#psi} [mm]");
	if(correctCtau) sprintf(lifetimeTitle,"l^{#chi} [mm]");

	// define variables necessary for J/Psi(Psi(2S)) mass,lifetime fit
	RooRealVar* chicMass =
		new RooRealVar("chicMass", "M^{#chi} [GeV]", onia::chimassMin, onia::chimassMax);
	RooRealVar* chicRap =
		new RooRealVar("chicRap", "y^{#chi}", -onia::chirap, onia::chirap);
	RooRealVar* chicPt =
		new RooRealVar("chicPt", "p^{#chi}_{T} [GeV]", 0. ,100.);
	RooRealVar* Jpsict =
		new RooRealVar("Jpsict", lifetimeTitle, onia::ctVarMin, onia::ctVarMax);
	RooRealVar* JpsictErr =
		new RooRealVar("JpsictErr", Form("Error on %s",lifetimeTitle), 0.0001, 10);

	// Set bins
	Jpsict->setBins(10000,"cache");
	Jpsict->setBins(100);
	chicMass->setBins(100);
	JpsictErr->setBins(100);

	// The list of data variables    
	RooArgList dataVars(*chicMass,*chicRap,*chicPt,*Jpsict,*JpsictErr);

	// construct dataset to contain events
	RooDataSet* fullData = new RooDataSet("fullData","The Full Data From the Input ROOT Trees",dataVars);

	int entries = tree->GetEntries();
	cout << "entries " << entries << endl;

	int numEntriesTotal=0;
	int numEntriesInAnalysis=0;
	int numEntriesNotInAnalysis=0;

	// loop through events in tree and save them to dataset
	for (int ientries = 0; ientries < entries; ientries++) {
		numEntriesTotal++;
		if (ientries%10000==0) std::cout << "event " << ientries << " of " << entries <<  std::endl;

		tree->GetEntry(ientries);

		double M =chic_rf->M();
		double y=chic->Rapidity();
		double pt=chic->Pt();



		if (M > chicMass->getMin() && M < chicMass->getMax()
				&& pt > chicPt->getMin() && pt < chicPt->getMax()
				&& y > chicRap->getMin() && y < chicRap->getMax()
				&& lifetime > Jpsict->getMin() && lifetime < Jpsict->getMax()
				&& lifetimeErr > JpsictErr->getMin() && lifetimeErr < JpsictErr->getMax()
			 ){

			chicPt      ->setVal(pt);
			chicRap     ->setVal(y);
			chicMass    ->setVal(M);
			Jpsict      ->setVal(lifetime);
			JpsictErr   ->setVal(lifetimeErr);

			fullData->add(dataVars);
			numEntriesInAnalysis++;
		}
		else{
			numEntriesNotInAnalysis++;
			//if (M < chicMass->getMin() || M > chicMass->getMax()) cout << "M " << M << endl;
			//if (pt < chicPt->getMin() || pt > chicPt->getMax()) cout << "pt " << pt << endl;
			//if (y < chicRap->getMin() || y > chicRap->getMax()) cout << "y " << y << endl;
			//if (lifetime < Jpsict->getMin() || lifetime > Jpsict->getMax()) cout << "lifetime " << lifetime << endl;
			//if (lifetimeErr < JpsictErr->getMin() || lifetimeErr > JpsictErr->getMax()) cout << "lifetimeErr " << lifetimeErr << endl;
			//cout << "M " << M << endl;
			//cout << "pt " << pt << endl;
			//cout << "y " << y << endl;
			//cout << "lifetime " << lifetime << endl;
			//cout << "lifetimeErr " << lifetimeErr << endl;
			//cout << " " << endl;

		}

	}//ientries

	cout << "entries entering all bins " << fullData->sumEntries() << endl;
	cout << "numEntriesTotal " << numEntriesTotal << endl;
	cout << "numEntriesInAnalysis " << numEntriesInAnalysis << endl;
	cout << "numEntriesNotInAnalysis " << numEntriesNotInAnalysis << endl;

	//------------------------------------------------------------------------------------------------------------------
	// Define workspace and import datasets

	////Get datasets binned in pT an y

	for(int iRap = 0; iRap <= onia::kNbRapForPTBins; iRap++){

		Double_t yMin;
		Double_t yMax;
		if(iRap==0){
			yMin = onia::rapForPTRange[0];
			yMax = onia::rapForPTRange[onia::kNbRapForPTBins];
		} else{
			yMin = onia::rapForPTRange[iRap-1];
			yMax = onia::rapForPTRange[iRap];
		}

		for(int iPT = 0; iPT <= onia::kNbPTBins[iRap]; iPT++){
			//for(int iPT = 0; iPT <= 0; iPT++)

			Double_t ptMin;
			Double_t ptMax;
			if(iPT==0){
				ptMin = onia::pTRange[iRap][0];
				ptMax = onia::pTRange[iRap][onia::kNbPTBins[0]];
			} else{
				ptMin = onia::pTRange[iRap][iPT-1];
				ptMax = onia::pTRange[iRap][iPT];
			}

			// output file name and workspace
			std::stringstream outfilename;
			outfilename << "tmpFiles/backupWorkSpace/ws_createWorkspace_Chi_rap" << iRap << "_pt" << iPT << ".root";
			RooWorkspace* ws = new RooWorkspace(workspacename.c_str());

			// define pt and y cuts on dataset
			std::stringstream cutString;
			cutString << "(chicPt >= " << ptMin << " && chicPt < "<< ptMax << ") && "
				<< "(TMath::Abs(chicRap) >= " << yMin << " && TMath::Abs(chicRap) < " << yMax << ")";

			cout << "cutString: " << cutString.str().c_str() << endl;

			// get the dataset for the fit
			RooDataSet* binData = (RooDataSet*)fullData->reduce(cutString.str().c_str());
			std::stringstream name;
			name << "data_rap" << iRap << "_pt" << iPT;
			binData->SetNameTitle(name.str().c_str(), "Data For Fitting");    

			double chicMeanPt = binData->mean(*chicPt);
		    RooRealVar var_chicMeanPt("var_chicMeanPt","var_chicMeanPt",chicMeanPt); if(!ws->var("var_chicMeanPt")) ws->import(var_chicMeanPt); else ws->var("var_chicMeanPt")->setVal(chicMeanPt);
			cout << "numEvents = " << binData->sumEntries() << endl;
			cout << "chicMeanPt = " << chicMeanPt << endl;

			// Import variables to workspace
			ws->import(*binData);
			ws->writeToFile(outfilename.str().c_str());
		}//iPT
	}//iRap

	////---------------------------------------------------------------
	////--Integrating rapidity and pt bins, in +/- 3*sigma mass window
	////---------------------------------------------------------------
	if(drawRapPt2D){
		double yMin = onia::rapForPTRange[0];
		double yMax = 1.6;//onia::rapForPTRange[onia::kNbRapForPTBins];
		double ptMin =  onia::pTRange[0][0];
		double ptMax =  onia::pTRange[0][onia::kNbPTBins[0]];

		std::stringstream cutRapPt;
		cutRapPt << "(chicPt > " << ptMin << " && chicPt < "<< ptMax << ") && "
			<< "(TMath::Abs(chicRap) > " << yMin << " && TMath::Abs(chicRap) < " << yMax << ")";
		cout<<"cutRapPt: "<<cutRapPt.str().c_str()<<endl;

		RooDataSet* rapPtData = (RooDataSet*)fullData->reduce(cutRapPt.str().c_str());
		std::stringstream nameRapPt;
		nameRapPt << "data_rap0_pt0";
		rapPtData->SetNameTitle(nameRapPt.str().c_str(), "Data For full rap and pt");

		// output file name and workspace
		std::stringstream outfilename;
		outfilename << "tmpFiles/backupWorkSpace/ws_createWorkspace_Chi_rap0_pt0.root";
		RooWorkspace* ws_RapPt = new RooWorkspace(workspacename.c_str());
		//Import variables to workspace
		ws_RapPt->import(*rapPtData);
		ws_RapPt->writeToFile(outfilename.str().c_str());

		TH2D* rapPt;
		TH1D* rap1p2;
		double MassMin;
		double MassMax;

		rap1p2 = new TH1D("rap1p2","rap1p2",30,1.2, 1.8); 
		if(nState==4){
			rapPt = new TH2D( "rapPt", "rapPt", 52,-1.3,1.3,144,0,72);
			MassMin=3.011;//massPsi1S-onia::nSigMass*sigma1S;
			MassMax=3.174;//massPsi1S+onia::nSigMass*sigma1S;
			// sigma  27.2 MeV
			// mean 3.093 GeV
		}
		if(nState==5){
			rapPt = new TH2D( "rapPt", "rapPt", 64,-1.6,1.6,144,0,72); //  rap<1.5
			//rapPt = new TH2D( "rapPt", "rapPt", 52,-1.3,1.3,144,0,72); //  rap<1.2
			MassMin=3.576;//massPsi2S-onia::nSigMass*sigma2S;
			MassMax=3.786;//massPsi2S+onia::nSigMass*sigma2S;
			// sigma 34.9 MeV // pT > 7
			// sigma 34.3 MeV // pT > 10
			// mean 3.681 GeV
		}

		cout<<"Plotting rap-Pt for Psi"<<nState-3<<"S"<<endl;
		cout<<"MassMin for rap-Pt plot = "<<MassMin<<endl;
		cout<<"MassMax for rap-Pt plot = "<<MassMax<<endl;

		TTree *rapPtTree = (TTree*)rapPtData->tree();
		std::stringstream cutMass;
		cutMass<<"(chicMass > " << MassMin << " && chicMass < "<< MassMax << ")";
		//following two methods can only be used in root_v30, 34 does not work
		rapPtTree->Draw("chicPt:chicRap>>rapPt",cutMass.str().c_str(),"colz");
		cout<<"debug"<<endl;
		rapPtTree->Draw("TMath::Abs(chicRap)>>rap1p2",cutMass.str().c_str());

		TCanvas* c2 = new TCanvas("c2","c2",1200,1500);
		rapPt->SetYTitle("p_{T}(#mu#mu) [GeV]");
		rapPt->SetXTitle("y(#mu#mu)");
		gStyle->SetPalette(1);
		gPad->SetFillColor(kWhite);
		rapPt->SetTitle(0);
		rapPt->SetStats(0);
		gPad->SetLeftMargin(0.15);
		gPad->SetRightMargin(0.17);
		rapPt->GetYaxis()->SetTitleOffset(1.5);

		rapPt->Draw("colz");

		TLine* rapPtLine;

		for(int iRap=0;iRap<onia::kNbRapForPTBins+1;iRap++){
			rapPtLine= new TLine( -onia::rapForPTRange[iRap], onia::pTRange[0][0], -onia::rapForPTRange[iRap], onia::pTRange[0][onia::kNbPTBins[iRap]] );
			rapPtLine->SetLineWidth( 2 );
			rapPtLine->SetLineStyle( 1 );
			rapPtLine->SetLineColor( kWhite );
			rapPtLine->Draw();
			rapPtLine= new TLine( onia::rapForPTRange[iRap], onia::pTRange[0][0], onia::rapForPTRange[iRap], onia::pTRange[0][onia::kNbPTBins[iRap]] );
			rapPtLine->SetLineWidth( 2 );
			rapPtLine->SetLineStyle( 1 );
			rapPtLine->SetLineColor( kWhite );
			rapPtLine->Draw();

			int pTBegin = 0;
			if(nState==5) pTBegin = 1;
			for(int iPt=pTBegin;iPt<onia::kNbPTBins[iRap+1]+1;iPt++){
				rapPtLine= new TLine( -onia::rapForPTRange[onia::kNbRapForPTBins], onia::pTRange[0][iPt], onia::rapForPTRange[onia::kNbRapForPTBins], onia::pTRange[0][iPt] );
				rapPtLine->SetLineWidth( 2 );
				rapPtLine->SetLineStyle( 1 );
				rapPtLine->SetLineColor( kWhite );
				rapPtLine->Draw();
			}
		}

		char savename[200];
		sprintf(savename,"Fit/rapPt_Chi.pdf");
		c2->SaveAs(savename);

		TCanvas* c3 = new TCanvas("c3","c3",1500,1200);
		rap1p2->SetYTitle("Events");
		rap1p2->SetXTitle("y(#mu#mu)");
		rap1p2->SetTitle(0);
		rap1p2->SetStats(0);
		rap1p2->GetYaxis()->SetTitleOffset(1.2);
		rap1p2->Draw();
		sprintf(savename,"Fit/rap_Chi_1p2.pdf");
		c3->SaveAs(savename);
	}

	f->Close();
}
