#include "../interface/rootIncludes.inc"
#include "../interface/commonVar_Psi1S.h"
#include "../interface/effsAndCuts_Psi1S.h"

#include <string>
#include <iostream>
#include <sstream>
using namespace std;
using namespace onia;

void overlapTnPTree(){
	TGaxis::SetMaxDigits(3);
	int FidCuts = 11;

	char filename[500];
	sprintf(filename, "/afs/ihep.ac.cn/users/z/zhangll/data/data/Onia2MuMu/TTree/Onia2MuMu_v30_10May2012/TTree_Onia2MuMu_v30_PromptRecoAB_10May2012_Jpsi.root");
	TFile *f = TFile::Open(filename);
	TTree *fChain = (TTree*)f->Get("data");

	TLorentzVector  *onia, *muPos, *muNeg;
	int eventNb;
	int runNb;
	int lumiBlock;
	int nPriVtx;

	double Jpsict;
	double JpsictErr;
	double JpsiMassErr;
	double JpsiVprob;
	double JpsiDistM1;
	double JpsiDphiM1;
	double JpsiDrM1;
	double JpsiDistM2;
	double JpsiDphiM2;
	double JpsiDrM2;

	int HLT_Dimuon0_Jpsi_v1;
	int HLT_Dimuon0_Jpsi_v2;
	int HLT_Dimuon0_Jpsi_v3;
	int HLT_Dimuon0_Jpsi_v5;
	int HLT_Dimuon0_Jpsi_v6;
	int HLT_Dimuon0_Jpsi_v9;
	int HLT_Dimuon10_Jpsi_Barrel_v1;
	int HLT_Dimuon10_Jpsi_Barrel_v2;
	int HLT_Dimuon10_Jpsi_Barrel_v3;
	int HLT_Dimuon10_Jpsi_Barrel_v5;
	int HLT_Dimuon10_Jpsi_Barrel_v6;
	int HLT_Dimuon10_Jpsi_Barrel_v9;
	int HLT_Dimuon13_Jpsi_Barrel_v1;
	int HLT_Dimuon13_Jpsi_Barrel_v4;

	int HLT_Mu3_Track3_Jpsi_v4; 
	int HLT_Mu7_Track5_Jpsi_v1;
	int HLT_Mu7_Track7_Jpsi_v1;
	int HLT_Mu3_Track3_Jpsi_v5;
	int HLT_Mu5_Track2_Jpsi_v1;
	int HLT_Mu7_Track5_Jpsi_v2;
	int HLT_Mu7_Track7_Jpsi_v2;
	int HLT_Mu5_Track2_Jpsi_v2;
	int HLT_Mu7_Track7_Jpsi_v3;
	int HLT_Mu5_Track2_Jpsi_v4;
	int HLT_Mu7_Track7_Jpsi_v5;
	int HLT_Mu5_Track2_Jpsi_v5;
	int HLT_Mu7_Track7_Jpsi_v6;
	int HLT_Mu5_Track2_Jpsi_v6;
	int HLT_Mu7_Track7_Jpsi_v7;
	int HLT_Mu5_Track2_Jpsi_v8;
	int HLT_Mu7_Track7_Jpsi_v9;
	int HLT_Mu5_Track2_Jpsi_v9;
	int HLT_Mu7_Track7_Jpsi_v10;
	int HLT_Mu5_Track2_Jpsi_v12;
	int HLT_Mu7_Track7_Jpsi_v13;

	int HLT_Mu5_L2Mu2_v1;
	int HLT_Mu5_L2Mu2_Jpsi_v1;
	int HLT_Mu5_L2Mu2_v2;
	int HLT_Mu5_L2Mu2_Jpsi_v2;
	int HLT_Mu5_L2Mu2_Jpsi_v3;
	int HLT_Mu5_L2Mu2_Jpsi_v4;
	int HLT_Mu5_L2Mu2_Jpsi_v5;
	int HLT_Mu5_L2Mu2_Jpsi_v6;
	int HLT_Mu5_L2Mu2_Jpsi_v8;
	int HLT_Mu5_L2Mu2_Jpsi_v9;
	int HLT_Mu5_L2Mu2_Jpsi_v12;

	fChain->SetBranchAddress("JpsiP", &onia);
	fChain->SetBranchAddress("muNegP", &muNeg);
	fChain->SetBranchAddress("muPosP", &muPos);

	fChain->SetBranchAddress("eventNb", &eventNb);
	fChain->SetBranchAddress("runNb", &runNb);
	fChain->SetBranchAddress("lumiBlock", &lumiBlock);
	fChain->SetBranchAddress("nPriVtx", &nPriVtx);

	fChain->SetBranchAddress("Jpsict", &Jpsict);
	fChain->SetBranchAddress("JpsictErr", &JpsictErr);
	fChain->SetBranchAddress("JpsiMassErr", &JpsiMassErr);
	fChain->SetBranchAddress("JpsiVprob", &JpsiVprob);
	fChain->SetBranchAddress("JpsiDistM1", &JpsiDistM1);
	fChain->SetBranchAddress("JpsiDphiM1", &JpsiDphiM1);
	fChain->SetBranchAddress("JpsiDrM1", &JpsiDrM1);
	fChain->SetBranchAddress("JpsiDistM2", &JpsiDistM2);
	fChain->SetBranchAddress("JpsiDphiM2", &JpsiDphiM2);
	fChain->SetBranchAddress("JpsiDrM2", &JpsiDrM2);

	fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_v1", &HLT_Dimuon0_Jpsi_v1);
	fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_v2", &HLT_Dimuon0_Jpsi_v2);
	fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_v3", &HLT_Dimuon0_Jpsi_v3);
	fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_v5", &HLT_Dimuon0_Jpsi_v5);
	fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_v6", &HLT_Dimuon0_Jpsi_v6);
	fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_v9", &HLT_Dimuon0_Jpsi_v9);
	fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v1", &HLT_Dimuon10_Jpsi_Barrel_v1);
	fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v2", &HLT_Dimuon10_Jpsi_Barrel_v2);
	fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v3", &HLT_Dimuon10_Jpsi_Barrel_v3);
	fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v5", &HLT_Dimuon10_Jpsi_Barrel_v5);
	fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v6", &HLT_Dimuon10_Jpsi_Barrel_v6);
	fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v9", &HLT_Dimuon10_Jpsi_Barrel_v9);
	fChain->SetBranchAddress("HLT_Dimuon13_Jpsi_Barrel_v1", &HLT_Dimuon13_Jpsi_Barrel_v1);
	fChain->SetBranchAddress("HLT_Dimuon13_Jpsi_Barrel_v4", &HLT_Dimuon13_Jpsi_Barrel_v4);

	fChain->SetBranchAddress("HLT_Mu3_Track3_Jpsi_v4", &HLT_Mu3_Track3_Jpsi_v4);
	fChain->SetBranchAddress("HLT_Mu7_Track5_Jpsi_v1", &HLT_Mu7_Track5_Jpsi_v1);
	fChain->SetBranchAddress("HLT_Mu7_Track7_Jpsi_v1", &HLT_Mu7_Track7_Jpsi_v1);
	fChain->SetBranchAddress("HLT_Mu3_Track3_Jpsi_v5", &HLT_Mu3_Track3_Jpsi_v5);
	fChain->SetBranchAddress("HLT_Mu5_Track2_Jpsi_v1", &HLT_Mu5_Track2_Jpsi_v1);
	fChain->SetBranchAddress("HLT_Mu7_Track5_Jpsi_v2", &HLT_Mu7_Track5_Jpsi_v2);
	fChain->SetBranchAddress("HLT_Mu7_Track7_Jpsi_v2", &HLT_Mu7_Track7_Jpsi_v2);
	fChain->SetBranchAddress("HLT_Mu5_Track2_Jpsi_v2", &HLT_Mu5_Track2_Jpsi_v2);
	fChain->SetBranchAddress("HLT_Mu7_Track7_Jpsi_v3", &HLT_Mu7_Track7_Jpsi_v3);
	fChain->SetBranchAddress("HLT_Mu5_Track2_Jpsi_v4", &HLT_Mu5_Track2_Jpsi_v4);
	fChain->SetBranchAddress("HLT_Mu7_Track7_Jpsi_v5", &HLT_Mu7_Track7_Jpsi_v5);
	fChain->SetBranchAddress("HLT_Mu5_Track2_Jpsi_v5", &HLT_Mu5_Track2_Jpsi_v5);
	fChain->SetBranchAddress("HLT_Mu7_Track7_Jpsi_v6", &HLT_Mu7_Track7_Jpsi_v6);
	fChain->SetBranchAddress("HLT_Mu5_Track2_Jpsi_v6", &HLT_Mu5_Track2_Jpsi_v6);
	fChain->SetBranchAddress("HLT_Mu7_Track7_Jpsi_v7", &HLT_Mu7_Track7_Jpsi_v7);
	fChain->SetBranchAddress("HLT_Mu5_Track2_Jpsi_v8", &HLT_Mu5_Track2_Jpsi_v8);
	fChain->SetBranchAddress("HLT_Mu7_Track7_Jpsi_v9", &HLT_Mu7_Track7_Jpsi_v9);
	fChain->SetBranchAddress("HLT_Mu5_Track2_Jpsi_v9", &HLT_Mu5_Track2_Jpsi_v9);
	fChain->SetBranchAddress("HLT_Mu7_Track7_Jpsi_v10", &HLT_Mu7_Track7_Jpsi_v10);
	fChain->SetBranchAddress("HLT_Mu5_Track2_Jpsi_v12", &HLT_Mu5_Track2_Jpsi_v12);
	fChain->SetBranchAddress("HLT_Mu7_Track7_Jpsi_v13", &HLT_Mu7_Track7_Jpsi_v13);

	fChain->SetBranchAddress("HLT_Mu5_L2Mu2_v1", &HLT_Mu5_L2Mu2_v1);
	fChain->SetBranchAddress("HLT_Mu5_L2Mu2_Jpsi_v1", &HLT_Mu5_L2Mu2_Jpsi_v1);
	fChain->SetBranchAddress("HLT_Mu5_L2Mu2_v2", &HLT_Mu5_L2Mu2_v2);
	fChain->SetBranchAddress("HLT_Mu5_L2Mu2_Jpsi_v2", &HLT_Mu5_L2Mu2_Jpsi_v2);
	fChain->SetBranchAddress("HLT_Mu5_L2Mu2_Jpsi_v3", &HLT_Mu5_L2Mu2_Jpsi_v3);
	fChain->SetBranchAddress("HLT_Mu5_L2Mu2_Jpsi_v4", &HLT_Mu5_L2Mu2_Jpsi_v4);
	fChain->SetBranchAddress("HLT_Mu5_L2Mu2_Jpsi_v5", &HLT_Mu5_L2Mu2_Jpsi_v5);
	fChain->SetBranchAddress("HLT_Mu5_L2Mu2_Jpsi_v6", &HLT_Mu5_L2Mu2_Jpsi_v6);
	fChain->SetBranchAddress("HLT_Mu5_L2Mu2_Jpsi_v8", &HLT_Mu5_L2Mu2_Jpsi_v8);
	fChain->SetBranchAddress("HLT_Mu5_L2Mu2_Jpsi_v9", &HLT_Mu5_L2Mu2_Jpsi_v9);
	fChain->SetBranchAddress("HLT_Mu5_L2Mu2_Jpsi_v12", &HLT_Mu5_L2Mu2_Jpsi_v12);


	Long64_t nentries = fChain->GetEntries();
	Long64_t cutAtRecEvent = nentries;
	Long64_t countTreeTot = 0;
	Long64_t countMuTrackTot = 0;
	Long64_t countMuL2MuTot  = 0;
	Long64_t nb = 0;

	Long64_t countTree[onia::kNbRapForPTBins][kNbPTMaxBins],
					 countMuTrack[onia::kNbRapForPTBins][kNbPTMaxBins],
					 countMuL2Mu[onia::kNbRapForPTBins][kNbPTMaxBins];

	for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){
		for(int iPT = 0; iPT < onia::kNbPTMaxBins; iPT++){
			countTree[iRap][iPT] = 0;
			countMuTrack[iRap][iPT] = 0;
			countMuL2Mu[iRap][iPT] = 0;
		}
	}

	std::cout << "number of entries = " << nentries << std::endl;

	for (Long64_t jentry=0; jentry<nentries;jentry++) {

		if(jentry % 100000 == 0) std::cout << "event " << jentry << std::endl;

		//Long64_t ientry = LoadTree(jentry);
		//if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);

		if(onia->Pt() > 990.) continue;
		if(JpsiVprob < 0.01) continue;

		int trigDecision = -99;

		if( HLT_Dimuon10_Jpsi_Barrel_v1 == 1 ||
				HLT_Dimuon10_Jpsi_Barrel_v2 == 1 ||
				HLT_Dimuon10_Jpsi_Barrel_v3 == 1 ||
				HLT_Dimuon10_Jpsi_Barrel_v5 == 1 ||
				HLT_Dimuon10_Jpsi_Barrel_v6 == 1 ||
				HLT_Dimuon10_Jpsi_Barrel_v9 == 1)  
			trigDecision = 1;

		if(trigDecision != 1) continue;

		if( (HLT_Dimuon10_Jpsi_Barrel_v1 == 1 && onia->Pt() < 10) ||
				(HLT_Dimuon10_Jpsi_Barrel_v2 == 1 && onia->Pt() < 10) ||
				(HLT_Dimuon10_Jpsi_Barrel_v3 == 1 && onia->Pt() < 10) ||
				(HLT_Dimuon10_Jpsi_Barrel_v5 == 1 && onia->Pt() < 10) ||
				(HLT_Dimuon10_Jpsi_Barrel_v6 == 1 && onia->Pt() < 10) ||
				(HLT_Dimuon10_Jpsi_Barrel_v9 == 1 && onia->Pt() < 10))
			continue;

		double onia_mass = onia->M();
		double onia_pt = onia->Pt();
		double onia_P = onia->P();
		double onia_eta = onia->PseudoRapidity();
		double onia_rap = onia->Rapidity();
		double onia_phi = onia->Phi();
		double onia_mT = sqrt(onia_mass*onia_mass + onia_pt*onia_pt);

		if(TMath::Abs(onia_rap) > onia::rap) continue;

		double etaMuPos = muPos->PseudoRapidity();
		double etaMuNeg = muNeg->PseudoRapidity();
		double pTMuPos = muPos->Pt();
		double pTMuNeg = muNeg->Pt();
		double pMuPos = muPos->P();
		double pMuNeg = muNeg->P();
		double deltaPhi = muNeg->Phi() - muPos->Phi();
		if(deltaPhi > TMath::Pi()) deltaPhi -= 2.*TMath::Pi();
		else if(deltaPhi < -TMath::Pi()) deltaPhi += 2.*TMath::Pi();

		if(deltaPhi < 0.)  continue; //cowboy rejection

		// select mass window
		double massMin = 0, massMax = 0;
		massMin = onia::massMin;
		massMax = onia::massMax;
		if(onia_mass < massMin || onia_mass > massMax) continue;

		bool muonsInAcc = kFALSE;
		if(isMuonInAcceptance(FidCuts-1, pTMuPos, etaMuPos) && isMuonInAcceptance(FidCuts-1, pTMuNeg, etaMuNeg)){
			muonsInAcc = kTRUE;
		}

		if(!muonsInAcc) continue;

		if( ! ( onia_pt > onia::pTRange[0][2] && onia_pt < onia::pTRange[0][kNbPTMaxBins] && 
					TMath::Abs(onia_rap) < onia::rapForPTRange[kNbRapForPTBins] ) )  continue;

		countTreeTot++;

		for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){
			for(int iPT = 0; iPT < onia::kNbPTMaxBins; iPT++){
				// events for different pT and y bins
				if(onia_pt > onia::pTRange[0][iPT] && onia_pt < onia::pTRange[0][iPT+1] && 
						TMath::Abs(onia_rap) > onia::rapForPTRange[iRap] && TMath::Abs(onia_rap) < onia::rapForPTRange[iRap+1]){
					countTree[iRap][iPT]++;
				}
			}
		}

		//// event in T&P efficiency triggers

		if( HLT_Mu3_Track3_Jpsi_v4 != 0 ||  
				HLT_Mu7_Track5_Jpsi_v1 != 0 || 
				HLT_Mu7_Track7_Jpsi_v1 != 0 || 
				HLT_Mu3_Track3_Jpsi_v5 != 0 || 
				HLT_Mu5_Track2_Jpsi_v1 != 0 || 
				HLT_Mu7_Track5_Jpsi_v2 != 0 || 
				HLT_Mu7_Track7_Jpsi_v2 != 0 || 
				HLT_Mu5_Track2_Jpsi_v2 != 0 || 
				HLT_Mu7_Track7_Jpsi_v3 != 0 || 
				HLT_Mu5_Track2_Jpsi_v4 != 0 || 
				HLT_Mu7_Track7_Jpsi_v5 != 0 || 
				HLT_Mu5_Track2_Jpsi_v5 != 0 || 
				HLT_Mu7_Track7_Jpsi_v6 != 0 || 
				HLT_Mu5_Track2_Jpsi_v6 != 0 || 
				HLT_Mu7_Track7_Jpsi_v7 != 0 || 
				HLT_Mu5_Track2_Jpsi_v8 != 0 || 
				HLT_Mu7_Track7_Jpsi_v9 != 0 || 
				HLT_Mu5_Track2_Jpsi_v9 != 0 || 
				HLT_Mu7_Track7_Jpsi_v10 != 0 || 
				HLT_Mu5_Track2_Jpsi_v12 != 0 ||
				HLT_Mu7_Track7_Jpsi_v13 != 0 ) {

			countMuTrackTot ++ ;
			for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){
				for(int iPT = 0; iPT < onia::kNbPTMaxBins; iPT++){
					// events for different pT and y bins
					if(onia_pt > onia::pTRange[0][iPT] && onia_pt < onia::pTRange[0][iPT+1] && 
							TMath::Abs(onia_rap) > onia::rapForPTRange[iRap] && TMath::Abs(onia_rap) < onia::rapForPTRange[iRap+1]){
						countMuTrack[iRap][iPT]++;
					}
				}
			}


		}


		//if( HLT_Mu5_L2Mu2_v1 == 1 ||
		//		HLT_Mu5_L2Mu2_Jpsi_v1 == 1 || 
		//		HLT_Mu5_L2Mu2_v2 == 1 ||
		//		HLT_Mu5_L2Mu2_Jpsi_v2 == 1 ||  
		//		HLT_Mu5_L2Mu2_Jpsi_v3 == 1 || 
		//		HLT_Mu5_L2Mu2_Jpsi_v4 == 1 || 
		//		HLT_Mu5_L2Mu2_Jpsi_v5 == 1 || 
		//		HLT_Mu5_L2Mu2_Jpsi_v6 == 1 || 
		//		HLT_Mu5_L2Mu2_Jpsi_v8 == 1 || 
		//		HLT_Mu5_L2Mu2_Jpsi_v9 == 1 || 
		//		HLT_Mu5_L2Mu2_Jpsi_v12 == 1 ){
		if( HLT_Mu5_L2Mu2_v1 != 0 ||
				HLT_Mu5_L2Mu2_Jpsi_v1 != 0 || 
				HLT_Mu5_L2Mu2_v2 != 0 ||
				HLT_Mu5_L2Mu2_Jpsi_v2 != 0 ||  
				HLT_Mu5_L2Mu2_Jpsi_v3 != 0 || 
				HLT_Mu5_L2Mu2_Jpsi_v4 != 0 || 
				HLT_Mu5_L2Mu2_Jpsi_v5 != 0 || 
				HLT_Mu5_L2Mu2_Jpsi_v6 != 0 || 
				HLT_Mu5_L2Mu2_Jpsi_v8 != 0 || 
				HLT_Mu5_L2Mu2_Jpsi_v9 != 0 || 
				HLT_Mu5_L2Mu2_Jpsi_v12 != 0 ){

			countMuL2MuTot ++ ;

			for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){
				for(int iPT = 0; iPT < onia::kNbPTMaxBins; iPT++){
					// events for different pT and y bins
					if(onia_pt > onia::pTRange[0][iPT] && onia_pt < onia::pTRange[0][iPT+1] && 
							TMath::Abs(onia_rap) > onia::rapForPTRange[iRap] && TMath::Abs(onia_rap) < onia::rapForPTRange[iRap+1]){
						countMuL2Mu[iRap][iPT]++;
					}
				}
			}


		}

	}



	cout << "countTreeTot:    " << countTreeTot << endl; 
	cout << "countMuTrackTot: " << countMuTrackTot << endl; 
	cout << "countMuL2MuTot: " << countMuL2MuTot << endl; 
	cout << "ratio(MuTrack):  " << (double)countMuTrackTot / countTreeTot << endl;
	cout << "ratio(MuL2Mu):  " << (double)countMuL2MuTot / countTreeTot << endl;

	for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){
		  for(int iPT = 0; iPT < onia::kNbPTMaxBins; iPT++){
				if(iRap==0 && iPT==0) cout<< "countTree countMuTrack ratio " <<endl;
				cout<< countTree[iRap][iPT] << " " <<  countMuTrack[iRap][iPT] << " " << 
					(double)countMuTrack[iRap][iPT] / countTree[iRap][iPT]  << endl;
			}
	}

	for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){
		  for(int iPT = 0; iPT < onia::kNbPTMaxBins; iPT++){
				if(iRap==0 && iPT==0) cout<< "countTree countMuL2Mu ratio " <<endl;
				cout<< countTree[iRap][iPT] << " " <<  countMuL2Mu[iRap][iPT] << " " <<
					(double)countMuL2Mu[iRap][iPT] / countTree[iRap][iPT]  << endl;
			}
	}


}
