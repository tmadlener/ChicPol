#define muonTests_cxx
#include "muonTests.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>



const int kNbRapForPTBins = 2;
double rapForPTRange[kNbRapForPTBins+1] = {0., 0.6, 1.2};
const int kNbPTMaxBins = 5;
const int kNbPTBins[kNbRapForPTBins+1] = {kNbPTMaxBins, kNbPTMaxBins, kNbPTMaxBins};
double pTRange[kNbRapForPTBins+1][kNbPTMaxBins+1] = {
  {10., 15., 20., 25., 30., 50.},//all rapidities
  {10., 15., 20., 25., 30., 50.},//mid-rapidity
  {10., 15., 20., 25., 30., 50.}};//forward rapidity


void muonTests::Loop()
{

  TCanvas *c1 = new TCanvas("c", "c");

  // TFile* outfile = new TFile("muonTests_output.root", "RECREATE");

  TH2F *pt1vspt2[kNbRapForPTBins+1][kNbPTMaxBins+1];
  TH2F *pt1vseta1[kNbRapForPTBins+1][kNbPTMaxBins+1];
  TH2F *pt2vseta2[kNbRapForPTBins+1][kNbPTMaxBins+1];
  TH2F* dPhivsdPt[kNbRapForPTBins+1][kNbPTMaxBins+1];
  TH2F* dPhivsdEta[kNbRapForPTBins+1][kNbPTMaxBins+1];
  for (int rap=0;rap<kNbRapForPTBins+1;rap++) {
    for (int pt=0;pt<kNbPTMaxBins+1;pt++){

      pt1vspt2[rap][pt] = new TH2F("ptNvsptP","negMu pT vs posMu pT",1000,0,40,1000,0,40);
      pt1vseta1[rap][pt] = new TH2F("ptNvsetaN","negMu pT vs negMu eta",1000,0,40,320,-2.4,2.4);
      pt2vseta2[rap][pt] = new TH2F("ptPvsetaP","posMu pT vs posMu eta",1000,0,40,320,-2.4,2.4);
      dPhivsdPt[rap][pt] = new TH2F("dPhiVSdPt", "delta phi vs delta pT", 400, -3.3, 3.3, 1000, 0, 40);
      dPhivsdEta[rap][pt] = new TH2F("dPhiVSdEta", "delta phi vs delta eta", 400, -3.3, 3.3, 500, -5, 5);
    }
  }

  bool rejectCowboys;

  int trigDecision = -99;

  //   In a ROOT session, you can do:
  //      Root > .L muonTests.C
  //      Root > muonTests t
  //      Root > t.GetEntry(12); // Fill t data members with entry number 12
  //      Root > t.Show();       // Show values of entry 12
  //      Root > t.Show(16);     // Read and show values of entry 16
  //      Root > t.Loop();       // Loop on all entries
  //

  //     This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  //by  b_branchname->GetEntry(ientry); //read only this branch


  if (fChain == 0) return;


  gStyle->SetOptStat(0);

  Long64_t nentries = fChain->GetEntriesFast();

  std::cout << "number of entries = " << nentries << std::endl;

  nentries = 1000000;

  Long64_t nbytes = 0, nb = 0;

  for(int rejCows = 0; rejCows < 2; rejCows++)
    {
      if (rejCows == 0) rejectCowboys=false;
      if (rejCows == 1) rejectCowboys=true;

      for (Long64_t jentry=0; jentry<nentries;jentry++) {

        Long64_t ientry = LoadTree(jentry);

        if (ientry < 0) break;

        nb = fChain->GetEntry(jentry);   nbytes += nb;

	// if( HLT_Dimuon8_Jpsi_v3 == 1  ||
        //     HLT_Dimuon8_Jpsi_v4 == 1  ||
        //     HLT_Dimuon8_Jpsi_v5 == 1  ||
        //     HLT_Dimuon8_Jpsi_v6 == 1  ||
        //     HLT_Dimuon8_Jpsi_v7 == 1  ||
        //     HLT_Dimuon10_Jpsi_v3 == 1 ||
        //     HLT_Dimuon10_Jpsi_v4 == 1 ||
        //     HLT_Dimuon10_Jpsi_v5 == 1 ||
        //     HLT_Dimuon10_Jpsi_v6 == 1 )
        // replace 9 possible comparisons with 9 cheap additions and one comparison
        if (HLT_Dimuon8_Jpsi_v3 + HLT_Dimuon8_Jpsi_v4 + HLT_Dimuon8_Jpsi_v5 + HLT_Dimuon8_Jpsi_v6 + HLT_Dimuon8_Jpsi_v7
            + HLT_Dimuon10_Jpsi_v3 + HLT_Dimuon10_Jpsi_v4 + HLT_Dimuon10_Jpsi_v5 + HLT_Dimuon10_Jpsi_v6)
          trigDecision = 1;

	if(trigDecision != 1 ) continue;
	if(onia->Pt() > 990.) continue;


	double deltaPhi = muNeg->Phi() - muPos->Phi();
	if(deltaPhi > TMath::Pi()) deltaPhi -= 2.*TMath::Pi();
	else if(deltaPhi < -TMath::Pi()) deltaPhi += 2.*TMath::Pi();

        // fill the delta phi histograms before rejecting anything
        // also do this only for one run (in this case the second)
        if (rejCows){
          dPhivsdPt[0][0]->Fill(deltaPhi, muNeg->Pt() - muPos->Pt());
          dPhivsdEta[0][0]->Fill(deltaPhi, muNeg->Rapidity() - muPos->Rapidity());
          for (int iRap = 1; iRap < kNbRapForPTBins + 1; ++iRap) {
            for (int iPt = 1; iPt < kNbPTMaxBins + 1; ++iPt) {
              double pT = onia->Pt();
              double absY = TMath::Abs(onia->Rapidity());
              if (pT > pTRange[iRap-1][iPt-1] && pT < pTRange[iRap][iPt] &&
                  absY > rapForPTRange[iRap-1] && absY < rapForPTRange[iRap]) {
                dPhivsdPt[iRap][iPt]->Fill(deltaPhi, muNeg->Pt() - muPos->Pt());
                dPhivsdEta[iRap][iPt]->Fill(deltaPhi, muNeg->Rapidity() - muPos->Rapidity());
              }
            }
          }
        }

	if(rejectCowboys){
          if(deltaPhi < 0.)  continue;
	}

	if(!rejectCowboys){
          if(deltaPhi > 0.)  continue;
	}


        pt1vspt2[0][0]->Fill(muPos->Pt(),muNeg->Pt());
        pt1vseta1[0][0]->Fill(muNeg->Pt(),muNeg->Rapidity());
        pt2vseta2[0][0]->Fill(muPos->Pt(),muPos->Rapidity());
        for(int rapbins = 1; rapbins < kNbRapForPTBins+1; rapbins++)
          {
            for(int ptbins = 1; ptbins < kNbPTMaxBins+1; ptbins++)
              {
                if(onia->Pt()>pTRange[rapbins-1][ptbins-1] && onia->Pt()<pTRange[rapbins][ptbins])
                  {
                    if(TMath::Abs(onia->Rapidity())>rapForPTRange[rapbins-1] && TMath::Abs(onia->Rapidity())<rapForPTRange[rapbins]){
                      pt1vspt2[rapbins][ptbins]->Fill(muPos->Pt(),muNeg->Pt());
                      pt1vseta1[rapbins][ptbins]->Fill(muNeg->Pt(),muNeg->Rapidity());
                      pt2vseta2[rapbins][ptbins]->Fill(muPos->Pt(),muPos->Rapidity());
                    }
                  }

              }
          }

      }

      char pdfname[200];
      // save the delta phi histos, again only necessary once
      if (rejCows) {
        for (int i = 0; i < kNbRapForPTBins + 1; ++i) {
          for (int j = 0; j < kNbPTMaxBins + 1; ++j) {
            dPhivsdPt[i][j]->GetXaxis()->SetTitle("#phi(#mu^{-}) - #phi(#mu^{-})");
            dPhivsdPt[i][j]->GetYaxis()->SetTitle("p_{T}(#mu^{-}) - p_{T}(#mu^{+})");
            dPhivsdEta[i][j]->GetXaxis()->SetTitle("#phi(#mu^{-}) - #phi(#mu^{-})");
            dPhivsdEta[i][j]->GetYaxis()->SetTitle("y(#mu^{-}) - y(#mu^{+})");

            dPhivsdPt[i][j]->Draw("colz");
            if (!i && !j) {
              c1->SaveAs("global_dPhi_vs_dPt.pdf");
            }
            if (i && j) { // do not draw rap or pt summarized plots
              sprintf(pdfname, "dPhi_vs_dPt_rap%d_pt%d.pdf", i, j);
              c1->SaveAs(pdfname);
            }

            dPhivsdEta[i][j]->Draw("colz");
            if (!i && !j) {
              c1->SaveAs("global_dPhi_vs_dEta.pdf");
            }
            if (i && j) {
              sprintf(pdfname, "dPhi_vs_dEta_rap%d_pt%d.pdf", i, j);
              c1->SaveAs(pdfname);
            }
          }
        }
      }

      pt1vspt2[0][0]->GetXaxis()->SetTitle("posMu p_{T}");
      pt1vspt2[0][0]->GetYaxis()->SetTitle("negMu p_{T}");
      pt1vspt2[0][0]->Draw("colz");
      if(rejectCowboys) {c1->SaveAs("global_pTN_vs_pTP_Seagulls.pdf");}
      else {c1->SaveAs("global_pTN_vs_pTP_Cowboys.pdf");}


      pt1vseta1[0][0]->GetXaxis()->SetTitle("negMu p_{T}");
      pt1vseta1[0][0]->GetYaxis()->SetTitle("negMu eta");
      pt1vseta1[0][0]->Draw("colz");
      if(rejectCowboys) {c1->SaveAs("global_pTN_vs_etaN_Seagulls.pdf");}
      else {c1->SaveAs("global_pTN_vs_etaN_Cowboys.pdf");}


      pt2vseta2[0][0]->GetXaxis()->SetTitle("posMu p_{T}");
      pt2vseta2[0][0]->GetYaxis()->SetTitle("posMu eta");
      pt2vseta2[0][0]->Draw("colz");
      if(rejectCowboys) {c1->SaveAs("global_pTP_vs_etaP_Seagulls.pdf");}
      else {c1->SaveAs("global_pTP_vs_etaP_Cowboys.pdf");}

      for(int rapbins = 1; rapbins < kNbRapForPTBins+1; rapbins++)
        {
          for(int ptbins = 1; ptbins < kNbPTMaxBins+1; ptbins++)
            {
              pt1vspt2[rapbins][ptbins]->GetXaxis()->SetTitle("posMu p_{T}");
              pt1vspt2[rapbins][ptbins]->GetYaxis()->SetTitle("negMu p_{T}");
              pt1vspt2[rapbins][ptbins]->Draw("colz");
              if(rejectCowboys) {sprintf(pdfname,"pTN_vs_pTP_Seagulls_rap%d_pt%d.pdf",rapbins,ptbins);}
              else {sprintf(pdfname,"pTN_vs_pTP_Cowboys_rap%d_pt%d.pdf",rapbins,ptbins);}
              c1->SaveAs(pdfname);


              pt1vseta1[rapbins][ptbins]->GetXaxis()->SetTitle("negMu p_{T}");
              pt1vseta1[rapbins][ptbins]->GetYaxis()->SetTitle("negMu eta");
              pt1vseta1[rapbins][ptbins]->Draw("colz");
              if(rejectCowboys) {sprintf(pdfname,"pTN_vs_etaN_Seagulls_rap%d_pt%d.pdf",rapbins,ptbins);}
              else {sprintf(pdfname,"pTN_vs_etaN_Cowboys_rap%d_pt%d.pdf",rapbins,ptbins);}
              c1->SaveAs(pdfname);

              pt2vseta2[rapbins][ptbins]->GetXaxis()->SetTitle("posMu p_{T}");
              pt2vseta2[rapbins][ptbins]->GetYaxis()->SetTitle("posMu eta");
              pt2vseta2[rapbins][ptbins]->Draw("colz");
              if(rejectCowboys) {sprintf(pdfname,"pTP_vs_etaP_Seagulls_rap%d_pt%d.pdf",rapbins,ptbins);}
              else {sprintf(pdfname,"pTP_vs_etaP_Cowboys_rap%d_pt%d.pdf",rapbins,ptbins);}
              c1->SaveAs(pdfname);

            }
        }

    }
}
