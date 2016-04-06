/**
 * PlotCosThetaDistribution.cc
 *
 *  Created on: Dec 8, 2011
 *      Author: valentinknuenz
 *  Modified on: Jan 22, 2013 linlinzhang
 *  Modified on: Sept 16, 2015 ilse
 */

#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include "rootIncludes.inc"
#include "commonVar.h"
#include "clarg_parsing.h"


// binning
int const pT_binsFOR1Dhists=onia::kNbPTMaxBins;
int const Rap_binsFOR1Dhists=onia::kNbRapForPTBins;
// extremes and binning of lambda_gen extraction histos
const double l_min = -1;
const double l_max =  1;
const double l_step_1D = 0.02;

// cosTheta - phi - Distributions
TH2D* CosThPhDist[onia::kNbFrames][onia::kNbRapForPTBins+1][onia::kNbPTMaxBins+1][2];
TH1D* CosThDist[onia::kNbFrames][onia::kNbRapForPTBins+1][onia::kNbPTMaxBins+1][2];
TH1D* PhiDist[onia::kNbFrames][onia::kNbRapForPTBins+1][onia::kNbPTMaxBins+1][2];

// extraction histos
TH1D* h_costh2[onia::kNbFrames][onia::kNbRapForPTBins+1][onia::kNbPTMaxBins+1][2];
TH1D* h_cos2ph[onia::kNbFrames][onia::kNbRapForPTBins+1][onia::kNbPTMaxBins+1][2];
TH1D* h_sin2thcosph[onia::kNbFrames][onia::kNbRapForPTBins+1][onia::kNbPTMaxBins+1][2];

// results
TGraphAsymmErrors* graph_lamth[onia::kNbFrames][onia::kNbRapForPTBins+1][2];
TGraphAsymmErrors* graph_lamph[onia::kNbFrames][onia::kNbRapForPTBins+1][2];
TGraphAsymmErrors* graph_lamtp[onia::kNbFrames][onia::kNbRapForPTBins+1][2];
TGraphAsymmErrors* graph_lamthstar[onia::kNbFrames][onia::kNbRapForPTBins+1][2];
TGraphAsymmErrors* graph_lamphstar[onia::kNbFrames][onia::kNbRapForPTBins+1][2];
TGraphAsymmErrors* graph_lamtilde[onia::kNbFrames][onia::kNbRapForPTBins+1][2];

void LoadHistos(int iRapBin, int iPTBin, int nState, const std::string& dataPath);
void PlotHistos(int iRapBin, int iPTBin, int iFrame, int nState, const std::string& DataPath);

//===========================

int main(int argc, char* argv[]){
  // set default values
  int nState = 999;
  std::string dataPath;

  for( int i=1;i < argc; i++ ) {
    std::string arg = argv[i];
    fromSplit("nState", arg, nState);
    fromSplit("DataPath", arg, dataPath);
  }

  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);

  std::cout << dataPath << std::endl;

  // initialization of all graphs
  for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){
    for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
      for (int i = 0; i < 2; i++) {
        graph_lamth[iFrame][iRap][i] = new TGraphAsymmErrors();
        graph_lamph[iFrame][iRap][i] = new TGraphAsymmErrors();
        graph_lamtp[iFrame][iRap][i] = new TGraphAsymmErrors();
        graph_lamthstar[iFrame][iRap][i] = new TGraphAsymmErrors();
        graph_lamphstar[iFrame][iRap][i] = new TGraphAsymmErrors();
        graph_lamtilde[iFrame][iRap][i] = new TGraphAsymmErrors();
      }

    }

    for(int iPT = 0; iPT < onia::kNbPTBins[iRap+1]; iPT++){
      for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
        std::stringstream name, name1;
        if(nState == 6){
          name << dataPath.c_str() << "/tmpFiles/data_chic1_rap" << iRap+1 << "_pT" << iPT+1 << ".root";
          name1 << dataPath.c_str() << "/tmpFiles/data_chic2_rap" << iRap+1 << "_pT" << iPT+1 << ".root";
        } else {
          name << dataPath.c_str() << "/tmpFiles/data_Psi" << nState-3 << "_rap" << iRap+1 << "_pT" << iPT+1 << ".root";
          name1 << dataPath.c_str() << "/tmpFiles/data_Psi" << nState-3 << "_rap" << iRap+1 << "_pT" << iPT+1 << ".root";
        }

        TFile *fIn = new TFile(name.str().c_str());
        TFile *fIn1 = new TFile(name1.str().c_str());
        TH2D *hNPBG_cosThetaPhi_PHX = (TH2D*)fIn->Get("hNPBG_cosThetaPhi_PHX");
        TH2D *hNPBG1_cosThetaPhi_PHX = (TH2D*)fIn1->Get("hNPBG_cosThetaPhi_PHX");
        int nBinsCT1 = hNPBG_cosThetaPhi_PHX->GetNbinsX();
        int nBinsPH1 = hNPBG_cosThetaPhi_PHX->GetNbinsY();
        int nBinsCT2 = hNPBG1_cosThetaPhi_PHX->GetNbinsX();
        int nBinsPH2 = hNPBG1_cosThetaPhi_PHX->GetNbinsY();

        CosThPhDist[iFrame][iRap][iPT][0] = new TH2D("","",nBinsCT1, onia::cosTMin, onia::cosTMax, nBinsPH1, onia::phiPolMin, onia::phiPolMax);
        CosThPhDist[iFrame][iRap][iPT][1] = new TH2D("","",nBinsCT2, onia::cosTMin, onia::cosTMax, nBinsPH2, onia::phiPolMin, onia::phiPolMax);
        CosThDist[iFrame][iRap][iPT][0] = new TH1D("","",nBinsCT1, onia::cosTMin, onia::cosTMax);
        CosThDist[iFrame][iRap][iPT][1] = new TH1D("","",nBinsCT2, onia::cosTMin, onia::cosTMax);
        PhiDist[iFrame][iRap][iPT][0] = new TH1D("","",nBinsPH1, onia::phiPolMin, onia::phiPolMax);
        PhiDist[iFrame][iRap][iPT][1] = new TH1D("","",nBinsPH2, onia::phiPolMin, onia::phiPolMax);

        for (int i = 0; i < 2; i++){
          h_costh2[iFrame][iRap][iPT][i] = new TH1D( "", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
          h_cos2ph[iFrame][iRap][iPT][i] = new TH1D( "", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
          h_sin2thcosph[iFrame][iRap][iPT][i] = new TH1D( "", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
        }

      } //iFrame

      LoadHistos(iRap, iPT,nState,dataPath);

      for(int iFrame = 0; iFrame < 3; iFrame++){
        PlotHistos(iRap, iPT, iFrame,nState,dataPath);
      }

    } //iPT
  } //iRap

  std::stringstream filename, filename1, angFile, angFile1;
  if(nState == 6){
    filename << dataPath.c_str() << "/Figures/TGraphResults_gen_chic1.root";
    filename1 << dataPath.c_str() << "/Figures/TGraphResults_gen_chic2.root";
    angFile << dataPath.c_str() << "/Figures/AngDistHist_chic1.root";
    angFile1 << dataPath.c_str() << "/Figures/AngDistHist_chic2.root";
  } else {
    filename << dataPath.c_str() << "/Figures/TGraphResults_gen_Psi" << nState-3 << ".root";
    angFile << dataPath.c_str() << "/Figures/AngDistHist.root";
  }

  TFile* results = new TFile(filename.str().c_str(), "RECREATE");
  TFile* AngDistHistFile = new TFile(angFile.str().c_str(), "RECREATE", "AngDistHistFile");
  TFile *results1 = new TFile(filename1.str().c_str(), "RECREATE");
  TFile *AngDistHistFile1 = new TFile(angFile1.str().c_str(), "RECREATE", "AngDistHistFile");

  for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){
    for(int iFrame = 0; iFrame < 3; iFrame++){
      for(int iPT = 0; iPT < onia::kNbPTBins[iRap+1]; iPT++){

        // write histograms to file
        AngDistHistFile->cd();
        std::stringstream hct, hphi, hctphi, hct2, hc2p, hs2tcp;
        hct << "Proj_" << onia::frameLabel[iFrame] << "_costh_rap" << iRap+1 << "_pT" << iPT+1;
        CosThDist[iFrame][iRap][iPT][0]->SetName(hct.str().c_str());
        CosThDist[iFrame][iRap][iPT][0]->Write();
        hphi << "Proj_" << onia::frameLabel[iFrame] << "_phi_rap" << iRap+1 << "_pT" << iPT+1;
        PhiDist[iFrame][iRap][iPT][0]->SetName(hphi.str().c_str());
        PhiDist[iFrame][iRap][iPT][0]->Write();
        hctphi << "CosthPhi_" << onia::frameLabel[iFrame] << "_rap" << iRap+1 << "_pT" << iPT+1;
        CosThPhDist[iFrame][iRap][iPT][0]->SetName(hctphi.str().c_str());
        CosThPhDist[iFrame][iRap][iPT][0]->Write();
        hct2 << "costh2_" << onia::frameLabel[iFrame] << "_rap" << iRap+1 << "_pT" << iPT+1;
        h_costh2[iFrame][iRap][iPT][0]->SetName(hct2.str().c_str());
        h_costh2[iFrame][iRap][iPT][0]->Write();
        hc2p << "cos2ph_" << onia::frameLabel[iFrame] << "_rap" << iRap+1 << "_pT" << iPT+1;
        h_cos2ph[iFrame][iRap][iPT][0]->SetName(hc2p.str().c_str());
        h_cos2ph[iFrame][iRap][iPT][0]->Write();
        hs2tcp << "sin2thcosph_" << onia::frameLabel[iFrame] << "_rap" << iRap+1 << "_pT" << iPT+1;
        h_sin2thcosph[iFrame][iRap][iPT][0]->SetName(hs2tcp.str().c_str());
        h_sin2thcosph[iFrame][iRap][iPT][0]->Write();

        // calculate the polarization parameters
        double costh2 = h_costh2[iFrame][iRap][iPT][0]->GetMean();
        double lamth = (1. - 3. * costh2 ) / ( costh2 - 3./5. );
        double cos2ph = h_cos2ph[iFrame][iRap][iPT][0]->GetMean();
        double lamph = cos2ph * (3. + lamth);
        double sin2thcosph = h_sin2thcosph[iFrame][iRap][iPT][0]->GetMean();
        double lamtp = sin2thcosph * 5./4. * (3. + lamth);
        double lamtilde = (lamth + 3*lamph)/(1-lamph);
        double lamthstar = -9;
        double lamphstar = -9;

        // calculate mean pt of bin
        double meanpt = onia::pTRange[iRap+1][iPT] + (onia::pTRange[iRap+1][iPT+1] - onia::pTRange[iRap+1][iPT])/2;

        // save to TGraph
        graph_lamth[iFrame][iRap][0]->SetPoint(iPT, meanpt, lamth);
        graph_lamph[iFrame][iRap][0]->SetPoint(iPT, meanpt, lamph);
        graph_lamtp[iFrame][iRap][0]->SetPoint(iPT, meanpt, lamtp);
        graph_lamtilde[iFrame][iRap][0]->SetPoint(iPT, meanpt, lamtilde);
        graph_lamthstar[iFrame][iRap][0]->SetPoint(iPT, meanpt, lamthstar);
        graph_lamphstar[iFrame][iRap][0]->SetPoint(iPT, meanpt, lamphstar);

        // chic2
        if(nState == 6){
          AngDistHistFile1->cd();
          CosThDist[iFrame][iRap][iPT][1]->SetName(hct.str().c_str());
          CosThDist[iFrame][iRap][iPT][1]->Write();
          PhiDist[iFrame][iRap][iPT][1]->SetName(hphi.str().c_str());
          PhiDist[iFrame][iRap][iPT][1]->Write();
          CosThPhDist[iFrame][iRap][iPT][1]->SetName(hctphi.str().c_str());
          CosThPhDist[iFrame][iRap][iPT][1]->Write();
          h_costh2[iFrame][iRap][iPT][1]->SetName(hct2.str().c_str());
          h_costh2[iFrame][iRap][iPT][1]->Write();
          h_cos2ph[iFrame][iRap][iPT][1]->SetName(hc2p.str().c_str());
          h_cos2ph[iFrame][iRap][iPT][1]->Write();
          h_sin2thcosph[iFrame][iRap][iPT][1]->SetName(hs2tcp.str().c_str());
          h_sin2thcosph[iFrame][iRap][iPT][1]->Write();

          costh2 = h_costh2[iFrame][iRap][iPT][1]->GetMean();
          lamth = (1. - 3. * costh2 ) / ( costh2 - 3./5. );
          cos2ph = h_cos2ph[iFrame][iRap][iPT][1]->GetMean();
          lamph = cos2ph * (3. + lamth);
          sin2thcosph = h_sin2thcosph[iFrame][iRap][iPT][1]->GetMean();
          lamtp = sin2thcosph * 5./4. * (3. + lamth);
          lamtilde = (lamth + 3*lamph)/(1-lamph);
          lamthstar = -9;
          lamphstar = -9;
          meanpt = onia::pTRange[iRap+1][iPT] + (onia::pTRange[iRap+1][iPT+1] - onia::pTRange[iRap+1][iPT])/2;

          graph_lamth[iFrame][iRap][1]->SetPoint(iPT, meanpt, lamth);
          graph_lamph[iFrame][iRap][1]->SetPoint(iPT, meanpt, lamph);
          graph_lamtp[iFrame][iRap][1]->SetPoint(iPT, meanpt, lamtp);
          graph_lamtilde[iFrame][iRap][1]->SetPoint(iPT, meanpt, lamtilde);
          graph_lamthstar[iFrame][iRap][1]->SetPoint(iPT, meanpt, lamthstar);
          graph_lamphstar[iFrame][iRap][1]->SetPoint(iPT, meanpt, lamphstar);
        } // if nState

      } // iPT

      results->cd();
      std::stringstream lth, lph, ltp, ltilde, lthstar, lphstar;
      if(iFrame == 2){
        lth << "lth_PX_rap" << iRap+1;
        lph << "lph_PX_rap" << iRap+1;
        ltp << "ltp_PX_rap" << iRap+1;
        ltilde << "ltilde_PX_rap" << iRap+1;
        lthstar << "lthstar_PX_rap" << iRap+1;
        lphstar << "lphstar_PX_rap" << iRap+1;
      } else {
        lth << "lth_" << onia::frameLabel[iFrame] << "_rap" << iRap+1;
        lph << "lph_" << onia::frameLabel[iFrame] << "_rap" << iRap+1;
        ltp << "ltp_" << onia::frameLabel[iFrame] << "_rap" << iRap+1;
        ltilde << "ltilde_" << onia::frameLabel[iFrame] << "_rap" << iRap+1;
        lthstar << "lthstar_" << onia::frameLabel[iFrame] << "_rap" << iRap+1;
        lphstar << "lphstar_" << onia::frameLabel[iFrame] << "_rap" << iRap+1;
      }
      graph_lamth[iFrame][iRap][0]->SetName(lth.str().c_str());
      graph_lamth[iFrame][iRap][0]->Write();
      graph_lamph[iFrame][iRap][0]->SetName(lph.str().c_str());
      graph_lamph[iFrame][iRap][0]->Write();
      graph_lamtp[iFrame][iRap][0]->SetName(ltp.str().c_str());
      graph_lamtp[iFrame][iRap][0]->Write();
      graph_lamtilde[iFrame][iRap][0]->SetName(ltilde.str().c_str());
      graph_lamtilde[iFrame][iRap][0]->Write();
      graph_lamthstar[iFrame][iRap][0]->SetName(lthstar.str().c_str());
      graph_lamthstar[iFrame][iRap][0]->Write();
      graph_lamphstar[iFrame][iRap][0]->SetName(lphstar.str().c_str());
      graph_lamphstar[iFrame][iRap][0]->Write();

      // chic2
      if(nState == 6){
        results1->cd();
        graph_lamth[iFrame][iRap][1]->SetName(lth.str().c_str());
        graph_lamth[iFrame][iRap][1]->Write();
        graph_lamph[iFrame][iRap][1]->SetName(lph.str().c_str());
        graph_lamph[iFrame][iRap][1]->Write();
        graph_lamtp[iFrame][iRap][1]->SetName(ltp.str().c_str());
        graph_lamtp[iFrame][iRap][1]->Write();
        graph_lamtilde[iFrame][iRap][1]->SetName(ltilde.str().c_str());
        graph_lamtilde[iFrame][iRap][1]->Write();
        graph_lamthstar[iFrame][iRap][1]->SetName(lthstar.str().c_str());
        graph_lamthstar[iFrame][iRap][1]->Write();
        graph_lamphstar[iFrame][iRap][1]->SetName(lphstar.str().c_str());
        graph_lamphstar[iFrame][iRap][1]->Write();
      }
    } // iFrame
  } // iRap

  results->Write();
  results->Close();
  AngDistHistFile->Write();
  AngDistHistFile->Close();
  if(nState == 6){
    results1->Write();
    AngDistHistFile1->Write();
  }
  results1->Close();
  AngDistHistFile1->Close();
  return 0;
}
//===========================
void PlotHistos(int iRapBin, int iPTBin, int iFrame,int nState, const std::string& dataPath){

  TGaxis::SetMaxDigits(3);

  double lvalue = 0.28, tvalue = 0.92;
  double left=lvalue, top=tvalue, textSize=0.035;
  TLatex *latex=new TLatex();
  latex->SetTextFont(42);
  latex->SetNDC(kTRUE);
  latex->SetTextSize(textSize);
  double step=textSize*1.3;

  TCanvas *c1 = new TCanvas("", "", 500, 500);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gPad->SetFillColor(kWhite);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);

  std::stringstream XTitle, YTitle;
  XTitle << "cos#theta_" << onia::frameLabel[iFrame];
  YTitle << "#phi_" << onia::frameLabel[iFrame] << " [deg]";
  double yOffset=1.4;

  CosThPhDist[iFrame][iRapBin][iPTBin][0]->GetYaxis()->SetTitleOffset(yOffset);
  //gPad->SetLeftMargin(0.125);
  CosThPhDist[iFrame][iRapBin][iPTBin][0]->SetStats(0);
  CosThPhDist[iFrame][iRapBin][iPTBin][0]->GetYaxis()->SetTitle(YTitle.str().c_str());
  CosThPhDist[iFrame][iRapBin][iPTBin][0]->GetXaxis()->SetTitle(XTitle.str().c_str());
  CosThPhDist[iFrame][iRapBin][iPTBin][0]->Draw("colz");
  if(iRapBin==0)
    latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                   onia::rapForPTRange[iRapBin+1],
                                   onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
  else
    latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                   onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                   onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
  std::stringstream distname;
  if(nState == 6)
    distname << dataPath.c_str() << "/Figures/cosThetaPhi_chic1_" << onia::frameLabel[iFrame] << "_rap" << iRapBin+1 << "_pT" << iPTBin+1 << ".pdf";
  else
    distname << dataPath.c_str() << "/Figures/cosThetaPhi_" << onia::frameLabel[iFrame] << "_rap" << iRapBin+1 << "_pT" << iPTBin+1 << ".pdf";
  c1->SaveAs(distname.str().c_str());

  // chic2
  if(nState == 6){
    CosThPhDist[iFrame][iRapBin][iPTBin][1]->GetYaxis()->SetTitleOffset(yOffset);
    CosThPhDist[iFrame][iRapBin][iPTBin][1]->SetStats(0);
    CosThPhDist[iFrame][iRapBin][iPTBin][1]->GetYaxis()->SetTitle(YTitle.str().c_str());
    CosThPhDist[iFrame][iRapBin][iPTBin][1]->GetXaxis()->SetTitle(XTitle.str().c_str());
    CosThPhDist[iFrame][iRapBin][iPTBin][1]->Draw("colz");
    if(iRapBin==0)
      latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                     onia::rapForPTRange[iRapBin+1],
                                     onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
      latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                     onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                     onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));

    std::stringstream distname1;
    distname1 << dataPath.c_str() << "/Figures/cosThetaPhi_chic2_" << onia::frameLabel[iFrame] << "_rap" << iRapBin+1 << "_pT" << iPTBin+1 << ".pdf";
    c1->SaveAs(distname1.str().c_str());
  }
}

//===========================
void LoadHistos(int iRapBin, int iPTBin, int nState, const std::string& dataPath){

  std::stringstream name, name1;
  if(nState == 6){
    name << dataPath.c_str() << "/tmpFiles/data_chic1_rap" << iRapBin+1 << "_pT" << iPTBin+1 << ".root";
    name1 << dataPath.c_str() << "/tmpFiles/data_chic2_rap" << iRapBin+1 << "_pT" << iPTBin+1 << ".root";
  } else {
    name << dataPath.c_str() << "/tmpFiles/data_Psi" << nState-3 << "_rap" << iRapBin+1 << "_pT" << iPTBin+1 << ".root";
    name1 << dataPath.c_str() << "/tmpFiles/data_Psi" << nState-3 << "_rap" << iRapBin+1 << "_pT" << iPTBin+1 << ".root";
  }

  TFile *fIn = new TFile(name.str().c_str());
  TTree* selectedData = (TTree*)fIn->Get("selectedData");
  TFile *fIn1 = new TFile(name1.str().c_str());
  TTree* selectedData1 = (TTree*)fIn1->Get("selectedData");

  TLorentzVector *lepP = new TLorentzVector();
  TLorentzVector *lepN = new TLorentzVector();
  selectedData->SetBranchAddress("lepP", &lepP);
  selectedData->SetBranchAddress("lepN", &lepN);
  TLorentzVector *lepP1 = new TLorentzVector();
  TLorentzVector *lepN1 = new TLorentzVector();
  selectedData1->SetBranchAddress("lepP", &lepP1);
  selectedData1->SetBranchAddress("lepN", &lepN1);

  // const double pbeam_ = 7000.; // exact number irrelevant as long as pbeam >> Mprot
  // const double Mprot_ = 0.9382720;
  // const double Mlepton_ = 0.10566;  // (muon)
  // const double Ebeam_ = sqrt( pbeam_*pbeam_ + Mprot_*Mprot_ );
  // TLorentzVector beam1_LAB_( 0., 0., pbeam_, Ebeam_ );
  // TLorentzVector beam2_LAB_( 0., 0., -pbeam_, Ebeam_ );

  const double gPI_ = TMath::Pi();
  
  std::cout << "chic1 has " << selectedData->GetEntries() << " events." << std::endl;
  
  for(int i=0;i<selectedData->GetEntries();i++){

    selectedData->GetEvent( i );

    double lepP_pT  = lepP->Pt();
    double lepN_pT  = lepN->Pt();

    double lepP_eta = lepP->PseudoRapidity();
    double lepN_eta = lepN->PseudoRapidity();

    // dilepton 4-vector:
    TLorentzVector dilepton = *lepP + *lepN;
    double pT   = dilepton.Pt();
    double rap  = dilepton.Rapidity();
    double mass = dilepton.M();

    TVector3 lab_to_dilep = -dilepton.BoostVector();

    TLorentzVector beam1_DILEP = onia::beam1_LAB;
    beam1_DILEP.Boost(lab_to_dilep);         // beam1 in the dilepton rest frame
    TLorentzVector beam2_DILEP = onia::beam2_LAB;
    beam2_DILEP.Boost(lab_to_dilep);         // beam2 in the dilepton rest frame

    TVector3 beam1_direction     = beam1_DILEP.Vect().Unit();
    TVector3 beam2_direction     = beam2_DILEP.Vect().Unit();
    TVector3 dilep_direction     = dilepton.Vect().Unit();
    TVector3 beam1_beam2_bisect  = ( beam1_direction - beam2_direction ).Unit();

    // all polarization frames have the same Y axis = the normal to the plane formed by
    // the directions of the colliding hadrons:
    TVector3 Yaxis = ( beam1_direction.Cross( beam2_direction ) ).Unit();

    // flip of y axis with rapidity:
    if ( rap < 0. ) Yaxis = - Yaxis;

    TVector3 perpendicular_to_beam = ( beam1_beam2_bisect.Cross( Yaxis ) ).Unit();

    // positive lepton in the dilepton rest frame:
    TLorentzVector lepton_DILEP = *lepP;
    lepton_DILEP.Boost(lab_to_dilep);

    // CS frame angles:
    TVector3 newZaxis = beam1_beam2_bisect;
    TVector3 newYaxis = Yaxis;
    TVector3 newXaxis = newYaxis.Cross( newZaxis );

    TRotation rotation;
    rotation.SetToIdentity();
    rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation.Invert();  // transforms coordinates from the "xyz" frame to the new frame
    TVector3 lepton_DILEP_rotated = lepton_DILEP.Vect();
    lepton_DILEP_rotated.Transform(rotation);

    double costh_CS = lepton_DILEP_rotated.CosTheta();
    double phi_CS   = lepton_DILEP_rotated.Phi() * 180. / gPI_;
    double phith_CS;
    if ( costh_CS < 0. ) phith_CS = phi_CS - 135.;
    if ( costh_CS > 0. ) phith_CS = phi_CS - 45.;
    if ( phith_CS < -180. ) phith_CS = 360. + phith_CS;

    // HELICITY frame angles:
    newZaxis = dilep_direction;
    newYaxis = Yaxis;
    newXaxis = newYaxis.Cross( newZaxis );

    rotation.SetToIdentity();
    rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation.Invert();
    lepton_DILEP_rotated = lepton_DILEP.Vect();
    lepton_DILEP_rotated.Transform(rotation);

    double costh_HX = lepton_DILEP_rotated.CosTheta();
    double phi_HX   = lepton_DILEP_rotated.Phi() * 180. / gPI_;
    double phith_HX;
    if ( costh_HX < 0. ) phith_HX = phi_HX - 135.;
    if ( costh_HX > 0. ) phith_HX = phi_HX - 45.;
    if ( phith_HX < -180. ) phith_HX = 360. + phith_HX;

    // PERPENDICULAR HELICITY frame angles:
    newZaxis = perpendicular_to_beam;
    newYaxis = Yaxis;
    newXaxis = newYaxis.Cross( newZaxis );

    rotation.SetToIdentity();
    rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation.Invert();
    lepton_DILEP_rotated = lepton_DILEP.Vect();
    lepton_DILEP_rotated.Transform(rotation);
    double costh_PX = lepton_DILEP_rotated.CosTheta();
    double phi_PX   = lepton_DILEP_rotated.Phi() * 180. / gPI_;
    double phith_PX;
    if ( costh_PX < 0. ) phith_PX = phi_PX - 135.;
    if ( costh_PX > 0. ) phith_PX = phi_PX - 45.;
    if ( phith_PX < -180. ) phith_PX = 360. + phith_PX;

    double cosalpha = sqrt( 1. - pow(costh_PX, 2.) ) * sin( lepton_DILEP_rotated.Phi() );

    // Filling Histograms of costh2, cos2ph and sin2thcosph for the extraction of the actual generated polarization
    // CS frame:
    double costh2 = pow(costh_CS,2.);
    double Phi = phi_CS/180. * gPI_;
    double cos2ph = cos(2.*Phi);
    double sin2thcosph= sin(2.*acos(costh_CS))*cos(Phi);
    h_costh2[0][iRapBin][iPTBin][0]->Fill( costh2 );
    h_cos2ph[0][iRapBin][iPTBin][0]->Fill( cos2ph );
    h_sin2thcosph[0][iRapBin][iPTBin][0]->Fill( sin2thcosph );

    // HX frame:
    costh2 = pow(costh_HX,2.);
    Phi = phi_HX/180. * gPI_;
    cos2ph = cos(2.*Phi);
    sin2thcosph= sin(2.*acos(costh_HX))*cos(Phi);
    h_costh2[1][iRapBin][iPTBin][0]->Fill( costh2 );
    h_cos2ph[1][iRapBin][iPTBin][0]->Fill( cos2ph );
    h_sin2thcosph[1][iRapBin][iPTBin][0]->Fill( sin2thcosph );

    // PHX frame:
    costh2 = pow(costh_PX,2.);
    Phi = phi_PX/180. * gPI_;
    cos2ph = cos(2.*Phi);
    sin2thcosph= sin(2.*acos(costh_PX))*cos(Phi);
    h_costh2[2][iRapBin][iPTBin][0]->Fill( costh2 );
    h_cos2ph[2][iRapBin][iPTBin][0]->Fill( cos2ph );
    h_sin2thcosph[2][iRapBin][iPTBin][0]->Fill( sin2thcosph );

    // fill cosTheta-phi-distributions in different frames
    CosThPhDist[0][iRapBin][iPTBin][0]->Fill(costh_CS,phi_CS);
    CosThPhDist[1][iRapBin][iPTBin][0]->Fill(costh_HX,phi_HX);
    CosThPhDist[2][iRapBin][iPTBin][0]->Fill(costh_PX,phi_PX);

    CosThDist[0][iRapBin][iPTBin][0]->Fill(costh_CS);
    PhiDist[0][iRapBin][iPTBin][0]->Fill(phi_CS);
    CosThDist[1][iRapBin][iPTBin][0]->Fill(costh_HX);
    PhiDist[1][iRapBin][iPTBin][0]->Fill(phi_HX);
    CosThDist[2][iRapBin][iPTBin][0]->Fill(costh_PX);
    PhiDist[2][iRapBin][iPTBin][0]->Fill(phi_PX);

  } // loop through events

    // loop for chic2 events
  if(nState == 6){
    std::cout << "chic2 has " << selectedData1->GetEntries() << " events." << std::endl;
    for(int i=0;i<selectedData1->GetEntries();i++){

      selectedData1->GetEvent( i );

      double lepP_pT  = lepP1->Pt();
      double lepN_pT  = lepN1->Pt();
      double lepP_eta = lepP1->PseudoRapidity();
      double lepN_eta = lepN1->PseudoRapidity();

      // dilepton 4-vector:
      TLorentzVector dilepton = *lepP1 + *lepN1;
      double pT   = dilepton.Pt();
      double rap  = dilepton.Rapidity();
      double mass = dilepton.M();

      TVector3 lab_to_dilep = -dilepton.BoostVector();
      TLorentzVector beam1_DILEP = onia::beam1_LAB;
      beam1_DILEP.Boost(lab_to_dilep);         // beam1 in the dilepton rest frame
      TLorentzVector beam2_DILEP = onia::beam2_LAB;
      beam2_DILEP.Boost(lab_to_dilep);         // beam2 in the dilepton rest frame
      TVector3 beam1_direction     = beam1_DILEP.Vect().Unit();
      TVector3 beam2_direction     = beam2_DILEP.Vect().Unit();
      TVector3 dilep_direction     = dilepton.Vect().Unit();
      TVector3 beam1_beam2_bisect  = ( beam1_direction - beam2_direction ).Unit();

      // all polarization frames have the same Y axis = the normal to the plane formed by
      // the directions of the colliding hadrons:
      TVector3 Yaxis = ( beam1_direction.Cross( beam2_direction ) ).Unit();
      // flip of y axis with rapidity:
      if ( rap < 0. ) Yaxis = - Yaxis;
      TVector3 perpendicular_to_beam = ( beam1_beam2_bisect.Cross( Yaxis ) ).Unit();
      // positive lepton in the dilepton rest frame:
      TLorentzVector lepton_DILEP = *lepP1;
      lepton_DILEP.Boost(lab_to_dilep);
      // CS frame angles:
      TVector3 newZaxis = beam1_beam2_bisect;
      TVector3 newYaxis = Yaxis;
      TVector3 newXaxis = newYaxis.Cross( newZaxis );
      TRotation rotation;
      rotation.SetToIdentity();
      rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
      rotation.Invert();  // transforms coordinates from the "xyz" frame to the new frame
      TVector3 lepton_DILEP_rotated = lepton_DILEP.Vect();
      lepton_DILEP_rotated.Transform(rotation);
      double costh_CS = lepton_DILEP_rotated.CosTheta();
      double phi_CS   = lepton_DILEP_rotated.Phi() * 180. / gPI_;
      double phith_CS;
      if ( costh_CS < 0. ) phith_CS = phi_CS - 135.;
      if ( costh_CS > 0. ) phith_CS = phi_CS - 45.;
      if ( phith_CS < -180. ) phith_CS = 360. + phith_CS;

      // HELICITY frame angles:
      newZaxis = dilep_direction;
      newYaxis = Yaxis;
      newXaxis = newYaxis.Cross( newZaxis );
      rotation.SetToIdentity();
      rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
      rotation.Invert();
      lepton_DILEP_rotated = lepton_DILEP.Vect();
      lepton_DILEP_rotated.Transform(rotation);
      
      double costh_HX = lepton_DILEP_rotated.CosTheta();
      double phi_HX   = lepton_DILEP_rotated.Phi() * 180. / gPI_;
      double phith_HX;
      if ( costh_HX < 0. ) phith_HX = phi_HX - 135.;
      if ( costh_HX > 0. ) phith_HX = phi_HX - 45.;
      if ( phith_HX < -180. ) phith_HX = 360. + phith_HX;
      
      // PERPENDICULAR HELICITY frame angles:
      newZaxis = perpendicular_to_beam;
      newYaxis = Yaxis;
      newXaxis = newYaxis.Cross( newZaxis );
      rotation.SetToIdentity();
      rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
      rotation.Invert();
      lepton_DILEP_rotated = lepton_DILEP.Vect();
      lepton_DILEP_rotated.Transform(rotation);
      double costh_PX = lepton_DILEP_rotated.CosTheta();
      double phi_PX   = lepton_DILEP_rotated.Phi() * 180. / gPI_;
      double phith_PX;
      if ( costh_PX < 0. ) phith_PX = phi_PX - 135.;
      if ( costh_PX > 0. ) phith_PX = phi_PX - 45.;
      if ( phith_PX < -180. ) phith_PX = 360. + phith_PX;
      
      double cosalpha = sqrt( 1. - pow(costh_PX, 2.) ) * sin( lepton_DILEP_rotated.Phi() );

      // Filling Histograms of costh2, cos2ph and sin2thcosph for the extraction of the actual generated polarization
      // CS frame:
      double costh2 = pow(costh_CS,2.);
      double Phi = phi_CS/180. * gPI_;
      double cos2ph = cos(2.*Phi);
      double sin2thcosph= sin(2.*acos(costh_CS))*cos(Phi);
      h_costh2[0][iRapBin][iPTBin][1]->Fill( costh2 );
      h_cos2ph[0][iRapBin][iPTBin][1]->Fill( cos2ph );
      h_sin2thcosph[0][iRapBin][iPTBin][1]->Fill( sin2thcosph );

      // HX frame:
      costh2 = pow(costh_HX,2.);
      Phi = phi_HX/180. * gPI_;
      cos2ph = cos(2.*Phi);
      sin2thcosph= sin(2.*acos(costh_HX))*cos(Phi);
      h_costh2[1][iRapBin][iPTBin][1]->Fill( costh2 );
      h_cos2ph[1][iRapBin][iPTBin][1]->Fill( cos2ph );
      h_sin2thcosph[1][iRapBin][iPTBin][1]->Fill( sin2thcosph );

      // PHX frame:
      costh2 = pow(costh_PX,2.);
      Phi = phi_PX/180. * gPI_;
      cos2ph = cos(2.*Phi);
      sin2thcosph= sin(2.*acos(costh_PX))*cos(Phi);
      h_costh2[2][iRapBin][iPTBin][1]->Fill( costh2 );
      h_cos2ph[2][iRapBin][iPTBin][1]->Fill( cos2ph );
      h_sin2thcosph[2][iRapBin][iPTBin][1]->Fill( sin2thcosph );

      // fill cosTheta-phi-distributions in different frames
      CosThPhDist[0][iRapBin][iPTBin][1]->Fill(costh_CS,phi_CS);
      CosThPhDist[1][iRapBin][iPTBin][1]->Fill(costh_HX,phi_HX);
      CosThPhDist[2][iRapBin][iPTBin][1]->Fill(costh_PX,phi_PX);

      CosThDist[0][iRapBin][iPTBin][1]->Fill(costh_CS);
      PhiDist[0][iRapBin][iPTBin][1]->Fill(phi_CS);
      CosThDist[1][iRapBin][iPTBin][1]->Fill(costh_HX);
      PhiDist[1][iRapBin][iPTBin][1]->Fill(phi_HX);
      CosThDist[2][iRapBin][iPTBin][1]->Fill(costh_PX);
      PhiDist[2][iRapBin][iPTBin][1]->Fill(phi_PX);

    } // loop through events
  } // nState == 6

  std::cout << "histograms initialzed and filled" << std::endl;
}
