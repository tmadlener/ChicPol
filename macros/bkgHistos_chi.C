//#include <iostream>
//#include <string>
//#include <vector>
//#include <sstream>

//using namespace RooFit;
//using namespace std;

//int order(int n);
//int findEvenNum(double number);
//TH2D *subtract2D(TH2D* hist1, TH2D* hist2);
//TH3D *subtract3D(TH3D* hist1, TH3D* hist2);
//TH2D* ReSetBin(TH2D* hist, int nBinX, int nBinY, const std::stringstream& name, const std::stringstream& title);

//---------------------------------------------------------------------------------------------------------------
void bkgHistos_chi(const std::string infilename, int rapBin, int ptBin, bool folding, bool MC, bool PolLSB, bool PolRSB, bool PolNP, int FracLSB, bool normApproach, bool subtractNP,
                   bool useRefittedChic){

  const std::string
    datafilename = "tmpFiles/selEvents_data.root",
    treename = "selectedData",
    wsname = "ws_masslifetime";

  // input
  TFile *datafile = TFile::Open(datafilename.c_str());
  if(!datafile){
    std::cout << "Inputfile missing" << std::endl;
    return;
  }
  TTree *intree = (TTree *)datafile->Get(treename.c_str());
  TLorentzVector *lepP = 0, *lepN = 0, *jpsi = 0, *chic = 0;

  TFile *fitfile = TFile::Open(infilename.c_str());
  if(!fitfile){
    std::cout << "fitfile is missing" << std::endl;
    return;
  }
  RooWorkspace *ws = (RooWorkspace*)fitfile->Get(wsname.c_str());
  if(!ws){
    std::cout << "workspace not found" << std::endl;
    return;
  }

  gStyle->SetPadRightMargin(0.2);
  gROOT->SetStyle("Plain");
  gStyle->SetTitleBorderSize(0);

  // create output
  std::stringstream outfilename1, outfilename2;
  outfilename1 << "tmpFiles/data_chic1_rap" << rapBin << "_pT" << ptBin << ".root";
  outfilename2 << "tmpFiles/data_chic2_rap" << rapBin << "_pT" << ptBin << ".root";
  TFile *output1 = TFile::Open(outfilename1.str().c_str(), "RECREATE");
  TTree *outtree1 = intree->CloneTree(0);

  TFile *output2 = TFile::Open(outfilename2.str().c_str(), "RECREATE");
  TTree *outtree2 = intree->CloneTree(0);

  // background histos
  TH2D *hBG_cosThetaPhiL[onia::kNbFrames];
  TH2D *hBG_cosThetaPhiR[onia::kNbFrames];
  TH2D *hBG2_cosThetaPhiR[onia::kNbFrames];
  TH2D *hBG1_cosThetaPhi[onia::kNbFrames];
  TH2D *hBG2_cosThetaPhi[onia::kNbFrames];
  TH2D *hBGinNP1_cosThetaPhi[onia::kNbFrames];
  TH2D *hBGinNP2_cosThetaPhi[onia::kNbFrames];
  TH2D *hBGinNP_cosThetaPhiL[onia::kNbFrames];
  TH2D *hBGinNP_cosThetaPhiR[onia::kNbFrames];
  TH2D *hBGinNP2_cosThetaPhiR[onia::kNbFrames];
  TH2D *hNPBG1_cosThetaPhi[onia::kNbFrames];
  TH2D *hNPBG2_cosThetaPhi[onia::kNbFrames];
  TH2D *hNPS1_cosThetaPhi[onia::kNbFrames];
  TH2D *hNPS2_cosThetaPhi[onia::kNbFrames];
  TH2D *hTBG1_cosThetaPhi[onia::kNbFrames];
  TH2D *hTBG2_cosThetaPhi[onia::kNbFrames];
  TH2D *hSR_cosThetaPhiL[onia::kNbFrames];
  TH2D *hSR_cosThetaPhiR[onia::kNbFrames];
  TH2D *hSR1_cosThetaPhi[onia::kNbFrames];
  TH2D *hSR2_cosThetaPhi[onia::kNbFrames];

  for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
    //book the 2D (cosTheta, phi) histos for the L and R mass sideband
    std::stringstream nameL, nameR, nameNP, nameBGinNPL, nameBGinNPR, nameSR, title;
    nameL << "hBG_cosThetaPhi_" << onia::frameLabel[iFrame] << "_L";
    nameR << "hBG_cosThetaPhi_" << onia::frameLabel[iFrame] << "_R";
    nameNP << "hNPBG_cosThetaPhi_" << onia::frameLabel[iFrame];
    nameBGinNPL << "hBGinNP_cosThetaPhi_" << onia::frameLabel[iFrame] << "_L";
    nameBGinNPR << "hBGinNP_cosThetaPhi_" << onia::frameLabel[iFrame] << "_R";
    nameSR << "hSR_cosThetaPhi_" << onia::frameLabel[iFrame];
    title << ";cos#theta_{"<< onia::frameLabel[iFrame] << "};#phi_{" << onia::frameLabel[iFrame] << "} [deg]";
    hBG_cosThetaPhiL[iFrame] = new TH2D(nameL.str().c_str(), title.str().c_str(),
                                        onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
    hBG_cosThetaPhiL[iFrame]->Sumw2();
    hBG_cosThetaPhiR[iFrame] = new TH2D(nameR.str().c_str(), title.str().c_str(),
                                        onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
    hBG_cosThetaPhiR[iFrame]->Sumw2();

    hNPBG1_cosThetaPhi[iFrame] = new TH2D(nameNP.str().c_str(), title.str().c_str(),
                                          onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
    hNPBG1_cosThetaPhi[iFrame]->Sumw2();
    hNPBG2_cosThetaPhi[iFrame] = new TH2D(nameNP.str().c_str(), title.str().c_str(),
                                          onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
    hNPBG2_cosThetaPhi[iFrame]->Sumw2();

    hBGinNP_cosThetaPhiL[iFrame] = new TH2D(nameBGinNPL.str().c_str(), title.str().c_str(),
                                            onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
    hBGinNP_cosThetaPhiL[iFrame]->Sumw2();

    hBGinNP_cosThetaPhiR[iFrame] = new TH2D(nameBGinNPR.str().c_str(), title.str().c_str(),
                                            onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
    hBGinNP_cosThetaPhiR[iFrame]->Sumw2();

    hSR1_cosThetaPhi[iFrame] = new TH2D(nameSR.str().c_str(), title.str().c_str(),
                                        onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
    hSR1_cosThetaPhi[iFrame]->Sumw2();
    hSR2_cosThetaPhi[iFrame] = new TH2D(nameSR.str().c_str(), title.str().c_str(),
                                        onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
    hSR2_cosThetaPhi[iFrame]->Sumw2();
  } // iFrame

    // mean pT and y histos (background-subtracted)
  int nBins = 100;
  TH1D* pT_L   = new TH1D( "pTLSB", "pTLSB", nBins, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin]);
  TH1D* pT_R   = new TH1D( "pTRSB", "pTRSB", nBins, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin]);
  TH1D* pT_highct_L   = new TH1D( "pT_highct_LSB", "pT_highct_LSB", nBins, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin]);
  TH1D* pT_highct_R   = new TH1D( "pT_highcta_RSB", "pT_highct_RSB", nBins, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin]);
  TH1D* pT_NP1   = new TH1D( "pTNP", "pTNP", nBins, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin]);
  TH1D* pT_NP2   = new TH1D( "pTNP", "pTNP", nBins, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin]);
  TH1D* pT_PSR1   = new TH1D( "pTPSR", "pTPSR", nBins, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin]);
  TH1D* pT_PSR2   = new TH1D( "pTPSR", "pTPSR", nBins, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin]);
  TH1D* rap_L   = new TH1D( "rapLSB", "rapLSB", nBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);
  TH1D* rap_R   = new TH1D( "rapRSB", "rapRSB", nBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);
  TH1D* rap_highct_L   = new TH1D( "rap_highct_LSB", "rap_highct_LSB", nBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);
  TH1D* rap_highct_R   = new TH1D( "rap_highct_RSB", "rap_highct_RSB", nBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);
  TH1D* rap_NP1   = new TH1D( "rapNP", "rapNP",nBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);
  TH1D* rap_NP2   = new TH1D( "rapNP", "rapNP",nBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);
  TH1D* rap_PSR1   = new TH1D( "rapPSR", "rapPSR", nBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);
  TH1D* rap_PSR2   = new TH1D( "rapPSR", "rapPSR", nBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);

  //------------------------------------------------------------------------------------------------
  // store pT and y borders
  TVectorD* pTBorder = new TVectorD(1, 2, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin], "END");
  TVectorD* yBorder = new TVectorD(1, 2, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin], "END");
  pTBorder->Print();
  yBorder->Print();
  output1->cd();
  pTBorder->Write();
  yBorder->Write();
  output2->cd();
  pTBorder->Write();
  yBorder->Write();

  // set branches
  intree->SetBranchAddress("lepP", &lepP);
  intree->SetBranchAddress("lepN", &lepN);
  intree->SetBranchAddress("jpsi", &jpsi);
  // tmadlener, 15.04.2016:
  // differentiat between use of refitted chic and non-refitted chic.
  // If refitted chic is used, everything is obtained from it. If non-refitted chic is used mQ is used as mass
  // other necessary changes are done after intree->GetEntry(). (See lines 578, 996, 1563 or search for "tmadlener")
  double mQ = 0;
  if (useRefittedChic) {
    intree->SetBranchAddress("chic_rf", &chic);
  } else {
    intree->SetBranchAddress("mQ", &mQ);
    intree->SetBranchAddress("chic", &chic);
  }


  double jpsict = 0;
  intree->SetBranchAddress("Jpsict", &jpsict);

  //---------------------------------------------------------------------------------------------
  // INPUT FROM ROOWORKSPACE
  // pre-calculated fractions
  double fracLSB1 = ((RooRealVar*)ws->var("var_fLSBChic1"))->getValV();
  double fracLSB2 = ((RooRealVar*)ws->var("var_fLSBChic2"))->getValV();
  if(FracLSB!=-1){
    fracLSB1 = double(FracLSB)/100.;
    fracLSB2 = double(FracLSB)/100.;
  }
  std::cout << "fLSB = " << fracLSB1 << " for chic1 and " << fracLSB2 << " for chic2" << std::endl;
  double fBGsig1 = ((RooRealVar*)ws->var("var_fracBackgroundInPRSR1"))->getValV(); //comb. background in PRSR1
  double fBGsig2 = ((RooRealVar*)ws->var("var_fracBackgroundInPRSR2"))->getValV(); //comb. background in PRSR2
  double fBGinNP1 = ((RooRealVar*)ws->var("var_fracBackgroundInNPSR1"))->getValV(); //comb. background in NPSR1
  double fBGinNP2 = ((RooRealVar*)ws->var("var_fracBackgroundInNPSR2"))->getValV(); //comb. background in NPSR2
  double fNPB1 = ((RooRealVar*)ws->var("var_fracNPChic1InPRSR1"))->getValV(); // non prompt background in PRSR1
  double fNPB2 = ((RooRealVar*)ws->var("var_fracNPChic2InPRSR2"))->getValV(); // non prompt background in PRSR2
  double fTBGsig1 = fBGsig1;  // total background fraction
  double fTBGsig2 = fBGsig2;  // total background fraction
  double fP1 = 1 - fBGsig1 - fNPB1;
  double fP2 = 1 - fBGsig2 - fNPB2;
  double fSRinPLSB = ((RooRealVar*)ws->var("var_fracPRChic1InPRLSB"))->getValV(); // prompt chic1 contamination in LSB
  double fSRinPRSB = ((RooRealVar*)ws->var("var_fracPRChic2InPRRSB"))->getValV(); // prompt chic2 contamination in RSB
  double nBkg = ((RooRealVar*)ws->var("var_nBackground"))->getValV(); // total number of background events
  double nChic = ((RooRealVar*)ws->var("var_nChic"))->getValV(); // total number of chic events
  double nChic1 = ((RooRealVar*)ws->var("var_nChic1"))->getValV(); // total number of chic1 events
  double nChic2 = ((RooRealVar*)ws->var("var_nChic2"))->getValV(); // total number of chic2 events
  double nPRChic1InPRSR1 = ((RooRealVar*)ws->var("var_nPRChic1InPRSR1"))->getValV(); // number of prompt chic1 events in PRSR1
  double nPRChic2InPRSR2 = ((RooRealVar*)ws->var("var_nPRChic2InPRSR2"))->getValV(); // number of prompt chic2 events in PRSR2
  double fTotInSR1 = ((RooRealVar*)ws->var("var_fTotInSR1"))->getValV(); // fraction of total events in SR1
  double fBGinSR1 = ((RooRealVar*)ws->var("var_fracBackgroundInSR1"))->getValV();; // background fraction in SR1
  double nSR1 = (nBkg+nChic)*fTotInSR1; // total number of events in SR1
  double fTotInSR2 = ((RooRealVar*)ws->var("var_fTotInSR2"))->getValV(); // fraction of total events in SR2
  double fBGinSR2 = ((RooRealVar*)ws->var("var_fracBackgroundInSR2"))->getValV();; // background fraction in SR2
  double nSR2 = (nBkg+nChic)*fTotInSR2; // total number of events in SR2
  double fPRSR1InSR1 = ((RooRealVar*)ws->var("var_fPRSR1InSR1"))->getValV();; // fraction of PRSR1 in SR1
  double nPRSR1 = nSR1*fPRSR1InSR1; // total number of events in PRSR1
  double fNPSR1InSR1 = ((RooRealVar*)ws->var("var_fNPSR1InSR1"))->getValV();; // fraction of NPSR1 in SR1
  double nNPSR1 = nSR1*fNPSR1InSR1; // total number of events in NPSR1
  double fNPChic1InPRSR1 = ((RooRealVar*)ws->var("var_fracNPChic1InPRSR1"))->getValV();; // non prompt chic1 fraction in PRSR1
  double fPRChic1InNPSR1 = ((RooRealVar*)ws->var("var_fracPRChic1InNPSR1"))->getValV();; // prompt chic1 fraction in NPSR1
  double fNPChic1InNPSR1 = ((RooRealVar*)ws->var("var_fracNPChic1InNPSR1"))->getValV();; // non prompt chic1 fraction in NPSR1
  double fPRSR2InSR2 = ((RooRealVar*)ws->var("var_fPRSR2InSR2"))->getValV();; // fraction of PRSR2 in SR2
  double nPRSR2 = nSR2*fPRSR2InSR2; // total number of events in PRSR2
  double fNPSR2InSR2 = ((RooRealVar*)ws->var("var_fNPSR2InSR2"))->getValV();; // fraction of NPSR2 in SR2
  double nNPSR2 = nSR2*fNPSR2InSR2; // total number of events in NPSR2
  double fNPChic2InPRSR2 = ((RooRealVar*)ws->var("var_fracNPChic2InPRSR2"))->getValV();; // non prompt chic2 fraction in PRSR2
  double fPRChic2InNPSR2 = ((RooRealVar*)ws->var("var_fracPRChic2InNPSR2"))->getValV();; // prompt chic2 fraction in NPSR2
  double fNPChic2InNPSR2 = ((RooRealVar*)ws->var("var_fracNPChic2InNPSR2"))->getValV();; // non prompt chic2 fraction in NPSR2

  // errrors on fractions
  double fBGerr1 = ((RooRealVar*)ws->var("var_fracBackgroundInPRSR1"))->getError();
  double fBGerr2 = ((RooRealVar*)ws->var("var_fracBackgroundInPRSR2"))->getError();
  double fNPerr1 = ((RooRealVar*)ws->var("var_fracBackgroundInNPSR1"))->getError();
  double fNPerr2 = ((RooRealVar*)ws->var("var_fracBackgroundInNPSR2"))->getError();
  double fBGinNPerr1 = ((RooRealVar*)ws->var("var_fracNPChic1InPRSR1"))->getError();
  double fBGinNPerr2 = ((RooRealVar*)ws->var("var_fracNPChic1InPRSR2"))->getError();
  double fTBGerr1 = fBGerr1;
  double fTBGerr2 = fBGerr2;
  double fPerr1 = TMath::Sqrt(TMath::Power(fNPerr1, 2) + TMath::Power(fBGerr1, 2));
  double fPerr2 = TMath::Sqrt(TMath::Power(fNPerr2, 2) + TMath::Power(fBGerr2, 2));

  // in case of measuring non prompt polarization
  if(PolNP){
    fTBGsig1 = fBGinNP1;
    fTBGsig2 = fBGinNP2;
    fTBGerr1 = fBGinNPerr1;
    fTBGerr2 = fBGinNPerr2;
  }
  // for MC set fractions to 0.001 except prompt fraction
  if(MC || PolLSB || PolRSB){
    fBGinNP1 = 0.001;
    fBGinNP2 = 0.001;
    fNPB1    = 0.001;
    fNPB2    = 0.001;
    fBGsig1  = 0.001;
    fBGsig2  = 0.001;
    fP1 = 0.001;
    fP2 = 0.002;
    fBGerr1 = 0;
    fBGerr2 = 0;
    fTBGerr1 = 0;
    fTBGerr2 = 0;
    fBGinNPerr1 = 0;
    fBGinNPerr2 = 0;
    fNPerr1 = 0;
    fNPerr2 = 0;
    fPerr1 = 0;
    fPerr2 = 0;
    fTBGsig1 = 0.001;
    fTBGsig2 = 0.001;
    if(PolLSB){
      fTBGsig1 = fSRinPLSB;
      fTBGsig2 = fSRinPLSB;
    }else if(PolRSB){
      fTBGsig1 = fSRinPRSB;
      fTBGsig2 = fSRinPRSB;
    }
  }
  // if NP should also be subtracted
  if(subtractNP){
    fTBGsig1 =  fNPB1 + fBGsig1;
    fTBGsig2 =  fNPB1 + fBGsig2;
    fTBGerr1 = TMath::Sqrt(TMath::Power(fNPerr1, 2) + TMath::Power(fBGerr1, 2));
    fTBGerr2 = TMath::Sqrt(TMath::Power(fNPerr2, 2) + TMath::Power(fBGerr2, 2));
  }

  // variables
  RooRealVar *m = ws->var("chicMass");
  RooRealVar *mjpsi = ws->var("JpsiMass");
  RooRealVar *ct = ws->var("Jpsict");
  RooRealVar *poly1 = ws->var("BK_p1");
  RooRealVar *poly2 = ws->var("BK_p2");
  RooRealVar *CBmass1 = ws->var("CBmass1");
  RooRealVar *CBmass2 = ws->var("CBmass2");
  RooRealVar *CBsigma1 = ws->var("CBsigma1");
  RooRealVar *CBsigma2 = ws->var("CBsigma2");
  RooRealVar *CBalpha1 = ws->var("CBalpha1");
  RooRealVar *CBalpha2 = ws->var("CBalpha2");
  RooRealVar *CBn1 = ws->var("CBn");
  RooRealVar *CBn2 = ws->var("CBn2");
  RooRealVar *lambdaJpsi = ws->var("bkgLambda_jpsi");
  RooRealVar *CBmassJpsi = ws->var("CBmass_jpsi");
  RooRealVar *CBsigmaJpsi = ws->var("CBsigma_jpsi");

  // pdf
  RooAbsPdf *bkgMass = (RooAbsPdf*)ws->pdf("M_background");
  RooAbsPdf *signalMass1 = (RooAbsPdf*)ws->pdf("M_chic1");
  RooAbsPdf *signalMass2 = (RooAbsPdf*)ws->pdf("M_chic2");
  RooAbsPdf *bkgMassJpsi = (RooAbsPdf*)ws->pdf("bkgMassShape_jpsi");
  RooAbsPdf *jpsiMass = (RooAbsPdf*)ws->pdf("gaussMassShape_jpsi");

  // load snapshot with all results
  std::stringstream masssnapshotname;
  masssnapshotname << "m_snapshot_rap" << rapBin << "_pt" << ptBin;
  ws->loadSnapshot(masssnapshotname.str().c_str());

  // functions to draw random mass values
  TF1* funcBG = (TF1*)bkgMass->asTF(*m, RooArgList(*poly1, *poly2), *m);
  TF1* funcSig1 = (TF1*)signalMass1->asTF(*m, RooArgList(*CBmass1, *CBsigma1, *CBalpha1, *CBn1), *m);
  TF1* funcSig2 = (TF1*)signalMass2->asTF(*m, RooArgList(*CBmass2, *CBsigma2, *CBalpha2, *CBn2), *m);
  TF1* funcBGJpsi = (TF1*)bkgMassJpsi->asTF(*mjpsi, RooArgList(*lambdaJpsi), *mjpsi);
  TF1* funcSigJpsi = (TF1*)bkgMassJpsi->asTF(*mjpsi, RooArgList(*CBmassJpsi, *CBsigmaJpsi), *mjpsi);

  //-------------------------------------------------------------------------------------------------
  // regions
  double massMinL = onia::massChiSBMin;
  double massMaxR = onia::massChiSBMax;
  double massMaxL = ((RooRealVar*)ws->var("var_lsbMaxMass"))->getValV();
  double massMinR = ((RooRealVar*)ws->var("var_rsbMinMass"))->getValV();
  double massMinSR1 = ((RooRealVar*)ws->var("var_sig1MinMass"))->getValV();
  double massMinSR2 = ((RooRealVar*)ws->var("var_sig2MinMass"))->getValV();
  double massMaxSR1 = ((RooRealVar*)ws->var("var_sig1MaxMass"))->getValV();
  double massMaxSR2 = ((RooRealVar*)ws->var("var_sig2MaxMass"))->getValV();
  double PRmin = ((RooRealVar*)ws->var("var_PRMin"))->getValV();
  double PRmax = ((RooRealVar*)ws->var("var_PRMax"))->getValV();
  double NPmin = ((RooRealVar*)ws->var("var_NPMin"))->getValV();
  double NPmax = ((RooRealVar*)ws->var("var_NPMax"))->getValV();

  std::cout << "-------------------------------------------------------------\n" <<
    "left  sideband: mass window " << massMinL  << " < M < " << massMaxL  << " GeV\n" <<
    "right sideband: mass window " << massMinR  << " < M < " << massMaxR  << " GeV\n" <<
    "signal	 region chic1: mass window " << massMinSR1  << " < M < " << massMaxSR1	<< " GeV\n" <<
    "signal	 region chic2: mass window " << massMinSR2  << " < M < " << massMaxSR2	<< " GeV\n" <<
    "-------------------------------------------------------------\n" << std::endl;

  std::stringstream namedata;
  namedata << "jpsi_data_rap" << rapBin << "_pt" << ptBin;
  RooAbsData* data = ws->data(namedata.str().c_str());

  std::stringstream cutSR1, cutSR2, cutLSB, cutRSB;
  cutSR1 << "chicMass > " << massMinSR1 << " && chicMass < " << massMaxSR1;
  cutSR2 << "chicMass > " << massMinSR2 << " && chicMass < " << massMaxSR2;
  cutLSB << "chicMass > " << massMinL << " && chicMass < " << massMaxL;
  cutRSB << "chicMass > " << massMinR << " && chicMass < " << massMaxR;

  std::stringstream binNameSR1, binNameSR2, binNameLSB, binNameRSB;
  binNameSR1  << "data_rap" << rapBin << "_pt" << ptBin << "_SR1";
  binNameSR2  << "data_rap" << rapBin << "_pt" << ptBin << "_SR2";
  binNameLSB << "data_rap" << rapBin << "_pt" << ptBin << "_LSB";
  binNameRSB << "data_rap" << rapBin << "_pt" << ptBin << "_RSB";

  RooAbsData* dataSR1 = data->reduce(Cut(cutSR1.str().c_str()));
  RooAbsData* dataSR2 = data->reduce(Cut(cutSR2.str().c_str()));
  RooAbsData* dataLSB = data->reduce(Cut(cutLSB.str().c_str()));
  RooAbsData* dataRSB = data->reduce(Cut(cutRSB.str().c_str()));

  std::stringstream cutPR, cutNP;
  cutPR << "Jpsict > " << PRmin << " && Jpsict < " << PRmax;
  cutNP << "Jpsict > " << NPmin << " && Jpsict < " << NPmax;

  RooAbsData* dataPRSR1 = dataSR1->reduce(Cut(cutPR.str().c_str()));
  RooAbsData* dataPRSR2 = dataSR2->reduce(Cut(cutPR.str().c_str()));
  RooAbsData* dataPRLSB = dataLSB->reduce(Cut(cutPR.str().c_str()));
  RooAbsData* dataPRRSB = dataRSB->reduce(Cut(cutPR.str().c_str()));
  RooAbsData* dataNPSR1 = dataSR1->reduce(Cut(cutNP.str().c_str()));
  RooAbsData* dataNPSR2 = dataSR2->reduce(Cut(cutNP.str().c_str()));
  RooAbsData* dataNPLSB = dataLSB->reduce(Cut(cutNP.str().c_str()));
  RooAbsData* dataNPRSB = dataRSB->reduce(Cut(cutNP.str().c_str()));

  //------------------------------------------------------------------------------------------------
  // chic1
  output1->cd();
  // fill histogram with combinatorial background fraction
  std::string bkgname1 = ";;fraction of comb. BG in PRSR1";
  TH1D* hFracBG1 = new TH1D("comb_background_fraction", bkgname1.c_str(), 1, 0., 1.);
  hFracBG1->SetBinContent(1, fBGsig1);
  hFracBG1->SetBinError(1, fBGerr1);
  hFracBG1->Write();
  // fill histogram with non prompt background fraction
  std::string NPbkgname1 =  ";;fraction of non prompt BG in PRSR1";
  TH1D* hFracNPBG1 = new TH1D("nonprompt_background_fraction", NPbkgname1.c_str(), 1, 0., 1.);
  hFracNPBG1->SetBinContent(1, fNPB1);
  hFracNPBG1->SetBinError(1, fNPerr1);
  hFracNPBG1->Write();
  // fill histogram with total background fraction
  std::string tbkgname1 = ";;fraction of total BG in PRSR1";
  TH1D* hFracTBG1 = new TH1D("background_fraction", tbkgname1.c_str(), 1, 0., 1.);
  /*if(scaleFracBg){
    std::stringstream filename;
    filename << polDataPath.c_str() << "/results_Psi" << nState-3 << "S_rap" << rapBin << "_pT" << ptBin << ".root";
    TFile *resultFile = new TFile(filename.str().c_str89,"R");
    TH1D* SubtractedBG_test=(TH1D*)resultFile->Get("SubtractedBG_test");
    double ratio = SubtractedBG_test -> GetMean();
    std::cout << "fBG_sub / fBG : " << ratio << std::endl;
    fTBGsig = fTBGsig / ratio ;
    resultFile->Close();
    }
    output->cd();*/
  hFracTBG1->SetBinContent(1, fTBGsig1);
  hFracTBG1->SetBinError(1, fTBGerr1);
  hFracTBG1->Write();
  // fill histogram with prompt fraction
  std::string Pname1 =  ";;fraction of prompt events in PRSR1 (1-fNP-fBkg)";
  TH1D* hFracP1 = new TH1D("prompt_fraction", Pname1.c_str(), 1, 0., 1.);
  hFracP1->SetBinContent(1, fP1);
  hFracP1->SetBinError(1, fPerr1);
  hFracP1->Write();
  // fill histogram with fraction of LSB
  TH1D* hFracLSB1 = new TH1D("fraction_LSB", ";;f_{LSB}", 1, 0., 1.);
  hFracLSB1->SetBinContent(1, fracLSB1);
  hFracLSB1->Write();
  // fill histogram with events in signal region
  TH1D* hEvtSR1 = new TH1D("events_SR", ";;chic1 events in SR1", 3, 0., 3.);
  // fill histogram only for data
  if(!MC){
    hEvtSR1->SetBinContent(1, nPRChic1InPRSR1+fPRChic1InNPSR1*nNPSR1);
    hEvtSR1->SetBinContent(2, fNPChic1InPRSR1*nPRSR1+fNPChic1InNPSR1*nNPSR1);
    hEvtSR1->SetBinContent(3, fBGinSR1*nSR1);
    hEvtSR1->Write();
  }
  // fill histogram with events in prompt signal region
  TH1D* hEvtPSR1 = new TH1D("events_promptSR", ";;chic1 events in PRSR1", 3, 0., 3.);
  // fill histogram only for data
  if(!MC){
    hEvtPSR1->SetBinContent(1, nPRChic1InPRSR1);
    hEvtPSR1->SetBinContent(2, fNPChic1InPRSR1*nPRSR1);
    hEvtPSR1->SetBinContent(3, fBGsig1*nPRSR1);
    hEvtPSR1->Write();
  }
  //fill histogram with events in non-prompt signal region
  TH1D* hEvtNPSR1 = new TH1D("events_nonpromptSR", ";;chic1 events in NPSR1", 3, 0., 3.);
  // fill histogram only for data
  if(!MC){
    hEvtNPSR1->SetBinContent(1, fPRChic1InNPSR1*nNPSR1);
    hEvtNPSR1->SetBinContent(2, fNPChic1InNPSR1*nNPSR1);
    hEvtNPSR1->SetBinContent(3, fBGinNP1*nNPSR1);
    hEvtNPSR1->Write();
  }

  // chic2
  output2->cd();
  std::string bkgname2 = ";;fraction of comb. BG in PRSR2";
  TH1D* hFracBG2 = new TH1D("comb_background_fraction", bkgname2.c_str(), 1, 0., 1.);
  hFracBG2->SetBinContent(1, fBGsig2);
  hFracBG2->SetBinError(1, fBGerr2);
  hFracBG2->Write();
  std::string NPbkgname2 =  ";;fraction of non prompt BG in PRSR2";
  TH1D* hFracNPBG2 = new TH1D("nonprompt_background_fraction", NPbkgname2.c_str(), 1, 0., 1.);
  hFracNPBG2->SetBinContent(1, fNPB2);
  hFracNPBG2->SetBinError(1, fNPerr2);
  hFracNPBG2->Write();
  std::string tbkgname2 = ";;fraction of total BG in PRSR2";
  TH1D* hFracTBG2 = new TH1D("background_fraction", tbkgname2.c_str(), 1, 0., 1.);
  hFracTBG2->SetBinContent(1, fTBGsig2);
  hFracTBG2->SetBinError(1, fTBGerr2);
  hFracTBG2->Write();
  std::string Pname2 =  ";;fraction of prompt events in PRSR2 (1-fNP-fBkg)";
  TH1D* hFracP2 = new TH1D("prompt_fraction", Pname2.c_str(), 1, 0., 1.);
  hFracP2->SetBinContent(1, fP2);
  hFracP2->SetBinError(1, fPerr2);
  hFracP2->Write();
  TH1D* hFracLSB2 = new TH1D("fraction_LSB", ";;f_{LSB}", 1, 0., 1.);
  hFracLSB2->SetBinContent(1, fracLSB2);
  hFracLSB2->Write();
  TH1D* hEvtSR2 = new TH1D("events_SR", ";;chic2 events in SR2", 3, 0., 3.);
  if(!MC){
    hEvtSR2->SetBinContent(1, nPRChic2InPRSR2 + fPRChic2InNPSR2*nNPSR2);
    hEvtSR2->SetBinContent(2, fNPChic2InPRSR2*nPRSR2 + fNPChic2InNPSR2*nNPSR2);
    hEvtSR2->SetBinContent(3, fBGinSR2*nSR2);
    hEvtSR2->Write();
  }
  TH1D* hEvtPSR2 = new TH1D("events_promptSR", ";;chic2 events in PRSR2", 3, 0., 3.);
  if(!MC){
    hEvtPSR2->SetBinContent(1, nPRChic2InPRSR2);
    hEvtPSR2->SetBinContent(2, fNPChic2InPRSR2*nPRSR2);
    hEvtPSR2->SetBinContent(3, fBGsig2*nPRSR2);
    hEvtPSR2->Write();
  }
  TH1D* hEvtNPSR2 = new TH1D("events_nonpromptSR", ";;chic2 events in NPSR2", 3, 0., 3.);
  if(!MC){
    hEvtNPSR2->SetBinContent(1, fPRChic2InNPSR2*nNPSR2);
    hEvtNPSR2->SetBinContent(2, fNPChic2InNPSR2*nNPSR2);
    hEvtNPSR2->SetBinContent(3, fBGinNP2*nNPSR2);
    hEvtNPSR2->Write();
  }

  //------------------------------------------------------------------------------------------------
  // build the 3D (pT, |y|, M) histos for the L and R mass sideband
  TH3D* hBG1_pTRapMass_L = new TH3D("hBG_pTRapMass_L", ";p_{T} [GeV/c]; |y|; M [GeV]",
                                    7, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin],
                                    onia::kNbRapForPTBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin],
                                    7, massMinSR1, massMaxSR1); // signal mass window!
  hBG1_pTRapMass_L->Sumw2();

  TH3D* hBG2_pTRapMass_L = new TH3D("hBG_pTRapMass_L", ";p_{T} [GeV/c]; |y|; M [GeV]",
                                    7, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin],
                                    onia::kNbRapForPTBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin],
                                    7, massMinSR2, massMaxSR2); // signal mass window!
  hBG2_pTRapMass_L->Sumw2();

  TH3D* hBG1_pTRapMass_R = new TH3D("hBG_pTRapMass_R", ";p_{T} [GeV/c]; |y|; M [GeV]",
                                    7, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin],
                                    onia::kNbRapForPTBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin],
                                    7, massMinSR1, massMaxSR1); // signal mass window!
  hBG1_pTRapMass_R->Sumw2();

  TH3D* hBG2_pTRapMass_R = new TH3D("hBG_pTRapMass_R", ";p_{T} [GeV/c]; |y|; M [GeV]",
                                    7, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin],
                                    onia::kNbRapForPTBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin],
                                    7, massMinSR2, massMaxSR2); // signal mass window!
  hBG2_pTRapMass_R->Sumw2();

  TH3D* hBG1_pTRapMass_highct_L = new TH3D("hBG_pTRapMass_highct_L", ";p_{T} [GeV/c]; |y|; M [GeV]",
                                           7, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin],
                                           onia::kNbRapForPTBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin],
                                           7, massMinSR1, massMaxSR1);
  hBG1_pTRapMass_highct_L->Sumw2();

  TH3D* hBG2_pTRapMass_highct_L = new TH3D("hBG_pTRapMass_highct_L", ";p_{T} [GeV/c]; |y|; M [GeV]",
                                           7, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin],
                                           onia::kNbRapForPTBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin],
                                           7, massMinSR2, massMaxSR2);
  hBG2_pTRapMass_highct_L->Sumw2();

  TH3D* hBG1_pTRapMass_highct_R = new TH3D("hBG_pTRapMass_highct_R", ";p_{T} [GeV/c]; |y|; M [GeV]",
                                           7, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin],
                                           onia::kNbRapForPTBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin],
                                           7, massMinSR1, massMaxSR1);
  hBG1_pTRapMass_highct_R->Sumw2();

  TH3D* hBG2_pTRapMass_highct_R = new TH3D("hBG_pTRapMass_highct_R", ";p_{T} [GeV/c]; |y|; M [GeV]",
                                           7, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin],
                                           onia::kNbRapForPTBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin],
                                           7, massMinSR2, massMaxSR2);
  hBG2_pTRapMass_highct_R->Sumw2();

  TH3D* hNP1_pTRapMass = new TH3D("hNP_pTRapMass_NP", ";p_{T} [GeV/c]; |y|; M [GeV]",
                                  7, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin],
                                  onia::kNbRapForPTBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin],
                                  7, massMinSR1, massMaxSR1);
  hNP1_pTRapMass->Sumw2();

  TH3D* hNP2_pTRapMass = new TH3D("hNP_pTRapMass_NP", ";p_{T} [GeV/c]; |y|; M [GeV]",
                                  7, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin],
                                  onia::kNbRapForPTBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin],
                                  7, massMinSR2, massMaxSR2);
  hNP2_pTRapMass->Sumw2();

  TH3D* hSR1_pTRapMass = new TH3D("hSR_pTRapMass", ";p_{T} [GeV/c]; |y|; M [GeV]",
                                  7, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin],
                                  onia::kNbRapForPTBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin],
                                  7, massMinSR1, massMaxSR1);
  hSR1_pTRapMass->Sumw2();

  TH3D* hSR2_pTRapMass = new TH3D("hSR_pTRapMass", ";p_{T} [GeV/c]; |y|; M [GeV]",
                                  7, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin],
                                  onia::kNbRapForPTBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin],
                                  7, massMinSR2, massMaxSR2);
  hSR2_pTRapMass->Sumw2();

  //------------------------------------------------------------------------------------------------
  // loop through tree, fill background histos and save data in pt and y bins
  int index = -1;
  int MCindex = -1;
  int n = intree->GetEntries();
  int count1=0;
  int count2=0;
  int MCevents1 = 0;
  int MCevents2 = 0;

  for(int i = 0; i < n; i++){

    long iEntry = intree->LoadTree(i);
    intree->GetEntry(iEntry);
    if(i % 100000 == 0) {std::cout << "entry " << i << " out of " << n << std::endl;}

    // tmadlener, 15.04.2016:
    // have to define tha mass here depending on which chic should be used
    double usedChicMass = 0;
    if (useRefittedChic) {
      usedChicMass = chic->M(); // use mass from refitted chic (which just happens to be in the chic variable)
    } else {
      usedChicMass = mQ; // use mQ as mass for the non-refitted chic
    }

    // ------------------------- TLorentzVecotrs -------------------------
    if(chic->Pt() >= onia::pTRange[rapBin-1][ptBin-1] &&
       chic->Pt() < onia::pTRange[rapBin-1][ptBin] &&
       TMath::Abs(chic->Rapidity()) >= onia::rapForPTRange[rapBin-1] &&
       TMath::Abs(chic->Rapidity()) < onia::rapForPTRange[rapBin]){

      //store TLorentzVectors of the two muons in the given pT and rap cell
      // left sideband
      if(PolLSB){
        if(usedChicMass > massMinL && usedChicMass < massMaxL && jpsict > PRmin && jpsict < PRmax){
          outtree1->Fill();
          outtree2->Fill();
          pT_PSR1->Fill(chic->Pt());
          rap_PSR1->Fill(TMath::Abs(chic->Rapidity()));
          pT_PSR2->Fill(chic->Pt());
          rap_PSR2->Fill(TMath::Abs(chic->Rapidity()));
        }
      } // PolLSB
      // right sideband
      else if(PolRSB){
        if(usedChicMass > massMinR && usedChicMass < massMaxR && jpsict > PRmin && jpsict < PRmax){
          outtree1->Fill();
          outtree2->Fill();
          pT_PSR1->Fill(chic->Pt());
          rap_PSR1->Fill(TMath::Abs(chic->Rapidity()));
          pT_PSR2->Fill(chic->Pt());
          rap_PSR2->Fill(TMath::Abs(chic->Rapidity()));
        }
      } // PolRSB
      // for non prompt data
      else if(PolNP){
        if(usedChicMass > massMinSR1 && usedChicMass < massMaxSR1 && jpsict > NPmin && jpsict < NPmax){
          outtree1->Fill();
          pT_PSR1->Fill(chic->Pt());
          rap_PSR1->Fill(TMath::Abs(chic->Rapidity()));
        }
        else if(usedChicMass > massMinSR2 && usedChicMass < massMaxSR2 && jpsict > NPmin && jpsict < NPmax){
          outtree2->Fill();
          pT_PSR2->Fill(chic->Pt());
          rap_PSR2->Fill(TMath::Abs(chic->Rapidity()));
        }
      } // PolNP
      // for prompt data and MC
      //store only events from signal region
      else if(MC){
        if(usedChicMass > massMinSR1 && usedChicMass < massMaxSR1){
          outtree1->Fill();
          pT_PSR1->Fill(chic->Pt());
          rap_PSR1->Fill(TMath::Abs(chic->Rapidity()));
          MCevents1++;
        }
        else if(usedChicMass > massMinSR2 && usedChicMass < massMaxSR2){
          outtree2->Fill();
          pT_PSR2->Fill(chic->Pt());
          rap_PSR2->Fill(TMath::Abs(chic->Rapidity()));
          MCevents2++;
        }
      }
      else{
        if(usedChicMass > massMinSR1 && usedChicMass < massMaxSR1 && jpsict > PRmin && jpsict < PRmax){
          outtree1->Fill();
          pT_PSR1->Fill(chic->Pt());
          rap_PSR1->Fill(TMath::Abs(chic->Rapidity()));
          count1++;
        }
        else if(usedChicMass > massMinSR2 && usedChicMass < massMaxSR2 && jpsict > PRmin && jpsict < PRmax){
          outtree2->Fill();
          pT_PSR2->Fill(chic->Pt());
          rap_PSR2->Fill(TMath::Abs(chic->Rapidity()));
          count2++;
        }
      } // else

      //---------------------------- mass histograms for background model ------------------------------
      // for MC: fill mass histograms with random mass from signal region
      // fill rapidity and pT histograms with all events
      if(MC){
        pT_L->Fill(chic->Pt());
        pT_R->Fill(chic->Pt());
        pT_highct_L->Fill(chic->Pt());
        pT_highct_R->Fill(chic->Pt());
        rap_L->Fill(TMath::Abs(chic->Rapidity()));
        rap_R->Fill(TMath::Abs(chic->Rapidity()));
        rap_highct_L->Fill(TMath::Abs(chic->Rapidity()));
        rap_highct_R->Fill(TMath::Abs(chic->Rapidity()));

        if(usedChicMass > massMinSR1 && usedChicMass < massMaxSR1){
          MCindex = 0; // events in SR1 get index 0
          hNP1_pTRapMass->Fill(chic->Pt(), TMath::Abs(chic->Rapidity()), gRandom->Uniform(massMinSR1, massMaxSR1));
          hSR1_pTRapMass->Fill(chic->Pt(), TMath::Abs(chic->Rapidity()), gRandom->Uniform(massMinSR1, massMaxSR1));
          pT_NP1->Fill(chic->Pt());
          rap_NP1->Fill(TMath::Abs(chic->Rapidity()));
        }
        else if(usedChicMass > massMinSR2 && usedChicMass < massMaxSR2){
          MCindex = 1; // events in SR2 get index 1
          hNP2_pTRapMass->Fill(chic->Pt(), TMath::Abs(chic->Rapidity()), gRandom->Uniform(massMinSR2, massMaxSR2));
          hSR2_pTRapMass->Fill(chic->Pt(), TMath::Abs(chic->Rapidity()), gRandom->Uniform(massMinSR2, massMaxSR2));
          pT_NP2->Fill(chic->Pt());
          rap_NP2->Fill(TMath::Abs(chic->Rapidity()));
        }
        // split events up to fill in left and right background histograms
        // not all events are filled to be able to put together background histogram
        if(i%3==0){
          hBG1_pTRapMass_L->Fill(chic->Pt(), TMath::Abs(chic->Rapidity()), gRandom->Uniform(massMinSR1, massMaxSR1));
          hBG2_pTRapMass_L->Fill(chic->Pt(), TMath::Abs(chic->Rapidity()), gRandom->Uniform(massMinSR2, massMaxSR2));
          hBG1_pTRapMass_highct_L->Fill(chic->Pt(), TMath::Abs(chic->Rapidity()), gRandom->Uniform(massMinSR1, massMaxSR1));
          hBG2_pTRapMass_highct_L->Fill(chic->Pt(), TMath::Abs(chic->Rapidity()), gRandom->Uniform(massMinSR2, massMaxSR2));
        }
        else if(i%5==0){
          hBG1_pTRapMass_R->Fill(chic->Pt(), TMath::Abs(chic->Rapidity()), gRandom->Uniform(massMinSR1, massMaxSR1));
          hBG2_pTRapMass_R->Fill(chic->Pt(), TMath::Abs(chic->Rapidity()), gRandom->Uniform(massMinSR2, massMaxSR2));
          hBG1_pTRapMass_highct_R->Fill(chic->Pt(), TMath::Abs(chic->Rapidity()), gRandom->Uniform(massMinSR1, massMaxSR1));
          hBG2_pTRapMass_highct_R->Fill(chic->Pt(), TMath::Abs(chic->Rapidity()), gRandom->Uniform(massMinSR2, massMaxSR2));
        }
      } // if (MC)

      else{
        // store cosTheta and phi distributions of the background
        // events gets index 0 if it is in the prompt and 1 if its in the non prompt left sideband
        // events with index 2 and 3 are from the prompt and non prompt chic1 signal region
        // events with index 4 and 5 are from the prompt and non prompt chic2 signal region
        // events with index 6 and 7 are from the prompt and non prompt right sideband
        if(usedChicMass > massMinL && usedChicMass < massMaxL && jpsict > PRmin && jpsict < PRmax){
          index = 0;
          hBG1_pTRapMass_L->Fill(chic->Pt(), TMath::Abs(chic->Rapidity()), funcBG->GetRandom(massMinSR1, massMaxSR1));
          hBG2_pTRapMass_L->Fill(chic->Pt(), TMath::Abs(chic->Rapidity()), funcBG->GetRandom(massMinSR2, massMaxSR2));
          pT_L->Fill(chic->Pt());
          rap_L->Fill(TMath::Abs(chic->Rapidity()));
        }
        else if(usedChicMass > massMinL && usedChicMass < massMaxL && jpsict > NPmin && jpsict < NPmax){
          index = 1;
          hBG1_pTRapMass_highct_L->Fill(chic->Pt(), TMath::Abs(chic->Rapidity()), funcBG->GetRandom(massMinSR1, massMaxSR1));
          hBG2_pTRapMass_highct_L->Fill(chic->Pt(), TMath::Abs(chic->Rapidity()), funcBG->GetRandom(massMinSR2, massMaxSR2));
          pT_highct_L->Fill(chic->Pt());
          rap_highct_L->Fill(TMath::Abs(chic->Rapidity()));
        }
        else if(usedChicMass > massMinSR1 && usedChicMass < massMaxSR1 && jpsict > PRmin && jpsict < PRmax){
          index = 2;
          hSR1_pTRapMass->Fill(chic->Pt(), TMath::Abs(chic->Rapidity()), funcSig1->GetRandom(massMinSR1, massMaxSR1));
        }
        else if(usedChicMass > massMinSR1 && usedChicMass < massMaxSR1 && jpsict > NPmin && jpsict < NPmax){
          index = 3;
          hNP1_pTRapMass->Fill(chic->Pt(), TMath::Abs(chic->Rapidity()), funcSig1->GetRandom(massMinSR1, massMaxSR1));
          pT_NP1->Fill(chic->Pt());
          rap_NP1->Fill(TMath::Abs(chic->Rapidity()));
        }
        else if(usedChicMass > massMinSR2 && usedChicMass < massMaxSR2 && jpsict > PRmin && jpsict < PRmax){
          index = 4;
          hSR2_pTRapMass->Fill(chic->Pt(), TMath::Abs(chic->Rapidity()), funcSig2->GetRandom(massMinSR2, massMaxSR2));
        }
        else if(usedChicMass > massMinSR2 && usedChicMass < massMaxSR2 && jpsict > NPmin && jpsict < NPmax){
          index = 5;
          hNP2_pTRapMass->Fill(chic->Pt(), TMath::Abs(chic->Rapidity()), funcSig2->GetRandom(massMinSR2, massMaxSR2));
          pT_NP2->Fill(chic->Pt());
          rap_NP2->Fill(TMath::Abs(chic->Rapidity()));
        }
        else if(usedChicMass > massMinR && usedChicMass < massMaxR && jpsict > PRmin && jpsict < PRmax){
          index = 6;
          hBG1_pTRapMass_R->Fill(chic->Pt(), TMath::Abs(chic->Rapidity()), funcBG->GetRandom(massMinSR1, massMaxSR1));
          hBG2_pTRapMass_R->Fill(chic->Pt(), TMath::Abs(chic->Rapidity()), funcBG->GetRandom(massMinSR2, massMaxSR2));
          pT_R->Fill(chic->Pt());
          rap_R->Fill(TMath::Abs(chic->Rapidity()));
        }
        else if(usedChicMass > massMinR && usedChicMass < massMaxR && jpsict > NPmin && jpsict < NPmax){
          index = 7;
          hBG1_pTRapMass_highct_R->Fill(chic->Pt(), TMath::Abs(chic->Rapidity()), funcBG->GetRandom(massMinSR1, massMaxSR1));
          hBG2_pTRapMass_highct_R->Fill(chic->Pt(), TMath::Abs(chic->Rapidity()), funcBG->GetRandom(massMinSR2, massMaxSR2));
          pT_highct_R->Fill(chic->Pt());
          rap_highct_R->Fill(TMath::Abs(chic->Rapidity()));
        }
        else continue;
      }// else (filling data histograms)

      ///////////////////////
      calcPol(*lepP, *lepN);
      ///////////////////////

      for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){

        // folding in phi
        double phiFolded = thisPhi[iFrame];
        double thetaAdjusted = thisCosTh[iFrame];
        if(thisPhi[iFrame] >= -90. && thisPhi[iFrame] < 0.) phiFolded *= -1;
        else if(thisPhi[iFrame] >= 90 && thisPhi[iFrame] < 180){
          phiFolded = 180. - thisPhi[iFrame];
          thetaAdjusted *= -1;
        }
        else if(thisPhi[iFrame] >= -180. && thisPhi[iFrame] < -90.){
          phiFolded = 180. + thisPhi[iFrame];
          thetaAdjusted *= -1;
        }

        // if bool folding is true, folding is applied to all background histograms
        if(folding){
          thisPhi[iFrame] = phiFolded;
          thisCosTh[iFrame] = thetaAdjusted;
        }

        // filling histograms
        if(MC){
          if(MCindex == 0){
            hNPBG1_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
            hSR1_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
          } else if(MCindex == 1){
            hNPBG2_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
            hSR2_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
          }
          if(i%3==0){
            hBG_cosThetaPhiL[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
            hBGinNP_cosThetaPhiL[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
          }
          else if(i%5==0){
            hBG_cosThetaPhiR[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
            hBGinNP_cosThetaPhiR[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
          }
        }
        // for data: fill histograms according to the different regions
        else{
          if(index == 0) hBG_cosThetaPhiL[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
          else if(index == 1) hBGinNP_cosThetaPhiL[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
          else if(index == 2) hSR1_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
          else if(index == 3) hNPBG1_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
          else if(index == 4) hSR2_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
          else if(index == 5) hNPBG2_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
          else if(index == 6) hBG_cosThetaPhiR[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
          else if(index == 7) hBGinNP_cosThetaPhiR[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
        }
      } // iFrame
    } // if(onia...)
  } // i

  std::cout << "total events in PRSR: " << count1 << " chic1, " << count2 << " chic2" << std::endl;

  //fill histograms with number of events for MC
  if(MC){
    // number of events in signal region
    hEvtSR1->SetBinContent(1, fP1*MCevents1);
    hEvtSR2->SetBinContent(1, fP2*MCevents2);
    // number of events in prompt signal region (same as in signal region)
    hEvtPSR1->SetBinContent(1, fP1*MCevents1);
    hEvtPSR2->SetBinContent(1, fP2*MCevents2);
    output1->cd();
    hEvtSR1->Write();
    hEvtPSR1->Write();
    output2->cd();
    hEvtSR2->Write();
    hEvtPSR2->Write();
  }

  //-----------------------------------------------------------------------------------------------------------------
  //---------------- binning algorithm
  int nBinsPhi = 16, nBinsCosth = 160;
  int totalBins = 0, filledBins = 0;
  for(int binCosth = 0; binCosth < hBG_cosThetaPhiR[2]->GetNbinsX(); binCosth++){
    for(int binPhi = 0; binPhi < hBG_cosThetaPhiR[2]->GetNbinsY(); binPhi++){
      totalBins++;
      //use NP histo (same physical coverage, but better filled, no holes -> better estimate of coverage)
      int binContent = hNPBG1_cosThetaPhi[2]->GetBinContent(binCosth+1,binPhi+1);
      if(binContent>0) filledBins++;
    }
  }

  double coverage = 2*(double)filledBins/(double)totalBins;
  nBinsCosth = 16*2/coverage;
  // find 2^n closest to nBinsCosth, but above the actual number
  nBinsCosth = findEvenNum((double)nBinsCosth);
  // set maximum binning to 64
  if(nBinsCosth > 64) nBinsCosth = 64;

  std::cout << "------------------------------------------------" << "\n"
            << "Starting binning algorithm" << "\n"
            << "filled bins: " << filledBins << "\n"
            << "total bins: " << totalBins << "\n"
            << "bin coverage: " << coverage << "\n"
            << "starting point for binning in phi: " << nBinsPhi << "\n"
            << "starting point for binning in cosTheta: " << nBinsCosth << "\n"
            << "------------------------------------------------" <<std::endl;

  // calculate the integral of the lowstatBG histo (calculate all integrals of the 2 PR sideband regions, and use the one with the smallest integral)
  int IntBG = hBG_cosThetaPhiL[2]->Integral();
  if(IntBG > hBG_cosThetaPhiR[2]->Integral()){
    IntBG = hBG_cosThetaPhiR[2]->Integral();
    std::cout << "right low ct integral is smaller" << std::endl;
  }

  int nBinsPhiBG = nBinsPhi,
    nBinsCosthBG = nBinsCosth;

  // calculate average events per-bin cell for background histo
  double Naverage = (double)IntBG/((double)nBinsPhi*nBinsCosth*coverage/2.);
  std::cout << "average cell coverage: " << Naverage << std::endl;

  // if average events per bin is bigger than 10, no rebinning is needed
  if(Naverage > 10){
    std::cout << "Rebinning is not necessary in this case." << "\n"
              << "Ending binning algorithm." << "\n"
              << "------------------------------------------------" << std::endl;
  }
  // otherwise rebin
  else{
    std::cout << "------------------------------------------------" << "\n"
              << "old cosTheta binning: " << nBinsCosth << "\n"
              << "old phi binning: " << nBinsPhi << std::endl;

    //set nBinsPhi to the lowest 2^n, such that nBinsPhi > nBinsCosth*coverage/2
    nBinsPhi = findEvenNum(nBinsCosth*coverage/2.);

    std::cout << "closest 2^n number to cosTheta bins: " << nBinsCosth << "\n"
              << "lowest 2^n number so that phi bins > cosTheta bins * coverage/2: " << nBinsPhi << "\n"
              << "------------------------------------------------" << std::endl;

    // set minimum binning
    int nBinsPhiMin = 8,
      nBinsCosthMin = 8;
    if(folding) nBinsPhiMin = 16;

    //BG
    nBinsCosthBG = nBinsCosth;
    nBinsPhiBG = nBinsPhi;
    double NaverageBG = 0.;

    for(int i = 0; i < 500; i++){

      std::cout << "looping for correct binning in background histogram" << std::endl;
      // If the mimimum number of bins for both phi and costh are reached, stop the loop
      if(nBinsPhiBG/2 < nBinsPhiMin && nBinsCosthBG/2 < nBinsCosthMin) break;

      //Change the binning, first in phi, then in costh:
      if(nBinsPhiBG/2 >= nBinsPhiMin) nBinsPhiBG = nBinsPhiBG/2;  //This ensures a mimimum number of bins in phi, e.g. 4
      NaverageBG = (double)IntBG/((double)nBinsPhiBG*nBinsCosthBG*coverage/2.);
      std::cout << "average bin content per cell after " << i << " phi rebinning: " << NaverageBG << std::endl;
      if(NaverageBG > 10) break;

      if(nBinsCosthBG/2 >= nBinsCosthMin) nBinsCosthBG = nBinsCosthBG/2; //This ensures a mimimum number of bins in costh, e.g. 4
      NaverageBG = (double)IntBG/((double)nBinsPhiBG*nBinsCosthBG*coverage/2.);
      std::cout << "average bin content per cell after " << i << " cosTheta rebinning: " << NaverageBG << std::endl;
      if(NaverageBG > 10) break;
    }
    std::cout << "average bin content per cell exceeds 10: " << NaverageBG << "\n"
              << "phi bins = " << nBinsPhiBG << ", cosTheta bins = " << nBinsCosthBG << "\n"
              << "------------------------------------------------" << std::endl;

  } // else Naverage < 10

  std::cout << "final binning for background histogram: " << "\n"
            << "phi bins: " << nBinsPhiBG << "\n"
            << "cosTheta bins: " << nBinsCosthBG << "\n"
            << "final binning for non prompt histogram" << "\n"
            << "phi bins: " << nBinsPhiBG << "\n"
            << "cosTheta bins: " << nBinsCosthBG << "\n"
            << "------------------------------------------------" << std::endl;

  //loop again with new binning
  for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
    //book the 2D (cosTheta, phi) histos for the L and R mass sideband
    std::stringstream nameL, nameR, nameNP, nameBGinNPL, nameBGinNPR, nameSR, title;
    nameL << "hBG_cosThetaPhi_" << onia::frameLabel[iFrame] << "_L";
    nameR << "hBG_cosThetaPhi_" << onia::frameLabel[iFrame] << "_R";
    nameNP << "hNPBG_cosThetaPhi_" << onia::frameLabel[iFrame];
    nameBGinNPL << "hBGinNP_cosThetaPhi_" << onia::frameLabel[iFrame] << "_L";
    nameBGinNPR << "hBGinNP_cosThetaPhi_" << onia::frameLabel[iFrame] << "_R";
    nameSR << "hSR_cosThetaPhi_" << onia::frameLabel[iFrame];
    title << ";cos#theta_{"<< onia::frameLabel[iFrame] << "};#phi_{" << onia::frameLabel[iFrame] << "} [deg]";

    delete hBG_cosThetaPhiL[iFrame];
    hBG_cosThetaPhiL[iFrame] = new TH2D(nameL.str().c_str(), title.str().c_str(),
                                        nBinsCosthBG, onia::cosTMin, onia::cosTMax, nBinsPhiBG, onia::phiPolMin, onia::phiPolMax);
    hBG_cosThetaPhiL[iFrame]->Sumw2();
    delete hBG_cosThetaPhiR[iFrame];
    hBG_cosThetaPhiR[iFrame] = new TH2D(nameR.str().c_str(), title.str().c_str(),
                                        nBinsCosthBG, onia::cosTMin, onia::cosTMax, nBinsPhiBG, onia::phiPolMin, onia::phiPolMax);
    hBG_cosThetaPhiR[iFrame]->Sumw2();

    delete hBGinNP_cosThetaPhiL[iFrame];
    hBGinNP_cosThetaPhiL[iFrame] = new TH2D(nameBGinNPL.str().c_str(), title.str().c_str(),
                                            nBinsCosthBG, onia::cosTMin, onia::cosTMax, nBinsPhiBG, onia::phiPolMin, onia::phiPolMax);
    hBGinNP_cosThetaPhiL[iFrame]->Sumw2();
    delete hBGinNP_cosThetaPhiR[iFrame];
    hBGinNP_cosThetaPhiR[iFrame] = new TH2D(nameBGinNPR.str().c_str(), title.str().c_str(),
                                            nBinsCosthBG, onia::cosTMin, onia::cosTMax, nBinsPhiBG, onia::phiPolMin, onia::phiPolMax);
    hBGinNP_cosThetaPhiR[iFrame]->Sumw2();

    delete hNPBG1_cosThetaPhi[iFrame];
    hNPBG1_cosThetaPhi[iFrame] = new TH2D(nameNP.str().c_str(), title.str().c_str(),
                                          nBinsCosthBG, onia::cosTMin, onia::cosTMax, nBinsPhiBG, onia::phiPolMin, onia::phiPolMax);
    hNPBG1_cosThetaPhi[iFrame]->Sumw2();

    delete hNPBG2_cosThetaPhi[iFrame];
    hNPBG2_cosThetaPhi[iFrame] = new TH2D(nameNP.str().c_str(), title.str().c_str(),
                                          nBinsCosthBG, onia::cosTMin, onia::cosTMax, nBinsPhiBG, onia::phiPolMin, onia::phiPolMax);
    hNPBG2_cosThetaPhi[iFrame]->Sumw2();

    delete hSR1_cosThetaPhi[iFrame];
    hSR1_cosThetaPhi[iFrame] = new TH2D(nameSR.str().c_str(), title.str().c_str(),
                                        nBinsCosthBG, onia::cosTMin, onia::cosTMax, nBinsPhiBG, onia::phiPolMin, onia::phiPolMax);
    hSR1_cosThetaPhi[iFrame]->Sumw2();

    delete hSR2_cosThetaPhi[iFrame];
    hSR2_cosThetaPhi[iFrame] = new TH2D(nameSR.str().c_str(), title.str().c_str(),
                                        nBinsCosthBG, onia::cosTMin, onia::cosTMax, nBinsPhiBG, onia::phiPolMin, onia::phiPolMax);
    hSR2_cosThetaPhi[iFrame]->Sumw2();

  } // iFrame

  for(int i = 0; i < n; i++){

    long iEntry = intree->LoadTree(i);
    intree->GetEntry(iEntry);

    // tmadlener, 15.04.2016:
    // have to define tha mass here depending on which chic should be used
    double usedChicMass = 0;
    if (useRefittedChic) {
      usedChicMass = chic->M(); // use mass from refitted chic (which just happens to be in the chic variable)
    } else {
      usedChicMass = mQ; // use mQ as mass for the non-refitted chic
    }

    if(i % 100000 == 0) std::cout << "entry " << i << " out of " << n << std::endl;

    if(chic->Pt() >= onia::pTRange[rapBin-1][ptBin-1] &&
       chic->Pt() < onia::pTRange[rapBin-1][ptBin] &&
       TMath::Abs(chic->Rapidity()) >= onia::rapForPTRange[rapBin-1] &&
       TMath::Abs(chic->Rapidity()) < onia::rapForPTRange[rapBin]){

      if(!MC){
        if(usedChicMass > massMinL && usedChicMass < massMaxL && jpsict > PRmin && jpsict < PRmax) index = 0;
        else if(usedChicMass > massMinL && usedChicMass < massMaxL && jpsict > NPmin && jpsict < NPmax) index = 1;
        else if(usedChicMass > massMinSR1 && usedChicMass < massMaxSR1 && jpsict > PRmin && jpsict < PRmax) index = 2;
        else if(usedChicMass > massMinSR1 && usedChicMass < massMaxSR1 && jpsict > NPmin && jpsict < NPmax) index = 3;
        else if(usedChicMass > massMinSR2 && usedChicMass < massMaxSR2 && jpsict > PRmin && jpsict < PRmax) index = 4;
        else if(usedChicMass > massMinSR2 && usedChicMass < massMaxSR2 && jpsict > NPmin && jpsict < NPmax) index = 5;
        else if(usedChicMass > massMinR && usedChicMass < massMaxR && jpsict > PRmin && jpsict < PRmax) index = 6;
        else if(usedChicMass > massMinR && usedChicMass < massMaxR && jpsict > NPmin && jpsict < NPmax) index = 7;
        else continue;
      }
      else{
        if(usedChicMass > massMinSR1 && usedChicMass < massMaxSR1) MCindex = 0;
        else if(usedChicMass > massMinSR2 && usedChicMass < massMaxSR2) MCindex = 1;
      }

      ////////////////////////
      calcPol(*lepP, *lepN);
      ////////////////////////

      for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){

        // folding in phi
        double phiFolded = thisPhi[iFrame];
        double thetaAdjusted = thisCosTh[iFrame];
        if(thisPhi[iFrame] > -90. && thisPhi[iFrame] < 0.) phiFolded *= -1;
        else if(thisPhi[iFrame] > 90 && thisPhi[iFrame] < 180){
          phiFolded = 180. - thisPhi[iFrame];
          thetaAdjusted *= -1;
        }
        else if(thisPhi[iFrame] > -180. && thisPhi[iFrame] < -90.){
          phiFolded = 180. + thisPhi[iFrame];
          thetaAdjusted *= -1;
        }

        // if folding is true, apply folding in phi
        if(folding){
          thisPhi[iFrame] = phiFolded;
          thisCosTh[iFrame] = thetaAdjusted;
        }

        // filling histograms
        if(MC){
          if(MCindex == 0){
            hNPBG1_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
            hSR1_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
          } else if(MCindex == 1){
            hNPBG2_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
            hSR2_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
          }
          if(i%3==0){
            hBG_cosThetaPhiL[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
            hBGinNP_cosThetaPhiL[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
          }else if(i%5==0){
            hBG_cosThetaPhiR[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
            hBGinNP_cosThetaPhiR[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
          }
        }
        // for data: fill histograms according to the different regions
        else{
          if(index == 0) hBG_cosThetaPhiL[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
          else if(index == 1) hBGinNP_cosThetaPhiL[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
          else if(index == 2) hSR1_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
          else if(index == 3) hNPBG1_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
          else if(index == 4) hSR2_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
          else if(index == 5) hNPBG2_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
          else if(index == 6) hBG_cosThetaPhiR[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
          else if(index == 7) hBGinNP_cosThetaPhiR[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
        }
      } // iFrame
    } // if(onia...)
  } // i
  //------loop finished

  //---------------- end- binning algorithm

  //----------------------------------------------------------------------------------------------------
  // 3D (pT, |y|, M) histos
  // write background histos to file
  output1->cd();
  hBG1_pTRapMass_L->Write();
  hBG1_pTRapMass_R->Write();
  hSR1_pTRapMass->Write();
  hBG1_pTRapMass_highct_L->Write();
  hBG1_pTRapMass_highct_R->Write();
  hNP1_pTRapMass->Write();
  output2->cd();
  hBG2_pTRapMass_L->Write();
  hBG2_pTRapMass_R->Write();
  hSR2_pTRapMass->Write();
  hBG2_pTRapMass_highct_L->Write();
  hBG2_pTRapMass_highct_R->Write();
  hNP2_pTRapMass->Write();

  // add left and right sideband of low ctau region to combinatorial background histogram
  // chic1
  hBG1_pTRapMass_L->Scale(fracLSB1/(1.*hBG1_pTRapMass_L->Integral()));
  hBG1_pTRapMass_R->Scale((1.-fracLSB1)/(1.*hBG1_pTRapMass_R->Integral()));
  std::string namepTrapMasslowct = "comb_background_pTrapMass";
  TH3D* hBG1_pTRapMass_lowct = (TH3D*) hBG1_pTRapMass_L->Clone(namepTrapMasslowct.c_str());
  hBG1_pTRapMass_lowct->Add(hBG1_pTRapMass_R);
  output1->cd();
  hBG1_pTRapMass_lowct->Write();
  // chic2
  hBG2_pTRapMass_L->Scale(fracLSB2/(1.*hBG2_pTRapMass_L->Integral()));
  hBG2_pTRapMass_R->Scale((1.-fracLSB2)/(1.*hBG2_pTRapMass_R->Integral()));
  TH3D* hBG2_pTRapMass_lowct = (TH3D*) hBG2_pTRapMass_L->Clone(namepTrapMasslowct.c_str());
  hBG2_pTRapMass_lowct->Add(hBG2_pTRapMass_R);
  output2->cd();
  hBG2_pTRapMass_lowct->Write();

  // add left and right sideband of high ctau region to combinatorial background histogram in high ctau region
  // chic1
  hBG1_pTRapMass_highct_L->Scale(fracLSB1/(1.*hBG1_pTRapMass_highct_L->Integral()));
  hBG1_pTRapMass_highct_R->Scale((1.-fracLSB1)/(1.*hBG1_pTRapMass_highct_R->Integral()));
  std::string namepTrapMasshighct = "comb_background_highct_pTrapMass";
  TH3D* hBG1_pTRapMass_highct = (TH3D*) hBG1_pTRapMass_highct_L->Clone(namepTrapMasshighct.c_str());
  hBG1_pTRapMass_highct->Add(hBG1_pTRapMass_highct_R);
  output1->cd();
  hBG1_pTRapMass_highct->Write();
  // chic2
  hBG2_pTRapMass_highct_L->Scale(fracLSB2/(1.*hBG2_pTRapMass_highct_L->Integral()));
  hBG2_pTRapMass_highct_R->Scale((1.-fracLSB2)/(1.*hBG2_pTRapMass_highct_R->Integral()));
  TH3D* hBG2_pTRapMass_highct = (TH3D*) hBG2_pTRapMass_highct_L->Clone(namepTrapMasshighct.c_str());
  hBG2_pTRapMass_highct->Add(hBG2_pTRapMass_highct_R);
  output2->cd();
  hBG2_pTRapMass_highct->Write();

  // create non prompt background histogram in signal region: (hNP)_norm - fBGinNP * (hBG_highct)_norm
  // chic1
  hNP1_pTRapMass->Scale(1./(1.*hNP1_pTRapMass->Integral()));
  hBG1_pTRapMass_highct->Scale(fBGinNP1/(1.*hBG1_pTRapMass_highct->Integral()));
  std::string namepTrapMassNPS = "NPS_highct_pTrapMass";
  TH3D* hNPS1_pTRapMass = (TH3D*) hNP1_pTRapMass->Clone(namepTrapMassNPS.c_str());
  hNPS1_pTRapMass = subtract3D(hNPS1_pTRapMass, hBG1_pTRapMass_highct);
  output1->cd();
  hNPS1_pTRapMass->Write();
  // chic2
  hNP2_pTRapMass->Scale(1./(1.*hNP2_pTRapMass->Integral()));
  hBG2_pTRapMass_highct->Scale(fBGinNP2/(1.*hBG2_pTRapMass_highct->Integral()));
  TH3D* hNPS2_pTRapMass = (TH3D*) hNP2_pTRapMass->Clone(namepTrapMassNPS.c_str());
  hNPS2_pTRapMass = subtract3D(hNPS2_pTRapMass, hBG2_pTRapMass_highct);
  output2->cd();
  hNPS2_pTRapMass->Write();

  // create total background
  std::string namepTrapMass = "background_pTrapMass";
  TH3D* hBG1_pTRapMass = new TH3D();
  TH3D* hBG2_pTRapMass = new TH3D();
  if(PolLSB){
    // for polarization of LSB: total background = signal contamination of prompt chic1
    hSR1_pTRapMass->Scale(1./(1.*hSR1_pTRapMass->Integral()));
    hBG1_pTRapMass = (TH3D*) hSR1_pTRapMass->Clone(namepTrapMass.c_str());
    hSR2_pTRapMass->Scale(1./(1.*hSR1_pTRapMass->Integral()));
    hBG2_pTRapMass = (TH3D*) hSR1_pTRapMass->Clone(namepTrapMass.c_str());
  } else if(PolRSB){
    // for polarization of RSB: total background = signal contamination of prompt chic2
    hSR1_pTRapMass->Scale(1./(1.*hSR2_pTRapMass->Integral()));
    hBG1_pTRapMass = (TH3D*) hSR2_pTRapMass->Clone(namepTrapMass.c_str());
    hSR2_pTRapMass->Scale(1./(1.*hSR2_pTRapMass->Integral()));
    hBG2_pTRapMass = (TH3D*) hSR2_pTRapMass->Clone(namepTrapMass.c_str());
  } else if(PolNP){
    // for non prompt polarization: total background = high ct background
    hBG1_pTRapMass = (TH3D*) hBG1_pTRapMass_highct->Clone(namepTrapMass.c_str());
    hBG2_pTRapMass = (TH3D*) hBG2_pTRapMass_highct->Clone(namepTrapMass.c_str());
  } else if(subtractNP){
    // add low ct background and non prompt background
    hNPS1_pTRapMass->Scale(fNPB1/(1.*hNPS1_pTRapMass->Integral()));
    hBG1_pTRapMass = (TH3D*) hNPS1_pTRapMass->Clone(namepTrapMass.c_str());
    hBG1_pTRapMass_lowct->Scale(fBGsig1/(1.*hBG1_pTRapMass_lowct->Integral()));
    hBG1_pTRapMass->Add(hBG1_pTRapMass_lowct);
    hNPS2_pTRapMass->Scale(fNPB2/(1.*hNPS2_pTRapMass->Integral()));
    hBG2_pTRapMass = (TH3D*) hNPS2_pTRapMass->Clone(namepTrapMass.c_str());
    hBG2_pTRapMass_lowct->Scale(fBGsig2/(1.*hBG2_pTRapMass_lowct->Integral()));
    hBG2_pTRapMass->Add(hBG2_pTRapMass_lowct);
  } else {
    hBG1_pTRapMass = (TH3D*) hBG1_pTRapMass_lowct->Clone(namepTrapMass.c_str());
    hBG2_pTRapMass = (TH3D*) hBG2_pTRapMass_lowct->Clone(namepTrapMass.c_str());
  }
  output1->cd();
  hBG1_pTRapMass->Write();
  output2->cd();
  hBG2_pTRapMass->Write();

  // mean pT and y histos
  double meanPT1 = 0, meanPT2 = 0;
  double meanY1 = 0, meanY2 = 0;
  // prompt chic1 contamination in LSB for polarization of LSB
  if(PolLSB){
    pT_PSR1->Scale(1./(1.*pT_PSR1->Integral()));
    TH1D* pT_SRL = (TH1D*) pT_PSR1->Clone();
    pT_SRL->Scale(fSRinPLSB/(1.*pT_SRL->Integral()));
    pT_L->Scale(1./(1.*pT_L->Integral()));
    pT_L->Add(pT_SRL, -1.);
    meanPT1 = pT_L->GetMean();
    meanPT2 = pT_L->GetMean();

    rap_PSR1->Scale(1./(1.*rap_PSR1->Integral()));
    TH1D* rap_SRL = (TH1D*) rap_PSR1->Clone();
    rap_SRL->Scale(fSRinPLSB/(1.*rap_SRL->Integral()));
    rap_L->Scale(1./(1.*rap_L->Integral()));
    rap_L->Add(rap_SRL, -1.);
    meanY1 = rap_L->GetMean();
    meanY2 = rap_L->GetMean();
  }
  // prompt chic2 contamination in RSB for polarization of RSB
  else if(PolRSB){
    pT_PSR2->Scale(1./(1.*pT_PSR2->Integral()));
    TH1D* pT_SRR = (TH1D*) pT_PSR2->Clone();
    pT_SRR->Scale(fSRinPRSB/(1.*pT_SRR->Integral()));
    pT_R->Scale(1./(1.*pT_R->Integral()));
    pT_R->Add(pT_SRR, -1.);
    meanPT1 = pT_R->GetMean();
    meanPT2 = pT_R->GetMean();

    TH1D* rap_SRR = (TH1D*) rap_PSR2->Clone();
    rap_SRR->Scale(fSRinPRSB/(1.*rap_SRR->Integral()));
    rap_R->Scale(1./(1.*rap_R->Integral()));
    rap_R->Add(rap_SRR, -1.);
    meanY1 = rap_R->GetMean();
    meanY2 = rap_R->GetMean();
  } else {
    // continuum background using different fLSB for the two signal regions
    TH1D* pT_L2 = (TH1D*) pT_L->Clone();
    pT_L->Scale(fracLSB1/(1.*pT_L->Integral()));
    pT_L2->Scale(fracLSB2/(1.*pT_L2->Integral()));
    TH1D* pT_R2 = (TH1D*) pT_R->Clone();
    pT_R->Scale((1.-fracLSB1)/(1.*pT_R->Integral()));
    pT_R2->Scale((1.-fracLSB2)/(1.*pT_R2->Integral()));
    pT_L->Add(pT_R);
    pT_L2->Add(pT_R2);

    TH1D* rap_L2 = (TH1D*) rap_L->Clone();
    rap_L->Scale(fracLSB1/(1.*rap_L->Integral()));
    rap_L2->Scale(fracLSB2/(1.*rap_L2->Integral()));
    TH1D* rap_R2 = (TH1D*) rap_R->Clone();
    rap_R->Scale((1.-fracLSB1)/(1.*rap_R->Integral()));
    rap_R2->Scale((1.-fracLSB2)/(1.*rap_R2->Integral()));
    rap_L->Add(rap_R);
    rap_L2->Add(rap_R2);

    // non prompt continuum background in NPSB
    TH1D* pT_highct_L2 = (TH1D*) pT_highct_L->Clone();
    pT_highct_L->Scale(fracLSB1/(1.*pT_highct_L->Integral()));
    pT_highct_L2->Scale(fracLSB2/(1.*pT_highct_L2->Integral()));
    TH1D* pT_highct_R2 = (TH1D*) pT_highct_R->Clone();
    pT_highct_R->Scale((1.-fracLSB1)/(1.*pT_highct_R->Integral()));
    pT_highct_R2->Scale((1.-fracLSB2)/(1.*pT_highct_R2->Integral()));
    pT_highct_L->Add(pT_highct_R);
    pT_highct_L2->Add(pT_highct_R2);

    TH1D* rap_highct_L2 = (TH1D*) rap_highct_L->Clone();
    rap_highct_L->Scale(fracLSB1/(1.*rap_highct_L->Integral()));
    rap_highct_L2->Scale(fracLSB2/(1.*rap_highct_L2->Integral()));
    TH1D* rap_highct_R2 = (TH1D*) rap_highct_R->Clone();
    rap_highct_R->Scale((1.-fracLSB1)/(1.*rap_highct_R->Integral()));
    rap_highct_R2->Scale((1.-fracLSB2)/(1.*rap_highct_R2->Integral()));
    rap_highct_L->Add(rap_highct_R);
    rap_highct_L2->Add(rap_highct_R2);

    // non prompt background
    pT_NP1->Scale(1./(1.*pT_NP1->Integral()));
    pT_NP2->Scale(1./(1.*pT_NP2->Integral()));
    pT_highct_L->Scale(fBGinNP1/(1.*pT_highct_L->Integral()));
    pT_highct_L2->Scale(fBGinNP2/(1.*pT_highct_L2->Integral()));
    pT_NP1->Add(pT_highct_L, -1.);
    pT_NP2->Add(pT_highct_L2, -1.);

    rap_NP1->Scale(1./(1.*rap_NP1->Integral()));
    rap_NP2->Scale(1./(1.*rap_NP2->Integral()));
    rap_highct_L->Scale(fBGinNP1/(1.*rap_highct_L->Integral()));
    rap_highct_L2->Scale(fBGinNP2/(1.*rap_highct_L2->Integral()));
    rap_NP1->Add(rap_highct_L, -1.);
    rap_NP2->Add(rap_highct_L2, -1.);

    if(PolNP){
      meanPT1 = pT_NP1->GetMean();
      meanPT2 = pT_NP2->GetMean();
      meanY1 = rap_NP1->GetMean();
      meanY2 = rap_NP2->GetMean();
    }

    pT_L->Scale(fBGsig1/(1.*pT_L->Integral()));
    pT_L2->Scale(fBGsig2/(1.*pT_L2->Integral()));
    if(subtractNP){
      pT_NP1->Scale(fNPB1/(1.*pT_NP1->Integral()));
      pT_L->Add(pT_NP1);
      pT_NP2->Scale(fNPB2/(1.*pT_NP2->Integral()));
      pT_L2->Add(pT_NP2);
    }
    pT_PSR1->Add(pT_L, -1.);
    pT_PSR2->Add(pT_L2, -1.);
    meanPT1 = pT_PSR1->GetMean();
    meanPT2 = pT_PSR2->GetMean();

    rap_L->Scale(fBGsig1/(1.*rap_L->Integral()));
    rap_L2->Scale(fBGsig2/(1.*rap_L2->Integral()));
    if(subtractNP){
      rap_NP1->Scale(fNPB1/(1.*rap_NP1->Integral()));
      rap_L->Add(rap_NP1);
      rap_NP2->Scale(fNPB2/(1.*rap_NP2->Integral()));
      rap_L2->Add(rap_NP2);
    }
    rap_PSR1->Add(rap_L, -1.);
    rap_PSR2->Add(rap_L2, -1.);
    meanY1 = rap_PSR1->GetMean();
    meanY2 = rap_PSR2->GetMean();
  }

  std::stringstream meanPTname;
  meanPTname << ";;mean p_{T}";
  TH1D* h_meanPT1 = new TH1D("mean_pT", meanPTname.str().c_str(), 1, 0., 1.);
  h_meanPT1->SetBinContent(1, meanPT1);
  output1->cd();
  h_meanPT1->Write();
  TH1D* h_meanPT2 = new TH1D("mean_pT", meanPTname.str().c_str(), 1, 0., 1.);
  h_meanPT2->SetBinContent(1, meanPT2);
  output2->cd();
  h_meanPT2->Write();

  std::stringstream meanYname;
  meanYname << ";;mean |y|";
  TH1D* h_meanY1 = new TH1D("mean_y", meanYname.str().c_str(), 1, 0., 1.);
  h_meanY1->SetBinContent(1, meanY1);
  output1->cd();
  h_meanY1->Write();
  TH1D* h_meanY2 = new TH1D("mean_y", meanYname.str().c_str(), 1, 0., 1.);
  h_meanY2->SetBinContent(1, meanY2);
  output2->cd();
  h_meanY2->Write();

  // cosTheta and phi
  for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){

    output1->cd();
    hBG_cosThetaPhiL[iFrame]->Write();
    hBG_cosThetaPhiR[iFrame]->Write();
    hBGinNP_cosThetaPhiL[iFrame]->Write();
    hBGinNP_cosThetaPhiR[iFrame]->Write();
    hNPBG1_cosThetaPhi[iFrame]->Write();
    hSR1_cosThetaPhi[iFrame]->Write();

    output2->cd();
    hBG_cosThetaPhiL[iFrame]->Write();
    hBG_cosThetaPhiR[iFrame]->Write();
    hBGinNP_cosThetaPhiL[iFrame]->Write();
    hBGinNP_cosThetaPhiR[iFrame]->Write();
    hNPBG2_cosThetaPhi[iFrame]->Write();
    hSR2_cosThetaPhi[iFrame]->Write();

    std::stringstream nameTBG, nameTBGfolded, nameTBGunfolded;
    nameTBG << "background_costhphi" << onia::frameLabel[iFrame];
    nameTBGfolded << "background_folded_costhphi" << onia::frameLabel[iFrame];
    nameTBGunfolded << "background_unfolded_costhphi" << onia::frameLabel[iFrame];

    // combinatorial background in signal region (prompt region)
    // combination of left and right sideband
    std::stringstream name;
    name << "comb_background_costhphi" << onia::frameLabel[iFrame];
    hBG1_cosThetaPhi[iFrame] = (TH2D *) hBG_cosThetaPhiL[iFrame]->Clone(name.str().c_str());
    hBG2_cosThetaPhi[iFrame] = (TH2D*) hBG_cosThetaPhiL[iFrame]->Clone(name.str().c_str());
    hBG1_cosThetaPhi[iFrame]->Scale(fracLSB1/(1.*hBG1_cosThetaPhi[iFrame]->Integral()));
    hBG2_cosThetaPhi[iFrame]->Scale(fracLSB2/(1.*hBG2_cosThetaPhi[iFrame]->Integral()));
    hBG2_cosThetaPhiR[iFrame] = (TH2D *) hBG_cosThetaPhiR[iFrame]->Clone();
    hBG_cosThetaPhiR[iFrame]->Scale((1.-fracLSB1)/(1.*hBG_cosThetaPhiR[iFrame]->Integral()));
    hBG2_cosThetaPhiR[iFrame]->Scale((1.-fracLSB2)/(1.*hBG2_cosThetaPhiR[iFrame]->Integral()));
    hBG1_cosThetaPhi[iFrame]->Add(hBG_cosThetaPhiR[iFrame]);
    hBG2_cosThetaPhi[iFrame]->Add(hBG2_cosThetaPhiR[iFrame]);
    output1->cd();
    hBG1_cosThetaPhi[iFrame]->Write();
    output2->cd();
    hBG2_cosThetaPhi[iFrame]->Write();
    // total background = combinatorial background
    hTBG1_cosThetaPhi[iFrame] = (TH2D *) hBG1_cosThetaPhi[iFrame]->Clone(nameTBGfolded.str().c_str());
    hTBG2_cosThetaPhi[iFrame] = (TH2D *) hBG2_cosThetaPhi[iFrame]->Clone(nameTBGfolded.str().c_str());

    // combinatorial background in high ctau region
    // combination of left and right sideband in high ctau region
    std::stringstream nameBGinNP;
    nameBGinNP << "background_NPR_costhphi" << onia::frameLabel[iFrame];
    hBGinNP1_cosThetaPhi[iFrame] = (TH2D *) hBGinNP_cosThetaPhiL[iFrame]->Clone(nameBGinNP.str().c_str());
    hBGinNP2_cosThetaPhi[iFrame] = (TH2D *) hBGinNP_cosThetaPhiL[iFrame]->Clone(nameBGinNP.str().c_str());
    hBGinNP1_cosThetaPhi[iFrame]->Scale(fracLSB1/(1.*hBGinNP1_cosThetaPhi[iFrame]->Integral()));
    hBGinNP2_cosThetaPhi[iFrame]->Scale(fracLSB2/(1.*hBGinNP2_cosThetaPhi[iFrame]->Integral()));
    hBGinNP2_cosThetaPhiR[iFrame] = (TH2D *) hBGinNP_cosThetaPhiR[iFrame]->Clone();
    hBGinNP_cosThetaPhiR[iFrame]->Scale((1.-fracLSB1)/(1.*hBGinNP_cosThetaPhiR[iFrame]->Integral()));
    hBGinNP2_cosThetaPhiR[iFrame]->Scale((1.-fracLSB2)/(1.*hBGinNP2_cosThetaPhiR[iFrame]->Integral()));
    hBGinNP1_cosThetaPhi[iFrame]->Add(hBGinNP_cosThetaPhiR[iFrame]);
    hBGinNP2_cosThetaPhi[iFrame]->Add(hBGinNP2_cosThetaPhiR[iFrame]);
    output1->cd();
    hBGinNP1_cosThetaPhi[iFrame]->Write();
    output2->cd();
    hBGinNP2_cosThetaPhi[iFrame]->Write();
    // for non prompt polarization: only use high ctau background
    if(PolNP){
      hTBG1_cosThetaPhi[iFrame] = (TH2D *) hBGinNP1_cosThetaPhi[iFrame]->Clone(nameTBGfolded.str().c_str());
      hTBG2_cosThetaPhi[iFrame] = (TH2D *) hBGinNP2_cosThetaPhi[iFrame]->Clone(nameTBGfolded.str().c_str());
    }
    // non prompt background in high ctau region
    // (hNPBG_cosThetaPhi)_norm - fBGinNP * (hBGinNP_cosThetaPhi)_norm
    hNPBG1_cosThetaPhi[iFrame]->Scale(1./(1.*hNPBG1_cosThetaPhi[iFrame]->Integral()));
    hNPBG2_cosThetaPhi[iFrame]->Scale(1./(1.*hNPBG2_cosThetaPhi[iFrame]->Integral()));
    hBGinNP1_cosThetaPhi[iFrame]->Scale(fBGinNP1/(1.*hBGinNP1_cosThetaPhi[iFrame]->Integral()));
    hBGinNP2_cosThetaPhi[iFrame]->Scale(fBGinNP2/(1.*hBGinNP2_cosThetaPhi[iFrame]->Integral()));
    std::stringstream nameNPS;
    nameNPS << "background_NPSR_costhphi" << onia::frameLabel[iFrame];
    hNPS1_cosThetaPhi[iFrame] = (TH2D *) hNPBG1_cosThetaPhi[iFrame]->Clone(nameNPS.str().c_str());
    hNPS2_cosThetaPhi[iFrame] = (TH2D *) hNPBG2_cosThetaPhi[iFrame]->Clone(nameNPS.str().c_str());
    hNPS1_cosThetaPhi[iFrame] = subtract2D(hNPS1_cosThetaPhi[iFrame], hBGinNP1_cosThetaPhi[iFrame]);
    hNPS2_cosThetaPhi[iFrame] = subtract2D(hNPS2_cosThetaPhi[iFrame], hBGinNP2_cosThetaPhi[iFrame]);
    output1->cd();
    hNPS1_cosThetaPhi[iFrame]->Write();
    output2->cd();
    hNPS2_cosThetaPhi[iFrame]->Write();

    // total background
    // polarization of LSB and RSB
    // bkg = signal contamination in left and right sideband
    if(PolLSB){
      hTBG1_cosThetaPhi[iFrame] = (TH2D *) hSR1_cosThetaPhi[iFrame]->Clone(nameTBGfolded.str().c_str());
      hTBG2_cosThetaPhi[iFrame] = (TH2D *) hSR1_cosThetaPhi[iFrame]->Clone(nameTBGfolded.str().c_str());
    } else if(PolRSB){
      hTBG1_cosThetaPhi[iFrame] = (TH2D *) hSR1_cosThetaPhi[iFrame]->Clone(nameTBGfolded.str().c_str());
      hTBG2_cosThetaPhi[iFrame] = (TH2D *) hSR1_cosThetaPhi[iFrame]->Clone(nameTBGfolded.str().c_str());
    }
    // subtract NP contribution
    else if(subtractNP){
      // fNPBG * (hNPS_cosThetaPhi)_norm + fBGsig * (hBG_cosThetaPhi)_norm
      hBG1_cosThetaPhi[iFrame]->Scale(fBGsig1/(1.*hBG1_cosThetaPhi[iFrame]->Integral()));
      hNPS1_cosThetaPhi[iFrame]->Scale(fNPB1/(1.*hNPS1_cosThetaPhi[iFrame]->Integral()));
      hTBG1_cosThetaPhi[iFrame] = (TH2D *) hNPS1_cosThetaPhi[iFrame]->Clone(nameTBGfolded.str().c_str());
      hTBG1_cosThetaPhi[iFrame]->Add(hBG1_cosThetaPhi[iFrame]);
      hBG2_cosThetaPhi[iFrame]->Scale(fBGsig2/(1.*hBG2_cosThetaPhi[iFrame]->Integral()));
      hNPS2_cosThetaPhi[iFrame]->Scale(fNPB2/(1.*hNPS2_cosThetaPhi[iFrame]->Integral()));
      hTBG2_cosThetaPhi[iFrame] = (TH2D *) hNPS2_cosThetaPhi[iFrame]->Clone(nameTBGfolded.str().c_str());
      hTBG2_cosThetaPhi[iFrame]->Add(hBG2_cosThetaPhi[iFrame]);
    }

    // write folded histogram to file
    output1->cd();
    hTBG1_cosThetaPhi[iFrame]->Write();
    output2->cd();
    hTBG2_cosThetaPhi[iFrame]->Write();

    // get binning of total background histogram
    int nx = hTBG1_cosThetaPhi[iFrame]->GetXaxis()->GetNbins();
    int ny = hTBG1_cosThetaPhi[iFrame]->GetYaxis()->GetNbins();
    int yPhi = ny/4;
    int xCosTheta = nx/2;

    // unfold the total background histogram
    if(folding){

      if(iFrame == 0){
        std::cout << "---------------------------------------------------" << "\n"
                  << "Total background histogram" << "\n"
                  << "number of cosTheta bins: " << nx << "\n"
                  << "number of phi bins: " << ny << "\n"
                  << "phi bins " << 2*yPhi+1 << " to " << 3*yPhi << " are filled." << "\n"
                  << "---------------------------------------------------" << std::endl;
      }

      for (int j = 0; j <= nx; j++){
        for (int k = 2*yPhi+1; k <= 3*yPhi; k++){

          double c1 = hTBG1_cosThetaPhi[iFrame]->GetBinContent(j,k);
          double e1 = hTBG1_cosThetaPhi[iFrame]->GetBinError(j,k);
          double c2 = hTBG2_cosThetaPhi[iFrame]->GetBinContent(j,k);
          double e2 = hTBG2_cosThetaPhi[iFrame]->GetBinError(j,k);

          // flip in cosTheta
          double l = nx + 1 - j;

          // set bin content and error of phiFolded in the other 3 (not yet filled) phi regions
          // 90 - 180: flip phi (upwards), flip cosTheta
          hTBG1_cosThetaPhi[iFrame]->SetBinContent(l,6*yPhi+1-k,c1);
          hTBG1_cosThetaPhi[iFrame]->SetBinError(l,6*yPhi+1-k,e1);
          hTBG2_cosThetaPhi[iFrame]->SetBinContent(l,6*yPhi+1-k,c2);
          hTBG2_cosThetaPhi[iFrame]->SetBinError(l,6*yPhi+1-k,e2);
          // 0 - -90: flip phi (downwards)
          hTBG1_cosThetaPhi[iFrame]->SetBinContent(j,ny+1-k,c1);
          hTBG1_cosThetaPhi[iFrame]->SetBinError(j,ny+1-k,e1);
          hTBG2_cosThetaPhi[iFrame]->SetBinContent(j,ny+1-k,c2);
          hTBG2_cosThetaPhi[iFrame]->SetBinError(j,ny+1-k,e2);
          // -90 - -180: flip cosTheta, shift phi
          hTBG1_cosThetaPhi[iFrame]->SetBinContent(l,k-2*yPhi,c1);
          hTBG1_cosThetaPhi[iFrame]->SetBinError(l,k-2*yPhi,e1);
          hTBG2_cosThetaPhi[iFrame]->SetBinContent(l,k-2*yPhi,c2);
          hTBG2_cosThetaPhi[iFrame]->SetBinError(l,k-2*yPhi,e2);

        }
      }
    } // folding
    // write unfolded and final histogram to file (same histogram twice because of normApproach relic)
    output1->cd();
    hTBG1_cosThetaPhi[iFrame]->SetName(nameTBGunfolded.str().c_str());
    hTBG1_cosThetaPhi[iFrame]->Write();
    hTBG1_cosThetaPhi[iFrame]->SetName(nameTBG.str().c_str());
    hTBG1_cosThetaPhi[iFrame]->Write();

    output2->cd();
    hTBG2_cosThetaPhi[iFrame]->SetName(nameTBGunfolded.str().c_str());
    hTBG2_cosThetaPhi[iFrame]->Write();
    hTBG2_cosThetaPhi[iFrame]->SetName(nameTBG.str().c_str());
    hTBG2_cosThetaPhi[iFrame]->Write();
  } // iFrame

  // in case of normalization approach
  if(normApproach == true){
    for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
      // set binning of background histogram to 16 x 16 if smaller than 16 x 16
      std::stringstream title, nameSR, nameTBG;
      title << ";cos#theta_{"<< onia::frameLabel[iFrame] << "};#phi_{" << onia::frameLabel[iFrame] << "} [deg]";
      nameSR << "hSR_rebinned_cosThetaPhi_" << onia::frameLabel[iFrame];
      nameTBG << "background_costhphi" << onia::frameLabel[iFrame];
      int nx = hTBG1_cosThetaPhi[iFrame]->GetXaxis()->GetNbins();
      int ny = hTBG1_cosThetaPhi[iFrame]->GetYaxis()->GetNbins();
      if(nx < 16) nx = 16;
      if(ny < 16) ny = 16;
      hTBG1_cosThetaPhi[iFrame] = ReSetBin(hTBG1_cosThetaPhi[iFrame], nx, ny, nameTBG, title);
      hTBG2_cosThetaPhi[iFrame] = ReSetBin(hTBG2_cosThetaPhi[iFrame], nx, ny, nameTBG, title);

      if(iFrame == 0){
        std::cout << "----------------------------------------------" << "\n"
                  << "Final binning of total background histogram: " << "\n"
                  << "cosTheta: " << nx << "\n"
                  << "phi: " << ny << "\n"
                  << "----------------------------------------------" << std::endl;
      }

      //------------ Normalization issue: set bins that are not filled in hSR to 0 in hTBG
      hSR1_cosThetaPhi[iFrame] = (TH2D*)hTBG1_cosThetaPhi[iFrame]->Clone(nameSR.str().c_str());
      hSR2_cosThetaPhi[iFrame] = (TH2D*)hTBG2_cosThetaPhi[iFrame]->Clone(nameSR.str().c_str());

      // delete contents of hSR and later fill it with events from prompt signal region
      for (int j = 0; j <= nx; j++){
        for (int k = 0; k <= ny; k++){
          hSR1_cosThetaPhi[iFrame]->SetBinContent(j,k,0);
          hSR2_cosThetaPhi[iFrame]->SetBinError(j,k,0);
        }
      }
    } // iFrame

    // loop through tree and fill hSR histogram
    std::cout << "Filling prompt signal region histogram" << std::endl;
    for(int i = 0; i < n; i++){

      long iEntry = intree->LoadTree(i);
      intree->GetEntry(iEntry);
      if(i % 100000 == 0) std::cout << "entry " << i << " out of " << n << std::endl;

      // tmadlener, 15.04.2016:
      // have to define tha mass here depending on which chic should be used
      double usedChicMass = 0;
      if (useRefittedChic) {
        usedChicMass = chic->M(); // use mass from refitted chic (which just happens to be in the chic variable)
      } else {
        usedChicMass = mQ; // use mQ as mass for the non-refitted chic
      }

      if(chic->Pt() >= onia::pTRange[rapBin-1][ptBin-1] &&
         chic->Pt() < onia::pTRange[rapBin-1][ptBin] &&
         TMath::Abs(chic->Rapidity()) >= onia::rapForPTRange[rapBin-1] &&
         TMath::Abs(chic->Rapidity()) < onia::rapForPTRange[rapBin]){

        if(PolLSB){
          if(usedChicMass > massMinL && usedChicMass < massMaxL && jpsict > PRmin && jpsict < PRmax){
            calcPol(*lepP, *lepN);
            for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
              hSR1_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
              hSR2_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
            }
          }
        } // PolLSB
        else if(PolRSB){
          if(usedChicMass > massMinR && usedChicMass < massMaxR && jpsict > PRmin && jpsict < PRmax){
            calcPol(*lepP, *lepN);
            for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
              hSR1_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
              hSR2_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
            }
          }
        } // PolRSB
        else if(PolNP){
          if(usedChicMass > massMinSR1 && usedChicMass < massMaxSR1 && jpsict > NPmin && jpsict < NPmax){
            calcPol(*lepP, *lepN);
            for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++)
              hSR1_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
          }
          else if(usedChicMass > massMinSR2 && usedChicMass < massMaxSR2 && jpsict > NPmin && jpsict < NPmax){
            calcPol(*lepP, *lepN);
            for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++)
              hSR2_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
          }
        } // PolNP
        else{
          if(usedChicMass > massMinSR1 && usedChicMass < massMaxSR1 && jpsict > PRmin && jpsict < PRmax){
            calcPol(*lepP, *lepN);
            for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++)
              hSR1_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
          }
          else if(usedChicMass > massMinSR2 && usedChicMass < massMaxSR2 && jpsict > PRmin && jpsict < PRmax){
            calcPol(*lepP, *lepN);
            for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++)
              hSR2_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
          }
        } // else
      } // if(onia...)
    } // for(int i ..)

    // set bins in hTBG_cosThetaPhi to 0 when bin is 0 in hSR_cosThetaPhi
    std::cout << "Setting bins in total background histogram to 0 when they are unfilled in prompt signal region histogram:" << std::endl;
    for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){

      output1->cd();
      hSR1_cosThetaPhi[iFrame]->Write();
      output2->cd();
      hSR2_cosThetaPhi[iFrame]->Write();

      int nx = hTBG1_cosThetaPhi[iFrame]->GetXaxis()->GetNbins();
      int ny = hTBG1_cosThetaPhi[iFrame]->GetYaxis()->GetNbins();
      int zeroBins1 = 0, zeroBins2 = 0;

      for (int j = 0; j <= nx; j++){
        for (int k = 0; k <= ny; k++){
          double c1 = hSR1_cosThetaPhi[iFrame]->GetBinContent(j,k);
          double c2 = hTBG1_cosThetaPhi[iFrame]->GetBinContent(j,k);
          double e2 = hTBG1_cosThetaPhi[iFrame]->GetBinError(j,k);
          std::cout << c2 << " " << e2 << std::endl;
          if (c1 == 0 && c2 != 0){
            hTBG1_cosThetaPhi[iFrame]->SetBinContent(j,k,0);
            hTBG1_cosThetaPhi[iFrame]->SetBinError(j,k,0);
            zeroBins1++;
            std::cout << "chic1: bin " << j << ", " << k << " was set to 0" << std::endl;
          }
          double d1 = hSR2_cosThetaPhi[iFrame]->GetBinContent(j,k);
          double d2 = hTBG2_cosThetaPhi[iFrame]->GetBinContent(j,k);
          double f2 = hTBG2_cosThetaPhi[iFrame]->GetBinError(j,k);
          std::cout << d2 << " " << f2 << std::endl;
          if (d1 == 0 && d2 != 0){
            hTBG2_cosThetaPhi[iFrame]->SetBinContent(j,k,0);
            hTBG2_cosThetaPhi[iFrame]->SetBinError(j,k,0);
            zeroBins2++;
            std::cout << "chic2: bin " << j << ", " << k << " was set to 0" << std::endl;
          }
        } // k
      } // j
      std::cout << "In the " << onia::frameLabel[iFrame] << " frame, " << zeroBins1 << " (chic1) and " << zeroBins2 << " (chic2) of " << nx*ny << " bins were set to 0." << std::endl;

      std::stringstream nameTBG;
      nameTBG << "background_costhphi" << onia::frameLabel[iFrame];
      hTBG1_cosThetaPhi[iFrame]->SetName(nameTBG.str().c_str());
      output1->cd();
      hTBG1_cosThetaPhi[iFrame]->Write();
      hTBG2_cosThetaPhi[iFrame]->SetName(nameTBG.str().c_str());
      output2->cd();
      hTBG2_cosThetaPhi[iFrame]->Write();
    }// iFrame
  }// normApproach

  std::cout << "end of this file" << std::endl;
  output1->cd();
  outtree1->Write();
  output2->cd();
  outtree2->Write();
  datafile->Close();
  fitfile->Close();
} // void

//=================================================
// manual subtraction
/*TH3D *subtract3D(TH3D* hist1, TH3D* hist2){

  int nx = hist1->GetXaxis()->GetNbins();
  int ny = hist1->GetYaxis()->GetNbins();
  int nz = hist1->GetZaxis()->GetNbins();

  for (int j = 0; j <= nx; j++){
  for (int k = 0; k <= ny; k++){
  for(int l = 0; l <= nz; l++){

  double c1 = hist1->GetBinContent(j,k,l);
  if (c1 > 0) {
  double c2 = hist2->GetBinContent(j,k,l);
  double c3 = c1 - 1.*c2;
  double e1 = hist1->GetBinError(j,k,l);
  double e2 = hist2->GetBinError(j,k,l);
  double e3 = TMath::Sqrt(e1*e1 + e2*e2);
  if(c3 < 0){
  c3 = 0;
  e3 = 0;
  }
  hist1->SetBinContent(j,k,l,c3);
  hist1->SetBinError(j,k,l,e3);
  }

  } // j
  } // k
  } // l

  return hist1;
  }

  //=================================================
  TH2D *subtract2D(TH2D* hist1, TH2D* hist2){

  int nx = hist1->GetXaxis()->GetNbins();
  int ny = hist1->GetYaxis()->GetNbins();

  for (int j = 0; j <= nx; j++){
  for (int k = 0; k <= ny; k++){
  double c1 = hist1->GetBinContent(j,k);
  if (c1 > 0) {
  double c2 = hist2->GetBinContent(j,k);
  double c3 = c1 - 1.*c2;
  double e1 = hist1->GetBinError(j,k);
  double e2 = hist2->GetBinError(j,k);
  double e3 = TMath::Sqrt(e1*e1 + e2*e2);
  if(c3 < 0){
  c3 = 0;
  e3 = 0;
  }
  hist1->SetBinContent(j,k,c3);
  hist1->SetBinError(j,k, e3);
  }
  }
  }

  return hist1;
  }


  //=================================================
  TH2D* ReSetBin(TH2D* hist, int nBinX, int nBinY, const std::stringstream& name, const std::stringstream& title){
  TH2D *tempHist = (TH2D*)hist->Clone("temp_BG_cosThetaPhiL");
  delete hist;
  hist = new TH2D(name.str().c_str(), title.str().c_str(),
  nBinX, onia::cosTMin, onia::cosTMax, nBinY, onia::phiPolMin, onia::phiPolMax);
  hist->Sumw2();
  TAxis *Xold = tempHist->GetXaxis();
  TAxis *Yold = tempHist->GetYaxis();
  TAxis *Xnew = hist->GetXaxis();
  TAxis *Ynew = hist->GetYaxis();
  for(int binX = 1; binX <= Xnew->GetNbins(); binX++){
  for(int binY = 1; binY <= Ynew->GetNbins(); binY++){
  double centerX = Xnew->GetBinCenter(binX);
  double centerY = Ynew->GetBinCenter(binY);

  //find the corresponding bin and bin error
  double binCont=0.,binErr=0.;
  bool findBin=false;
  for(int BinX = 1; BinX <= Xold->GetNbins(); BinX++){
  for(int BinY = 1; BinY <= Yold->GetNbins(); BinY++){
  double lowX = Xold->GetBinLowEdge(BinX);
  double upX  = Xold->GetBinUpEdge(BinX);
  double lowY = Yold->GetBinLowEdge(BinY);
  double upY  = Yold->GetBinUpEdge(BinY);
  if(centerX > lowX && centerX < upX && centerY > lowY && centerY < upY){
  binCont = tempHist->GetBinContent(BinX,BinY);
  binErr = tempHist->GetBinError(BinX,BinY);
  findBin=true;
  }
  if(findBin) break;
  }//BinY
  }//BinX
  //done
  hist->SetBinContent(binX,binY,binCont);
  hist->SetBinError(binX,binY,binErr);
  }//binY
  }//binX

  return hist;
  }


  //=================================================
  int order(int n){
  int total=1;
  for(int i=0;i<n;i++)
  total=total*2;
  return total;
  }

  int findEvenNum(double number){
  int thisNum=0;
  for(int n=0;n<100;n++){
  if(number >= order(n) && number <= order(n+1)){
  thisNum=order(n+1);
  break;
  }
  }
  return thisNum;
  }
*/
