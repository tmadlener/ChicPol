#include <vector>
#include <string>
#include <sstream>
#include <utility>

#include "rootIncludes.inc"

// ============================== directly copied from bkgHistos_leptonBased
double rapSigma(double p0, double p1, double p2, double rap){
  return p0 + p1 * TMath::Abs(rap) + p2 * TMath::Power(TMath::Abs(rap),2);
}

/**
 * Get the value from the RooRealVar stored in the workspace by name.
 * small helper function for less typing effort. COULDDO: safety measures.
 */
double getVarVal(RooWorkspace* ws, const std::string& name) {
  // for development:
  std::cout << "getVarVal, " << name << ": " << static_cast<RooRealVar*>(ws->var(name.c_str())) << std::endl;

  return static_cast<RooRealVar*>(ws->var(name.c_str()))->getVal();
}

/**
 * Get the pdf stored in the workspace by name
 * small helper function for less typing effort. COULDDO: safety measures.
 */
RooAbsPdf* getPdf(RooWorkspace* ws, const std::string& name) {
  // for develpment:
  std::cout << "getPdf, " << name << ": " << static_cast<RooAbsPdf*>(ws->pdf(name.c_str())) << std::endl;

  return static_cast<RooAbsPdf*>(ws->pdf(name.c_str()));
}

// forward declarations:
void createHistograms(TH2D* hists[], const std::string& nameBase, /*const std::string& titleBase,*/
                      const std::string& suffix, const int nBinsX, const int nBinsY, bool replace = false);

void storeFactor(TFile* file, const std::string& name, const std::string& title, const double val, const double valErr);

std::vector<std::vector<double> > calcCosThetaPhiValues(const TLorentzVector& lepP, const TLorentzVector& lepN, bool folding);

double calcMeanPtRap(TH1D* pT_L, TH1D* pT_R, TH1D* pT_PSR, double fracLSB, double fBGsig);

void fillCosThPhiHistos(TH2D** promptH, TH2D** nonPromptH, const std::vector<std::vector<double> >& values);

void bkgHistos_MCclosure(const std::string& infilename, int rapBin, int ptBin, int nState, bool folding)
{
  // NOTE: expanding the acceptance 15 times (as compared to data) for keeping all MC events
  const double nSigMassMC = onia::nSigMass * 15; // keep this here at the moment to not let it be global!

  const std::string datafilename = "tmpFiles/selEvents_data.root";

  // input
  TFile *datafile = TFile::Open(datafilename.c_str());
  if (!datafile) {
    std::cout << "Inputfile missing" << std::endl;
    return;
  }

  const std::string treename = "selectedData";
  TTree *intree = (TTree *)datafile->Get(treename.c_str());
  TLorentzVector* lepP = NULL;
  TLorentzVector* lepN = NULL;
  TLorentzVector* jpsi = NULL;

  TFile* fitfile = TFile::Open(infilename.c_str());
  if (!fitfile) {
    std::cout << "fitfile is missing" << std::endl;
    return;
  }

  const std::string wsname = "ws_masslifetime";
  RooWorkspace* ws = static_cast<RooWorkspace*>(fitfile->Get(wsname.c_str()));
  if (!ws) {
    std::cout << "workspace not found" << std::endl;
    return;
  }

  gStyle->SetPadRightMargin(0.2);
  gROOT->SetStyle("Plain");
  gStyle->SetTitleBorderSize(0);

  // create output
  std::stringstream outfilename;
  outfilename << "tmpFiles/data_Psi" << nState-3 << "S_rap" << rapBin << "_pT" << ptBin << ".root";
  TFile *output = TFile::Open(outfilename.str().c_str(), "RECREATE");

  TTree *outtree = intree->CloneTree(0);

  // background histos
  TH2D *hBG_cosThetaPhiL[onia::kNbFrames];
  TH2D *hBG_cosThetaPhiR[onia::kNbFrames];
  // TH2D *hBG_cosThetaPhi[onia::kNbFrames];
  TH2D *hNPBG_cosThetaPhi[onia::kNbFrames];
  // TH2D *hNPS_cosThetaPhi[onia::kNbFrames];
  // TH2D *hBGinNP_cosThetaPhi[onia::kNbFrames];
  TH2D *hBGinNP_cosThetaPhiL[onia::kNbFrames];
  TH2D *hBGinNP_cosThetaPhiR[onia::kNbFrames];
  // TH2D *hTBG_cosThetaPhi[onia::kNbFrames];
  // TH2D *hSR_cosThetaPhiL[onia::kNbFrames];
  // TH2D *hSR_cosThetaPhiR[onia::kNbFrames];
  TH2D *hSR_cosThetaPhi[onia::kNbFrames];

  createHistograms(hBG_cosThetaPhiL, "hBG_cosThetaPhi_", "_L", onia::kNbBinsCosT, onia::kNbBinsPhiPol, false);
  createHistograms(hBG_cosThetaPhiR, "hBG_cosThetaPhi_", "_R", onia::kNbBinsCosT, onia::kNbBinsPhiPol, false);
  createHistograms(hNPBG_cosThetaPhi, "hNPBG_cosThetaPhi_", "", onia::kNbBinsCosT, onia::kNbBinsPhiPol, false);
  createHistograms(hBGinNP_cosThetaPhiL, "hBGinNP_cosThetaPhi_", "_L", onia::kNbBinsCosT, onia::kNbBinsPhiPol, false);
  createHistograms(hBGinNP_cosThetaPhiR, "hBGinNP_cosThetaPhi_", "_R", onia::kNbBinsCosT, onia::kNbBinsPhiPol, false);
  createHistograms(hSR_cosThetaPhi, "hSR_cosThetaPhi_", "", onia::kNbBinsCosT, onia::kNbBinsPhiPol, false);

  // mean pT and y histos (background-subtracted)
  int nBins = 100;
  TH1D* pT_L   = new TH1D( "pTLSB", "pTLSB", nBins, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin]);
  TH1D* pT_R   = new TH1D( "pTRSB", "pTRSB", nBins, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin]);
  // TH1D* pT_highct_L   = new TH1D( "pT_highct_LSB", "pT_highct_LSB", nBins, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin]);
  // TH1D* pT_highct_R   = new TH1D( "pT_highcta_RSB", "pT_highct_RSB", nBins, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin]);
  // TH1D* pT_NP   = new TH1D( "pTNP", "pTNP", nBins, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin]);
  TH1D* pT_PSR   = new TH1D( "pTPSR", "pTPSR", nBins, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin]);
  TH1D* rap_L   = new TH1D( "rapLSB", "rapLSB", nBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);
  TH1D* rap_R   = new TH1D( "rapRSB", "rapRSB", nBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);
  // TH1D* rap_highct_L   = new TH1D( "rap_highct_LSB", "rap_highct_LSB", nBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);
  // TH1D* rap_highct_R   = new TH1D( "rap_highct_RSB", "rap_highct_RSB", nBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);
  // TH1D* rap_NP   = new TH1D( "rapNP", "rapNP",nBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);
  TH1D* rap_PSR   = new TH1D( "rapPSR", "rapPSR", nBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);

  //------------------------------------------------------------------------------------------------
  // store pT and y borders
  TVectorD* pTBorder = new TVectorD(1, 2, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin], "END");
  TVectorD* yBorder = new TVectorD(1, 2, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin], "END");
  pTBorder->Print();
  yBorder->Print();
  output->cd();
  pTBorder->Write();
  yBorder->Write();

  // set branches
  intree->SetBranchAddress("lepP", &lepP);
  intree->SetBranchAddress("lepN", &lepN);
  intree->SetBranchAddress("jpsi", &jpsi);
  double jpsict = 0;
  intree->SetBranchAddress("Jpsict", &jpsict);

  // ================================================================================
  // HERE SOME WILD COPY AND PASTING IS HAPPENING AT THE MOMENT. HOPEFULLY THERE WILL BE TIME FOR CLEANUP
  // SOURCES ARE: bkgHistos.C and bkgHistos_leptonBased.C
  // AT THE MOMENT A LOT OF INTERMEDIATE STEPS ARE SKIPPED!!!
  // ================================================================================
  // variables
  RooRealVar* mjpsi = ws->var("JpsiMass");
  RooRealVar* lambdaJpsi = ws->var("bkgLambda_jpsi");
  RooRealVar *CBmassJpsi = ws->var("CBmass_jpsi");
  RooRealVar *CBalphaJpsi = ws->var("CBalpha_jpsi");
  RooRealVar *CBsigmaJpsi = ws->var("CBsigma_jpsi");
  RooRealVar *CBnJpsi = ws->var("CBn_jpsi");

  // pdf
  RooAbsPdf* bkgMassJpsi = getPdf(ws, "bkgMassShape_jpsi");
  RooAbsPdf* sigMassJpsi = getPdf(ws, "sigMassShape_jpsi");

  TF1* funcBGJpsi = static_cast<TF1*>(bkgMassJpsi->asTF(*mjpsi, RooArgList(*lambdaJpsi), *mjpsi));
  TF1* funcSigJpsi = static_cast<TF1*>(sigMassJpsi->asTF(*mjpsi, RooArgList(*CBmassJpsi, *CBsigmaJpsi, *CBalphaJpsi, *CBnJpsi), *mjpsi));

  // regions jpsi
  double jpsiMassMin = onia::massMin;
  double jpsiMassMax = onia::massMax;
  double jpsiMassPeak = getVarVal(ws, "JpsiMass");
  double p0 = getVarVal(ws, "CBsigma_p0_jpsi");
  double p1 = getVarVal(ws, "CBsigma_p1_jpsi");
  double p2 = getVarVal(ws, "CBsigma_p2_jpsi");
  double CBmass_p0 = getVarVal(ws, "CBmass_p0_jpsi");
  double jpsiMassSigma = rapSigma(p0, p1, p2, onia::rapForPTRange[onia::kNbRapForPTBins]);

  double jpsiMassBkgMaxL = CBmass_p0 - onia::nSigBkgLow * jpsiMassSigma;
  double jpsiMassSigMin = CBmass_p0 - onia::nSigMass * jpsiMassSigma;
  double jpsiMassSigMax = CBmass_p0 + onia::nSigMass * jpsiMassSigma;
  double jpsiMassBkgMinR = CBmass_p0 + onia::nSigBkgHigh * jpsiMassSigma;
  std::cout << "-------------------------------------------------------------\n" <<
    "signal region Jpsi: mass window " << jpsiMassSigMin  << " < M < " << jpsiMassSigMax  << " GeV" << std::endl;

  // no background in MC
  storeFactor(output, "background_fraction", ";;fraction of total BG in SR", 0.001, 0);

  double minPt = onia::pTRange[rapBin][ptBin-1];
  double maxPt = onia::pTRange[rapBin][ptBin];
  double minRap = onia::rapForPTRange[rapBin-1];
  double maxRap = onia::rapForPTRange[rapBin];

  // TODO: FINDING HIST BORDERS IN DIMUON KINEMATICS

  TH3D* hBG_pTrapMass_L = new TH3D("hBG_pTrapMass_L", ";p_{T} [GeV/c]; |y|; M [GeV]",
                                   7, minPt, maxPt,
                                   2, minRap, maxRap,
                                   7, jpsiMassSigMin, jpsiMassSigMax);
  hBG_pTrapMass_L->Sumw2();

  TH3D* hBG_pTrapMass_R = new TH3D("hBG_pTrapMass_R", ";p_{T} [GeV/c]; |y|; M [GeV]",
                                   7, minPt, maxPt,
                                   2, minRap, maxRap,
                                   7, jpsiMassSigMin, jpsiMassSigMax);
  hBG_pTrapMass_R->Sumw2();

  TH3D* hSR_pTrapMass = new TH3D("hSR_pTrapMass", ";p_{T} [GeV/c]; |y|; M [GeV]",
                                   7, minPt, maxPt,
                                   2, minRap, maxRap,
                                   7, jpsiMassSigMin, jpsiMassSigMax);
  hSR_pTrapMass->Sumw2();

  // some counter variables
  unsigned nAll = 0;
  unsigned nSR = 0;

  const int nEntries = intree->GetEntries();
  for (int iEntry = 0; iEntry < nEntries; ++iEntry) {
    long entryIdx = intree->LoadTree(iEntry); // NOTE: not entirely sure why we need this!
    intree->GetEntry(entryIdx);
    if(iEntry % 100000 == 0) {std::cout << "entry " << iEntry << " out of " << nEntries << std::endl;}

    double partPt = jpsi->Pt();
    double partAbsRap = TMath::Abs(jpsi->Rapidity());
    double partMass = jpsi->M();

    if (partPt >= onia::pTRange[rapBin][ptBin-1] && partPt < onia::pTRange[rapBin][ptBin] &&
        partAbsRap >= onia::rapForPTRange[rapBin-1] && partAbsRap < onia::rapForPTRange[rapBin]) {
      nAll++;

      if (TMath::Abs(jpsi->M() - CBmass_p0) < nSigMassMC * rapSigma(p0, p1, p2, jpsi->Rapidity())) {
        nSR++;

        // calculate these values here, they are needed in any case afterwards
        const std::vector<std::vector<double> > cosThPhiValues = calcCosThetaPhiValues(*lepP, *lepN, folding);

        if (partMass > jpsiMassSigMin && partMass < jpsiMassSigMax) { // signal region (mass)
          // store TLorentzVectors of the two muons in the given pT and rap cell
          outtree->Fill();
          pT_PSR->Fill(partPt);
          rap_PSR->Fill(partAbsRap);
          hSR_pTrapMass->Fill(partPt, partAbsRap, gRandom->Uniform(jpsiMassSigMin, jpsiMassSigMax));
          fillCosThPhiHistos(hSR_cosThetaPhi, hNPBG_cosThetaPhi, cosThPhiValues);
        } else if (partMass <= jpsiMassSigMin) { // LSB
          pT_L->Fill(partPt);
          rap_L->Fill(partAbsRap);
          hBG_pTrapMass_L->Fill(partPt, partAbsRap, funcBGJpsi->GetRandom(jpsiMassSigMin, jpsiMassSigMax));
          fillCosThPhiHistos(hBG_cosThetaPhiL, hBGinNP_cosThetaPhiL, cosThPhiValues);
        } else if (partMass >= jpsiMassSigMax) { // RSB
          pT_R->Fill(partPt);
          rap_R->Fill(partAbsRap);
          hBG_pTrapMass_R->Fill(partPt, partAbsRap, funcBGJpsi->GetRandom(jpsiMassSigMin, jpsiMassSigMax));
          fillCosThPhiHistos(hBG_cosThetaPhiR, hBGinNP_cosThetaPhiR, cosThPhiValues);
        }
      } // acceptance (rap depending mass region)
    } // pT and rap binning
  } // end loop over entries in input

  std::cout << "--------------------------- \n"
            << "nAll = " << nAll << ", nSR = " << nSR << std::endl;

  // ------------------------------ binning algorithm
  int nBinsPhi = 16, nBinsCosth = 160;
  int totalBins = 0, filledBins = 0;
  for(int binCosth = 0; binCosth < hBG_cosThetaPhiR[2]->GetNbinsX(); binCosth++){
    for(int binPhi = 0; binPhi < hBG_cosThetaPhiR[2]->GetNbinsY(); binPhi++){
      totalBins++;
      //use NP histo (same physical coverage, but better filled, no holes -> better estimate of coverage)
      int binContent = hNPBG_cosThetaPhi[2]->GetBinContent(binCosth+1,binPhi+1);
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
    std::cout << "right bg integral is smaller" << std::endl;
  }

  int nBinsPhiBG = nBinsPhi;
  int nBinsCosthBG = nBinsCosth;

  // calculate average events per-bin cell for background histo
  double Naverage = (double)IntBG/((double)nBinsPhi*nBinsCosth*coverage/2.);
  std::cout << "average cell coverage: " << Naverage << std::endl;

  // if average events per bin is bigger than 10, no rebinning is needed
  if(Naverage > 10){
    std::cout << "Rebinning is not necessary in this case." << "\n"
              << "Ending binning algorithm." << "\n"
              << "------------------------------------------------" << std::endl;
  } else { // otherwise rebin
    std::cout << "------------------------------------------------" << "\n"
              << "old cosTheta binning: " << nBinsCosth << "\n"
              << "old phi binning: " << nBinsPhi << std::endl;

    //set nBinsPhi to the lowest 2^n, such that nBinsPhi > nBinsCosth*coverage/2
    nBinsPhi = findEvenNum(nBinsCosth*coverage/2.);

    std::cout << "closest 2^n number to cosTheta bins: " << nBinsCosth << "\n"
              << "lowest 2^n number so that phi bins > cosTheta bins * coverage/2: " << nBinsPhi << "\n"
              << "------------------------------------------------" << std::endl;

    // set minimum binning
    int nBinsPhiMin = 8;
    int nBinsCosthMin = 8;
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

  // recreate the histograms with new binning
  createHistograms(hBG_cosThetaPhiL, "hBG_cosThetaPhi_", "_L", nBinsCosthBG, nBinsPhiBG, true);
  createHistograms(hBG_cosThetaPhiR, "hBG_cosThetaPhi_", "_R", nBinsCosthBG, nBinsPhiBG, true);
  createHistograms(hBGinNP_cosThetaPhiL, "hBGinNP_cosThetaPhi_", "_L", nBinsCosthBG, nBinsPhiBG, true);
  createHistograms(hBGinNP_cosThetaPhiR, "hBGinNP_cosThetaPhi_", "_R", nBinsCosthBG, nBinsPhiBG, true);
  createHistograms(hSR_cosThetaPhi, "hSR_cosThetaPhi_", "", nBinsCosthBG, nBinsPhiBG, true);
  createHistograms(hNPBG_cosThetaPhi, "hNPBG_cosThetaPhi_", "", nBinsCosthBG, nBinsPhiBG, true);

  // loop again with new binning
  for (int iEntry = 0; iEntry < nEntries; ++iEntry) {
    long entryIdx = intree->LoadTree(iEntry);
    intree->GetEntry(entryIdx);
    if(iEntry % 100000 == 0) std::cout << "entry " << iEntry << " out of " << nEntries << std::endl;

    double partPt = jpsi->Pt();
    double partAbsRap = TMath::Abs(jpsi->Rapidity());
    double partMass = jpsi->M();

    if (partPt >= onia::pTRange[rapBin][ptBin-1] && partPt < onia::pTRange[rapBin][ptBin] &&
        partAbsRap >= onia::rapForPTRange[rapBin-1] && partAbsRap < onia::rapForPTRange[rapBin] &&
        TMath::Abs(jpsi->M() - CBmass_p0) < nSigMassMC * rapSigma(p0, p1, p2, jpsi->Rapidity()) ) {

      const std::vector<std::vector<double> > cosThPhiValues = calcCosThetaPhiValues(*lepP, *lepN, folding);

      if (partMass > jpsiMassSigMin && partMass < jpsiMassSigMax) {
        fillCosThPhiHistos(hSR_cosThetaPhi, hNPBG_cosThetaPhi, cosThPhiValues);
      } else if (partMass <= jpsiMassSigMin) {
        fillCosThPhiHistos(hBG_cosThetaPhiL, hBGinNP_cosThetaPhiL, cosThPhiValues);
      } else if (partMass >= jpsiMassSigMax) {
        fillCosThPhiHistos(hBG_cosThetaPhiR, hBGinNP_cosThetaPhiR, cosThPhiValues);
      }
    }
  } // end loop entries

  // write background histos to file
  output->cd();
  hBG_pTrapMass_L->Write();
  hBG_pTrapMass_R->Write();
  hSR_pTrapMass->Write();

  std::cout << "---------------------------------\n"
            << "Build background histograms\n"
            << "pT-rap-mass" << std::endl;

  // add left and right sideband of low ctau region to combinatorial background histogram
  double fracLSB = getVarVal(ws, "var_fLSBpsi");
  hBG_pTrapMass_L->Scale(fracLSB / (1. * hBG_pTrapMass_L->Integral()));
  hBG_pTrapMass_R->Scale(1. - fracLSB / (1. * hBG_pTrapMass_R->Integral()));
  TH3D* hBG_pTrapMass_lowct = static_cast<TH3D*>(hBG_pTrapMass_L->Clone("comb_background_pTrapMass"));
  hBG_pTrapMass_lowct->Add(hBG_pTrapMass_R);
  output->cd();
  hBG_pTrapMass_lowct->Write();

  // only have signal polarization at the moment (no PolLSB, etc...)
  std::cout << "total background = combinatorial background" << std::endl;
  TH3D* hBG_pTrapMass = static_cast<TH3D*>(hBG_pTrapMass_lowct->Clone("background_pTrapMass"));
  output->cd(); // redundant!
  hBG_pTrapMass->Write();

  // mean pT and rap histos
  double fBGsig = 0.001; // following bkgHistos.C:480
  double meanPT = calcMeanPtRap(pT_L, pT_R, pT_PSR, fracLSB, fBGsig);
  storeFactor(output, "mean_pT", ";;mean p_{T}", meanPT, 0.);

  double meanRap = calcMeanPtRap(rap_L, rap_R, rap_PSR, fracLSB, fBGsig);
  storeFactor(output, "mean_y", ";;mean |y|", meanRap, 0.);

  // outtree->Write();
  output->Write();

} // bkgHistos_MCclosure

/** setup the histograms (i.e. put a valid histogram into each of the pointer of the array).
 * COULDDO: make this a function that returns a vector of TH2D*
 */
void createHistograms(TH2D* hists[], const std::string& nameBase, /*const std::string& titleBase,*/
                      const std::string& suffix, const int nBinsX, const int nBinsY, bool replace)
{
  for (int iFrame = 0; iFrame < onia::kNbFrames; ++iFrame) {
    std::stringstream name;
    name << nameBase << onia::frameLabel[iFrame] << suffix;
    std::stringstream title;
    title << ";cos#theta_{" << onia::frameLabel[iFrame] << "};#phi_{" << onia::frameLabel[iFrame] << "} [deg]";

    if (replace) delete hists[iFrame]; // delete currently stored histogram
    hists[iFrame] = new TH2D(name.str().c_str(), title.str().c_str(),
                             nBinsX, onia::cosTMin, onia::cosTMax,
                             nBinsY, onia::phiPolMin, onia::phiPolMax);
    hists[iFrame]->Sumw2();
  }
}

/** create and store a TH1D in the rootfile pointed to by file */
void storeFactor(TFile* file, const std::string& name, const std::string& title, const double val, const double valErr)
{
  file->cd();
  TH1D* h = new TH1D(name.c_str(), title.c_str(), 1, 0., 1.);
  h->SetBinContent(1, val);
  h->SetBinError(1, valErr);
  h->Write();
  delete h;
}

/** calculate the cosTheta and phi values to be stored in the histograms.
 * Returntype is an onia::kNbFrames x 2 array of doubles, index 0 is CosTh, index 1 is Phi for each frame
 */
std::vector<std::vector<double> > calcCosThetaPhiValues(const TLorentzVector& lepP, const TLorentzVector& lepN, bool folding)
{
  // NOTE: This writes to global variables, which are then used below
  // usage of these variables is restricted to scope of this function in this case
  // TODO: Refactor / Redesign
  ////////////////////////
  calcPol(lepP, lepN);
  ////////////////////////

  std::vector<std::vector<double> > returnVals;
  for (int iFrame = 0; iFrame < onia::kNbFrames; ++iFrame) {
    std::vector<double> frameVals;
    double phi = thisPhi[iFrame]; // only usage of global vars!
    double cosTh = thisCosTh[iFrame]; // only usage of global vars!

    if (folding) { // map all values into first (?) octant
      double phiFolded = phi;
      double thetaAdjusted = cosTh;
      if (phi > -90. && phi < 0.) {
        phiFolded *= -1;
      } else if (phi > 90. && phi < 180.) {
        phiFolded = 180. - phi;
        thetaAdjusted *= -1;
      } else if (phi > -180. && phi < -90.) {
        phiFolded = 180. + phi;
        thetaAdjusted *= -1;
      }
      // asign again to values that will be stored
      phi = phiFolded;
      cosTh = thetaAdjusted;
    }

    frameVals.push_back(cosTh);
    frameVals.push_back(phi);

    returnVals.push_back(frameVals);
  }

  return returnVals;
}

/** scale the passed histogram to the Integral of itself.
 * Mainly for less typing effort.
 */
inline void selfScale(TH1D* h, double f = 1.0)
{
  h->Scale(f / (1. * h->Integral()));
}

/** calculate the mean pT or rap that should be stored.
 * NOTE: modifies the passed histograms (scaling and adding them together)
 * TODO: check if this actually does what is expected. The code in bkgHistos.C and bkgHistos_leptonBased.C is
 * somewhat conflicting
 */
double calcMeanPtRap(TH1D* h_L, TH1D* h_R, TH1D* h_PSR, double fracLSB, double fBGsig)
{
  selfScale(h_L, fracLSB);
  selfScale(h_R, 1. - fracLSB);
  h_L->Add(h_R);
  selfScale(h_L, fBGsig);
  // tmadlener: why no "normalization" of h_PSR?
  h_PSR->Add(h_L, -1.);

  return h_PSR->GetMean();
}

/** fill the costheta phi histograms.*/
void fillCosThPhiHistos(TH2D** promptH, TH2D** nonPromptH, const std::vector<std::vector<double> >& values)
{
  for (int iFrame = 0; iFrame < onia::kNbFrames; ++iFrame) {
    promptH[iFrame]->Fill(values[iFrame][0], values[iFrame][1]);
    nonPromptH[iFrame]->Fill(values[iFrame][0], values[iFrame][1]);
  }
}
