#ifndef BKGHISTOPRODUCER_H__
#define BKGHISTOPRODUCER_H__

#include "commonVar.h"

#include "BkgHistoProducer_helper.h"
#include "bkgHistos_helper.h"

// all the root stuff
#include "rootIncludes.inc"

// stl
#include <string>
#include <iostream>
#include <vector>
#include <cmath>

/** Interface class for the templatized version of the actual class, to have a common interface. */
class IBkgHistoProducer {
public:
  /** initialize everything that is necessary. Call for every pt and rapidity bin separately.*/
  virtual void initialize(const std::string& infileName, const int rapBin, const int ptBin, bool MC,
                          const int FracLSB, bool useRefittedChic) = 0;
  /** do the actual work and fill the histograms. Call for every pt and rapidity bin separately.*/
  virtual void fillHistos(const int rapBin, const int ptBin, bool useRefittedChic, bool MC, bool PolLSB,
                          bool PolRSB, bool PolNP, bool folding) = 0;
  /** store the appropriate histograms to the output files. Call for every pt and rapidity bin separately. */
  virtual void storeHistos(bool PolLSB, bool PolRSB, bool PolNP, bool subtractNP, bool folding) = 0;
  virtual ~IBkgHistoProducer() = 0; /**< destructor. */
};

inline IBkgHistoProducer::~IBkgHistoProducer() { /* nothing to do here */ }

/** enum for spezialization of the BkgHistoProducer. Each State gets its own enum value.
 * As soon as c++11 gets available make this an enum class.
 */
enum /*class*/ StateT {
  Jpsi = 4,
  Psi2S = 5,
  Chic = 6,
  undefined = -1
};

StateT stateTfromInt(const int state)
{
  switch (state) {
  case 4: return Jpsi;
  case 5: return Psi2S;
  case 6: return Chic;
  default: return undefined;
  }
}

/** ostream operator overload for StateT enum. */
std::ostream& operator<<(std::ostream& os, const StateT& state)
{
  switch (state) {
  case Jpsi:
    os << "4 - Jpsi";
    break;
  case Psi2S:
    os << "5 - Psi2S";
    break;
  case Chic:
    os << "6 - Chic";
    break;
  default:
    os << "-1 - undefined";
  }
  return os;
}

/**
 * BkgHistoProducer class that does all the work, including filehandling, etc...
 * Only initialize and fillHistos are publicly available and will be specialized for every state, to acommodate
 * for the differences.
 */
template<StateT State>
class BkgHistoProducer : public IBkgHistoProducer {
public:
  /** Constructor that does all the work that can be done for initialization independently of the State. */
  BkgHistoProducer();
  /**
   * initialize the fit file and setup all variables needed internally for the given State. Call for every
   * pt and rapidity bin separately.
   */
  virtual void initialize(const std::string& infileName, const int rapBin, const int ptBin, bool MC,
                          const int FracLSB, bool useRefittedChic) /*override*/;
  /** fill all needed histograms and write them to file. Call for every pt and rapitdity bin separately. */
  virtual void fillHistos(const int rapBin, const int ptBin, bool useRefittedChic, bool MC, bool PolLSB,
                          bool PolRSB, bool PolNP, bool folding) /*override*/;
  /** store all the filled histograms, and combine them together to have the correct output histograms. */
  virtual void storeHistos(bool PolLSB, bool PolRSB, bool PolNP, bool subtractNP, bool folding) /*override*/;
  /** Destructor, do all the state independent clean-up work that is needed. */
  ~BkgHistoProducer();

private:
  TFile* m_dataFile; /**< pointer to the root file containing all data. */
  TTree* m_inTree; /**< pointer to the TTree in the dataFile containing all data. */
  TFile* m_fitFile; /**< pointer to the root file containing the fit data. */
  RooWorkspace* m_ws; /**< RooWorkspace containing the fit results. */
  BkgHistoRootVars m_inputVars; /**< input root variables container. */
  TLorentzVector* m_particle; /**< 4 momentum of the particle that is actually used in chic case. */
  std::vector<TFile*> m_outFiles; /**< vector containing the rootfiles created by this class. */
  std::vector<TTree*> m_outTrees; /**< vector containing the TTrees that are stored in the output files */
  std::vector<BkgHistoCosThetaHists> m_cosThetaPhi; /**< vector containing the TH2D histo container. */
  BkgHisto1DHists m_ptHists; /**< the pT TH1D histogram container. */
  BkgHisto1DHists m_rapHists; /**< the rap TH1D histogram container. */
  std::vector<BkgHistoPtRapMassHists> m_ptRapMass; /**< vector for the pTRapMass histos. */
  NamedVarStore<double> m_fitVars; /**< storage for the different double values, like the fit variables. */
  /** storage for the different distributions from which random numbers are drawn. */
  NamedVarStore<TF1*> m_randDists;

  /** Initialization that is common to every state, like setup of the workspace, etc... */
  void initCommon(const std::string& infilename);
  /**
   * store the pT and y borders in every output file that is currently open.
   * NOTE: call this only after the outpufiles have been opened!
   */
  void storePtRapBorders(const int rapBin, const int ptBin);

  /** add a new outputile to the m_outFiles vector. */
  void addOutputFile(const std::string& filename);

  /**
   * Setup all TH2D cosThPhi histograms that are filled in the fillHistos function.
   * NOTE: when called for the first time, the histos get created as well
   * NOTE: leaves the histograms that are present in BkgHistoCosThetaHists but not used in fillHistos untouched
   * @param nBCT, number of bins for cosTheta
   * @param nBP, number of bins for phi
   * TODO: refactor this into a function that needs no spezialization (should be possible by passing another
   * argument)!
   */
  void setupCosThPhiHists(const int nBCT, const int NBP);

  /**
   * Setup all the TH3Ds that are needed before storing the final results (i.e. intermediate histos).
   * Creates nHists BkgHistoPtRapMass containers in the m_ptRapMass vector.
   * NOTE: needs input from fit results, ergo call this after the fit results are present.
   * TODO: check for appropriate
   */
  void setupPtRapMassHists(const size_t nHists,
                           const int ptBins, const double minPt, const double maxPt,
                           const int rapBins, const double minRap, const double maxRap,
                           const int massBins, const double minMass, const double maxMass);

  /**
   * Initialize the Fit Variables store by retrieving all necessary values from the RooWorkspace or setting them
   * to the appropriate values for MC closure. Also defines som other variables that would else be "global";
   */
  void setupFitVariables(const int rapBin, const int ptBin, bool MC, const int FracLSB);

  /** Setup the TF1s that are needed for filling the mass histograms */
  void setupRandomDists(const int rapBin, const int ptBin);

  /** Store all values (+errors) in TH1Ds that are needed in the polFit framework (?). */
  void store1DFactors(bool MC);

  /** Determine the pT and rapidity borders in dimuon kinematics in case of chi kinematics. */
  void findRapPtBordersDiMuonKin(const int rapBin, const int ptBin);

  /**
   * Fill the TH1D histograms. Chic case.
   * NOTE: The bools RSB, LSB, NP and MC have to be mutually exclusive if one of them is set! (They are queried
   * in a chain of if-else clauses in this order).
   */
  void fill1DHists(bool RSB, bool LSB, bool NP, bool MC, const double pt, const double y, const double mass,
                   const BkgHistoRangeReport& eventRegion);

  /**
   * Fill the PtRapMass TH3Ds.
   * NOTE: mass can probably be dropped! (All the needed information is already in the eventRegion)
   */
  void fill3DHists(bool RSB, bool LSB, bool NP, bool MC, const double pt, const double absY, const double mass,
                   const BkgHistoRangeReport& eventRegion);

  /**
   * Fill the cosThetaPhi histograms.
   */
  void fill2DHists(const std::vector<std::vector<double> >& cosThPVals, const BkgHistoRangeReport& eventRegion);

  /** store the pTRapMass Histos. */
  void store3DHists(bool PolLSB, bool PolRSB, bool PolNP, bool subtractNP);

  /**
   * store the pT and |y| histograms.
   * NOTE: actually stores just mean pt and mean |y|.
   * TODO: Refactor
   */
  void store1DHists(bool PolLSB, bool PolRSB, bool PolNP, bool subtractNP);

  /**
   * store the cosTheta Phi histograms (and do some combinations)
   * TODO: Refactor
   */
  void store2DHists(bool PolLSB, bool PolRSB, bool PolNP, bool subtractNP, bool folding);

  typedef typename std::vector<TFile*>::iterator TFileIt; /**< private typedef for easier looping over files. */
  typedef typename std::vector<TTree*>::iterator TTreeIt; /**< private typedef for easier looping over TTrees. */
};

// ================================================================================
//                               IMPLEMENTATION
// ================================================================================
template<StateT State>
BkgHistoProducer<State>::BkgHistoProducer() : m_particle(NULL)
{
  std::cout << "---------- CONSTRUCTOR, STATE == " << State << std::endl;

  const std::string dataFileName = "tmpFiles/selEvents_data.root";
  m_dataFile = TFile::Open(dataFileName.c_str());
  if (!m_dataFile) {
    std::cerr << "Inputfile: \'" << dataFileName << "\' cannot be opened in BkgHistoProducer()" << std::endl;
  } else {
    std::cout << "opened file: \'" << dataFileName << "\'." << std::endl;
  }
  const std::string treeName = "selectedData";
  m_inTree = static_cast<TTree*>(m_dataFile->Get(treeName.c_str()));

  gStyle->SetPadRightMargin(0.2);
  gROOT->SetStyle("Plain");
  gStyle->SetTitleBorderSize(0);

  std::cout << "********** CONSTRUCTOR, STATE == " << State << std::endl;
}

template<StateT State>
BkgHistoProducer<State>::~BkgHistoProducer()
{
  std::cout << "---------- DESTRUCTOR, STATE == " << State << std::endl;
  // close the input files (the way I understand root is, that this also destroys these objects)
  if (m_dataFile) { m_dataFile->Close(); }
  if (m_fitFile) { m_fitFile->Close(); }

  // write and close all output files
  for (size_t iFile = 0; iFile < m_outFiles.size(); ++iFile) {
    if (m_outFiles[iFile]) {
      m_outFiles[iFile]->cd();
      m_outTrees[iFile]->Write();
      m_outFiles[iFile]->Write();
      // NOTE: due to some strange ROOT behavior this call also seems to free some objects that are "attatched" to
      // this file. At the moment I have no solution for this and simply commented all deletes that cause problems
      // in any of the used helper structs.
      m_outFiles[iFile]->Close();
      delete m_outFiles[iFile];
    }
  }

  // delete m_particle; // DO NOT DO THIS! The ressource this pointer is pointing too is freed elsewhere!

  std::cout << "********** DESTRUCTOR, STATE == " << State << std::endl;
}

// ================================================================================
//        SETUP COS THETA PHI HISTOGRAMS (declaration needed before initialize!)
// ================================================================================
// only chic needs spezialized version
template<>
void BkgHistoProducer<Chic>::setupCosThPhiHists(const int nBCT, const int nBP)
{
  std::cout << "---------- IN BkgHistoProducer<Chic>::setupCosThPhiHists" << std::endl;

  // chic needs two sets of cosThetaPhi histos (create only in first call)
  if (m_cosThetaPhi.empty()) {
    m_cosThetaPhi.push_back(BkgHistoCosThetaHists(onia::kNbFrames));
    m_cosThetaPhi.push_back(BkgHistoCosThetaHists(onia::kNbFrames));
  }

  const int ctMin = onia::cosTMin;
  const int ctMax = onia::cosTMax;
  const int pMin = onia::phiPolMin;
  const int pMax = onia::phiPolMax;
  const char** frameLbls = onia::frameLabel;

  for (size_t iS = 0; iS < m_cosThetaPhi.size(); ++iS) {
    createHists(m_cosThetaPhi[iS].hBG_L, "hBG_cosThetaPhi_", frameLbls, "_L", nBCT, nBP, ctMin, ctMax, pMin, pMax);
    createHists(m_cosThetaPhi[iS].hBG_R, "hBG_cosThetaPhi_", frameLbls, "_R", nBCT, nBP, ctMin, ctMax, pMin, pMax);
    createHists(m_cosThetaPhi[iS].hNPBG, "hNPBG_cosThetaPhi_", frameLbls, "", nBCT, nBP, ctMin, ctMax, pMin, pMax);
    createHists(m_cosThetaPhi[iS].hBGinNP_L, "hBGinNP_cosThetaPhi_", frameLbls, "_L", nBCT, nBP, ctMin, ctMax, pMin, pMax);
    createHists(m_cosThetaPhi[iS].hBGinNP_R, "hBGinNP_cosThetaPhi_", frameLbls, "_R", nBCT, nBP, ctMin, ctMax, pMin, pMax);
    createHists(m_cosThetaPhi[iS].hSR, "hSR_cosThetaPhi_", frameLbls, "", nBCT, nBP, ctMin, ctMax, pMin, pMax);
  }

  std::cout << "********** IN BkgHistoProducer<Chic>::setupCosThPhiHists" << std::endl;
}

// all others can be dealt with in the same way
template<StateT State>
void BkgHistoProducer<State>::setupCosThPhiHists(const int nBCT, const int nBP)
{
  std::cout << "---------- IN BkgHistoProducer<" << State << ">::setupCosThPhiHists" << std::endl;

  // create all needed histos (only in first call)
  if (m_cosThetaPhi.empty()) {
    m_cosThetaPhi.push_back(BkgHistoCosThetaHists(onia::kNbFrames));
  }

  const int ctMin = onia::cosTMin;
  const int ctMax = onia::cosTMax;
  const int pMin = onia::phiPolMin;
  const int pMax = onia::phiPolMax;
  const char** frameLbls = onia::frameLabel;

  // set the necessary hists up properly
  createHists(m_cosThetaPhi[0].hBG_L, "hBG_cosThetaPhi_", frameLbls, "_L", nBCT, nBP, ctMin, ctMax, pMin, pMax);
  createHists(m_cosThetaPhi[0].hBG_R, "hBG_cosThetaPhi_", frameLbls, "_R", nBCT, nBP, ctMin, ctMax, pMin, pMax);
  createHists(m_cosThetaPhi[0].hNPBG, "hNPBG_cosThetaPhi_", frameLbls, "", nBCT, nBP, ctMin, ctMax, pMin, pMax);
  createHists(m_cosThetaPhi[0].hBGinNP_L, "hBGinNP_cosThetaPhi_", frameLbls, "_L", nBCT, nBP, ctMin, ctMax, pMin, pMax);
  createHists(m_cosThetaPhi[0].hBGinNP_R, "hBGinNP_cosThetaPhi_", frameLbls, "_R", nBCT, nBP, ctMin, ctMax, pMin, pMax);
  createHists(m_cosThetaPhi[0].hSR, "hSR_cosThetaPhi_", frameLbls, "", nBCT, nBP, ctMin, ctMax, pMin, pMax);

  std::cout << "********** IN BkgHistoProducer<" << State << ">::setupCosThPhiHists" << std::endl;
}

// ================================================================================
//                            SETUP FIT VARIABLES
// ================================================================================
// this definetely needs spezialization
template<>
void BkgHistoProducer<Chic>::setupFitVariables(const int rapBin, const int ptBin, bool MC, const int FracLSB)
{
  // NOTE: This is not yet done and no MC closure initialization is done at the moment. This is just a mere
  // "translation" from 'bkgHistos_leptonBased.C'. There is probably some unecessary storage done here
  // TODO: implement MC closure
  // TODO: find a way to implement this more readable. At the moment this is just a very long list of variables

  m_fitVars.setFromWS(m_ws, "var_fLSBChic1", "fracLSB1");
  m_fitVars.setFromWS(m_ws, "var_fLSBChic2", "fracLSB2");
  m_fitVars.setFromWS(m_ws, "var_fLSBpsi", "fracLSBpsi");
  if (FracLSB >= 0) { // NOTE: no sanity check for FracLSB > 100 !
    m_fitVars.set("fracLSB1", (double)FracLSB / 100.);
    m_fitVars.set("fracLSB2", (double)FracLSB / 100.);
    m_fitVars.set("fracLSBpsi", (double)FracLSB / 100.);
  }

  m_fitVars.setFromWS(m_ws, "var_fracBackgroundInPRSR1", "fBGsig1"); // comb. background in PRSR1
  m_fitVars.setFromWS(m_ws, "var_fracBackgroundInPRSR2", "fBGsig2"); // comb. background in PRSR2
  m_fitVars.setFromWS(m_ws, "var_fracBackgroundInNPSR1", "fBGinNP1"); // comb. background in NPSR1
  m_fitVars.setFromWS(m_ws, "var_fracBackgroundInNPSR2", "fBGinNP2"); // comb. background in NPSR2
  m_fitVars.setFromWS(m_ws, "var_fracNPChic1InPRSR1", "fNPB1"); // non prompt background in PRSR1
  m_fitVars.setFromWS(m_ws, "var_fracNPChic2InPRSR2", "fNPB2"); // non prompt background in PRSR2
  m_fitVars.set("fTBGsig1", m_fitVars["fBGsig1"]); // total background fraction
  m_fitVars.set("fTBGsig2", m_fitVars["fBGsig2"]); // total background fraction
  m_fitVars.set("fP1", 1. - m_fitVars["fBGsig1"] - m_fitVars["fNPB1"]);
  m_fitVars.set("fP2", 1. - m_fitVars["fBGsig2"] - m_fitVars["fNPB2"]);

  m_fitVars.setFromWS(m_ws, "var_fracPRChic1InPRLSB", "fSRinPLSB"); // prompt chic1 contamination in LSB
  m_fitVars.setFromWS(m_ws, "var_fracPRChic2InPRRSB", "fSRinPRSB"); // prompt chic2 contamination in RSB
  m_fitVars.setFromWS(m_ws, "var_fractionCombBGofTotalBackground", "fComBgToTotBg");
  // fraction of combinatorial Jpsi + gamma events to total background
  m_fitVars.setFromWS(m_ws, "var_fractionJpsiBGofTotalBackground", "fJpsiBgToTotBg");

  // errors on fractions
  m_fitVars.set("fBGerr1", getVarError(m_ws, "var_fracBackgroundInPRSR1"));
  m_fitVars.set("fBGerr2", getVarError(m_ws, "var_fracBackgroundInPRSR2"));
  m_fitVars.set("fNPerr1", getVarError(m_ws, "var_fracBackgroundInNPSR1"));
  m_fitVars.set("fNPerr2", getVarError(m_ws, "var_fracBackgroundInNPSR2"));
  m_fitVars.set("fBGinNPerr1", getVarError(m_ws, "var_fracNPChic1InPRSR1"));
  m_fitVars.set("fBGinNPerr2", getVarError(m_ws, "var_fracNPChic1InPRSR2"));
  m_fitVars.set("fTBGerr1", m_fitVars["fBGerr1"]);
  m_fitVars.set("fTBGerr2", m_fitVars["fBGerr2"]);
  m_fitVars.set("fPerr1", TMath::Sqrt(TMath::Power(m_fitVars["fNPerr1"], 2) + TMath::Power(m_fitVars["fBGerr1"], 2)));
  m_fitVars.set("fPerr2", TMath::Sqrt(TMath::Power(m_fitVars["fNPerr2"], 2) + TMath::Power(m_fitVars["fBGerr2"], 2)));

  std::cout << "-------------------------------------------------------------\n"
            << "fraction of left sideband: " << m_fitVars["fracLSBpsi"] << " (Jpsi), " << m_fitVars["fracLSB1"]
            << " (chic1), " << m_fitVars["fracLSB2"] << " (chic2)\n" << "combinatorial background fraction: "
            << m_fitVars["fBGsig1"] << " (chic1), " << m_fitVars["fBGsig2"] << " (chic2)\n"
            << "nonprompt background fraction: " << m_fitVars["fNPB1"] << " (chic1), " << m_fitVars["fNPB2"]
            << " (chic2)\n" << "total background fraction: " << m_fitVars["fTBGsig1"] << " (chic1), "
            << m_fitVars["fTBGsig2"] << " (chic2)" << std::endl;

  // TODO: chic kinematic dependent pt/rap min/max values
  m_fitVars.set("jpsiMassMin", onia::massMin);
  m_fitVars.set("jpsiMassMax", onia::massMax);
  m_fitVars.setFromWS(m_ws, "JpsiMass", "jpsiMassPeak");
  m_fitVars.setFromWS(m_ws, "CBsigma_p0_jpsi", "p0");
  m_fitVars.setFromWS(m_ws, "CBsigma_p1_jpsi", "p1");
  m_fitVars.setFromWS(m_ws, "CBsigma_p2_jpsi", "p2");
  m_fitVars.set("jpsiMassSigma", rapSigma(m_fitVars["p0"], m_fitVars["p1"], m_fitVars["p2"],
                                          onia::rapForPTRange[onia::kNbRapForPTBins]));
  m_fitVars.setFromWS(m_ws, "CBmass_p0_jpsi", "CBmass_p0");
  m_fitVars.set("jpsiMassBkgMaxL", m_fitVars["CBmass_p0"] - onia::nSigBkgLow * m_fitVars["jpsiMassSigma"]);
  m_fitVars.set("jpsiMassSigMin", m_fitVars["CBmass_p0"] - onia::nSigMass * m_fitVars["jpsiMassSigma"]);
  m_fitVars.set("jpsiMassSigMax", m_fitVars["CBmass_p0"] + onia::nSigMass * m_fitVars["jpsiMassSigma"]);
  m_fitVars.set("jpsiMassBkgMinR", m_fitVars["CBmass_p0"] + onia::nSigBkgHigh * m_fitVars["jpsiMassSigma"]);

  std::cout << "-------------------------------------------------------------\n"
            << "signal region Jpsi: mass window " << m_fitVars["jpsiMassSigMin"] << " < M "
            << m_fitVars["jpsiMassSigMax"] << " GeV" << std::endl;


  m_fitVars.set("minPt", onia::pTRange[rapBin][ptBin-1]);
  m_fitVars.set("maxPt", onia::pTRange[rapBin][ptBin]);
  m_fitVars.set("minRap", onia::rapForPTRange[rapBin-1]);
  m_fitVars.set("maxRap", onia::rapForPTRange[rapBin]);

  m_fitVars.set("massMinL", onia::massChiSBMin);
  m_fitVars.set("massMaxR", onia::massChiSBMax);
  m_fitVars.setFromWS(m_ws, "var_lsbMaxMass", "massMaxL");
  m_fitVars.setFromWS(m_ws, "var_rsbMinMass", "massMinR");
  m_fitVars.setFromWS(m_ws, "var_sig1MinMass", "massMinSR1");
  m_fitVars.setFromWS(m_ws, "var_sig2MinMass", "massMinSR2");
  m_fitVars.setFromWS(m_ws, "var_sig1MaxMass", "massMaxSR1");
  m_fitVars.setFromWS(m_ws, "var_sig2MaxMass", "massMaxSR2");
  m_fitVars.setFromWS(m_ws, "var_PRMin", "PRmin");
  m_fitVars.setFromWS(m_ws, "var_PRMax", "PRmax");
  m_fitVars.setFromWS(m_ws, "var_NPMin", "NPmin");
  m_fitVars.setFromWS(m_ws, "var_NPMax", "NPmax");

  std::cout << "-------------------------------------------------------------\n" <<
    "left  sideband: mass window " << m_fitVars["massMinL"]  << " < M < " << m_fitVars["massMaxL"]  << " GeV\n" <<
    "right sideband: mass window " << m_fitVars["massMinR"]  << " < M < " << m_fitVars["massMaxR"]  << " GeV\n" <<
    "signal region chic1: mass window " << m_fitVars["massMinSR1"]  << " < M < " << m_fitVars["massMaxSR1"]   << " GeV\n" <<
    "signal region chic2: mass window " << m_fitVars["massMinSR2"]  << " < M < " << m_fitVars["massMaxSR2"]   << " GeV\n" <<
    "prompt region: " << m_fitVars["PRmin"] << " < ctau < " << m_fitVars["PRmax"] << "\n"
    "nonprompt region: " << m_fitVars["NPmin"] << " < ctau < " << m_fitVars["NPmax"] << "\n"
    "-------------------------------------------------------------\n" << std::endl;


}

// all other states (i.e. not chic) currently handled by the same function
template<StateT State>
void BkgHistoProducer<State>::setupFitVariables(const int rapBin, const int ptBin, bool MC, const int FracLSB)
{
  // NOTE: this contains only the variables for MC closure at the moment. (I.e. only the ones from the massfit
  // are present). The variables should also be present for data! Only others will not be here!
  // TODO: implement others (see bkgHistos.C)
  // COULDDO: store some repeatedly used variables in local variables to avoid a possibly costly look-up everytime
  m_fitVars.set("jpsiMassMin", onia::massMin);
  m_fitVars.set("jpsiMassMax", onia::massMax);
  m_fitVars.setFromWS(m_ws, "JpsiMass", "jpsiMassPeak");
  m_fitVars.setFromWS(m_ws, "CBsigma_p0_jpsi", "p0");
  m_fitVars.setFromWS(m_ws, "CBsigma_p1_jpsi", "p1");
  m_fitVars.setFromWS(m_ws, "CBsigma_p2_jpsi", "p2");
  m_fitVars.set("jpsiMassSigma", rapSigma(m_fitVars["p0"], m_fitVars["p1"], m_fitVars["p2"],
                                          onia::rapForPTRange[onia::kNbRapForPTBins]));
  m_fitVars.setFromWS(m_ws, "CBmass_p0_jpsi", "CBmass_p0");
  m_fitVars.set("jpsiMassBkgMaxL", m_fitVars["CBmass_p0"] - onia::nSigBkgLow * m_fitVars["jpsiMassSigma"]);
  m_fitVars.set("jpsiMassSigMin", m_fitVars["CBmass_p0"] - onia::nSigMass * m_fitVars["jpsiMassSigma"]);
  m_fitVars.set("jpsiMassSigMax", m_fitVars["CBmass_p0"] - onia::nSigMass * m_fitVars["jpsiMassSigma"]);
  m_fitVars.set("jpsiMassBkgMinR", m_fitVars["CBmass_p0"] - onia::nSigBkgHigh * m_fitVars["jpsiMassSigma"]);

  m_fitVars.set("minPt", onia::pTRange[rapBin][ptBin-1]);
  m_fitVars.set("maxPt", onia::pTRange[rapBin][ptBin]);
  m_fitVars.set("minRap", onia::rapForPTRange[rapBin-1]);
  m_fitVars.set("maxRap", onia::rapForPTRange[rapBin]);
}

// ================================================================================
//                               SETUP RANDOM DISTS
// ================================================================================
template<>
void BkgHistoProducer<Chic>::setupRandomDists(const int rapBin, const int ptBin)
{
  // variables
  RooRealVar* m = m_ws->var("chicMass");
  RooRealVar* mjpsi = m_ws->var("JpsiMass");
  RooRealVar *poly1 = m_ws->var("BK_p1");
  RooRealVar *poly2 = m_ws->var("BK_p2");
  RooRealVar *CBmass1 = m_ws->var("CBmass1");
  RooRealVar *CBmass2 = m_ws->var("CBmass2");
  RooRealVar *CBsigma1 = m_ws->var("CBsigma1");
  RooRealVar *CBsigma2 = m_ws->var("CBsigma2");
  RooRealVar *CBalpha1 = m_ws->var("CBalpha1");
  RooRealVar *CBalpha2 = m_ws->var("CBalpha2");
  RooRealVar *CBn1 = m_ws->var("CBn");
  RooRealVar *CBn2 = m_ws->var("CBn2");
  RooRealVar *lambdaJpsi = m_ws->var("bkgLambda_jpsi");
  RooRealVar *CBmassJpsi = m_ws->var("CBmass_jpsi");
  RooRealVar *CBalphaJpsi = m_ws->var("CBalpha_jpsi");
  RooRealVar *CBsigmaJpsi = m_ws->var("CBsigma_jpsi");
  RooRealVar *CBnJpsi = m_ws->var("CBn_jpsi");

  // pdf
  RooAbsPdf *bkgMass = getPdf(m_ws, "M_background");
  RooAbsPdf *signalMass1 = getPdf(m_ws, "M_chic1");
  RooAbsPdf *signalMass2 = getPdf(m_ws, "M_chic2");
  RooAbsPdf *bkgMassJpsi = getPdf(m_ws, "bkgMassShape_jpsi");
  RooAbsPdf *sigMassJpsi = getPdf(m_ws, "sigMassShape_jpsi");

  // load snapshot with all results
  std::stringstream masssnapshotname;
  masssnapshotname << "m_snapshot_rap" << rapBin << "_pt" << ptBin;
  m_ws->loadSnapshot(masssnapshotname.str().c_str());

  // function to draw random mass values
  m_randDists.set("funcBG", (TF1*)bkgMass->asTF(*m, RooArgList(*poly1, *poly2, *m)));
  m_randDists.set("funcSig1", (TF1*)signalMass1->asTF(*m, RooArgList(*CBmass1, *CBsigma1, *CBalpha1, *CBn1, *m)));
  m_randDists.set("funcSig2", (TF1*)signalMass2->asTF(*m, RooArgList(*CBmass2, *CBsigma2, *CBalpha2, *CBn2, *m)));
  m_randDists.set("funcBGJpsi", (TF1*)bkgMassJpsi->asTF(*mjpsi, RooArgList(*lambdaJpsi), *mjpsi));
  m_randDists.set("funcSigJpsi", (TF1*)sigMassJpsi->asTF(*mjpsi, RooArgList(*CBmassJpsi, *CBsigmaJpsi, *CBalphaJpsi, *CBnJpsi), *mjpsi));
}

template<StateT State>
void BkgHistoProducer<State>::setupRandomDists(const int rapBin, const int ptBin)
{
  // TODO
}

// ================================================================================
//                              STORE 1D FACTORS
// ================================================================================
template<>
void BkgHistoProducer<Chic>::store1DFactors(bool MC)
{
  std::cout << "---------- IN BkgHistoProducer<Chic>::store1DFactors()" << std::endl;

  // for looping over chic1 and chic2 (less duplicate code). Replace with std::array once available
  const std::string cr[] = {"1", "2"};
  for (size_t i = 0; i < m_outFiles.size(); ++i) {
    storeFactor(m_outFiles[i], "comb_background_fraction", ";;fraction of comb. BG in PRSR" + cr[i],
                m_fitVars["fBGsig" + cr[i]], m_fitVars["fBGerr" + cr[i]]); // combinatorial background fraction
    storeFactor(m_outFiles[i], "nonprompt_background_fraction", ";;fraction of non prompt BG in PRSR" + cr[i],
                m_fitVars["fNPB" + cr[i]], m_fitVars["fNPerr" + cr[i]]); // non prompt background fraction
    storeFactor(m_outFiles[i], "background_fraction", ";;fraction of total BG in PRSR" + cr[i],
                m_fitVars["fTBGsig" + cr[i]], m_fitVars["fTBGerr" + cr[i]]); // total background fraction
    storeFactor(m_outFiles[i], "prompt_fraction", ";;fraction of prompt events in PRSR" + cr[i] + " (1-fNP-fBkg)",
                m_fitVars["fP" + cr[i]], m_fitVars["fPerr" + cr[i]]); // prompt fraction
    storeFactor(m_outFiles[i], "fraction_LSB", ";;f_{LSB}", m_fitVars["fracLSB" + cr[i]], 0); // no error here
    if (!MC) {
      // get some variables from the workspace
      // non prompt chic1 fraction in NPSR1
      double fNPChicInNPSR = getVarVal(m_ws, "var_fracNPChic"+cr[i]+"InNPSR"+cr[i]);
       // prompt chic1 fraction in NPSR1
      double fPRChicInNPSR = getVarVal(m_ws, "var_fracPRChic"+cr[i]+"InNPSR"+cr[i]);
      double nTot = getVarVal(m_ws, "var_nBackground") + getVarVal(m_ws, "var_nChic"); // total event number
      double nSR = nTot * getVarVal(m_ws, "var_fTotInSR" + cr[i]); // total events in signal region
      double nNPSR = nSR * getVarVal(m_ws, "var_fNPSR" + cr[i] + "InSR" + cr[i]); // total number of events in NPSR
      double nPRSR = nSR * getVarVal(m_ws, "var_fPRSR" + cr[i] + "InSR" + cr[i]); // total number of events in PRSR
      double fBGinSR = getVarVal(m_ws, "var_fracBackgroundInSR" + cr[i]); // background fraction in SR
      // non prompt chicX fraction in PRSR1
      double fNPChicInPRSR = getVarVal(m_ws, "var_fracNPChic"+cr[i]+"InPRSR"+cr[i]);
      // number of prompt chicX events in PRSR1
      double nPRChicInPRSR = getVarVal(m_ws, "var_nPRChic"+cr[i]+"InPRSR"+cr[i]);
      double fBGinNP = m_fitVars["fBGinNP" + cr[i]]; // comb. background in NPSR1

      store3Factors(m_outFiles[i], "events_SR", ";;chic" + cr[i]+ " events in SR" + cr[i],
                    nPRChicInPRSR + fPRChicInNPSR * nNPSR,
                    fNPChicInPRSR * nPRSR + fNPChicInNPSR * nNPSR,
                    fBGinSR * nSR); // events in signal region

      store3Factors(m_outFiles[i], "events_promptSR", ";;chic" + cr[i] + " events in PRSR" + cr[i],
                    nPRChicInPRSR, fNPChicInPRSR * nPRSR, m_fitVars["fBGsig" + cr[i]] * nPRSR);

      store3Factors(m_outFiles[i], "events_nonpromptSR", ";;chic" + cr[i] + " events in NPSR" + cr[i],
                    fPRChicInNPSR * nNPSR, fNPChicInNPSR * nNPSR, fBGinNP * nNPSR);

    }
  }

  std::cout << "********** IN BkgHistoProducer<Chic>::store1DFactors()" << std::endl;
}

template<StateT State>
void BkgHistoProducer<State>::store1DFactors(bool MC)
{
  // TODO
}

// ================================================================================
//                   FIND RAP PT BORDERS DI MUON KIN
// ================================================================================
template<>
void BkgHistoProducer<Chic>::findRapPtBordersDiMuonKin(const int rapBin, const int ptBin)
{
  // set these values eplicitly (although they should in principle also be in the fitVars store)
  double minPt = onia::pTRange[rapBin][ptBin-1];
  double maxPt = onia::pTRange[rapBin][ptBin];
  double minRap = onia::rapForPTRange[rapBin-1];
  double maxRap = onia::rapForPTRange[rapBin];

  // retrieve some variables from store to avoid lookups in the loop body
  double p0 = m_fitVars["p0"];
  double p1 = m_fitVars["p1"];
  double p2 = m_fitVars["p2"];
  double CBmass_p0 = m_fitVars["CBmass_p0"];

  TH1D h_pt("", "", 1000, 0, 100);
  TH1D h_y("", "", 240, 0, 2.4);

  int nEntries = m_inTree->GetEntries();
  for (int i = 0; i < nEntries; ++i) {
    long iEntry = m_inTree->LoadTree(i);
    if (i % 100000 == 0) { std:: cout << "entry " << i << " out of " << nEntries << std::endl; }

    m_inTree->GetEntry(iEntry);
    double pt = m_particle->Pt();
    double absY = TMath::Abs(m_particle->Rapidity());
    if (pt >= minPt && pt < maxPt && absY >= minRap && absY < maxRap && // NOTE: assumption standard def of min/max
        TMath::Abs(m_inputVars.jpsi->M() - CBmass_p0) < onia::nSigMass * rapSigma(p0, p1, p2, absY)) {

      h_pt.Fill(m_inputVars.jpsi->Pt());
      h_y.Fill(TMath::Abs(m_inputVars.jpsi->Rapidity()));
    }
  }

  double pt1 = h_pt.GetBinLowEdge(h_pt.FindFirstBinAbove());
  minPt = std::floor(pt1 * 2) / 2;
  double pt2 = h_pt.GetBinLowEdge(h_pt.FindLastBinAbove() + 1);
  maxPt = std::ceil(pt2 * 2) / 2;
  double y1 = h_y.GetBinLowEdge(h_y.FindFirstBinAbove());
  minRap = std::floor(y1 * 10) / 10;
  double y2 = h_y.GetBinLowEdge(h_y.FindLastBinAbove() + 1);
  maxRap = std::ceil(y2 * 10) / 10;
  if(maxRap < onia::rapForPTRange[rapBin]) maxRap = onia::rapForPTRange[rapBin];

  std::cout << "first bins = " << pt1 << ", " << y1 << "; last bins = " << pt2 << ", " << y2 << std::endl;

  m_fitVars.set("minPt", minPt); m_fitVars.set("maxPt", maxPt);
  m_fitVars.set("minRap", minRap); m_fitVars.set("maxRap", maxRap);
}

template<StateT State>
void BkgHistoProducer<State>::findRapPtBordersDiMuonKin(const int, const int) { /* No-OP */ }

// ================================================================================
//                         INITIALIZE SPEZIALIZATION
// ================================================================================
template<>
void BkgHistoProducer<Jpsi>::initialize(const std::string& infileName, const int rapBin, const int ptBin, bool MC,
                                        const int FracLSB, bool)
{
  std::cout << "---------- INITIALIZE FOR JPSI" << std::endl;

  initCommon(infileName);
  std::stringstream outfilename;
  outfilename << "tmpFiles/data_Psi" << 1 << "S_rap" << rapBin << "_pT" << ptBin << ".root";
  addOutputFile(outfilename.str());
  storePtRapBorders(rapBin, ptBin);

  // get the fit variables
  setupFitVariables(rapBin, ptBin, MC, FracLSB);
  std::cout << "m_fitVars after setup: " << m_fitVars << std::endl;

  store1DFactors(MC);

  setupRandomDists(rapBin, ptBin);

  setupCosThPhiHists(onia::kNbBinsCosT, onia::kNbBinsPhiPol);
  m_ptHists.createHists("pT", 1, 100, m_fitVars["minPt"], m_fitVars["maxPt"]);
  m_rapHists.createHists("rap", 1, 100, m_fitVars["minRap"], m_fitVars["maxRap"]);
  setupPtRapMassHists(1, 7, m_fitVars["minPt"], m_fitVars["maxPt"], 2, m_fitVars["minRap"], m_fitVars["maxRap"],
                      7, m_fitVars["jpsiMassSigMin"], m_fitVars["jpsiMassSigMax"]);

  std::cout << "********** INITIALIZE FOR JPSI" << std::endl;
}

template<>
void BkgHistoProducer<Psi2S>::initialize(const std::string& infileName, const int rapBin, const int ptBin, bool MC,
                                         const int FracLSB, bool)
{
  std::cout << "---------- INITIALIZE FOR PSI2S" << std::endl;

  initCommon(infileName);

  std::stringstream outfilename;
  outfilename << "tmpFiles/data_Psi" << 2 << "S_rap" << rapBin << "_pT" << ptBin << ".root";
  addOutputFile(outfilename.str());
  storePtRapBorders(rapBin, ptBin);

  // get the fit variables
  setupFitVariables(rapBin, ptBin, MC, FracLSB);
  std::cout << "m_fitVars after setup: " << m_fitVars << std::endl;

  store1DFactors(MC);

  setupRandomDists(rapBin, ptBin);

  setupCosThPhiHists(onia::kNbBinsCosT, onia::kNbBinsPhiPol);
  m_ptHists.createHists("pT", 1, 100, onia::pTRange[rapBin][ptBin-1], onia::pTRange[rapBin][ptBin]);
  m_rapHists.createHists("rap", 1, 100, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);
  setupPtRapMassHists(1, 7, m_fitVars["minPt"], m_fitVars["maxPt"], 2, m_fitVars["minRap"], m_fitVars["maxRap"],
                      7, m_fitVars["jpsiMassSigMin"], m_fitVars["jpsiMassSigMax"]);

  std::cout << "********** INITIALIZE FOR PSI2S" << std::endl;
}

template<>
void BkgHistoProducer<Chic>::initialize(const std::string& infileName, const int rapBin, const int ptBin, bool MC,
                                        const int FracLSB, bool useRefittedChic)
{
  std::cout << "---------- INITIALIZE FOR CHIC" << std::endl;

  initCommon(infileName);
  // set Branch Addresses that are only present for chic
  m_inTree->SetBranchAddress("chic", &m_inputVars.chic);
  m_inTree->SetBranchAddress("chic_rf", &m_inputVars.chic_rf);
  m_inTree->SetBranchAddress("mQ", &m_inputVars.mQ);

  if (onia::KinParticleChi) {
    if (useRefittedChic) {
      m_inTree->SetBranchAddress("chic_rf", &m_particle);
      m_inputVars.chic_rf = m_particle;
    } else {
      m_inTree->SetBranchAddress("chic", &m_particle);
      m_inputVars.chic = m_particle;
    }
  } else {
    m_inTree->SetBranchAddress("jpsi", &m_particle);
    m_inputVars.jpsi = m_particle;
    std::cout << "Jpsi kinematics used!" << std::endl;
  }

  // create the outpu files
  std::stringstream chic1file, chic2file;
  chic1file << "tmpFiles/data_chic1_rap" << rapBin << "_pT" << ptBin << ".root";
  chic2file << "tmpFiles/data_chic2_rap" << rapBin << "_pT" << ptBin << ".root";
  addOutputFile(chic1file.str());
  addOutputFile(chic2file.str());
  storePtRapBorders(rapBin, ptBin);

  // get the fit variables
  setupFitVariables(rapBin, ptBin, MC, FracLSB);
  std::cout << "m_fitVars after setup: " << m_fitVars << std::endl;

  store1DFactors(MC);

  setupRandomDists(rapBin, ptBin);
  std::cout << m_randDists << std::endl;

  if (onia::KinParticleChi) {
    findRapPtBordersDiMuonKin(rapBin, ptBin);
  }

  std::cout << "--------------------------------------------\n"
            << "pT binning of jpsi: " << m_fitVars["minPt"] << " - " << m_fitVars["maxPt"] << "\n"
            << "rapidity binning of jpsi: " << m_fitVars["minRap"] << " - " << m_fitVars["maxRap"] << std::endl;

  setupCosThPhiHists(onia::kNbBinsCosT, onia::kNbBinsPhiPol);
  m_ptHists.createHists("pT", 2, 100, onia::pTRange[rapBin][ptBin-1], onia::pTRange[rapBin][ptBin]);
  m_rapHists.createHists("rap", 2, 100, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);
  setupPtRapMassHists(2, 7, m_fitVars["minPt"], m_fitVars["maxPt"], 2, m_fitVars["minRap"], m_fitVars["maxRap"],
                      7, m_fitVars["jpsiMassSigMin"], m_fitVars["jpsiMassSigMax"]);


  std::cout << "********** INITIALIZE FOR CHIC" << std::endl;
}

// This case should never happen, since it should get caught at the factory creating this object already.
template<StateT State>
void BkgHistoProducer<State>::initialize(const std::string&, const int, const int, bool, const int, bool)
{
  std::cerr << "Calling BkgHistoProducer::initialize() with StateT (" << State << ")" << " which is not defined." << std::endl;
  return;
}

// ================================================================================
//                          FILL 1D HISTS
// ================================================================================
template<>
void BkgHistoProducer<Chic>::fill1DHists(bool RSB, bool LSB, bool NP, bool MC, const double pt, const double absY,
                                         const double mass, const BkgHistoRangeReport& eventRegion)
{
  // signal histos and outTree Filling
  if (LSB) { // left sideband
    if (eventRegion.isLSB() && eventRegion.isPR()) {
      for (size_t i = 0; i < m_outFiles.size(); ++i) {
        m_outTrees[i]->Fill();
        m_ptHists.h_PSR[i]->Fill(pt);
        m_rapHists.h_PSR[i]->Fill(absY);
      }
    }
  } else if (RSB) { // right sideband
    if (eventRegion.isRSB() && eventRegion.isPR()) {
      for (size_t i = 0; i < m_outFiles.size(); ++i) {
        m_outTrees[i]->Fill();
        m_ptHists.h_PSR[i]->Fill(pt);
        m_rapHists.h_PSR[i]->Fill(absY);
      }
    }
  } else if (NP) { // non prompt data
    if (eventRegion.isSR1() && eventRegion.isNP()) {
      m_outTrees[0]->Fill();
      m_ptHists.h_PSR[0]->Fill(pt);
      m_rapHists.h_PSR[0]->Fill(absY);
    } else if (eventRegion.isSR2() && eventRegion.isNP()) {
      m_outTrees[1]->Fill();
      m_ptHists.h_PSR[1]->Fill(pt);
      m_rapHists.h_PSR[1]->Fill(absY);
    }
  } else if (MC) { // only prompt data in MC
    if (eventRegion.isSR1()) {
      m_outTrees[0]->Fill();
      m_ptHists.h_PSR[0]->Fill(pt);
      m_rapHists.h_PSR[0]->Fill(absY);
    } else if (eventRegion.isSR2()) {
      m_outTrees[1]->Fill();
      m_ptHists.h_PSR[1]->Fill(pt);
      m_rapHists.h_PSR[1]->Fill(absY);
    }
  } else { // prompt data in signal region
    if (eventRegion.isSR1() && eventRegion.isPR()) {
      m_outTrees[0]->Fill();
      m_ptHists.h_PSR[0]->Fill(pt);
      m_rapHists.h_PSR[0]->Fill(absY);
    } else if (eventRegion.isSR2() && eventRegion.isPR()) {
      m_outTrees[1]->Fill();
      m_ptHists.h_PSR[1]->Fill(pt);
      m_rapHists.h_PSR[1]->Fill(absY);
    }
  }

  // rest of the pT and rap histograms
  if (MC) { // fill rap and pT hists with all events
    m_ptHists.h_L->Fill(pt);
    m_ptHists.h_R->Fill(pt);
    m_ptHists.h_highct_L->Fill(pt);
    m_ptHists.h_highct_R->Fill(pt);
    m_rapHists.h_L->Fill(absY);
    m_rapHists.h_R->Fill(absY);
    m_rapHists.h_highct_L->Fill(absY);
    m_rapHists.h_highct_R->Fill(absY);
    if (eventRegion.isSR1()) {
      m_ptHists.h_NP[0]->Fill(pt);
      m_rapHists.h_NP[0]->Fill(absY);
    } else if (eventRegion.isSR2()) {
      m_ptHists.h_NP[1]->Fill(pt);
      m_rapHists.h_NP[1]->Fill(absY);
    }
  } else { // for data fill the histograms accordingly
    if (eventRegion.isLSB() && eventRegion.isPR()) {
      m_ptHists.h_L->Fill(pt);
      m_rapHists.h_L->Fill(absY);
    } else if (eventRegion.isLSB() && eventRegion.isNP()) {
      m_ptHists.h_highct_L->Fill(pt);
      m_rapHists.h_highct_L->Fill(absY);
    } else if (eventRegion.isSR1() && eventRegion.isNP()) {
      m_ptHists.h_NP[0]->Fill(pt);
      m_ptHists.h_NP[0]->Fill(absY);
    } else if (eventRegion.isSR2() && eventRegion.isNP()) {
      m_ptHists.h_NP[1]->Fill(pt);
      m_ptHists.h_NP[1]->Fill(absY);
    } else if (eventRegion.isRSB() && eventRegion.isPR()) {
      m_ptHists.h_R->Fill(pt);
      m_rapHists.h_R->Fill(absY);
    } else if (eventRegion.isRSB() && eventRegion.isNP()) {
      m_ptHists.h_highct_R->Fill(pt);
      m_rapHists.h_highct_R->Fill(absY);
    } else { // Since SR and Sideband regions are not handled together it is possible that this is reached!
      // std::cerr << "Not in any of the regions (Fill1DHists)" << std::endl; // should be caught earlier
    }
  }
}

// TODO for Jpsis
template<StateT State>
void BkgHistoProducer<State>::fill1DHists(bool, bool, bool, bool, const double, const double, const double,
                                          const BkgHistoRangeReport&)
{ /* No-OP */ }

// ================================================================================
//                      FILL 3D HISTS
// ================================================================================
template<>
void BkgHistoProducer<Chic>::fill3DHists(bool RSB, bool LSB, bool NP, bool MC, const double pt, const double absY,
                                         const double mass, const BkgHistoRangeReport& eventRegion)
{
  /*const*/ TF1* funcBGJpsi = m_randDists["funcBGJpsi"]; // avoid too many lookups
  const double jpsiMassSigMin = m_fitVars["jpsiMassSigMin"];
  const double jpsiMassSigMax = m_fitVars["jpsiMassSigMax"];

  // TODO: MC
  if (eventRegion.isLSB() && eventRegion.isPR()) { // prompt left sideband
    m_ptRapMass[0].hBG_L->Fill(pt, absY, funcBGJpsi->GetRandom(jpsiMassSigMin, jpsiMassSigMax));
    m_ptRapMass[1].hBG_L->Fill(pt, absY, funcBGJpsi->GetRandom(jpsiMassSigMin, jpsiMassSigMax));
  } else if (eventRegion.isLSB() && eventRegion.isNP()) { // non prompt LSB
    m_ptRapMass[0].h_highct_L->Fill(pt, absY, funcBGJpsi->GetRandom(jpsiMassSigMin, jpsiMassSigMax));
    m_ptRapMass[1].h_highct_L->Fill(pt, absY, funcBGJpsi->GetRandom(jpsiMassSigMin, jpsiMassSigMax));
  } else if (eventRegion.isSR1() && eventRegion.isPR()) { // prompt SR1
    m_ptRapMass[0].hSR->Fill(pt, absY, funcBGJpsi->GetRandom(jpsiMassSigMin, jpsiMassSigMax));
  } else if (eventRegion.isSR1() && eventRegion.isNP()) { // non prompt SR1
    m_ptRapMass[0].hNP->Fill(pt, absY, funcBGJpsi->GetRandom(jpsiMassSigMin, jpsiMassSigMax));
  } else if (eventRegion.isSR2() && eventRegion.isPR()) { // prompt SR2
    m_ptRapMass[1].hSR->Fill(pt, absY, funcBGJpsi->GetRandom(jpsiMassSigMin, jpsiMassSigMax));
  } else if (eventRegion.isSR2() && eventRegion.isNP()) { // non prompt SR2
    m_ptRapMass[1].hNP->Fill(pt, absY, funcBGJpsi->GetRandom(jpsiMassSigMin, jpsiMassSigMax));
  } else if (eventRegion.isRSB() && eventRegion.isPR()) { // prompt right sideband
    m_ptRapMass[0].hBG_R->Fill(pt, absY, funcBGJpsi->GetRandom(jpsiMassSigMin, jpsiMassSigMax));
    m_ptRapMass[1].hBG_R->Fill(pt, absY, funcBGJpsi->GetRandom(jpsiMassSigMin, jpsiMassSigMax));
  } else if (eventRegion.isRSB() && eventRegion.isNP()) { // non prompt RSB
    m_ptRapMass[0].h_highct_R->Fill(pt, absY, funcBGJpsi->GetRandom(jpsiMassSigMin, jpsiMassSigMax));
    m_ptRapMass[1].h_highct_R->Fill(pt, absY, funcBGJpsi->GetRandom(jpsiMassSigMin, jpsiMassSigMax));
  } else {
    std::cerr << "Not in any of the regions (Fill3DHists)" << std::endl; // should be caught earlier
  }


}

template<StateT State>
void BkgHistoProducer<State>::fill3DHists(bool RSB, bool LSB, bool NP, bool MC, const double pt,
                                          const double absY, const double mass,
                                          const BkgHistoRangeReport& eventRegion)
{
  // TODO
}

// ================================================================================
//                              FILL 2D HISTS
// ================================================================================
template<>
void BkgHistoProducer<Chic>::fill2DHists(const std::vector<std::vector<double> >& cosThPVals,
                                         const BkgHistoRangeReport& eventRegion)
{
  for (int iFrame = 0; iFrame < onia::kNbFrames; ++iFrame) {
    if (eventRegion.isLSB() && eventRegion.isPR()) { // prompt left sideband
      m_cosThetaPhi[0].hBG_L[iFrame]->Fill(cosThPVals[iFrame][0], cosThPVals[iFrame][1]);
      m_cosThetaPhi[1].hBG_L[iFrame]->Fill(cosThPVals[iFrame][0], cosThPVals[iFrame][1]);
    } else if (eventRegion.isLSB() && eventRegion.isNP()) { // non prompt left sideband
      m_cosThetaPhi[0].hBGinNP_L[iFrame]->Fill(cosThPVals[iFrame][0], cosThPVals[iFrame][1]);
      m_cosThetaPhi[1].hBGinNP_L[iFrame]->Fill(cosThPVals[iFrame][0], cosThPVals[iFrame][1]);
    } else if (eventRegion.isSR1() && eventRegion.isPR()) { // prompt SR1
      m_cosThetaPhi[0].hSR[iFrame]->Fill(cosThPVals[iFrame][0], cosThPVals[iFrame][1]);
    } else if (eventRegion.isSR1() && eventRegion.isNP()) { // non prompt SR1
      m_cosThetaPhi[0].hNPBG[iFrame]->Fill(cosThPVals[iFrame][0], cosThPVals[iFrame][1]);
    } else if (eventRegion.isSR2() && eventRegion.isPR()) { // prompt SR2
      m_cosThetaPhi[1].hSR[iFrame]->Fill(cosThPVals[iFrame][0], cosThPVals[iFrame][1]);
    } else if (eventRegion.isSR2() && eventRegion.isNP()) { // non prompt SR2
      m_cosThetaPhi[1].hNPBG[iFrame]->Fill(cosThPVals[iFrame][0], cosThPVals[iFrame][1]);
    } else if (eventRegion.isRSB() && eventRegion.isPR()) { // prompt RSB
      m_cosThetaPhi[0].hBG_R[iFrame]->Fill(cosThPVals[iFrame][0], cosThPVals[iFrame][1]);
      m_cosThetaPhi[1].hBG_R[iFrame]->Fill(cosThPVals[iFrame][0], cosThPVals[iFrame][1]);
    } else if (eventRegion.isRSB() && eventRegion.isNP()) { // non prompt RSB
      m_cosThetaPhi[0].hBGinNP_R[iFrame]->Fill(cosThPVals[iFrame][0], cosThPVals[iFrame][1]);
      m_cosThetaPhi[1].hBGinNP_R[iFrame]->Fill(cosThPVals[iFrame][0], cosThPVals[iFrame][1]);
    } else {
      std::cerr << "Not in any region (fill2DHists)" << std::endl;
    }
  }
}

template<StateT State>
void BkgHistoProducer<State>::fill2DHists(const std::vector<std::vector<double> >& cosThPVals,
                                          const BkgHistoRangeReport& eventRegion)
{
  // TODO
}

// ================================================================================
//                               FILLHISTOS SPEZIALIZATION
// ================================================================================
template<>
void BkgHistoProducer<Jpsi>::fillHistos(const int rapBin, const int ptBin, bool useRefittedChic, bool MC,
                                        bool PolLSB, bool PolRSB, bool PolNP, bool folding)
{
  std::cout << "---------- FILLHISTOS FOR JPSI" << std::endl;
  // TODO
  std::cout << "********** FILLHISTOS FOR JPSI" << std::endl;
}

template<>
void BkgHistoProducer<Psi2S>::fillHistos(const int rapBin, const int ptBin, bool useRefittedChic, bool MC,
                                         bool PolLSB, bool PolRSB, bool PolNP, bool folding)
{
  std::cout << "---------- FILLHISTOS FOR PSI2S" << std::endl;
  // TODO
  std::cout << "********** FILLHISTOS FOR PSI2S" << std::endl;
}

template<>
void BkgHistoProducer<Chic>::fillHistos(const int rapBin, const int ptBin, bool useRefittedChic, bool MC,
                                        bool PolLSB, bool PolRSB, bool PolNP, bool folding)
{
  std::cout << "---------- FILLHISTOS FOR CHIC" << std::endl;
  typedef std::vector<std::vector<double> > cosTPT; // typedef for less typing

  const int nEntries = m_inTree->GetEntries(); //  100;
  unsigned int nAll = 0;
  unsigned int nSR = 0;
  unsigned int noValidRegion = 0;
  unsigned int inSR1 = 0;
  unsigned int inSR2 = 0;

  // cache some fit variables to avoid lookups in the loop body
  double p0 = m_fitVars["p0"], p1 = m_fitVars["p1"], p2 = m_fitVars["p2"];
  double CBmass_p0 = m_fitVars["CBmass_p0"];
  double nSigMass = MC ? onia::nSigMass * 15 : onia::nSigMass; // accept everything in MC-closure
  const BkgHistoRange prRange(m_fitVars["PRmin"], m_fitVars["PRmax"]);
  const BkgHistoRange npRange(m_fitVars["NPmin"], m_fitVars["NPmax"]);
  const BkgHistoRange lsbRange(m_fitVars["massMinL"], m_fitVars["massMaxL"]);
  const BkgHistoRange rsbRange(m_fitVars["massMinR"], m_fitVars["massMaxR"]);
  const BkgHistoRange sr1Range(m_fitVars["massMinSR1"], m_fitVars["massMaxSR1"]);
  const BkgHistoRange sr2Range(m_fitVars["massMinSR2"], m_fitVars["massMaxSR2"]);

  for (int i = 0; i < nEntries; ++i) {
    long iEntry = m_inTree->LoadTree(i);
    m_inTree->GetEntry(iEntry);
    if (i % 100000 == 0) { std::cout << "entry " << i << " out of " << nEntries << std::endl; }

    double usedChicMass = useRefittedChic ? m_particle->M() : m_inputVars.mQ;
    double pt = m_particle->Pt();
    double absY = TMath::Abs(m_particle->Rapidity());

    if (pt >= onia::pTRange[rapBin][ptBin-1] && pt < onia::pTRange[rapBin][ptBin] &&
        absY >= onia::rapForPTRange[rapBin-1] && absY < onia::rapForPTRange[rapBin]) {
      nAll++;

      if (TMath::Abs(m_inputVars.jpsi->M() - CBmass_p0) <
          nSigMass * rapSigma(p0, p1, p2, m_inputVars.jpsi->Rapidity())) {
        nSR++;

        // std::cout << "----------------------------------------" << std::endl;
        // classify the event only once and pass the results down to the filling functions
        // COULDDO: check already here if the event is in any of the desired regions
        BkgHistoRangeReport eventReg(sr1Range.accept(usedChicMass), sr2Range.accept(usedChicMass),
                                     npRange.accept(m_inputVars.jpsict), prRange.accept(m_inputVars.jpsict),
                                     lsbRange.accept(usedChicMass), rsbRange.accept(usedChicMass));

        if (!eventReg.isValidChicEvent()) {
          // std::cerr << "Event is not in any valid range!" << std::endl;
          noValidRegion++;
          continue;
        }
        if (eventReg.isSR1() && eventReg.isPR()) inSR1++;
        if (eventReg.isSR2() && eventReg.isPR()) inSR2++;

        // std::cout << "mass: " << usedChicMass << ", ct: " << m_inputVars.jpsict << std::endl;
        // std::cout << "SR1: " << sr1Range << ", SR2: " << sr2Range << ", LSB: " << lsbRange << ", RSB: "
        //           << rsbRange << ", PR: " << prRange << ", NP: " << npRange << std::endl;

        fill1DHists(PolRSB, PolLSB, PolNP, MC, pt, absY, usedChicMass, eventReg);
        fill3DHists(PolRSB, PolLSB, PolNP, MC, m_inputVars.jpsi->Pt(), absY, usedChicMass, eventReg);

        cosTPT cosThPhiVals = calcCosThetaPhiValues(*m_inputVars.lepP, *m_inputVars.lepN, folding);
        fill2DHists(cosThPhiVals, eventReg);
        // std::cout << "----------------------------------------" << std::endl;
      }
    }
  }
  std::cout << "------------------------------" << std::endl
            << "nAll = " << nAll << ", nSR = " << nSR << ", noValidRegion = " << noValidRegion << std::endl
            << "total events in PRSR: " << inSR1 << " chic1, " << inSR2 << " chic2" << std::endl;

  // rebinning
  std::pair<int, int> cosThPhiBins = determineBinning(*m_cosThetaPhi[0].hBG_R[2], *m_cosThetaPhi[0].hBG_L[2],
                                                      *m_cosThetaPhi[0].hNPBG[2], folding);

  std::cout << "final binning for background histogram:" << std::endl
            << "phi bins: " << cosThPhiBins.first << std::endl
            << "cosTheta bins: " << cosThPhiBins.second << std::endl
            << "------------------------------------------------------------" << std::endl;

  // re-run the filling for the cosTheta Phi histograms with the correct binning
  // COULDO: refactor the whole looping somehow, to remove this code duplicacy!
  // COULDO: check if this re-run is necessary, or if the binning is already correct.
  // COULDO: set the binning to the maximum number of bins and then use TH inbuilt Rebin() method to avoid
  // a second filling
  setupCosThPhiHists(cosThPhiBins.second, cosThPhiBins.first); // set them up

  for (int i = 0; i < nEntries; ++i) {
    long iEntry = m_inTree->LoadTree(i);
    m_inTree->GetEntry(iEntry);
    if (i % 100000 == 0) { std::cout << "entry " << i << " out of " << nEntries << std::endl; }

    double usedChicMass = useRefittedChic ? m_particle->M() : m_inputVars.mQ;
    double pt = m_particle->Pt();
    double absY = TMath::Abs(m_particle->Rapidity());

    if (pt >= onia::pTRange[rapBin][ptBin-1] && pt < onia::pTRange[rapBin][ptBin] &&
        absY >= onia::rapForPTRange[rapBin-1] && absY < onia::rapForPTRange[rapBin] &&
        TMath::Abs(m_inputVars.jpsi->M() - CBmass_p0) < nSigMass * rapSigma(p0,p1,p2,m_inputVars.jpsi->Eta())) {

      BkgHistoRangeReport eventReg(sr1Range.accept(usedChicMass), sr2Range.accept(usedChicMass),
                                   npRange.accept(m_inputVars.jpsict), prRange.accept(m_inputVars.jpsict),
                                   lsbRange.accept(usedChicMass), rsbRange.accept(usedChicMass));

      if (eventReg.isValidChicEvent()) {
        cosTPT cosThPhiVals = calcCosThetaPhiValues(*m_inputVars.lepP, *m_inputVars.lepN, folding);
        fill2DHists(cosThPhiVals, eventReg);
      }
    }
  }

  std::cout << "********** FILLHISTOS FOR CHIC" << std::endl;
}

// This case should never happen, since it should get caught at the factory creating this object already.
template<StateT State>
void BkgHistoProducer<State>::fillHistos(const int, const int, bool, bool, bool, bool, bool, bool)
{
  std::cerr << "Calling BkgHistoProducer::fillHistos() with StateT (" << State << ")" << " which is not defined." << std::endl;
  return;
}

// ================================================================================
//                              STORE 3D HISTS
// ================================================================================
template<>
void BkgHistoProducer<Chic>::store3DHists(bool PolLSB, bool PolRSB, bool PolNP, bool subtractNP)
{
  std::cout << "---------- IN BkgHistoProducer<Chic>::store3DHists" << std::endl;
  double fracLSB[] = { m_fitVars["fracLSB1"], m_fitVars["fracLSB2"] }; // std::array when available!
  double fBGinNP[] = { m_fitVars["fBGinNP1"], m_fitVars["fBGinNP2"] }; // std::array when available!
  double fNPB[] = { m_fitVars["fNPB1"], m_fitVars["fNPB2"] }; // std::array when available
  double fBGsig[] = { m_fitVars["fBGsig1"], m_fitVars["fBGsig2"] }; // std::array when available

  const std::string bgName = "background_pTrapMass";

  for (size_t iF = 0; iF < m_outFiles.size(); ++iF) {
    m_ptRapMass[iF].storeToFile(m_outFiles[iF]); // pt, |y|, M hists (store the starting points)

    TH3D* hBG_pTRapMass_lowct = addScaled(m_ptRapMass[iF].hBG_L, m_ptRapMass[iF].hBG_R, fracLSB[iF],
                                          "comb_background_pTrapMass", m_outFiles[iF]);
    TH3D* hBG_pTRapMass_highct = addScaled(m_ptRapMass[iF].h_highct_L, m_ptRapMass[iF].h_highct_R, fracLSB[iF],
                                           "comb_background_highct_pTrapMass", m_outFiles[iF]);
    TH3D* hNP_pTRapMass = subtractScaled(m_ptRapMass[iF].hNP, hBG_pTRapMass_highct, fBGinNP[iF],
                                         "NPS_highct_pTrapMass", m_outFiles[iF]);

    // TODO: there is probably a lot of potential for saving some lines and some instructions with some thinking
    TH3D* hBG_pTRapMass = new TH3D();
    if (PolLSB) { // LSB polarization: total background = signal contamination of prompt chic1;
      m_ptRapMass[iF].hSR->Scale(1. / (1. * m_ptRapMass[0].hSR->Integral()));
      hBG_pTRapMass = static_cast<TH3D*>(m_ptRapMass[0].hSR->Clone(bgName.c_str()));
    } else if (PolRSB) { // RSB polarization: total background = signal contamination of prompt chic2
      m_ptRapMass[iF].hSR->Scale(1. / (1. * m_ptRapMass[1].hSR->Integral()));
      hBG_pTRapMass = static_cast<TH3D*>(m_ptRapMass[1].hSR->Clone(bgName.c_str()));
    } else if (PolNP) { // non prompt polarization: total background = high ct background
      hBG_pTRapMass = static_cast<TH3D*>(hBG_pTRapMass_highct->Clone(bgName.c_str()));
    } else if (subtractNP) { // add low ct and non prompt background
      selfScale(hNP_pTRapMass, fNPB[iF]);
      hBG_pTRapMass = static_cast<TH3D*>(hNP_pTRapMass->Clone(bgName.c_str()));
      selfScale(hBG_pTRapMass_lowct, fBGsig[iF]);
      hBG_pTRapMass->Add(hBG_pTRapMass_lowct);
    } else {
      std::cout << "total background = combinatorial background" << std::endl;
      hBG_pTRapMass = static_cast<TH3D*>(hBG_pTRapMass_lowct->Clone(bgName.c_str()));
    }

    m_outFiles[iF]->cd();
    hBG_pTRapMass->Write();
  }

  std::cout << "********** IN BkgHistoProducer<Chic>::store3DHists" << std::endl;
}

template<StateT State>
void BkgHistoProducer<State>::store3DHists(bool, bool, bool, bool)
{
  // TODO
}

// ================================================================================
//                               STORE 1D HISTS
// ================================================================================
template<>
void BkgHistoProducer<Chic>::store1DHists(bool PolLSB, bool PolRSB, bool PolNP, bool subtractNP)
{
  // COULDDO: use the values directly
  const double fSRinPLSB = m_fitVars["fSRinPLSB"];
  const double fSRinPRSB = m_fitVars["fSRinPRSB"];
  const double fracLSB[] = { m_fitVars["fracLSB1"], m_fitVars["fracLSB1"] }; // std::array when available
  const double fBGinNP[] = { m_fitVars["fBGinNP1"], m_fitVars["fBGinNP1"] }; // std::array when available
  const double fBGsig[] = { m_fitVars["fBGsig1"], m_fitVars["fBGsig2"] };
  const double fNPB[] = { m_fitVars["fNPB1"], m_fitVars["fNPB1"] };

  // TODO: subtractScaledMean() is a really confusing function that modifies both input hists. -> Refactor!

  double meanPT[] = { -1., -1. }; // std::array when available
  double meanY[] = { -1., -1. }; // std::array when available
  if (PolLSB) { // prompt chic1 contamination in LSB for polarization of LSB
    double tmp = subtractScaledMean(m_ptHists.h_PSR[0], m_ptHists.h_L, fSRinPLSB);
    meanPT[0] = tmp; meanPT[1] = tmp;
    tmp = subtractScaledMean(m_rapHists.h_PSR[0], m_rapHists.h_L, fSRinPLSB);
    meanY[0] = tmp; meanY[1] = tmp;
  } else if (PolRSB) { // prompt chic2 contamination in RSB for polarization of RSB
    double tmp = subtractScaledMean(m_ptHists.h_PSR[1], m_ptHists.h_R, fSRinPRSB);
    meanPT[0] = tmp; meanPT[1] = tmp;
    tmp = subtractScaledMean(m_rapHists.h_PSR[1], m_ptHists.h_R, fSRinPRSB);
    meanY[0] = tmp; meanY[1] = tmp;
  } else {
    // continuum background using different fLSB for the two signal regions
    TH1D* pT_L2 = addSideBands(m_ptHists.h_L, m_ptHists.h_R, fracLSB[0], fracLSB[1]);
    TH1D* rap_L2 = addSideBands(m_rapHists.h_L, m_rapHists.h_R, fracLSB[0], fracLSB[1]);
    // non prompt continuum background in NPSB
    TH1D* pT_highct_L2 = addSideBands(m_ptHists.h_highct_L, m_ptHists.h_highct_R, fracLSB[0], fracLSB[1]);
    TH1D* rap_highct_L2 = addSideBands(m_rapHists.h_highct_L, m_rapHists.h_highct_R, fracLSB[0], fracLSB[1]);

    // non prompt background (NOTE: two statements per line!)
    selfScale(m_ptHists.h_NP[0]);                       selfScale(m_ptHists.h_NP[1]); // normalize non-prompts
    selfScale(m_ptHists.h_highct_L, fBGinNP[0]);        selfScale(pT_highct_L2, fBGinNP[1]);
    m_ptHists.h_NP[0]->Add(m_ptHists.h_highct_L, -1.);  m_ptHists.h_NP[1]->Add(pT_highct_L2, -1.);

    selfScale(m_rapHists.h_NP[0]);                        selfScale(m_rapHists.h_NP[1]); // normalize non-prompts
    selfScale(m_rapHists.h_highct_L, fBGinNP[0]);         selfScale(rap_highct_L2, fBGinNP[1]);
    m_rapHists.h_NP[0]->Add(m_rapHists.h_highct_L, -1.);  m_rapHists.h_NP[1]->Add(rap_highct_L2, -1.);

    if (PolNP) {
      meanPT[0] = m_ptHists.h_NP[0]->GetMean(); meanPT[1] = m_ptHists.h_NP[1]->GetMean();
      meanY[0] = m_rapHists.h_NP[0]->GetMean(); meanY[1] = m_rapHists.h_NP[1]->GetMean();
    }

    // pt
    selfScale(m_ptHists.h_L, fBGsig[0]);    selfScale(pT_L2, fBGsig[1]);
    if (subtractNP) {
      selfScale(m_ptHists.h_NP[0], fNPB[0]);    selfScale(m_ptHists.h_NP[1], fNPB[1]);
      m_ptHists.h_L->Add(m_ptHists.h_NP[0]);   pT_L2->Add(m_ptHists.h_NP[1]);
    }
    m_ptHists.h_PSR[0]->Add(m_ptHists.h_L, -1);  m_ptHists.h_PSR[1]->Add(pT_L2, -1);
    meanPT[0] = m_ptHists.h_PSR[0]->GetMean();   meanPT[1] = m_ptHists.h_PSR[1]->GetMean();

    // rap
    selfScale(m_rapHists.h_L, fBGsig[0]);    selfScale(rap_L2, fBGsig[1]);
    if (subtractNP) {
      selfScale(m_rapHists.h_NP[0], fNPB[0]);    selfScale(m_rapHists.h_NP[1], fNPB[1]);
      m_rapHists.h_L->Add(m_rapHists.h_NP[0]);   rap_L2->Add(m_rapHists.h_NP[1]);
    }
    m_rapHists.h_PSR[0]->Add(m_rapHists.h_L, -1);  m_rapHists.h_PSR[1]->Add(rap_L2, -1);
    meanY[0] = m_rapHists.h_PSR[0]->GetMean();   meanY[1] = m_rapHists.h_PSR[1]->GetMean();
  }

  for (size_t iF = 0; iF < m_outFiles.size(); ++iF) {
    storeFactor(m_outFiles[iF], "mean_pT", ";;mean p_{T}", meanPT[iF], 0);
    storeFactor(m_outFiles[iF], "mean_y", ";;mean |y|", meanY[iF], 0);
  }
}

template<StateT State>
void BkgHistoProducer<State>::store1DHists(bool PolLSB, bool PolRSB, bool PolNP, bool subtractNP)
{
  // TODO
}

// ================================================================================
//                               STORE 2D HISTS
// ================================================================================
template<>
void BkgHistoProducer<Chic>::store2DHists(bool PolLSB, bool PolRSB, bool PolNP, bool subtractNP, bool folding)
{
  std::cout << "---------- IN BkgHistoProducer<Chic>::store2DHists()" << std::endl;

  const std::vector<std::string> labels = charArrayToStrVec(onia::frameLabel, onia::kNbFrames);
  const double fracLSB[] = { m_fitVars["fracLSB1"], m_fitVars["fracLSB2"] };
  const double fBGinNP[] = { m_fitVars["fBGinNP1"], m_fitVars["fBGinNP2"] };
  const double fBGsig[] = { m_fitVars["fBGsig1"], m_fitVars["fBGsig2"] };
  const double fNPB[] = { m_fitVars["fNPB1"], m_fitVars["fNPB2"] };


  const std::string tbgFolded = "background_folded_costhphi";
  const std::string tbgUnfolded = "background_unfolded_costhphi";

  for (size_t iFr = 0; iFr < labels.size(); ++iFr) {
    for (size_t i = 0; i < 2; ++i) { // chic 1 and two (also size of m_outFiles, etc..)
      // combinatorial background in signal region (prompt region), combination of left and right sideband
      m_cosThetaPhi[i].hBG[iFr] = addScaled(m_cosThetaPhi[i].hBG_L[iFr], m_cosThetaPhi[i].hBG_R[iFr], fracLSB[i],
                                            "com_background_costhphi" + labels[iFr], m_outFiles[i]);

      // total background = combinatorial background
      m_cosThetaPhi[i].hTBG[iFr] = clone(m_cosThetaPhi[i].hBG[iFr], tbgFolded + labels[iFr]);
      if (!iFr && !i) { // output only once in the first run through
        std::cout << "background histogram output:" << std::endl
                  << "comb.: " << m_cosThetaPhi[i].hBG[iFr]->GetBinContent(3,12) << std::endl
                  << "tot.: " << m_cosThetaPhi[i].hTBG[iFr]->GetBinContent(3,12) << std::endl;
      }

      // combinatorial background in high ctau region, combination of left and right sideband in high ctau region
      m_cosThetaPhi[i].hBGinNP[iFr] = addScaled(m_cosThetaPhi[i].hBGinNP_L[iFr], m_cosThetaPhi[i].hBGinNP_R[iFr],
                                                fracLSB[i],"background_NPR_costhphi" + labels[iFr],m_outFiles[i]);
      if (PolNP) { // for non prompt polarization: only use high ctau background
        m_cosThetaPhi[i].hTBG[iFr] = clone(m_cosThetaPhi[0].hBGinNP[iFr], tbgFolded + labels[iFr]);
      }

      // non prompt background in high ctau region
      // (hNPBG)_norm - fBGinNP * (hBGinNP)_norm
      m_cosThetaPhi[i].hNPS[iFr] = subtractScaled(m_cosThetaPhi[i].hNPBG[iFr], m_cosThetaPhi[i].hBGinNP[iFr],
                                                  fBGinNP[i], "background_NPSR_costhphi" + labels[iFr],
                                                  m_outFiles[i]);

      // total background
      // polarization of LSB and RSB, bkg = signal contamination in left and right sideband
      if (PolLSB) { // NOTE: storing only SR1 (as in leptonBased.C)
        m_cosThetaPhi[i].hTBG[iFr] = clone(m_cosThetaPhi[0].hSR[iFr], tbgFolded + labels[iFr]);
      } else if (PolRSB) { // NOTE: storing only SR2 (NOT as in leptonBased.C!!!)
        m_cosThetaPhi[i].hTBG[iFr] = clone(m_cosThetaPhi[1].hSR[iFr], tbgFolded + labels[iFr]);
      } else if (subtractNP) { // fNPBG * (hNPS)_norm + fBGisg * (hBG)_norm
        m_cosThetaPhi[i].hTBG[iFr] = addWeighted(m_cosThetaPhi[i].hBG[iFr], m_cosThetaPhi[i].hNPS[iFr],
                                                 fBGsig[i], fNPB[i], tbgFolded + labels[iFr]); // don't store yet
      }
      if(!iFr && !i) { // once again only output once
        std::cout << "final tot.: " << m_cosThetaPhi[i].hTBG[iFr]->GetBinContent(3,12) << std::endl;
      }

      // write folded histograms to files
      m_outFiles[i]->cd();
      m_cosThetaPhi[i].hTBG[iFr]->Write();

      if (folding) {
        if (!iFr && !i) { // some output again
          int nx = m_cosThetaPhi[i].hTBG[iFr]->GetXaxis()->GetNbins();
          int ny = m_cosThetaPhi[i].hTBG[iFr]->GetYaxis()->GetNbins();

          std::cout << "------------------------------------------------------------" << std::endl
                    << "Total background histogram" << std::endl
                    << "number of cosTheta bins: " << nx << std::endl
                    << "number of phi bins: " << ny << std::endl
                    << "phi bins " << ny / 2 + 1 << " to " << 3 * ny / 4 << " are filled." << std::endl
                    << "------------------------------------------------------------" << std::endl;
        }
        unfold(m_cosThetaPhi[i].hTBG[iFr]);
      }

      // write unfolded and final histogram to file (same histogram twice becaus of normApproach relic)
      m_cosThetaPhi[i].hTBG[iFr]->SetName((tbgUnfolded + labels[iFr]).c_str());
      m_cosThetaPhi[i].hTBG[iFr]->Write();
    }
  }

  std::cout << "********** IN BkgHistoProducer<Chic>::store2DHists()" << std::endl;
}

template<StateT State>
void BkgHistoProducer<State>::store2DHists(bool PolLSB, bool PolRSB, bool PolNP, bool subtractNP, bool folding)
{
  // TODO
}


// ================================================================================
//                            STOREHISTOS SPEZIALIZATION
// ================================================================================
template<>
void BkgHistoProducer<Chic>::storeHistos(bool PolLSB, bool PolRSB, bool PolNP, bool subtractNP, bool folding)
{
  std::cout << "---------- IN BkgHistoProducer<Chic>::storeHistos()" << std::endl;

  std::cout << "---------------------------------\n"
            << "Build background histograms\n"
            << "pT-rap-mass" << std::endl;
  store3DHists(PolLSB, PolRSB, PolNP, subtractNP);

  std::cout << "pT and y" << std::endl;
  store1DHists(PolLSB, PolRSB, PolNP, subtractNP);

  std::cout << "cosTheta-phi" << std::endl;
  // store the currently filled cosThetaPhi histograms
  for (size_t iF = 0; iF < m_outFiles.size(); ++iF) {
    m_cosThetaPhi[iF].storeToFile(m_outFiles[iF]);
  }
  store2DHists(PolLSB, PolRSB, PolNP, subtractNP, folding);

  std::cout << "********** IN BkgHistoProducer<Chic>::storeHistos()" << std::endl;
}

template<StateT State>
void BkgHistoProducer<State>::storeHistos(bool PolLSB, bool PolRSB, bool PolNP, bool subtractNP, bool folding)
{
  std::cerr << "BkgHistoProducer::storeHistos() not yet implemented for State (" << State << ")" << std::endl;
  // temporary output for testing!!!
  for (size_t iF = 0; iF < m_outFiles.size(); ++iF) {
    m_cosThetaPhi[iF].storeToFile(m_outFiles[iF]);
  }
  // temporary saving for testing!!!
  m_ptHists.storeToFiles(m_outFiles);
  m_rapHists.storeToFiles(m_outFiles);

  std::cout << m_fitVars.printUsages() << std::endl;
  std::cout << m_randDists.printUsages() << std::endl;

  return;
}

// ================================================================================
//                          COMMON INITIALIZATION
// ================================================================================
template<StateT State>
void BkgHistoProducer<State>::initCommon(const std::string& infileName)
{
  m_fitFile = TFile::Open(infileName.c_str());
  if (!m_fitFile) {
    std::cerr << "Fitfile: \'" << infileName << "\' cannot be opened in BkgHistoProducer()" << std::endl;
  } else {
    std::cout << "Opened file: \'" << infileName << "\'." << std::endl;
  }
  const std::string wsname = "ws_masslifetime";
  m_ws = static_cast<RooWorkspace*>(m_fitFile->Get(wsname.c_str()));
  if (!m_ws) {
    std::cerr << "workspace \'" << wsname << "\' not found in Fitfile: \'" << infileName << "\'" << std::endl;
  }

  // set the branch adresses for the branches that are common in all rootfiles
  m_inTree->SetBranchAddress("lepP", &m_inputVars.lepP);
  m_inTree->SetBranchAddress("lepN", &m_inputVars.lepN);
  m_inTree->SetBranchAddress("jpsi", &m_inputVars.jpsi);
  m_inTree->SetBranchAddress("Jpsict", &m_inputVars.jpsict);

}

// ================================================================================
//                        STORE PT AND RAP BORDERS
// ================================================================================
template<StateT State>
void BkgHistoProducer<State>::storePtRapBorders(const int rapBin, const int ptBin)
{
  TVectorD* pTBorder = new TVectorD(1, 2, onia::pTRange[rapBin][ptBin - 1], onia::pTRange[rapBin][ptBin], "END");
  TVectorD* yBorder = new TVectorD(1, 2, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin], "END");

  for (TFileIt fIt = m_outFiles.begin(); fIt != m_outFiles.end(); ++fIt) {
    TFile* file = *fIt;
    if (file) {
      file->cd();
      pTBorder->Write();
      yBorder->Write();
    }
  }
}

// ================================================================================
//                           ADD OUTPUT FILE
// ================================================================================
template<StateT State>
inline void BkgHistoProducer<State>::addOutputFile(const std::string& filename)
{
  TFile* file = TFile::Open(filename.c_str(), "RECREATE");
  TTree* tree = m_inTree->CloneTree(0);

  m_outFiles.push_back(file);
  m_outTrees.push_back(tree);
}

// ================================================================================
//                          SETUP PT RAP MASS HISTS
// ================================================================================
template<StateT State>
void BkgHistoProducer<State>::setupPtRapMassHists(const size_t nHists,
                                                  const int ptBins, const double minPt, const double maxPt,
                                                  const int rapBins, const double minRap, const double maxRap,
                                                  const int massBins, const double minMass, const double maxMass)
{
  std::cout << "---------- IN BkgHistoProducer<" << State << ">::setupPtRapMassHists" << std::endl;
  for (size_t iH = 0; iH < nHists; ++iH) {
    m_ptRapMass.push_back(BkgHistoPtRapMassHists());
    m_ptRapMass[iH].hBG_L = new TH3D("hBG_pTRapMass_L", ";p_{T} [Gev/c]; |y|; M[GeV]",
                                     ptBins, minPt, maxPt, rapBins, minRap, maxRap, massBins, minMass, maxMass);
    m_ptRapMass[iH].hBG_L->Sumw2();

    m_ptRapMass[iH].hBG_R = new TH3D("hBG_pTRapMass_R", ";p_{T} [Gev/c]; |y|; M[GeV]",
                                     ptBins, minPt, maxPt, rapBins, minRap, maxRap, massBins, minMass, maxMass);
    m_ptRapMass[iH].hBG_R->Sumw2();

    m_ptRapMass[iH].h_highct_L = new TH3D("hBG_pTRapMass_highct_L", ";p_{T} [Gev/c]; |y|; M[GeV]",
                                          ptBins, minPt, maxPt, rapBins, minRap, maxRap, massBins, minMass, maxMass);
    m_ptRapMass[iH].h_highct_L->Sumw2();

    m_ptRapMass[iH].h_highct_R = new TH3D("hBG_pTRapMass_highct_R", ";p_{T} [Gev/c]; |y|; M[GeV]",
                                          ptBins, minPt, maxPt, rapBins, minRap, maxRap, massBins, minMass, maxMass);
    m_ptRapMass[iH].h_highct_R->Sumw2();

    m_ptRapMass[iH].hNP = new TH3D("hNP_pTRapMass_NP", ";p_{T} [Gev/c]; |y|; M[GeV]",
                                   ptBins, minPt, maxPt, rapBins, minRap, maxRap, massBins, minMass, maxMass);
    m_ptRapMass[iH].hNP->Sumw2();

    m_ptRapMass[iH].hNP = new TH3D("hNP_pTRapMass_NP", ";p_{T} [Gev/c]; |y|; M[GeV]",
                                   ptBins, minPt, maxPt, rapBins, minRap, maxRap, massBins, minMass, maxMass);
    m_ptRapMass[iH].hNP->Sumw2();

    m_ptRapMass[iH].hSR = new TH3D("hSR_pTRapMass", ";p_{T} [Gev/c]; |y|; M[GeV]",
                                   ptBins, minPt, maxPt, rapBins, minRap, maxRap, massBins, minMass, maxMass);
    m_ptRapMass[iH].hSR->Sumw2();
  }

  std::cout << "********** IN BkgHistoProducer<" << State << ">::setupPtRapMassHists" << std::endl;
}

#endif
