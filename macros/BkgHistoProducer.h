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

/** Interface class for the templatized version of the actual class, to have a common interface. */
class IBkgHistoProducer {
public:
  /** initialize everything that is necessary. Call for every pt and rapidity bin separately.*/
  virtual void initialize(const std::string& infileName, const int rapBin, const int ptBin, bool MC,
                          const int FracLSB) = 0;
  /** do the actual work and fill the histograms. Call for every pt and rapidity bin separately.*/
  virtual void fillHistos() = 0;
  /** store the appropriate histograms to the output files. Call for every pt and rapidity bin separately. */
  virtual void storeHistos() = 0;
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
                          const int FracLSB) /*override*/;
  /** fill all needed histograms and write them to file. Call for every pt and rapitdity bin separately. */
  virtual void fillHistos() /*override*/;
  /** store all the filled histograms, and combine them together to have the correct output histograms. */
  virtual void storeHistos() /*override*/;
  /** Destructor, do all the state independent clean-up work that is needed. */
  ~BkgHistoProducer();

private:
  TFile* m_dataFile; /**< pointer to the root file containing all data. */
  TTree* m_inTree; /**< pointer to the TTree in the dataFile containing all data. */
  TFile* m_fitFile; /**< pointer to the root file containing the fit data. */
  RooWorkspace* m_ws; /**< RooWorkspace containing the fit results. */
  BkgHistoRootVars m_inputVars; /**< input root variables container. */
  std::vector<TFile*> m_outFiles; /**< vector containing the rootfiles created by this class. */
  std::vector<TTree*> m_outTrees; /**< vector containing the TTrees that are stored in the output files */
  std::vector<BkgHistoCosThetaHists> m_cosThetaPhi; /**< vector containing the TH2D histo container. */
  BkgHisto1DHists m_ptHists; /**< the pT TH1D histogram container. */
  BkgHisto1DHists m_rapHists; /**< the rap TH1D histogram container. */
  std::vector<BkgHistoPtRapMassHists> m_ptRapMass; /**< vector for the pTRapMass histos. */
  NamedVarStore<double> m_fitVars; /**< storage for the different double values, like the fit variables. */
  /** storage for the different distributions from which random numbers are drawn. */
  NamedVarStore<TF1*> m_randomDists;

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
  void setupRandomDists();

  typedef typename std::vector<TFile*>::iterator TFileIt; /**< private typedef for easier looping over files. */
  typedef typename std::vector<TTree*>::iterator TTreeIt; /**< private typedef for easier looping over TTrees. */
};

// ================================================================================
//                               IMPLEMENTATION
// ================================================================================
template<StateT State>
BkgHistoProducer<State>::BkgHistoProducer()
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
  m_fitVars.set("fTBsig1", m_fitVars["fBGsig1"]); // total background fraction
  m_fitVars.set("fTBsig2", m_fitVars["fBGsig2"]); // total background fraction
  m_fitVars.set("fP1", 1. - m_fitVars["fBGsig1"] - m_fitVars["fNPB1"]);
  m_fitVars.set("fP2", 1. - m_fitVars["fBGsig2"] - m_fitVars["fNPB2"]);

  m_fitVars.setFromWS(m_ws, "var_fracPRChic1InPRLSB", "fSRinPLSB"); // prompt chic1 contamination in LSB
  m_fitVars.setFromWS(m_ws, "var_fracPRChic2InPRRSB", "fSRinPRSB"); // prompt chic2 contamination in RSB
  m_fitVars.setFromWS(m_ws, "var_nBackground", "nBkg"); // total number of background events
  m_fitVars.setFromWS(m_ws, "var_nChic", "nChic"); // total number of chic events
  m_fitVars.setFromWS(m_ws, "var_nChic1", "nChic1"); // total number of chic1 events
  m_fitVars.setFromWS(m_ws, "var_nChic2", "nChic2"); // total number of chic2 events
  m_fitVars.setFromWS(m_ws, "var_nPRChic1InPRSR1","nPRChic1InPRSR1"); // number of chic1 events in PRSR1
  m_fitVars.setFromWS(m_ws, "var_nPRChic2InPRSR2","nPRChic2InPRSR2"); // number of chic2 events in PRSR2

  m_fitVars.setFromWS(m_ws, "var_fTotInSR1", "fTotInSR1"); // fraction of total events in SR1
  m_fitVars.setFromWS(m_ws, "var_fracBackgroundInSR1", "fBGinSR1"); // background fraction in SR1
  m_fitVars.set("nSR1", (m_fitVars["nBkg"] + m_fitVars["nChic"]) * m_fitVars["fTotInSR1"]); // total number of events in SR1

  m_fitVars.setFromWS(m_ws, "var_fTotInSR2", "fTotInSR2"); // fraction of total events in SR2
  m_fitVars.setFromWS(m_ws, "var_fracBackgroundInSR2", "fBGinSR2"); // background fraction in SR2
  m_fitVars.set("nSR2", (m_fitVars["nBkg"] + m_fitVars["nChic"]) * m_fitVars["fTotInSR2"]); // total number of events in SR2

  m_fitVars.setFromWS(m_ws, "var_fPRSR1InSR1", "fPRSR1InSR1"); // fraction of PRSR1 in SR1
  m_fitVars.set("nPRSR1", m_fitVars["nSR1"] * m_fitVars["fPRSR1InSR1"]); // total number of events in PRSR1
  m_fitVars.setFromWS(m_ws, "var_fNPSR1InSR1", "fNPSR1InPRSR1"); // fraction of PRSR1 in SR1
  m_fitVars.set("nNPSR1", m_fitVars["nSR1"] * m_fitVars["fNPSR1InPRSR1"]); // total number of events in NPSR1
  m_fitVars.setFromWS(m_ws, "var_fracNPChic1InPRSR1", "fNPChic1InPRSR1"); // non prompt chic1 fraction in PRSR1
  m_fitVars.setFromWS(m_ws, "var_fracPRChic1InNPSR1", "fNPChic1InPRSR1"); // prompt chic1 fraction in NPSR1
  m_fitVars.setFromWS(m_ws, "var_fracNPChic1InNPSR1", "fNPChic1InNPSR1"); // non prompt chic1 fraction in NPSR1

  m_fitVars.setFromWS(m_ws, "var_fPRSR2InSR2", "fPRSR2InSR2"); // fraction of PRSR2 in SR2
  m_fitVars.set("nPRSR2", m_fitVars["nSR2"] * m_fitVars["fPRSR2InSR2"]); // total number of events in PRSR2
  m_fitVars.setFromWS(m_ws, "var_fNPSR2InSR2", "fNPSR2InPRSR2"); // fraction of PRSR2 in SR2
  m_fitVars.set("nNPSR2", m_fitVars["nSR2"] * m_fitVars["fNPSR2InPRSR2"]); // total number of events in NPSR2
  m_fitVars.setFromWS(m_ws, "var_fracNPChic2InPRSR2", "fNPChic2InPRSR2"); // non prompt chic2 fraction in PRSR2
  m_fitVars.setFromWS(m_ws, "var_fracPRChic2InNPSR2", "fNPChic2InPRSR2"); // prompt chic2 fraction in NPSR2
  m_fitVars.setFromWS(m_ws, "var_fracNPChic2InNPSR2", "fNPChic2InNPSR2"); // non prompt chic2 fraction in NPSR2

  // no differentiation between mumu+gamma and Jpsi+gamma events as mass distribution has no influence
  // fraction of combinatorial mumu+gamma events to total background
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
  m_fitVars.set("jpsiMassSigMax", m_fitVars["CBmass_p0"] - onia::nSigMass * m_fitVars["jpsiMassSigma"]);
  m_fitVars.set("jpsiMassBkgMinR", m_fitVars["CBmass_p0"] - onia::nSigBkgHigh * m_fitVars["jpsiMassSigma"]);

  m_fitVars.set("minPt", onia::pTRange[rapBin][ptBin-1]);
  m_fitVars.set("maxPt", onia::pTRange[rapBin][ptBin]);
  m_fitVars.set("minRap", onia::rapForPTRange[rapBin-1]);
  m_fitVars.set("maxRap", onia::rapForPTRange[rapBin]);
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
void BkgHistoProducer<Chic>::setupRandomDists()
{
  // TODO
}

template<StateT State>
void BkgHistoProducer<State>::setupRandomDists()
{
  // TODO
}

// ================================================================================
//                               INITIALIZE SPEZIALIZATION
// ================================================================================
template<>
void BkgHistoProducer<Jpsi>::initialize(const std::string& infileName, const int rapBin, const int ptBin, bool MC,
                                        const int FracLSB)
{
  std::cout << "---------- INITIALIZE FOR JPSI" << std::endl;

  initCommon(infileName);
  std::stringstream outfilename;
  outfilename << "tmpFiles/data_Psi" << 1 << "S_rap" << rapBin << "_pT" << ptBin << ".root";
  addOutputFile(outfilename.str());

  // get the fit variables
  setupFitVariables(rapBin, ptBin, MC, FracLSB);
  std::cout << "m_fitVars after setup: " << m_fitVars << std::endl;

  setupRandomDists();

  setupCosThPhiHists(onia::kNbBinsCosT, onia::kNbBinsPhiPol);
  m_ptHists.createHists("pT", 1, 100, m_fitVars["minPt"], m_fitVars["maxPt"]);
  m_rapHists.createHists("rap", 1, 100, m_fitVars["minRap"], m_fitVars["maxRap"]);
  setupPtRapMassHists(1, 7, m_fitVars["minPt"], m_fitVars["maxPt"], 2, m_fitVars["minRap"], m_fitVars["maxRap"],
                      7, m_fitVars["jpsiMassSigMin"], m_fitVars["jpsiMassSigMax"]);

  std::cout << "********** INITIALIZE FOR JPSI" << std::endl;
}

template<>
void BkgHistoProducer<Psi2S>::initialize(const std::string& infileName, const int rapBin, const int ptBin, bool MC,
                                         const int FracLSB)
{
  std::cout << "---------- INITIALIZE FOR PSI2S" << std::endl;

  initCommon(infileName);

  std::stringstream outfilename;
  outfilename << "tmpFiles/data_Psi" << 2 << "S_rap" << rapBin << "_pT" << ptBin << ".root";
  addOutputFile(outfilename.str());

  // get the fit variables
  setupFitVariables(rapBin, ptBin, MC, FracLSB);
  std::cout << "m_fitVars after setup: " << m_fitVars << std::endl;

  setupRandomDists();
  
  setupCosThPhiHists(onia::kNbBinsCosT, onia::kNbBinsPhiPol);
  m_ptHists.createHists("pT", 1, 100, onia::pTRange[rapBin][ptBin-1], onia::pTRange[rapBin][ptBin]);
  m_rapHists.createHists("rap", 1, 100, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);
  setupPtRapMassHists(1, 7, m_fitVars["minPt"], m_fitVars["maxPt"], 2, m_fitVars["minRap"], m_fitVars["maxRap"],
                      7, m_fitVars["jpsiMassSigMin"], m_fitVars["jpsiMassSigMax"]);

  std::cout << "********** INITIALIZE FOR PSI2S" << std::endl;
}

template<>
void BkgHistoProducer<Chic>::initialize(const std::string& infileName, const int rapBin, const int ptBin, bool MC,
                                        const int FracLSB)
{
  std::cout << "---------- INITIALIZE FOR CHIC" << std::endl;

  initCommon(infileName);
  // set Branch Addresses that are only present for chic
  m_inTree->SetBranchAddress("chic", &m_inputVars.chic);
  m_inTree->SetBranchAddress("chic_rf", &m_inputVars.chic_rf);
  m_inTree->SetBranchAddress("mQ", &m_inputVars.mQ);

  // create the outpu files
  std::stringstream chic1file, chic2file;
  chic1file << "tmpFiles/data_chic1_rap" << rapBin << "_pT" << ptBin << ".root";
  chic2file << "tmpFiles/data_chic2_rap" << rapBin << "_pT" << ptBin << ".root";
  addOutputFile(chic1file.str());
  addOutputFile(chic2file.str());

  // get the fit variables
  setupFitVariables(rapBin, ptBin, MC, FracLSB);
  std::cout << "m_fitVars after setup: " << m_fitVars << std::endl;

  setupRandomDists();
  
  setupCosThPhiHists(onia::kNbBinsCosT, onia::kNbBinsPhiPol);
  m_ptHists.createHists("pT", 2, 100, onia::pTRange[rapBin][ptBin-1], onia::pTRange[rapBin][ptBin]);
  m_rapHists.createHists("rap", 2, 100, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);
  setupPtRapMassHists(2, 7, m_fitVars["minPt"], m_fitVars["maxPt"], 2, m_fitVars["minRap"], m_fitVars["maxRap"],
                      7, m_fitVars["jpsiMassSigMin"], m_fitVars["jpsiMassSigMax"]);


  std::cout << "********** INITIALIZE FOR CHIC" << std::endl;
}

// This case should never happen, since it should get caught at the factory creating this object already.
template<StateT State>
void BkgHistoProducer<State>::initialize(const std::string&, const int, const int, bool, const int)
{
  std::cerr << "Calling BkgHistoProducer::initialize() with StateT (" << State << ")" << " which is not defined." << std::endl;
  return;
}

// ================================================================================
//                               FILLHISTOS SPEZIALIZATION
// ================================================================================
template<>
void BkgHistoProducer<Jpsi>::fillHistos()
{
  std::cout << "---------- FILLHISTOS FOR JPSI" << std::endl;
  // TODO
  std::cout << "********** FILLHISTOS FOR JPSI" << std::endl;
}

template<>
void BkgHistoProducer<Psi2S>::fillHistos()
{
  std::cout << "---------- FILLHISTOS FOR PSI2S" << std::endl;
  // TODO
  std::cout << "********** FILLHISTOS FOR PSI2S" << std::endl;
}

template<>
void BkgHistoProducer<Chic>::fillHistos()
{
  std::cout << "---------- FILLHISTOS FOR CHIC" << std::endl;
  // TODO
  std::cout << "********** FILLHISTOS FOR CHIC" << std::endl;
}

// This case should never happen, since it should get caught at the factory creating this object already.
template<StateT State>
void BkgHistoProducer<State>::fillHistos()
{
  std::cerr << "Calling BkgHistoProducer::fillHistos() with StateT (" << State << ")" << " which is not defined." << std::endl;
  return;
}

// ================================================================================
//                            STOREHISTOS SPEZIALIZATION
// ================================================================================
template<StateT State>
void BkgHistoProducer<State>::storeHistos()
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
