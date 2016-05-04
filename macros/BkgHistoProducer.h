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
  virtual void initialize(const std::string& infileName, const int rapBin, const int ptBin) = 0;
  /** do the actual work and fill the histograms. Call for every pt and rapidity bin separately.*/
  virtual void fillHistos() = 0;
  /** store the appropriate histograms to the output files. */
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

/** BkgHistoProducer class that does all the work, including filehandling, etc...
 * Only initialize and fillHistos are publicly available and will be specialized for every state, to acommodate
 * for the differences.
 */
template<StateT State>
class BkgHistoProducer : public IBkgHistoProducer {
public:
  /** Constructor that does all the work that can be done for initialization independently of the State. */
  BkgHistoProducer();
  /** initialize the fit file and setup all variables needed internally for the given State. Call for every
   * pt and rapidity bin separately.
   */
  virtual void initialize(const std::string& infileName, const int rapBin, const int ptBin) /*override*/;
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

  /** Initialization that is common to every state, like setup of the workspace, etc... */
  void initCommon(const std::string& infilename);
  /** store the pT and y borders in every output file that is currently open.
   * NOTE: call this only after the outpufiles have been opened!
   */
  void storePtRapBorders(const int rapBin, const int ptBin);

  /** add a new outputile to the m_outFiles vector. */
  void addOutputFile(const std::string& filename);

  /** setup all TH2D cosThPhi histograms that are filled in the fillHistos function.
   * NOTE: when called for the first time, the histos get created as well
   * NOTE: leaves the histograms that are present in BkgHistoCosThetaHists but not used in fillHistos untouched
   * @param nBCT, number of bins for cosTheta
   * @param nBP, number of bins for phi
   */
  void setupCosThPhiHists(const int nBCT, const int NBP);

  typedef typename std::vector<TFile*>::iterator TFileIt; /**< private typedef for easier looping over files. */
  typedef typename std::vector<TTree*>::iterator TTreeIt; /**< private typedef for easier looping over TTrees. */
};

// ================================================================================
//                               IMPLEMENTATION
// ================================================================================
template<StateT State>
BkgHistoProducer<State>::BkgHistoProducer()
{
  std::cout << "CONSTRUCTOR, STATE == " << State << std::endl;

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
}

template<StateT State>
BkgHistoProducer<State>::~BkgHistoProducer()
{
  std::cout << "DESTRUCTOR, STATE == " << State << std::endl;
  // close the input files (the way I understand root is, that this also destroys these objects)
  if (m_dataFile) { m_dataFile->Close(); }
  if (m_fitFile) { m_fitFile->Close(); }

  // write and close all output files
  for (size_t iFile = 0; iFile < m_outFiles.size(); ++iFile) {
    if (m_outFiles[iFile]) {
      m_outFiles[iFile]->cd();
      m_outTrees[iFile]->Write();
      m_outFiles[iFile]->Write();
      m_outFiles[iFile]->Close(); // should alse destroy the object
    }
  }
}

// ================================================================================
//        SETUP COS THETA PHI HISTOGRAMS (declaration needed before initialize!)
// ================================================================================
// only chic needs spezialized version
template<>
void BkgHistoProducer<Chic>::setupCosThPhiHists(const int nBCT, const int nBP)
{
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
}

// all others can be dealt with in the same way
template<StateT State>
void BkgHistoProducer<State>::setupCosThPhiHists(const int nBCT, const int nBP)
{
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
}

// ================================================================================
//                               INITIALIZE SPEZIALIZATION
// ================================================================================
template<>
void BkgHistoProducer<Jpsi>::initialize(const std::string& infileName, const int rapBin, const int ptBin)
{
  std::cout << "INITIALIZE FOR JPSI" << std::endl;

  initCommon(infileName);
  std::stringstream outfilename;
  outfilename << "tmpFiles/data_Psi" << 1 << "S_rap" << rapBin << "_pT" << ptBin << ".root";
  addOutputFile(outfilename.str());

  setupCosThPhiHists(onia::kNbBinsCosT, onia::kNbBinsPhiPol);
  // TODO
}

template<>
void BkgHistoProducer<Psi2S>::initialize(const std::string& infileName, const int rapBin, const int ptBin)
{
  std::cout << "INITIALIZE FOR PSI2S" << std::endl;

  initCommon(infileName);

  std::stringstream outfilename;
  outfilename << "tmpFiles/data_Psi" << 2 << "S_rap" << rapBin << "_pT" << ptBin << ".root";
  addOutputFile(outfilename.str());

  setupCosThPhiHists(onia::kNbBinsCosT, onia::kNbBinsPhiPol);
  // TODO
}

template<>
void BkgHistoProducer<Chic>::initialize(const std::string& infileName, const int rapBin, const int ptBin)
{
  std::cout << "INITIALIZE FOR CHIC" << std::endl;

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

  setupCosThPhiHists(onia::kNbBinsCosT, onia::kNbBinsPhiPol);

  // TODO
}

// This case should never happen, since it should get caught at the factory creating this object already.
template<StateT State>
void BkgHistoProducer<State>::initialize(const std::string&, const int, const int)
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
  std::cout << "FILLHISTOS FOR JPSI" << std::endl;
  // TODO
}

template<>
void BkgHistoProducer<Psi2S>::fillHistos()
{
  std::cout << "FILLHISTOS FOR PSI2S" << std::endl;
  // TODO
}

template<>
void BkgHistoProducer<Chic>::fillHistos()
{
  std::cout << "FILLHISTOS FOR CHIC" << std::endl;
  // TODO
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

#endif
