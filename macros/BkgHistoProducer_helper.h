#ifndef BKGHISTOPRODUCER_HELPER_H__
#define BKGHISTOPRODUCER_HELPER_H__

#include "bkgHistos_helper.h"

#include "rootIncludes.inc"

#include <vector>
#include <iostream> // this can probably removed after develpment
#include <string>
#include <map> // unordered_map, if c++11 is supported
#include <sstream>
#include <stdexcept> // remove after development/debugging


/** struct containing the variables that are stored in the data root file. */
struct BkgHistoRootVars {
  /** default constructor. Initialize all pointers to nullptr (as soon as c++11 is used), and double values to -1.*/
  BkgHistoRootVars() : lepP(NULL), lepN(NULL), jpsi(NULL), chic(NULL), chic_rf(NULL), jpsict(-1.), mQ(-1) {;}
  /** destructor. Delete all 4 vectors. */
  ~BkgHistoRootVars();
  TLorentzVector* lepP; /**< 4 momentum of positive muon. */
  TLorentzVector* lepN; /**< 4 momentum of negative muon. */
  TLorentzVector* jpsi; /**< 4 momentum of Jpsi. */
  TLorentzVector* chic; /**< 4 momentum of chic. */
  TLorentzVector* chic_rf; /**< 4 momentum of refitted chic. */
  double jpsict; /**< jpsi ct variable. */
  double mQ; /**< mQ variable. */
};

/**
 * struct containing all the TH2Ds for the different cosTheta Phi histograms. Holds intermediate histos as well
 * as output histograms.
 */
struct BkgHistoCosThetaHists {
  /** constructor. Creates the necessary number of histos (one for each frame) and initializes them to nullptr.*/
  BkgHistoCosThetaHists(const int nHists);
  /** destructor. Deletes all present histograms.*/
  ~BkgHistoCosThetaHists();
  /** store hists to root file. (all that are not NULL!) */
  void storeToFile(TFile* file);

  std::vector<TH2D*> hBG_L; /**< LSB histogram. */
  std::vector<TH2D*> hBG_R; /**< RSB histogram. */
  std::vector<TH2D*> hBG; /**< BG histogram. */
  std::vector<TH2D*> hNPBG; /**< NonPrompt BG histogram. */
  std::vector<TH2D*> hNPS; /**< NonPrompt signal histogram. */
  std::vector<TH2D*> hBGinNP_L; /**< BG in NP histogram for LSB. */
  std::vector<TH2D*> hBGinNP_R; /**< BG in NP histogram for RSB. */
  std::vector<TH2D*> hBGinNP; /**< BG in NP histogram. */
  std::vector<TH2D*> hTBG; /**< (?) histogram. */
  std::vector<TH2D*> hSR_L; /**< signal region LSB (?). */
  std::vector<TH2D*> hSR_R; /**< signal region RSB (?).*/
  std::vector<TH2D*> hSR; /**< signal region (?). */
};

/** struct containing all TH1Ds necessary to store either the pT or the rapidity data.  */
struct BkgHisto1DHists {
  /** constructor. Initializes everything to nullptrs. Leaves h_NP and h_PSR empty vectors. */
  BkgHisto1DHists();
  /** destructor. Deletes all present histograms. TODO: FIX! */
  ~BkgHisto1DHists();
  /** Creates all histograms with nBins and the passed min and max values.
   * Taking a prefix that is used for storing in root files and also the number of needed NP and PSR1 histograms.
   */
  void createHists(const std::string& prefix, const size_t nHists,
                   const int nBins, const double min, const double max);
  /** store the histograms to the files. the histograms not in vectors will be stored to each file. The histograms
   * in vectors will be stored to only one file (index).
   * NOTE: no check is done whether the vectors actually have the same size!
   */
  void storeToFiles(const std::vector<TFile*>& files);

  TH1D* h_L; /**< LSB histogram. */
  TH1D* h_R; /**< RSB histogram. */
  TH1D* h_highct_L; /**< LSB histogram for highct. */
  TH1D* h_highct_R; /**< RSB histogram for highct. */
  std::vector<TH1D*> h_NP; /**< non-prompt histograms. */
  std::vector<TH1D*> h_PSR; /**< prompt signal regions histograms. */
};

/** struct containing all TH3Ds for intermediate and final storage. */
struct BkgHistoPtRapMassHists {
  /** constructor. Initializes everything to nullptrs. */
  BkgHistoPtRapMassHists();
  /** destructor. Deletes all pointers. TODO: FIX! */
  ~BkgHistoPtRapMassHists();
  /** store all histograms in the passed file. */
  void storeToFile(TFile* file);

  TH3D* hBG_L; /**< LSB background */
  TH3D* hBG_R; /**< RSB background. */
  TH3D* h_highct_L; /**< LSB background highct.*/
  TH3D* h_highct_R; /**< RSB background highct.*/
  TH3D* hNP; /**< NP histo. */
  TH3D* hSR; /**< SR histo. */
};

/**
 * flexible class to hold differently named variables of the same type  with easy retrieval and set functions.
 * NOTE: makes use of a map internally, so consider storing variables that are needed frequently in a local
 * variable.
 * TODO: upgrade to unordered_map with c++11
 */
template<typename T>
class NamedVarStore {
public:
  /**
   * retrieve a value form the RooWorkspace and store it under the same name or another (if passed).
   * NOTE: if a value already exists with the name it is replaced (silently).
   */
  void setFromWS(RooWorkspace* ws, const std::string& wsName, const std::string& name = "");

  /**
   * store a value under the given name.
   * NOTE: if a value already exists with the name it is replaced (silently).
   */
  void set(const std::string& name, const T value);

  /**
   * get the value to the appropriate name.
   * NOTE: made somewhat save by using at(). Throws an exception if no entry with that name exists.
   * TODO: remove try-catch block after development!
   */
  T operator[](const std::string& name) const;

  /** dump the whole store onto a stream (for debugging purposes). */
  template <typename U>
  friend std::ostream& operator<<(std::ostream&, const NamedVarStore<U>&);

  /** debugging. Print how many times each variable has been retrieved from the store. */
  const std::string printUsages() const;

private:
  typedef typename std::map<std::string, T> MapT; /**< private typedef for the internally used map. */
  MapT m_store; /**< actual storing entity. */
  /** development entity for counting the times of usage for each variable. */
  mutable std::map<std::string, size_t> m_counter;
};

/** helper struct to represent simple (double valued) ranges. */
class BkgHistoRange {
public:
  BkgHistoRange(const double low, const double high) : m_low(low), m_high(high) {;} /**< ctor. */
  /** val is in range? */
  bool accept(const double val) const { return (val > m_low && val < m_high); }
  /** ostream operator */
  friend std::ostream& operator<<(std::ostream& os, const BkgHistoRange& range);
private:
  double m_low; /**< lower bound of range. */
  double m_high; /**< upper bound of range. */
};

std::ostream& operator<<(std::ostream& os, const BkgHistoRange& range)
{
  os << "[" << range.m_low << "," << range.m_high << "]";
  return os;
}

/**
 * class to store the results of the different range checkes of an event. Can be queried to get in which
 * region a given event is.
 */
class BkgHistoRangeReport {
public:
  /** ctor for chic case. */
  BkgHistoRangeReport(bool SR1, bool SR2, bool NP, bool PR, bool LSB, bool RSB) :
    m_NP(NP), m_PR(PR), m_LSB(LSB), m_RSB(RSB)
  { m_SRs.push_back(SR1); m_SRs.push_back(SR2); }
  bool isNP() const { return m_NP; } /**< is in non-prompt region. */
  bool isPR() const { return m_PR; } /**< is in prompt region. */
  bool isLSB() const { return m_LSB; } /**< is in LSB. */
  bool isRSB() const { return m_RSB; } /**< is in RSB. */
  bool isSR1() const { return m_SRs[0]; } /**< is in SR1. */
  bool isSR2() const { return m_SRs[1]; } /**< is in SR2. */
  /** check if the event can be categorized into the predefined regions. */
  bool isValidChicEvent() const { return ((m_NP || m_PR) && (m_LSB || m_RSB || m_SRs[0] || m_SRs[1])); }
private:
  std::vector<bool> m_SRs;
  bool m_NP;
  bool m_PR;
  bool m_LSB;
  bool m_RSB;
};

// ================================================================================
//                 IMPLEMENTATION BKG HISTO ROOTVARS
// ================================================================================
BkgHistoRootVars::~BkgHistoRootVars()
{
  std::cout << "---------- DESTRUCTOR OF BkgHistoRootVars" << std::endl;
  delete lepP;
  delete lepN;
  delete jpsi;
  delete chic;
  delete chic_rf;
  std::cout << "********** DESTRUCTOR OF BkgHistoRootVars" << std::endl;
}

// ================================================================================
//              IMPLEMENTATION BKG HISTO COSTHETA PHI HISTS
// ================================================================================
BkgHistoCosThetaHists::BkgHistoCosThetaHists(const int nHists)
{
  std::cout << "---------- CONSTRUCTOR OF BkgHistoCosThetaHists" << std::endl;
  // TODO: replace NULL with nullptr after migragtion to c++11
  for (int iH = 0; iH < nHists; ++iH) {
    hBG_L.push_back(NULL);
    hBG_R.push_back(NULL);
    hBG.push_back(NULL);
    hNPBG.push_back(NULL);
    hNPS.push_back(NULL);
    hBGinNP_L.push_back(NULL);
    hBGinNP_R.push_back(NULL);
    hBGinNP.push_back(NULL);
    hTBG.push_back(NULL);
    hSR_L.push_back(NULL);
    hSR_R.push_back(NULL);
    hSR.push_back(NULL);
  }
  std::cout << "********** CONSTRUCTOR OF BkgHistoCosThetaHists" << std::endl;
}

BkgHistoCosThetaHists::~BkgHistoCosThetaHists()
{
  std::cout << "---------- DESTRUCTOR OF BkgHistoCosThetaHists" << std::endl;
  // TODO: currently get a seg fault here. This is caused by the call to TFile::Close() prior to this
  // destructor, which apparently means for ROOT to delete TH2D that are "attatched" to the TFile

  // for (size_t iH = 0; iH < hBG_L.size(); ++iH) { // every vector should contain the same number of TH2Ds
  //   delete hBG_L[iH];
  //   delete hBG_R[iH];
  //   delete hBG[iH];
  //   delete hNPBG[iH];
  //   delete hNPS[iH];
  //   delete hBGinNP_L[iH];
  //   delete hBGinNP_R[iH];
  //   delete hBGinNP[iH];
  //   delete hTBG[iH];
  //   delete hSR_L[iH];
  //   delete hSR_R[iH];
  //   delete hSR[iH];
  // }
  std::cout << "********** DESTRUCTOR OF BkgHistoCosThetaHists" << std::endl;
}

void BkgHistoCosThetaHists::storeToFile(TFile* file)
{
  std::cout << "---------- STORE TO FILE OF BkgHistoCosThetaHists" << std::endl;
  file->cd();
  for (size_t iH = 0; iH < hBG_L.size(); ++iH) { // every vector should contain the same number of TH2Ds
    if (hBG_L[iH]) hBG_L[iH]->Write();
    if (hBG_R[iH]) hBG_R[iH]->Write();
    if (hBG[iH]) hBG[iH]->Write();
    if (hNPBG[iH]) hNPBG[iH]->Write();
    if (hNPS[iH]) hNPS[iH]->Write();
    if (hBGinNP_L[iH]) hBGinNP_L[iH]->Write();
    if (hBGinNP_R[iH]) hBGinNP_R[iH]->Write();
    if (hBGinNP[iH]) hBGinNP[iH]->Write();
    if (hTBG[iH]) hTBG[iH]->Write();
    if (hSR_L[iH]) hSR_L[iH]->Write();
    if (hSR_R[iH]) hSR_R[iH]->Write();
    if (hSR[iH]) hSR[iH]->Write();
  }
  std::cout << "********** STORE TO FILE OF BkgHistoCosThetaHists" << std::endl;
}

// ================================================================================
//                    IMPLEMENTATION BKG HISTO 1D HISTS
// ================================================================================
BkgHisto1DHists::BkgHisto1DHists()
{
  std::cout <<"---------- CONSTRUCTOR OF BkgHisto1DHists" << std::endl;
  // TODO: replace NULL with nullptr with c++11 availability
  h_L = NULL;
  h_R = NULL;
  h_highct_L = NULL;
  h_highct_R = NULL;

  std::cout <<"********** CONSTRUCTOR OF BkgHisto1DHists" << std::endl;
}

BkgHisto1DHists::~BkgHisto1DHists()
{
  std::cout <<"---------- DESTRUCTOR OF BkgHisto1DHists" << std::endl;

  // TODO: currently get a seg fault here. This is caused by the call to TFile::Close() prior to this
  // destructor, which apparently means for ROOT to delete TH1D that are "attatched" to the TFile
  // The question remains, why this works for TH2Ds (e.g.)? -> If properly set it doesn't!
  // delete h_L;
  // delete h_R;
  // delete h_highct_L;
  // delete h_highct_R;
  // for (size_t i = 0; i < h_NP.size(); ++i) delete h_NP[i];
  // for (size_t i = 0; i < h_PSR.size(); ++i) delete h_PSR[i];

  std::cout <<"********** DESTRUCTOR OF BkgHisto1DHists" << std::endl;
}

void BkgHisto1DHists::createHists(const std::string& prefix, const size_t nHists,
                                  const int nBins, const double min, const double max)
{
  std::cout <<"---------- CREATE HISTS OF BkgHisto1DHists" << std::endl;

  std::string temp = prefix + "LSB";
  h_L = new TH1D(temp.c_str(), temp.c_str(), nBins, min, max);
  temp = prefix + "RSB";
  h_R = new TH1D(temp.c_str(), temp.c_str(), nBins, min, max);
  temp = prefix + "_highct_LSB";
  h_highct_L = new TH1D(temp.c_str(), temp.c_str(), nBins, min, max);
  temp = prefix + "_highct_RSB";
  h_highct_R = new TH1D(temp.c_str(), temp.c_str(), nBins, min, max);

  for (size_t i = 0; i < nHists; ++i) {
    temp = prefix + "NP";
    h_NP.push_back(new TH1D(temp.c_str(), temp.c_str(), nBins, min, max));
    temp = prefix + "PSR";
    h_PSR.push_back(new TH1D(temp.c_str(), temp.c_str(), nBins, min, max));
  }

  std::cout <<"********** CREATE HISTS OF BkgHisto1DHists" << std::endl;
}

void BkgHisto1DHists::storeToFiles(const std::vector<TFile*>& files)
{
  std::cout << "---------- STORE TO FILES OF BkgHisto1DHists" << std::endl;

  for (size_t iFile = 0; iFile < files.size(); ++iFile) {
    TFile* file = files[iFile];

    file->cd();
    h_L->Write();
    h_R->Write();
    h_highct_L->Write();
    h_highct_R->Write();
    h_NP[iFile]->Write();
    h_PSR[iFile]->Write();
  }

  std::cout << "********** STORE TO FILES OF BkgHisto1DHists" << std::endl;
}

// ================================================================================
//               IMPLEMENTATION BKG HISTO PT RAP MASS HISTS
// ================================================================================
BkgHistoPtRapMassHists::BkgHistoPtRapMassHists()
{
  std::cout << "---------- CONSTRUCTOR OF BkgHistoPtRapMassHists" << std::endl;
  // TODO: replace NULL w/ nullptr for c++11
  hBG_L = NULL;
  hBG_R = NULL;
  h_highct_L = NULL;
  h_highct_R = NULL;
  hNP = NULL;
  hSR = NULL;

  std::cout << "********** CONSTRUCTOR OF BkgHistoPtRapMassHists" << std::endl;
}

BkgHistoPtRapMassHists::~BkgHistoPtRapMassHists()
{
  std::cout << "---------- DESTRUCTOR OF BkgHistoPtRapMassHists" << std::endl;
  // TODO: currently get a seg fault here. This is caused by the call to TFile::Close() prior to this
  // destructor, which apparently means for ROOT to delete TH3D that are "attatched" to the TFile
  // The question remains, why this works for TH2Ds (e.g.)?
  // delete hBG_L;
  // delete hBG_R;
  // delete h_highct_L;
  // delete h_highct_R;
  // delete hNP;
  // delete hSR;

  std::cout << "********** DESTRUCTOR OF BkgHistoPtRapMassHists" << std::endl;
}

void BkgHistoPtRapMassHists::storeToFile(TFile* file)
{
  std::cout << "---------- BkgHistoPtRapMassHists::storeToFile()" << std::endl;

  file->cd();
  hBG_L->Write();
  hBG_R->Write();
  h_highct_L->Write();
  h_highct_R->Write();
  hNP->Write();
  hSR->Write();

  std::cout << "********** BkgHistoPtRapMassHists::storeToFile()" << std::endl;
}

// ================================================================================
//                IMPLEMENTATION OF BKG HISTO FIT VARS STORE
// ================================================================================
template<>
inline void NamedVarStore<double>::setFromWS(RooWorkspace* ws, const std::string& wsName, const std::string& name)
{
  std::string storeName = name.empty() ? wsName : name; // determine name to be used for storage
  set(storeName, getVarVal(ws, wsName));
}

template<typename T>
inline void NamedVarStore<T>::setFromWS(RooWorkspace*, const std::string&, const std::string&)
{
  // NOP
  std::cerr << "NamedVarStore<T> only defined vor T = double" << std::endl;
}

template<typename T>
inline void NamedVarStore<T>::set(const std::string& name, const T value)
{
   m_store[name] = value;
   m_counter[name] = 0; // reset the counter
}

template<typename T>
inline T NamedVarStore<T>::operator[](const std::string& name) const
{
  // TODO: remove try-catch block after dev.
  try {
    // don't care if the variable exists or gets created here! If its not present all fails in the next line anyway
    m_counter[name]++;
    return m_store.at(name);
  } catch(std::out_of_range& ex) {
    std::cerr << ex.what() << ", no variable " << name << " stored!" << std::endl;
    throw; // rethrow as this is an error that cannot be handled in the calling function
  }
};

template<typename T>
std::ostream& operator<<(std::ostream& os, const NamedVarStore<T>& store)
{
  // os << "content of BkgHistoFitVarsStore:" << std::endl;
  // typedef for easier exchange of internal map
  typedef typename NamedVarStore<T>::MapT::const_iterator mapIt;

  for (mapIt it = store.m_store.begin(); it != store.m_store.end(); ++it) {
    os << it->first << " = " << it->second << ", ";
  }

  return os;
}

template<typename T>
const std::string NamedVarStore<T>::printUsages() const
{
  typedef std::map<std::string, size_t>::const_iterator mapIt;
  std::stringstream output;
  output << "usage counter: ";

  for (mapIt it = m_counter.begin(); it != m_counter.end(); ++it) {
    output << it->first << ": " << it->second << ", ";
  }

  return output.str();
}

#endif
