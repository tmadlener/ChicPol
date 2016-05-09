#ifndef BKGHISTOS_HELPER_H__
#define BKGHISTOS_HELPER_H__
// this header defines some functions which are commonly used throughout the bkgHistos

#include "rootIncludes.inc"

#include <climits> // for CHAR_BIT
#include <string>
#include <sstream>
#include <vector>

/** find the power of 2 which is larger then the passed number. */
int findEvenNum(double number)
{
  // define the maximum number of allowed left bit shifts before shifting into the sign bit
  // should evaluate to (by definition at least) 31 for 'normal' ints
  const unsigned short maxPosBits = sizeof(int) * CHAR_BIT - 1;
  for (int n = 0; n < maxPosBits; ++n) { // going further resutls in an integer overflow
    const int pow2n = (1 << n); // use bitshifting for calculating 2^n
    const int pow2np1 = pow2n*2; // == 2^(n+1) NOTE: this results in an overflow for n == maxPosBits - 1
    if (number >= pow2n && number <= pow2np1) { // The overflow is 'caught' here and the default return value is used
      return pow2np1;
    }
  }
  return (1 << (maxPosBits - 1)); // max. power of two storable in an int
}

/** (re)create all TH2Ds in the passed vector with the same nameBase, a (frame) label and a given suffix with a
 * given number
 * of bins in X and Y direction.
 * NOTE: if a histogram is present it will be replaced!
 * TODO: make the const char** a vector of strings
 * NOTE: this relies very much on the right initialization of everything so this is a likely breaking point!
 */
void createHists(std::vector<TH2D*> hists, const std::string& nameBase, /*const std::string& titleBase*/
                 const char* labels[], const std::string& suffix, const int nBinsX, const int nBinsY,
                 const double xMin, const double xMax, const double yMin, const double yMax)
{
  std::cout << "### CREATE HISTS: " << nameBase << ", " << suffix << std::endl; // dev / testing

  for (size_t iH = 0; iH < hists.size(); ++iH) {
    std::stringstream name;
    name << nameBase << labels[iH] << suffix;
    std::stringstream title;
    title << ";cos#theta_{" << labels[iH] << "};#phi_{" << labels[iH] << "} [deg]";

    if (hists[iH]) delete hists[iH]; // destroy currently stored object if present
    hists[iH] = new TH2D(name.str().c_str(), title.str().c_str(),
                         nBinsX, xMin, xMax, nBinsY, yMin, yMax);
    hists[iH]->Sumw2();
  }
}

/**
 * Get the value from the RooRealVar stored in the workspace by name.
 * small helper function for less typing effort. COULDDO: safety measures.
 */
inline double getVarVal(RooWorkspace* ws, const std::string& name)
{
  // for development:
  std::cout << "### getVarVal, " << name << ": " << static_cast<RooRealVar*>(ws->var(name.c_str())) << std::endl;

  return static_cast<RooRealVar*>(ws->var(name.c_str()))->getVal();
}

/**
 * Get the error from the RooRealVar stored in the workspace by name.
 * small helper function for less typing effort. COULDDO: safety measures.
 */
inline double getVarError(RooWorkspace* ws, const std::string& name)
{
  // for development:
  std::cout << "### getVarError, " << name << ": " << static_cast<RooRealVar*>(ws->var(name.c_str())) << std::endl;

  return static_cast<RooRealVar*>(ws->var(name.c_str()))->getError();
}

/**
 * Get the pdf stored in the workspace by name
 * small helper function for less typing effort. COULDDO: safety measures.
 */
inline RooAbsPdf* getPdf(RooWorkspace* ws, const std::string& name)
{
  // for develpment:
  std::cout << "### getPdf, " << name << ": " << static_cast<RooAbsPdf*>(ws->pdf(name.c_str())) << std::endl;

  return static_cast<RooAbsPdf*>(ws->pdf(name.c_str()));
}

/** calculate the rapidity dependent mass sigma: \f$\sigma_{m}(|y|) = p_{0} + p_{1} * |y| + p_{2} *|y|^2\f$ */
inline double rapSigma(double p0, double p1, double p2, double rap)
{
  return p0 + p1 * TMath::Abs(rap) + p2 * TMath::Power(TMath::Abs(rap),2);
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

/**
 * create and store a TH1D with 3 bins and store the three passed values in the bins and the histo in the passed
 * TFile.
 * COULDDO: generalize this to N factors
 */
void store3Factors(TFile* file, const std::string& name, const std::string& title,
                   const double d1, const double d2, const double d3)
{
  file->cd();
  TH1D* h = new TH1D(name.c_str(), title.c_str(), 3, 0, 3.);
  h->SetBinContent(1, d1);
  h->SetBinContent(2, d2);
  h->SetBinContent(3, d3);
  h->Write();
  delete h;
}

#endif
