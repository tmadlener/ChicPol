#ifndef BKGHISTOS_HELPER_H__
#define BKGHISTOS_HELPER_H__
// this header defines some functions which are commonly used throughout the bkgHistos

#include "rootIncludes.inc"
#include "commonVar.h" // included here because it is needed in calcPol.C but could break things there -> FIXME
#include "calcPol.C" // for calcPol

#include <climits> // for CHAR_BIT
#include <string>
#include <sstream>
#include <vector>
#include <algorithm> // for std::min
#include <iostream> // for std::getline

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

/** check if string starts with another string. */
inline bool startsWith(const std::string& input, const std::string& prefix)
{
  return input.substr(0, prefix.length()) == prefix;
}

/** remove all that is made up of leading characters defined as arguments (default to space and tab). */
std::string removeLeading(const std::string& str, const std::string& leading = " \t")
{
  const size_t pos = str.find_first_not_of(leading);
  if (pos == std::string::npos) return std::string(""); // only leading characters in input
  return str.substr(pos);
}

/** remove trailing characters (passed as arguments) */
std::string removeTrailing(const std::string& str, const std::string& trailing = " \t")
{
  const size_t pos = str.find_last_not_of(trailing);
  return str.substr(0, pos + 1);
}

/** trim string by removing all leading and trainling characters defined by argument chars. */
std::string trim(const std::string& str, const std::string& chars = " \t")
{
  return removeTrailing(removeLeading(str));
}

/** split string at delimiter delim and return vector of all substrings. If a token is empty it will be ignored. */
std::vector<std::string> splitString(const std::string& in, const char delim)
{
  std::vector<std::string> tokens;
  std::stringstream sstr(in);
  std::string tok;
  while(std::getline(sstr,tok,delim)) {
    if(!tok.empty()) {
      tokens.push_back(tok);
    }
  }

  return tokens;
}

/** (re)create all TH2Ds in the passed vector with the same nameBase, a (frame) label and a given suffix with a
 * given number
 * of bins in X and Y direction.
 * NOTE: if a histogram is present it will be replaced!
 * TODO: make the const char** a vector of strings
 * NOTE: this relies very much on the right initialization of everything so this is a likely breaking point!
 */
void createHists(std::vector<TH2D*>& hists, const std::string& nameBase, /*const std::string& titleBase*/
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
  RooRealVar* var = static_cast<RooRealVar*>(ws->var(name.c_str()));
  if (var) return var->getVal();
  // std::cout << name << " is a RooFormulaVar!" << std::endl;
  var = static_cast<RooRealVar*>(ws->function(name.c_str()));
  if (var) return var->getVal();
  std::cout << "Could not get " << name << " from workspace" << std::endl;
  return -99999;
}

/**
 * Get the error from the RooRealVar stored in the workspace by name.
 * small helper function for less typing effort. COULDDO: safety measures.
 */
inline double getVarError(RooWorkspace* ws, const std::string& name)
{
  // for development:
  // std::cout << "### getVarError, " << name << ": " << static_cast<RooRealVar*>(ws->var(name.c_str())) << std::endl;
  RooRealVar* var = static_cast<RooRealVar*>(ws->var(name.c_str()));
  if (var) return var->getError();
  var = static_cast<RooRealVar*>(ws->function(name.c_str()));
  if (var) return var->getError();
  std::cout << "Could not get " << name << " from workspace" << std::endl;

  return -99999;
}

/**
 * Get the pdf stored in the workspace by name
 * small helper function for less typing effort. COULDDO: safety measures.
 */
inline RooAbsPdf* getPdf(RooWorkspace* ws, const std::string& name)
{
  // for develpment:
  // std::cout << "### getPdf, " << name << ": " << static_cast<RooAbsPdf*>(ws->pdf(name.c_str())) << std::endl;

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

/**
 * calculate the cosTheta and phi values to be stored in the histograms.
 * Returntype is an onia::kNbFrames x 2 array of doubles, index 0 is CosTh, index 1 is Phi for each frame
 * COULDDO: define a nice(r) return type for this.
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

/**
 * Calculate the ratio of filled bins over total bins
 * hR, hL and hC have to have the samme binning! (at least the same bin numbers).
 * hC is used to estimate the coverage
 */
double calcCoverage(const TH2D& hR, const TH2D& hL, const TH2D& hC)
{
  unsigned int totBins = 0;
  unsigned int filledBins = 0;

  for (int iX = 0; iX < hR.GetNbinsX(); ++iX) {
    for (int iY = 0; iY < hR.GetNbinsY(); ++iY) {
      totBins++;
      if (hC.GetBinContent(iX + 1, iY + 1) > 0) filledBins++;
    }
  }

  std::cout << "------------------------------------------------------------" << std::endl
            << "filled bins: " << filledBins << std::endl
            << "total bins: " << totBins << std::endl;

  return (double) filledBins / (double) totBins;
}

/** mainly for less typing! */
inline bool minBinningReached(const int phi, const int phiMin, const int cosTh, const int cosThMin)
{
  return ((phi / 2 < phiMin) && (cosTh /2 < cosThMin));
}

/**
 * determine the binning for the cosThPhi histograms.
 * .first is Phi Binning, .second is CosTheta Binning
 */
std::pair<int, int> determineBinning(const TH2D& hR, const TH2D& hL, const TH2D& hC, bool folding)
{
  double coverage = 2 * calcCoverage(hR, hL, hC);
  int nBinsCosth = findEvenNum(16 * 2 / coverage); // find 2^n closest (but above) the passed number
  if (nBinsCosth > 64) nBinsCosth = 64; // maximum is 64 bins
  int nBinsPhi = 16;

  std::cout << "------------------------------------------------------------" << std::endl
            << "Calculating coverage and average number of bin contents" << std::endl
            << "bin coverage: " << coverage << std::endl
            << "starting point for binning in phi: " << nBinsPhi << std::endl
            << "starting point for binning in cosTheta: " << nBinsCosth << std::endl
            << "------------------------------------------------------------" << std::endl;

  // calculate the integral of the lowstatBG histo (check which of the histos has less events)
  int IntBG = std::min(hL.Integral(), hR.Integral());

  double Naverage = (double)IntBG/((double)nBinsPhi * nBinsCosth * coverage / 2.);
  std::cout << "average cell coverage: " << Naverage << std::endl;

  if (Naverage > 10) {
    std::cout << "Rebinning is not necessary in this case." << "\n"
              << "Ending binning algorithm." << "\n"
              << "------------------------------------------------" << std::endl;
  } else {
    std::cout << "------------------------------------------------" << "\n"
              << "old cosTheta binning: " << nBinsCosth << "\n"
              << "old phi binning: " << nBinsPhi << std::endl;

    nBinsPhi = findEvenNum(nBinsCosth * coverage / 2.);

    std::cout << "closest 2^n number to cosTheta bins: " << nBinsCosth << "\n"
              << "lowest 2^n number so that phi bins > cosTheta bins * coverage/2: " << nBinsPhi << "\n"
              << "------------------------------------------------" << std::endl;

    int nBinsPhiMin = folding ? 16 : 8;
    int nBinsCosthMin = 8;

    int iLoop = 0;
    while (Naverage <= 10. || !minBinningReached(nBinsPhi, nBinsPhiMin, nBinsCosth, nBinsCosthMin)) {
      iLoop++;
      std::cout << "looping for correct binning in background histogram" << std::endl;

      // Change binning of  in phi first, then in costh. Check inbetween and escape loop if necessary
      if (nBinsPhi / 2 >= nBinsPhiMin) nBinsPhi /= 2; // split only if minimum binning is ensured
      Naverage = (double)IntBG / ((double) nBinsPhi * nBinsCosth * coverage / 2.);
      std::cout << "average bin content per cell after " << iLoop << " phi rebinning: " << Naverage << std::endl;
      if (Naverage > 10) break;

      if (nBinsCosth / 2 >= nBinsCosthMin) nBinsCosth /= 2; // split only if minimum binning is ensured
      Naverage = (double)IntBG / ((double) nBinsPhi * nBinsCosth * coverage / 2.);
      std::cout << "avg bin content per cell after " << iLoop << " cosTheta rebinning:" << Naverage << std::endl;
    }

    std::cout << "average bin content per cell exceeds 10: " << Naverage << "\n"
              << "phi bins = " << nBinsPhi << ", cosTheta bins = " << nBinsCosth << "\n"
              << "------------------------------------------------" << std::endl;
  }

  return std::make_pair(nBinsPhi, nBinsCosth);
}

/**
 * Normalize the histogram to the passed value f.
 */
template<typename HistT>
inline void selfScale(HistT* h, const double f = 1.0)
{
  h->Scale(f / (1. * h->Integral()));
}

/**
 * Clone any TObject. Mainly for less typing effort.
 * TODO: increase usage.
 */
template<typename T>
inline T* clone(const T* obj, const std::string& name = "")
{
  return static_cast<T*>(obj->Clone(name.c_str()));
}

/**
 * add the two histograms with a scaling factor. Returns a new histogram (for which the caller takes ownership!).
 * If a file is passed as well, the returned histogram is also stored there prior to returning.
 *
 * Adding follows the following formula: \f$ h_{result} = f_{1} * h_{1} / N_{1} + (1 - f_{1}) * h_{2} / N_{2}) \f$,
 * where \f$ N_{i} \f$ is the integral of the \f$ i \f$-th histogram.
 *
 * NOTE: h1 and h2 get scaled by a call to this function!
 */
template<typename HistT>
HistT* addScaled(HistT* h1, HistT* h2, const double f1, const std::string& name, TFile* file = NULL)
{
  selfScale(h1, f1);
  selfScale(h2, 1 - f1);
  HistT* hr = static_cast<HistT*>(h1->Clone(name.c_str()));
  hr->Add(h2);

  if (file) {
    file->cd();
    hr->Write();
  }

  return hr;
}

/**
 * add the two histograms with two different weight factors w1 and w2.
 * If a file is passed the histogram is also stored there prior to returning.
 *
 * The returned histogram is calculated as: \f$ h_{result} = w_{1} * h_{1} / N_{1} + w_{2} * h_{2} / N_{2} \f$,
 * where \f$ N_{i} \f$ is the integral of the \f$ i \f$-th histogram.
 * NOTE: h1 and h2 get modified with a call to this function!
 */
template<typename HistT>
HistT* addWeighted(HistT* h1, HistT* h2, const double w1, const double w2, const std::string& name,
                   TFile* file = NULL)
{
  selfScale(h1, w1);
  selfScale(h2, w2);
  HistT* hr = static_cast<HistT*>(h2->Clone(name.c_str()));
  hr->Add(h1);

  if (file) {
    file->cd();
    hr->Write();
  }

  return hr;
}

/**
 * manually subract h2 from h1.
 * Result: h(i,j,k) = h1(i,j,k) - h2(i,j,k), if h1(i,j,k) > 0 and h1(i,j,k) - h2(i,j,k) > 0.
 *         h(i,j,k) = h1(i,j,k), if h1(i,j,k) <= 0.
 *         h(i,j,k) = 0, if h1(i,j,k) - h2(i,j,k) < 0
 * Errors are calculated as squared sum.
 * NOTE: h1 gets modified by a call to this!
 */
TH3D* subtractManual(TH3D* h1, const TH3D* h2)
{
  const int nx = h1->GetXaxis()->GetNbins();
  const int ny = h1->GetYaxis()->GetNbins();
  const int nz = h1->GetZaxis()->GetNbins();

  for (int j = 0; j <= nx; ++j) {
    for (int k = 0; k <= ny; ++k) {
      for (int l = 0; l <= nz; ++l) {
        const double c1 = h1->GetBinContent(j,k,l);
        if (c1 > 0) {
          const double c2 = h2->GetBinContent(j,k,l);
          const double e1 = h1->GetBinError(j,k,l);
          const double e2 = h2->GetBinError(j,k,l);
          double e3 = TMath::Sqrt(e1*e1 + e2*e2);
          double c3 = c1 - 1. * c2;
          if (c3 < 0 ) {
            c3 = 0;
            e3 = 0;
          }
          h1->SetBinContent(j,k,l, c3);
          h1->SetBinError(j,k,l, e3);
        }
      } // l
    } // k
  } // j

  return h1;
}

/**
 * manually subtract h2 from h1.
 * For a more detailed description see the TH3D version of this.
 * NOTE: modifies h1 within (also h1 is basically the return)
 */
TH2D* subtractManual(TH2D* h1, const TH2D* h2)
{
  int nx = h1->GetXaxis()->GetNbins();
  int ny = h1->GetYaxis()->GetNbins();

  for (int j = 0; j <= nx; ++j) {
    for (int k = 0; k <= ny; ++k) {
      const double c1 = h1->GetBinContent(j,k);
      if (c1 > 0) {
        const double c2 = h2->GetBinContent(j,k);
        const double e1 = h1->GetBinError(j,k);
        const double e2 = h2->GetBinError(j,k);
        double c3 = c1 - 1. * c2;
        double e3 = TMath::Sqrt(e1*e1 + e2*e2);
        if (c3 < 0) {
          c3 = 0;
          e3 = 0;
        }
        h1->SetBinContent(j,k, c3);
        h1->SetBinError(j,k, e3);
      }
    }
  }

  return h1;
}

/**
 * subtract h2 scaled by f from h1. Returns a new histogram (for which the caller takes ownership!).
 * If a file is passed as well, the returned histogram is also stored there prior to returning.
 *
 * Subtracting follows the following formula: \f$ h_{result} = h_{1} / N_{1} - f * h_{2} / N_{2} \f$,
 * where \f$ N_{i} \f$ is the integral of the \f$ i \f$-th histogram.
 * NOTE: h1 gets modified by a call to this!
 * NOTE: uses the subtractManual method to subtract h1 from h2!. (This has to be defined for the type!)
 */
template<typename HistT>
HistT* subtractScaled(HistT* h1, const HistT* h2, double f, const std::string& name, TFile* file = NULL)
{
  selfScale(h1);
  HistT* htmp = static_cast<HistT*>(h2->Clone());
  selfScale(htmp, f);
  HistT* hr = static_cast<HistT*>(h1->Clone(name.c_str()));
  hr = subtractManual(hr, htmp);
  if (file) {
    file->cd();
    hr->Write();
  }

  return hr;
}

/**
 * Subtracts h1 scaled by f from h2 (after normalization of both).
 * returns the mean of h2 after all operations.
 *
 * NOTE: normalizes h1 and modifies h2 beyond that!
 * TODO: this is a real beast, as it modifies all non-const inputs!
 */
double subtractScaledMean(TH1D* h1, TH1D* h2, const double f)
{
  selfScale(h1);
  TH1D* htmp = static_cast<TH1D*>(h1->Clone());
  selfScale(htmp, f);
  selfScale(h2);
  h2->Add(htmp, -1.);

  return h2->GetMean();
}

/**
 * Clones hL and hR, and adds them weighted.
 * hL after running: \f$ h_{L, out} = f_{1} * h_{L} / N_{L} + (1 - f_{1}) * h_R / N_{R} \f$,
 * returned hist: \f$ h_{result} = f_{2} * h_{L} / N_{L} + (1 - f_{2} * h_R / N_{R})\f$,
 * where \f$ N_{i} \f$ is the integral of the \f$ i\f$-th histogram.
 *
 * NOTE: modifies both passed hists.
 * TODO: this is another beast that is mainly present for less typing
 */
TH1D* addSideBands(TH1D* hL, TH1D* hR, const double f1, const double f2)
{
  TH1D* hL2 = static_cast<TH1D*>(hL->Clone());
  selfScale(hL, f1);
  selfScale(hL2, f2);

  TH1D* hR2 = static_cast<TH1D*>(hR->Clone());
  selfScale(hR, 1. - f1);
  selfScale(hR2, 1. - f2);

  hL->Add(hR);
  hL2->Add(hR2);

  return hL2;
}

/** convert an array of c-style strings (char*) to a vector of strings. */
std::vector<std::string> charArrayToStrVec(const char** charArr, const int arrSize)
{
  std::vector<std::string> retVec;
  for (int iL = 0; iL < arrSize; ++iL) {
    std::stringstream tmp;
    tmp << charArr[iL];
    retVec.push_back(tmp.str());
  }
  return retVec;
}

/** unfold the passed TH2D.*/
void unfold(TH2D* h)
{
  const int nx = h->GetXaxis()->GetNbins();
  const int ny = h->GetYaxis()->GetNbins();
  const int yPhi = ny / 4;

  for (int j = 0; j <= nx; ++j) {
    for (int k = 2 * yPhi + 1; k <= 3*yPhi; ++k) {
      const double c = h->GetBinContent(j,k);
      const double e = h->GetBinError(j,k);

      // flip in cosTheta
      const int l = nx + 1 - j;

      // set bin content and error of phiFolded in the other 3 (not yet filled) phi regions
      // 90 - 180: flip phi (upwards), flip cosTheta
      h->SetBinContent(l, 6 * yPhi + 1 - k, c);
      h->SetBinError(l, 6 * yPhi + 1 - k, e);
      // 0 - -90: flip phi (downwards)
      h->SetBinContent(j, ny + 1 - k, c);
      h->SetBinError(j, ny + 1 - k, e);
      // -90 - -180: flip cosTheta, shift phi
      h->SetBinContent(l, k - 2 * yPhi, c);
      h->SetBinError(l, k - 2 * yPhi, e);
    }
  }
}

#endif
