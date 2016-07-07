#include "rootIncludes.inc"
// #include "effsAndCuts_chi.h"

#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>


typedef TGraphAsymmErrors EffFuncType;

// some global vars
const double etaRange[] = {0., 0.2, 0.3, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};
const int bins = 1;

int getBin(const TGraphAsymmErrors* f, double pT)
{
  double x,y; f->GetPoint(0, x, y);
  if (pT < x - f->GetErrorXlow(0)) {
    std::cerr << "pT = " << pT << " is below validity of " << f->GetName() << std::endl;
    return -1;
  }

  for (int i = 0; i < f->GetN(); ++i) {
    f->GetPoint(i, x, y);
    if (pT > x - f->GetErrorXlow(i) && pT <= x + f->GetErrorXhigh(i)) {
      return i;
    }
  }

  std::cerr << "Reached upper bound for pT = " << pT << " in " << f->GetName() << std::endl;
  return f->GetN() - 1;
}

/**
 * linearly interpolate between (x0,y0) and (x1,y1) given x.
 * /f$ y = y_{0} + (y_{1} - y_{0}) \frac{x - x_{0}}{x_{1} - x_{0}}/f$
 */
inline double linInt(double x, double x0, double x1, double y0, double y1)
{
  return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}

/** calculate the error via a linear interpolation between the two points of TGraphAsymmErrors enclosing pT. */
double calculateError(TGraphAsymmErrors* f, double pT, std::ostream& os)
{
  int pTbin = getBin(f, pT);
  int extraBin; // bin needed to get the errors, etc. for interpolation

  double retVal = -1;

  // do not extrapolate below the threshold of the TGraph; be conservative use larger error of lowest bin
  if (pTbin < 0) /*return*/ retVal = std::max(f->GetErrorYhigh(0), f->GetErrorYlow(0));


  double x0, y0; f->GetPoint(pTbin, x0, y0);
  double x1, y1;
  if (pT > x0) { // this means that we need the bin to our right for interpolation
    extraBin = pTbin + 1;
    f->GetPoint(extraBin, x1, y1);
  } else { // need bin to left. -> Swap the points (x0,y0) and (x1,y1) because 0 index is for left point!
    extraBin = pTbin; pTbin--;
    x1 = x0; y1 = y0;
    f->GetPoint(pTbin, x0, y0);
  }

  // if we have no left sample point, we cannot interpolate! (Only happens when point is in the lower half of the lowest bin)
  if (pTbin < 0) /*return*/ retVal = std::max(f->GetErrorYhigh(0), f->GetErrorYlow(0));
  // if we have no right sample point, we cannot interpolate! (Only happens when point is in the upper half of the highest bin)
  if (extraBin >= f->GetN()) /*return*/ retVal = std::max(f->GetErrorYhigh(f->GetN() - 1), f->GetErrorYlow(f->GetN() - 1));

  double eh0 = f->GetErrorYhigh(pTbin);
  double el0 = f->GetErrorYlow(pTbin);
  double eh1 = f->GetErrorYhigh(extraBin);
  double el1 = f->GetErrorYlow(extraBin);

  // calculate the absolute values of the error bands and correct them for the central value
  double errHigh = linInt(pT, x0, x1, y0 + eh0, y1 + eh1);
  double errLow = linInt(pT, x0, x1, y0 - el0, y1 - el1);
  double eff = linInt(pT, x0, x1, y0, y1);

  if (pTbin >= 0 && extraBin < f->GetN()) { // just for dev!!
    /*return*/ retVal = std::max(errHigh - eff, eff - errLow); // want absolute values!
  }

  os << pT << " " << pTbin << " " << extraBin << " " << x0 << " " << y0 << " " << x1 << " " << y1 << " "
     << eh0 << " " << el0 << " " << eh1 << " " << el1 << " " << errHigh << " " << errLow << " " << eff << " "
     << errHigh - eff << " " << eff - errLow << " " << retVal;

  return retVal;
}

double evalParametrizedEff(double pT, double eta, EffFuncType* f, bool statVar, std::ostream& os)
{
  if (TMath::Abs(eta) > 1.8) return 0;
  double eff = f->Eval(pT);

  os << eff << " ";

  if (statVar) {
    double err = calculateError(f, pT, os);
    eff = gRandom->Gaus(eff, err);

    os << " " << eff << std::endl;
  }

  if (eff > 1.) eff = 1;
  if (eff < 0) eff = 0;
  return eff;
}

int main(int argc, char* argv[])
{
  TFile* inFile = TFile::Open(argv[1]);

  for (int ieta = 0; ieta < bins; ++ieta) {
    double pT = 0;
    std::stringstream ss;
    ss << "evalEff_" << "MC" << "_etaBin_" << ieta << ".out";
    std::ofstream ofs(ss.str());

    std::stringstream ess;
    ess << "gEffHybrid_" << "MC" << "_PT_AETA" << ieta;
    EffFuncType* func = dynamic_cast<EffFuncType*>(inFile->Get(ess.str().c_str()));

    while (pT < 150) {
      double eta = 1.0;
      evalParametrizedEff(pT, eta, func, true, ofs);
      pT += 0.1;
    }
    ofs.close();
  }

  return 0;
}
