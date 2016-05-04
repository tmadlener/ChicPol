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
  std::cout << "CREATE HISTS: " << nameBase << ", " << suffix << std::endl; // dev / testing

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

#endif
