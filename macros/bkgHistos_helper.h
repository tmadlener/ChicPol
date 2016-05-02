#ifndef BKGHISTOS_HELPER_H__
#define BKGHISTOS_HELPER_H__
// this header defines some functions which are commonly used throughout the bkgHistos

#include <climits> // for CHAR_BIT

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

#endif
