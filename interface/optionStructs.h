#ifndef OPTION_STRUCTS_H__
#define OPTION_STRUCTS_H__

/** struct to store the options for evaluating efficiencies stored as TGraphAsymmErrors. */
struct EffEvalOptions {
  /** default constructor. Initializes everything to the "normal" use case with no shifts and variations. */
  EffEvalOptions(bool stat = false, bool up = false, bool down = false, double sig = 1.0) :
    statVar(stat), shiftUp(up), shiftDown(down), nSigma(sig) {;}

  /**
   * vary the efficieny according to a normal distribution where the mean is the efficiency and the sigma is the
   * error times the nSigma below.
   */
  bool statVar;
  bool shiftUp; /**< shift the efficiency up by nSigma * error. */
  bool shiftDown; /**< shift the efficiency down by nSigma * error. */
  double nSigma; /**< nSigma controls how strong the variation/shift is. */
};

#endif
