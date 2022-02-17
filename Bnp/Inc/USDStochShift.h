#ifndef _USD_STOCH_SHIFT_H_
#define _USD_STOCH_SHIFT_H_

struct StochShift_Param {
  enum SS_MODEL { SS_EXACT, SS_APPROX } iModel;
  int iIntPoints;
};

char *StochShift_PriceCall(double in_dF, double in_dL, double in_dalpha,
                           double in_dsigma, double in_drho, double in_dK,
                           double in_dTau, struct StochShift_Param *in_Param,
                           double *out_dPrice);

#endif