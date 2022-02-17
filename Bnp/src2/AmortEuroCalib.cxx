#include "AmortEuroCalib.h"
#include "AmortEuroPrice.h"

double Equi_Strike(long lToday, const char *szYC, const long *plFixStart_Begin,
                   const long *plFixStart_End, const long *plFixEnd_Begin,
                   const long *plFixPay_Begin, const double *pdFixCoupon_Begin,
                   const double *pdFixNotional_Begin, SrtBasisCode eFixBasis,
                   const long *plFltStart_Begin, const long *plFltStart_End,
                   const long *plFltEnd_Begin, const long *plFltPay_Begin,
                   const double *pdMargin_Begin, const double *pdSpread_Begin,
                   const double *pdFltNotional_Begin, SrtBasisCode eFltBasis,
                   double dEqui_Notional) {
  double dLevel_Flt, dLevel_Mult_Margin, dLevel_Mult_Spread;
  double dLevel_Fix, dLevel_Mult_Coupon;

  const double dFltPV = FltLeg(
      lToday, szYC, plFltStart_Begin, plFltStart_End, plFltEnd_Begin,
      plFltPay_Begin, pdMargin_Begin, pdSpread_Begin, pdFltNotional_Begin,
      eFltBasis, &dLevel_Flt, &dLevel_Mult_Margin, &dLevel_Mult_Spread);

  const double dFixPV =
      FixLeg(lToday, szYC, plFixStart_Begin, plFixStart_End, plFixEnd_Begin,
             plFixPay_Begin, pdFixCoupon_Begin, pdFixNotional_Begin, eFixBasis,
             &dLevel_Fix, &dLevel_Mult_Coupon);

  // floating size
  const int nFltSize = plFltStart_End - plFltStart_Begin;
  // floating leg first df
  const double dDF_Front = swp_f_df(lToday, plFltStart_Begin[0], szYC);
  // floating leg last df
  const double dDF_Back = swp_f_df(lToday, plFltPay_Begin[nFltSize - 1], szYC);
  // floating leg of a regular swap (const notional 1.)
  const double dFltPV_Regular = dDF_Front - dDF_Back;

  // match cash flow to return equivalent strike
  double dResult = (dFixPV - dFltPV) / dEqui_Notional;
  dResult += dFltPV_Regular;
  return dResult / dLevel_Fix;
}

double Equi_Notional(long lToday, const char *szYC,
                     const long *plFltStart_Begin, const long *plFltStart_End,
                     const long *plFltEnd_Begin, const long *plFltPay_Begin,
                     const double *pdMargin_Begin, const double *pdSpread_Begin,
                     const double *pdFltNotional_Begin, SrtBasisCode eBasis) {
  const int nSize = plFltStart_End - plFltStart_Begin;
  const double dDF_Front = swp_f_df(lToday, plFltStart_Begin[0], szYC);
  const double dDF_Back = swp_f_df(lToday, plFltPay_Begin[nSize - 1], szYC);
  double dLevel_Flt, dLevel_Mult_Margin, dLevel_Mult_Spread;
  const double dFltPV = FltLeg(
      lToday, szYC, plFltStart_Begin, plFltStart_End, plFltEnd_Begin,
      plFltPay_Begin, pdMargin_Begin, pdSpread_Begin, pdFltNotional_Begin,
      eBasis, &dLevel_Flt, &dLevel_Mult_Margin, &dLevel_Mult_Spread);

  _ASSERTE((dDF_Front - dDF_Back + dLevel_Mult_Spread + dLevel_Mult_Spread) >
           0.);

  return dFltPV /
         (dDF_Front - dDF_Back + dLevel_Mult_Margin + dLevel_Mult_Spread);
}
