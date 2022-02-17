// prevent multiple inclusions
#pragma once

//////////////////////
//	warnings
#pragma warning(disable : 4786) //"identifier was truncated to '255' characters
                                //in the debug information"

// NB: force warnings for unused arguments and local parameters
#pragma warning(1 : 4100 4101)

#include "AmortMidatADIUtils.h"

double Equi_Strike(long lToday, const char *szYC, const long *plFixStart_Begin,
                   const long *plFixStart_End, const long *plFixEnd_Begin,
                   const long *plFixPay_Begin, const double *pdFixCoupon_Begin,
                   const double *pdFixNotional_Begin, SrtBasisCode eFixBasis,
                   const long *plFltStart_Begin, const long *plFltStart_End,
                   const long *plFltEnd_Begin, const long *plFltPay_Begin,
                   const double *pdMargin_Begin, const double *pdSpread_Begin,
                   const double *pdFltNotional_Begin, SrtBasisCode eFltBasis,
                   double dEqui_Notional);

double Equi_Notional(long lToday, const char *szYC, const long *plStart_Begin,
                     const long *plStart_End, const long *plEnd_Begin,
                     const long *plPay_Begin, const double *pdMargin_Begin,
                     const double *pdSpread_Begin,
                     const double *pdNotional_Begin, SrtBasisCode eBasis);
