///////////////////////////////////////////////////////////////////////////////
// File:             common/src/finance/exoticoption/curran.cpp                     
// Purpose:          Price Asian option using Curran's formula
// Author:           ITO 33 Canada
// Created:          April 27, 2005
// RCS-ID:           $Id: curran.cpp,v 1.5 2006/07/21 21:55:40 dave Exp $                                                   
// Copyright         (C) 2005 - Trilemma LLP                                                 
///////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/finance/exoticoption/curran.cpp 

    @brief Curran Formula to compute discretely observed
           Asian option using some approxilation. This formula
           is very often used on the street for simple cases.
 */


// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/beforestd.h"
#include <cmath>
#include <vector>
#include "ito33/afterstd.h"
#include "ito33/useexception.h"

#include "ito33/numeric/normaldistribution.h"

#include "ito33/finance/exoticoption/curran.h"
#include "ito33/finance/optionerror.h"

extern const ito33::Error ITO33_BAD_PARAM;

extern const ito33::finance::Error ITO33_NEG_SPOT;

extern const ito33::finance::OptionError ITO33_OPTION_NEGATIVE_STRIKE;

using namespace ito33::numeric;

namespace ito33
{

namespace finance
{


void CheckParam(double dSpot, 
                  double dStrike, 
                  double dMaturityTime,
                  size_t nNbAveraging,
                  double dStartTime, 
                  double dVolatility, 
                  double dInterestRate, 
                  double dDividend,
                  double dA,
                  size_t nNbIntoAveraging)
{

  CHECK_COND( dSpot > 0, ITO33_NEG_SPOT );

  CHECK_COND( dStrike > 0, ITO33_OPTION_NEGATIVE_STRIKE );

  if ( dMaturityTime < 0.0)
    throw EXCEPTION_MSG
          (
           ITO33_BAD_PARAM,
           TRANS("Curran: negative maturity time.")
          );
   
  if ( dVolatility < 0.0 || dVolatility > 5.)
    throw EXCEPTION_MSG
          (
           ITO33_BAD_PARAM,
           TRANS("Curran: negative volatility or above 500%.")
          );
 
  if ( dDividend < 0.0)
    throw EXCEPTION_MSG
          (
           ITO33_BAD_PARAM,
           TRANS("Curran: negative dividend.")
          );

  if ( dA < 0.0 )
     throw EXCEPTION_MSG
           (
            ITO33_BAD_PARAM,
            TRANS("Curran: negative current average.")
           );

  if ( nNbIntoAveraging > nNbAveraging )
    throw EXCEPTION_MSG
          (
           ITO33_BAD_PARAM,
           TRANS("Curran: number of averaging point is inconsistent.")
          );

  if ( dStartTime < 0.0 )
    throw EXCEPTION_MSG
          (
           ITO33_BAD_PARAM,
           TRANS("Curran: negative start time.")
          );

  if ( dInterestRate < -1.0  || dInterestRate > 1.0)
    throw EXCEPTION_MSG
          (
           ITO33_BAD_PARAM,
           TRANS("Curran: interest rate too large or too small.")
          );

} //CheckParam

double PriceCurran(int nIsCall,
                   double dSpot, 
                   double dStrike, 
                   double dMaturityTime,
                   size_t nNbAveraging,
                   double dStartTime, 
                   double dVolatility, 
                   double dInterestRate, 
                   double dDividend,
                   double dA,
                   size_t nNbIntoAveraging)
{
  
  if ( nNbIntoAveraging > 0 )
    dStrike = double(nNbAveraging)/double(nNbAveraging - nNbIntoAveraging) * dStrike
    - double(nNbIntoAveraging)/double(nNbAveraging - nNbIntoAveraging) * dA;

  double dDeltaT       = (dMaturityTime - dStartTime)/(nNbAveraging - 1.0);

  double dVolSquared = dVolatility*dVolatility;

  double dMu = log(dSpot) + (dInterestRate - dDividend - dVolSquared / 2.) *
    ( dStartTime + double(nNbAveraging - 1) * dDeltaT / 2.);

  double dVolXSquared = dVolSquared*( dStartTime + 
    dDeltaT * double(nNbAveraging - 1)
    * double(2 * nNbAveraging-1) / 6. / double(nNbAveraging) );

  double dVolX = sqrt(dVolXSquared);

  double dTemp_K = 0.0;
  size_t nIdx;
  std::vector<double> pdTI(nNbAveraging+1);
  std::vector<double> pdMuI(nNbAveraging+1);
  std::vector<double> pdVolISquared(nNbAveraging+1);
  std::vector<double> pdVolXI(nNbAveraging+1);
  std::vector<double> pdVolXISquared(nNbAveraging+1);

 for (  nIdx = 1; nIdx <= nNbAveraging; nIdx++)
 {
   pdTI[nIdx] = dStartTime + (nIdx-1)*dDeltaT;

  
   pdMuI[nIdx]= log(dSpot) 
    + ( dInterestRate - dDividend - dVolSquared/2. )*pdTI[nIdx];

   pdVolISquared[nIdx] = dVolSquared*( dStartTime  + (nIdx-1)*dDeltaT );

   pdVolXI[nIdx] = dVolSquared*( dStartTime 
     + dDeltaT*( (nIdx-1) - nIdx*(nIdx-1)/2./nNbAveraging) );

   pdVolXISquared[nIdx] = pdVolXI[nIdx]*pdVolXI[nIdx];

   dTemp_K += 
    exp( pdMuI[nIdx] + (pdVolXI[nIdx]/dVolXSquared)*(log(dStrike) - dMu) 
    + ( pdVolISquared[nIdx] - pdVolXISquared[nIdx]/dVolXSquared )/2. );

 } //end for loop

 double Khat = 2.*dStrike - dTemp_K/double(nNbAveraging);
 double dTemp_p = 0.0;

 for (nIdx = 1; nIdx <= nNbAveraging; nIdx++)
 { 

   dTemp_p += exp( pdMuI[nIdx] + pdVolISquared[nIdx]/2. ) *
           CumulatedNormal( nIsCall*( (dMu - log(Khat))/dVolX + pdVolXI[nIdx]/dVolX ) );
 }
 
 dTemp_p = dTemp_p/double(nNbAveraging) - dStrike*CumulatedNormal( nIsCall*(dMu - log(Khat))/dVolX );
 
 dTemp_p = exp(-dInterestRate*dMaturityTime)*nIsCall*dTemp_p
   *double(nNbAveraging - nNbIntoAveraging)/double(nNbAveraging);

 return dTemp_p;
}

double CurranCall(double dSpot, 
                  double dStrike, 
                  double dMaturityTime,
                  size_t nNbAveraging, 
                  double dStartTime, 
                  double dVolatility, 
                  double dInterestRate, 
                  double dDividend,
                  double dA,
                  size_t nNbIntoAveraging)
{
  
  CheckParam(dSpot, dStrike, dMaturityTime, nNbAveraging, dStartTime, 
             dVolatility, dInterestRate, dDividend, dA, nNbIntoAveraging);

   int nIsCall = 1; //solving a call option here
   
   return PriceCurran(nIsCall, dSpot, dStrike, dMaturityTime, nNbAveraging, dStartTime, 
               dVolatility,  dInterestRate,  dDividend, dA, nNbIntoAveraging);
} // CurranCall

double CurranPut(double dSpot, 
                  double dStrike, 
                  double dMaturityTime,
                  size_t nNbAveraging,
                  double dStartTime, 
                  double dVolatility, 
                  double dInterestRate, 
                  double dDividend,      
                  double dA,
                  size_t nNbIntoAveraging)
{
  CheckParam(dSpot, dStrike, dMaturityTime, nNbAveraging, dStartTime, 
             dVolatility, dInterestRate, dDividend, dA, nNbIntoAveraging);

   int nIsCall = -1; //solving a put option here
   
   return PriceCurran(nIsCall, dSpot, dStrike, dMaturityTime, nNbAveraging, dStartTime, 
               dVolatility,  dInterestRate,  dDividend, dA, nNbIntoAveraging);
 
} //CurranPut


} // namespace finance

} // namespace ito33
