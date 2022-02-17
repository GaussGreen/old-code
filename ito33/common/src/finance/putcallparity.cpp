///////////////////////////////////////////////////////////////////////////////
// File:             ito33/finance/putcallparity.cpp
// Purpose:          Functions supporting put-call parity calculations
// Author:           ITO 33                                           
// Created:          2005/05/12
// RCS-ID:           $Id: putcallparity.cpp,v 1.5 2006/08/19 23:06:55 wang Exp $                                                  
// Copyright         (C) 2005 - Trilemma LLP
///////////////////////////////////////////////////////////////////////////////

/**

  @file ito33/finance/putcallparity.cpp

*/
#include "ito33/beforestd.h"
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/date.h"
#include "ito33/sharedptr.h"
#include "ito33/dateutils.h"
#include "ito33/useexception.h"

#include "ito33/finance/option.h"
#include "ito33/finance/yieldcurve.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/putcallparity.h"
#include "ito33/finance/modeloutput.h"

extern const ito33::Error ITO33_BAD_PARAM;

namespace ito33
{

namespace finance
{


/**
  Check Put call parity

  Put = Call - Call(K = 0) + K exp( -r * T )

  @param option option contract
  @param tm Theoretical model
  @param dTolerance tolerance of put call parity check


  @return true/false
*/
bool IsPutCallParityValid(const Option &option, const TheoreticalModel& tm, 
                          double dTolerance)
{

  if ( option.GetExerciseType() == ExerciseType_American)
    throw EXCEPTION_MSG
          (
           ITO33_BAD_PARAM,
           TRANS("Put Call parity does not apply for American option.")
          );

  shared_ptr<SessionData> pSessionData = option.GetSessionData();
  Date maturityDate                   = option.GetMaturityDate();
  double dStrike                      = option.GetStrike();

  shared_ptr<Option> pCallOption
     ( 
      new Option(dStrike, maturityDate, Option_Call, ExerciseType_European) 
     );
  pCallOption->SetSessionData( pSessionData );
  
  shared_ptr<Option> pCallOptionStrikeZero
    ( 
      new Option(1.e-12, maturityDate, Option_Call, ExerciseType_European) 
    );
  pCallOptionStrikeZero->SetSessionData( pSessionData );
  
  shared_ptr<Option> pPutOption
    ( 
      new Option(dStrike, maturityDate, Option_Put, ExerciseType_European) 
    );
  pPutOption->SetSessionData( pSessionData );

  // Extract the relevent data from the session
  Date valuationDate                = pSessionData->GetValuationDate();
  shared_ptr<YieldCurve> pYieldCurve = pSessionData->GetYieldCurve();

  // Get the yield rate
  std::vector<double> pdDates(2);
  pdDates[0] = GetDoubleFrom(valuationDate);
  pdDates[1] = GetDoubleFrom(maturityDate);
     
  std::vector<double> pdValues(2);
   
  pYieldCurve->GetDiscountFactor( &pdDates[0], &pdValues[0], 2);
  double dYieldRate = pdValues[1] / pdValues[0];

  //computation part
  double dCallPrice           = tm.Compute( *pCallOption )->GetPrice();
  double dCallPriceStrikeZero = tm.Compute(*pCallOptionStrikeZero)->GetPrice();
  double dPutPrice            = tm.Compute( *pPutOption )->GetPrice();

  // P = C - C( K = 0 ) + K exp(-r*T);
  double dPutCallParity = dCallPrice - 
    dCallPriceStrikeZero + dStrike*dYieldRate - dPutPrice;

  if ( fabs( dPutCallParity ) < dTolerance ) //penny accurate
    return true;

  return false;
}

} // namespace finance

} // namespace ito33
