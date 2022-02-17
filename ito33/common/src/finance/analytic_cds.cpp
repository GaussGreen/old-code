///////////////////////////////////////////////////////////////////////////////
// File:             common/src/finance/analytic_cds.cpp
// Purpose:          Analytic solution for cds 
// Author:           ITO 33                                           
// Created:          14/02/2005
// RCS-ID:           $Id: analytic_cds.cpp,v 1.9 2006/08/21 14:26:44 wang Exp $                                                  
// Copyright         (c) 2005 - 2006 Trilemma LLP     
///////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include <cmath>

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"
#include "ito33/dateutils.h"

#include "ito33/finance/error.h"
#include "ito33/finance/cdslike.h"
#include "ito33/finance/cashflowstream_uniform.h"

extern const ito33::finance::Error ITO33_MATURITYBEFOREVALUATION;

namespace ito33
{

namespace finance
{

double AnalyticCDSPrice(const CDSLike& cds, double dRate, double dLambda, 
                        Date valuationDate)
{
  // The variable names follow the equation in Philippe's IHG document
  
  // Some common terms and data extraction from the financial cds class
  shared_ptr<CashFlowStream> pSpreadStream = cds.GetSpreadStream();
  double dRecovery = cds.GetRecoveryRate();
  double dRatePlusLambda = dRate + dLambda;
  
  // Find the first spread date after the valuation date
  // Also find the first spread date before the pricing date (could be
  // the issue date of the cds)
  CashFlowStream::const_iterator pFirstSpread = pSpreadStream->begin();
 
  Date initialDate = cds.GetIssueDate();

  while (   pFirstSpread != pSpreadStream->end() 
         && pFirstSpread->first < valuationDate )
  {
    initialDate = pFirstSpread->first;
    ++pFirstSpread;
  }

  CHECK_COND(pFirstSpread != pSpreadStream->end(),
             ITO33_MATURITYBEFOREVALUATION);

  double dInitialTime = GetDoubleFrom(initialDate);
  double dFirstTime = GetDoubleFrom( pFirstSpread->first );
  double dValuationTime = GetDoubleFrom(valuationDate);
  double dLastTime = GetDoubleFrom(pSpreadStream->GetLastPaymentDate());

  // Compute the first term
  double dTerm1  = ( 1.0 - dRecovery ) * dLambda / dRatePlusLambda 
            * ( 1.0 - exp(- dRatePlusLambda * (dLastTime - dValuationTime)) );

  // Compute the second term
  double dDelta = dFirstTime - dInitialTime;

  double dDeltaTmp = dValuationTime - dInitialTime;
  
  double dExpTmp = exp( - dRatePlusLambda * dDeltaTmp );

  double dTerm2 = pFirstSpread->second * dLambda
                / (dRatePlusLambda * dRatePlusLambda * dDelta)
                * exp( dRatePlusLambda * dDeltaTmp );

  dTerm2 *= (  1. - dExpTmp - dRatePlusLambda * dDeltaTmp * dExpTmp);

  // compute the first sum
  CashFlowStream::const_iterator iter;
  double dSum1 = 0;
  for ( iter = pFirstSpread; iter != pSpreadStream->end(); ++iter)
  {
    double dt = GetDoubleFrom( iter->first ) - dValuationTime;

    dSum1 += iter->second * exp( - dRatePlusLambda * dt);
  }
 
  // compute the second sum
  double
    dTmp,
    dTmpOld = dInitialTime;

  double dSum2 = 0;

  for (iter = pFirstSpread; iter != pSpreadStream->end(); ++iter)
  {
    dTmp = GetDoubleFrom( iter->first );
    dDelta = dTmp - dTmpOld;

    double dExpTmp = exp( - dRatePlusLambda * dDelta );
    double 
      dTerm = iter->second * dLambda //* dExpTmp1 
            / (dRatePlusLambda * dRatePlusLambda * dDelta)
            * exp( dRatePlusLambda * (dValuationTime - dTmpOld) );

    dSum2 += dTerm 
           * ( 1. - dExpTmp - dRatePlusLambda * dDelta * dExpTmp);

    dTmpOld = dTmp;
  }

  // Add the terms together to get the analytic price
   double dPrice = dTerm1 - dSum1 + dTerm2 - dSum2;

  //double dPrice = dTerm1 - dSum1 + dRecovery * (dTerm2 - dSum2);

  return dPrice;

}

} // namespace finance

} // namespace ito33
