///////////////////////////////////////////////////////////////////////////////
// File:         ito33/finance/analytic_cds.h  
// Purpose:      Analytical solution for cds
// Author:       ITO 33
// Created:      14/02/2005
// RCS-ID:       $Id: analytic_cds.h,v 1.6 2006/08/21 14:26:14 wang Exp $                                                    
// Copyright     (c) 2005 - 2006 Trilemma LLP                                                
///////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/analytic_cds.h 

    @brief Analytical solution for computing cds.
 */

#ifndef _ITO33_FINANCE_ANALYTICALCDS_H_
#define _ITO33_FINANCE_ANALYTICALCDS_H_

#include "ito33/date.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL CDSLike;

/**
    Analytic CDS formula assuming constant interest rate and hazard rate.
    Note that dRate must be the continuous interest rate.

    @param cds The cds to be priced
    @param dRate The continuous, constant interest rate 
    @param dLambda The constant hazard rate
    @param valuationDate The date at which to compute the price

    @return The analytic price of the cds
 */
double AnalyticCDSPrice(const CDSLike& cds, double dRate, double dLambda, 
                        Date valuationDate);

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_ANALYTICALCDS_H_
