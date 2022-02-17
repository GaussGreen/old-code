/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/impliedparametercalculator.h
// Purpose:     global functions for computing implied spread/strike etc
// Created:     2006/06/01
// RCS-ID:      $Id: impliedparametercalculator.h,v 1.5 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/impliedparametercalculator.h
    @brief global functions for computing implied spread/strike etc
 */

#ifndef _ITO33_FINANCE_IMPLIEDPARAMETERCALCULATOR_H_
#define _ITO33_FINANCE_IMPLIEDPARAMETERCALCULATOR_H_

#include "ito33/common.h"

#ifdef __CPP2ANY__
#include "ito33/finance/theoreticalmodel.h"
#endif

namespace ito33 
{
  
namespace finance
{
  class ITO33_DLLDECL SessionData;
  class ITO33_DLLDECL TheoreticalModel;
  
  class ITO33_DLLDECL CDSLike;
  class ITO33_DLLDECL EDS;
  class ITO33_DLLDECL VarianceSwapTerms;
  class ITO33_DLLDECL VarianceSwapLike;

/**
    class contains global functions for computing implied spread/strike etc.
 */
class ITO33_DLLDECL ImpliedParameterCalculator
{
public:
  /**
      Computes the implied CDS spread.

      @param cds The CDS like that we want to compute the implied spread for.
      @param model The model used to compute the implied spread.

      @return the implied spread of the CDS like
   */
  static double   
  ComputeImpliedCDSSpread(const CDSLike& cds, 
                          const TheoreticalModel& model); 

  /**
      Computes the implied EDS spread.

      @param eds The EDS that we want to compute the implied spread for.
      @param model The model used to compute the implied spread.

      @return The implied spread of the EDS
   */
  static double 
  ComputeImpliedEDSSpread(const EDS& eds, const TheoreticalModel& model);

  /**
      Computes the implied volatility strike for a variance swap.

      @param varianceSwap The variance swap that we want to compute the
                          implied volatility strike.
      @param model The model used to compute the implied volatility strike.

      @return the implied volatility strike of the variance swap
   */
  static double 
  ComputeImpliedVolatilityStrike(const VarianceSwapLike& varianceSwap,
                                 const TheoreticalModel& model);

  /**
      Computes the implied volatility strike for a variance swap whose
      sampling is not yet started.

      @param pTerms The variance swap terms.
      @param pSessionData The session data of the variance swap
      @param model The model used to compute the implied volatility strike.

      @return the implied volatility strike of the variance swap
   */
  static double 
  ComputeImpliedVolatilityStrike(const shared_ptr<VarianceSwapTerms>& pTerms,
                                 const shared_ptr<SessionData>& pSessionData,
                                 const TheoreticalModel& model);
};

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_IMPLIEDPARAMETERCALCULATOR_H_
