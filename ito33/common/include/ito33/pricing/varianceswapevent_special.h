/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/varianceswapevent_special.h
// Purpose:     In the case of a variance swap with no cap
//              where the volatility and hazard rate are not
//              spot dependent it is then possible to solve
//              a variance swap using only two paths.
// Created:     June 8, 2006
// RCS-ID:      $Id: varianceswapevent_special.h,v 1.4 2006/07/25 20:09:42 dave Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file  ito33/pricing/varianceswapevent_special.h
    @brief variance swap special event only two paths are required.
 */

#ifndef _ITO33_PRICING_VARIANCESWAPEVENT_SPECIAL_H_
#define _ITO33_PRICING_VARIANCESWAPEVENT_SPECIAL_H_

#include "ito33/vector.h"

#include "ito33/finance/returntype.h"

#include "ito33/pricing/varianceswapevent.h"
 
namespace ito33
{

namespace pricing
{

class PathDepStructure;

class VarianceSwapEventSpecial : public VarianceSwapEvent
{

public:

  /** 
      Constructor.

      @param dTime time of the event
      @param nObservation observation number
      @param bIsLastEvent indicate if this is the last event
      @param bHasSimilarityReduction indicate if a similarity reduction 
             is possible
      @param varianceSwap the underlying variance swap from which we get
             corridors, payoff type, etc
   */
  VarianceSwapEventSpecial(double dTime, size_t nObservation, 
                    finance::ReturnType returnType, bool bIsLastEvent, 
                    bool bHasSimilarityReduction,
                    VarianceSwap& varianceSwap)
    : VarianceSwapEvent(dTime, nObservation, returnType, bIsLastEvent, 
       bHasSimilarityReduction, varianceSwap)
  {    
  }

  virtual ~VarianceSwapEventSpecial () {}
    

protected:

  virtual void DoEvent(PathDepStructure& pathDepStruct) const;

  double Interpolate(double dP, double dZ, 
          const std::vector<double>& pdPreviousGrid,
          const std::vector<double> &pdAvgSqrReturnGrid) const;
};

} // namespace pricing

} // namespace ito33 

#endif // #ifndef _ITO33_PRICING_VARIANCESWAPEVENT_SPECIAL_H_
