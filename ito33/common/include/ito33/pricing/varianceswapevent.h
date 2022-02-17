/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/varianceswapevent.h
// Purpose:     variance swap 3d event
// Created:     2006/02/21
// RCS-ID:      $Id: varianceswapevent.h,v 1.9 2006/08/03 21:25:18 dave Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file  ito33/pricing/variance swapevent.h
    @brief variance swap 3d event.
 */

#ifndef _ITO33_PRICING_VARIANCESWAPEVENT_H_
#define _ITO33_PRICING_VARIANCESWAPEVENT_H_

#include "ito33/vector.h"

#include "ito33/finance/returntype.h"

#include "ito33/pricing/pathdepevent.h"
#include "ito33/pricing/varianceswap.h"
 
namespace ito33
{

namespace pricing
{

class PathDepStructure;

class VarianceSwapEvent : public PathDepEvent
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
  VarianceSwapEvent(double dTime, size_t nObservation, 
                    finance::ReturnType returnType, bool bIsLastEvent,
                    bool bHasSimilarityReduction,
                    VarianceSwap& varianceSwap)
    : PathDepEvent(dTime, dTime, dTime), 
      m_nObservation(nObservation),
      m_bIsLastEvent(bIsLastEvent),
      m_returnType(returnType),
      m_bHasSimilarityReduction(bHasSimilarityReduction),
      m_bRecurse(false),
      m_varianceSwap(varianceSwap)
  {    
  }

  virtual ~VarianceSwapEvent () {}
    

protected:

  virtual void ApplyAtStartTime(PathDepStructure& ) const {}

  virtual void ApplyAtEndTime(PathDepStructure& ) const {}

  virtual void ApplyAtTime(PathDepStructure& pathDepStruct) const;

  virtual void DoEvent(PathDepStructure& pathDepStruct) const;
  
  /**
      Copy values from the diagonal S = P for the path to save.

      Used for forward starting variance swaps.  At the start of the sampling
      period (observation number is zero), the previous share price will
      equal the share price.  For example, if S = 20, then the previous share
      price should be 20, and the value on the path to save for S = 20
      should come from the path P = 20.  This means copying along the
      diagonal.

      For homogeneous models, this means the price will be spot independent
      before the sampling start date.

      @param pathDepStruct the path dependent structure containing prices, grids
   */
  void CopyFromDiagonal(PathDepStructure& pathDepStruct) const;

protected:

  /// Observation number
  size_t m_nObservation;

  /// Indicates is this is the last event
  bool m_bIsLastEvent;

  /// Indicates if a similarity reduction is possible
  bool m_bHasSimilarityReduction;

  /// How returns are calculated
  finance::ReturnType m_returnType;

  /// Reference to underlying variance swap contract
  VarianceSwap& m_varianceSwap;

  /// Indicates if a recursive call to ApplyAtTime has been made
  mutable bool m_bRecurse;

  /// Linear interpolate in both the Z and P directions
  double LinearInterpolate(double dS, double dP, 
    double dZ, const std::vector<double>& pdPreviousGrid,
    const std::vector<double> &pdAvgSqrReturnGrid, 
    PathDepStructure& pathDepStruct) const;

  /// Quadratic interpolate in Z, linear in P
  double QuadraticInterpolate(double dS, double dP, 
    double dZ, double dRSquared,
    const std::vector<double>& pdPreviousGrid,
    const std::vector<double> &pdAvgSqrReturnGrid, 
    PathDepStructure& pathDepStruct) const;

  /// Upstream quadratic interpolation in Z, diagonal interpolation in P
  double QuadZDiagPInterpolate(
    double dP, double dZ, double dRSquared,
    const std::vector<double>& pdPreviousGrid,
    const std::vector<double> &pdAvgSqrReturnGrid) const;

  /// Store the similarity interpolant for each path
  /// This can be done since we have a similarity interpolation
  /// of degree 0 and S is equal to P.
  mutable std::vector<double> m_pdSimilarityValues;

private:

  NO_COPY_CLASS(VarianceSwapEvent);

};

} // namespace pricing

} // namespace ito33 

#endif // #ifndef _ITO33_PRICING_VARIANCESWAPEVENT_H_
