/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/volatilitytimeonly.h
// Purpose:     time only volatility class
// Author:      ZHANG Yunzhi, David
// Created:     2004/05/19
// RCS-ID:      $Id: volatilitytimeonly.h,v 1.17 2006/08/20 09:47:41 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
  @file ito33/ihg/volatilitytimeonly.h
  @brief time only volatility class
  */


#ifndef _ITO33_IHG_VOLATILITYTIMEONLY_H_
#define _ITO33_IHG_VOLATILITYTIMEONLY_H_

#include "ito33/ihg/volatilitywithtimecomponent.h"

namespace ito33
{


namespace ihg
{

  class VolatilityVisitor;

/// Time only volatility class
class ITO33_IHG_DLLDECL VolatilityTimeOnly : public VolatilityWithTimeComponent
{
public:
  
  /**
     @internal
     @brief Dummy ctor used for calibration

     @noexport
   */
  VolatilityTimeOnly() : VolatilityWithTimeComponent() { }

  /**
     @internal 
     @brief Ctor constructs from two arrays.
     
     @param pDates the times that split the tike into nNbTimes intervals
     @param pdValues the values at each intervals
     @param nNbTimes the number of time intervals

     @noexport
   */
  VolatilityTimeOnly(const Date* pDates, const double* pdValues, 
                     size_t nNbTimes);

  /**
     Ctor constructs from two vectors.

     @param dates the times that split the tike into nNbTimes intervals
     @param values the values at each intervals
   */
  VolatilityTimeOnly(const std::vector<Date>& dates,
                     const std::vector<double>& values);

  // Default dtor is ok.

  void GetVols(double dTime,
               const double *pdS,
               double *pdVols,
               size_t nNbS) const;

  shared_ptr<Volatility> Perturb(double dShift) const;

  void Dump(ito33::XML::Tag& tagParent) const;

  void Visit(VolatilityVisitor& visitor) const;

  bool IsTimeOnly() const 
  { 
    return true;
  }

  virtual 
    void GetModelParameters(finance::ModelParametersConsumer& visitor) const;

}; // class VolatilityTimeOnly


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_VOLATILITYTIMEONLY_H_
