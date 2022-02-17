/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/hazardratepower.h
// Purpose:     spot power hazard rate class
// Author:      Wang
// Created:     2004/06/11
// RCS-ID:      $Id: hazardratepower.h,v 1.23 2005/12/02 15:52:57 zhang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/ihg/hazardratepower.h
   @brief Spot power (space only) hazard rate class
 */


#ifndef _ITO33_IHG_HAZARDRATEPOWER_H_
#define _ITO33_IHG_HAZARDRATEPOWER_H_

#include "ito33/ihg/hazardrate.h"

namespace ito33
{

namespace ihg
{

 class HazardRateVisitor;

/** 
    Class for a hazard rate as a power function of the spot.

    \f$ \lambda (S,t) = \alpha (S_0 / S)^\beta \f$ 

    This is a typical parametrization of a space only hazard rate.
 */
class ITO33_IHG_DLLDECL HazardRatePower : public HazardRate
{
public:

  /**
     constructor

     @param dAlpha alpha value
     @param dBeta beta
     @param dS0 s0
   */
  HazardRatePower(double dAlpha, double dBeta, double dS0);

  // Default dtor is ok
  
  /**
     Gets the value of alpha.

     @return alpha
   */
  double GetAlpha() const { return m_dAlpha; }

  /**
     Gets the value of beta.

     @return beta
   */
  double GetBeta() const { return m_dBeta; }

  /**
     Gets the value of S0.

     @return S0
   */
  double GetS0() const { return m_dS0; }  
  
  ///@name noexport functions  
  //@{

  bool IsTimeOnly() const
  {
    if (m_dBeta == 0)
      return true;
    else
      return false;
  }

  void GetHazardRates(const double *pdS, double *pdValues, 
                      size_t nNumber) const;

  void GetHazardRates(double dTime, const double *pdS, double *pdValues, 
                      size_t nNbS) const;

  virtual void Dump(XML::Tag& tagParent) const;  
  
  virtual 
    void GetModelParameters(finance::ModelParametersConsumer& visitor) const;

  virtual void Visit(HazardRateVisitor& visitor) const;
  
  //@}


private:

  double m_dAlpha;

  double m_dBeta;
  
  double m_dS0;

}; // class HazardRatePower


} // namespace ihg

} // namespace ito33


#endif // #ifndef _ITO33_IHG_HAZARDRATEPOWER_H_

