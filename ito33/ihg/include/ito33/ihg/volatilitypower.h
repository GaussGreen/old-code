/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/volatilitypower.h
// Purpose:     implementation of  power volatility sigma(S)=alpha x (S0/S)^B
// Author:      Ito33
// Created:     2004/12/13
// RCS-ID:      $Id: volatilitypower.h,v 1.13 2006/07/20 16:51:29 zhang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/volatilitypower.h
 */

#ifndef _ITO33_IHG_VOLATILITYPOWER_H_
#define _ITO33_IHG_VOLATILITYPOWER_H_

#include "ito33/ihg/volatility.h"

namespace ito33
{

namespace ihg
{

 class VolatilityVisitor;

/**
    Power volatility class.
 */
class ITO33_IHG_DLLDECL VolatilityPower : public Volatility
{
public:

  /**
     Ctor constructs a power volatility

	   sigma(S) = alpha x (S0/S)^B

     @param dAlpha the coefficient (or vol value at dS0)
	   @param dBeta the power
	   @param dS0 stock price at which vol is dAlpha (typically the current spot)
   */
  VolatilityPower(double dAlpha, double dBeta, double dS0);

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


  void GetVols(double dTime,
               const double *pdS, double *pdVols, size_t nNbS) const;

  /**
     Gets the squares of the volatility at a given time for an array of spot

     It's more efficient than the default implementation in the base class.

     @param dTime the time to get the vol values at
     @param pdS Array of spot.
     @param pdVolsSquared Array of the squared values of (computed) volatilities
     @param nNbS The number of spots.
   */
  void GetVolsSquared(double dTime,
                      const double *pdS, 
                      double *pdVolsSquared, 
                      size_t nNbS) const;

  ///@name noexport functions  
  //@{
  void Dump(ito33::XML::Tag& tagParent) const;
  
  void Visit(VolatilityVisitor& visitor) const;


  bool IsTimeOnly() const
  {
    if (m_dBeta == 0)
      return true;
    else
      return false;
  }

  virtual
    void GetModelParameters(finance::ModelParametersConsumer& visitor) const;
  
  //@}

private:
 
  double m_dAlpha;

  double m_dBeta;

  double m_dS0;

}; // class VolatilityPower


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_VOLATILITYPOWER_H_

