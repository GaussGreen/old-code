/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/volatilitytanh.h
// Purpose:     Volatility as a tanh function of spot
// Created:     2005/02/04
// RCS-ID:      $Id: volatilitytanh.h,v 1.13 2006/01/10 17:25:07 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/volatilitytanh.h
    @brief Volatility as a tanh function of spot.
 */

#ifndef _ITO33_IHG_VOLATILITYTANH_H_
#define _ITO33_IHG_VOLATILITYTANH_H_

#include "ito33/ihg/volatility.h"

namespace ito33
{

namespace ihg
{


  class VolatilityVisitor;

/// Volatility as a tanh function of spot.
class ITO33_IHG_DLLDECL VolatilityTanh : public Volatility
{
public:

  /**
     Ctor constructs a tanh volatility.

	   sigma(S) = (dRight - dLeft) * 0.5 * ( 1. + tanh (scale * (S - S0)) )
              + dLeft

     @param dLeft The limit of the volatility when S -> - infinite
	   @param dRight The limit of the volatility when S -> + infinite
     @param dScale The scale on the tanh function
	   @param dS0 The shift on the tanh function
   */
  VolatilityTanh(double dLeft, double dRight, double dScale, double dS0);

  // Default dtor is ok
  
  /**
     Gets the left limit of the volatility.

     @return Left limit of the volatility.
   */
  double GetLeft() const { return m_dLeft; }

  /**
     Gets the right limit of the volatility.

     @return Right limit of the volatility.
   */
  double GetRight() const { return m_dRight; }

  /**
     Gets the scale on the tanh function.

     @return Scale on the tanh function.
   */
  double GetScale() const { return m_dScale; }

  /**
     Gets the value of S0.

     @return S0
   */
  double GetS0() const { return m_dS0; }  

  void 
  GetVols(double dTime, const double *pdS, double *pdVols, size_t nNbS) const;

  ///@name noexport functions  
  //@{

  void Dump(ito33::XML::Tag& tagParent) const;
  
  void Visit(VolatilityVisitor& visitor) const;

  virtual 
    void GetModelParameters(finance::ModelParametersConsumer& visitor) const;
  
  //@}


private:
 
  double m_dLeft;

  double m_dRight;

  double m_dScale;

  double m_dS0;

}; // class VolatilityTanh


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_VOLATILITYTANH_H_
