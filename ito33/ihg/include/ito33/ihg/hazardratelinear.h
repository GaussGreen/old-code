/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/hazardratelinear.h
// Purpose:     linear (space only) hazard rate class
// Author:      David
// Created:     03/12/19
// RCS-ID:      $Id: hazardratelinear.h,v 1.23 2005/08/23 12:00:15 zhang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
  @file ito33/ihg/hazardratelinear.h
  @brief linear (space only) hazard rate class
  */


#ifndef _ITO33_IHG_HAZARDRATELINEAR_H_
#define _ITO33_IHG_HAZARDRATELINEAR_H_

#include "ito33/useexception.h"

#include "ito33/ihg/hazardrate.h"

extern const ito33::Error ITO33_BAD_PARAM, ITO33_BAD_DATA;

namespace ito33
{

namespace ihg
{

class HazardRateVisitor;

/**
    @internal 
    @brief Linear (space only) hazard rate class 

    The user passes in the value of the hazard rate at S=0, as well
    as the value at one other point.  These two points define the 
    line.  In areas where the value is negative, the hazard rate
    is set to zero (as well as the derivative).

    This file is meant for testing, and does not neccessarily reflect a 
    market hazard rate function.

    @noexport
 */
class HazardRateLinear : public HazardRate
{
public:

  /**
     constructor by the value at 0 and the value at a given spot S0 

     @param dValue0 the value of the hazard rate at S=0
     @param dS0 spot different from 0
     @param dValueS0 the value of the hazard rate at S=S0
   */
  HazardRateLinear(double dValue0, double dS0, double dValueS0)
  { 
    if(dS0 == 0)
      throw EXCEPTION_MSG( ITO33_BAD_PARAM,
                  TRANS("Failed to define linear hazard rate: s0 is zero."));
    // Convert to y(x) = slope * x + b form
    m_dB = dValue0;
    m_dSlope = (dValueS0 - dValue0) / dS0;

    CheckSlope();
  }

  /* Override constructor since only
      dB and dSlope are available for the 
      the dump to XML */
  /**
    Trivial constructor

    @param dB b in y(x) = slope * x + b form
    @param dSlope slope in y(x) = slope * x + b form
    */
  HazardRateLinear(double dB,double dSlope)
  {
   m_dB     = dB;
   m_dSlope = dSlope;
   
    CheckSlope();
  }

  /// dtor
  virtual ~HazardRateLinear() { }
 
  virtual bool IsTimeOnly() const
  {
    if(m_dSlope == 0)
      return true;
    else
      return false;
  }  

  void GetHazardRates(double dTime, const double *pdS, double *pdValues, 
                      size_t nNumber) const;
  
  virtual void Dump(XML::Tag& tagParent) const;

  virtual void Visit(HazardRateVisitor& visitor) const;


protected:

  /**
    check if m_dSlope is less than 0
    */
  void CheckSlope()
  {
    if(m_dSlope > 0)
      throw EXCEPTION_MSG( ITO33_BAD_DATA,
            TRANS("Bad hazard rate: increasing function."));
  }
  
  /// b in y(x) = slope * x + b form
  double m_dB;

  /// slope in y(x) = slope * x + b form
  double m_dSlope;

}; // class HazardRateLinear


} // namespace ihg

} // namespace ito33


#endif // #ifndef _ITO33_IHG_HAZARDRATELINEAR_H_

