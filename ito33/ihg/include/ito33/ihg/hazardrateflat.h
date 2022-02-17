/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/hazardrateflat.h
// Purpose:     flat (constant) hazard rate class
// Author:      David
// Created:     2004/05/05
// RCS-ID:      $Id: hazardrateflat.h,v 1.24 2006/08/20 09:47:41 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/ihg/hazardrateflat.h
   @brief flat (constant) hazard rate class
 */


#ifndef _ITO33_IHG_HAZARDRATEFLAT_H_
#define _ITO33_IHG_HAZARDRATEFLAT_H_

#include "ito33/ihg/hazardratetimeonly.h"

namespace ito33
{

namespace ihg
{

class HazardRateVisitor;

/** 
   Flat (constant) hazard rate class.

   Serves the same purpose as a flat volatility.
 */
class HazardRateFlat : public HazardRateTimeOnly
{

public:

  /**
     Ctor constructs a flat hazard rate from a given value.

     Throws an exception if flat value out of accpetable bound.
     
     @param dValue the value of the flat hazard rate
   */
  HazardRateFlat(double dValue);

  // Default dtor is ok
   
  /**
      Gets the value of this hazard rate.

      @return the constant value of this hazard rate
   */  
  double GetValue() const { return m_dValue; }

  virtual double GetValueAtTime(double dTime) const;

  virtual shared_ptr<HazardRate> Perturb(double dShift);

  void GetHazardRates(double dTime, const double *pdS, double *pdValues, 
                      size_t nNumber) const;

  virtual void Dump(ito33::XML::Tag& tagParent) const;

  virtual void Visit(HazardRateVisitor& visitor) const;

private:

  double m_dValue;

}; // class HazardRateFlat


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_HAZARDRATEFLAT_H_
