/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/yieldcurve_flat.h
// Purpose:     Flat Zero Rate Yield Curve class
// Author:      Vadim Zeitlin, ZHANG Yunzhi
// Created:     12.02.04
// RCS-ID:      $Id: yieldcurve_flat.h,v 1.20 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/yieldcurve_flat.h
    @brief The declaration of YieldCurveFlat class.

    A class for the simplest possible yield curve with constant flat rate.
 */

#ifndef _ITO33_FINANCE_YIELDCURVE_FLAT_H_
#define _ITO33_FINANCE_YIELDCURVE_FLAT_H_

#include "ito33/sharedptr.h"
#include "ito33/finance/yieldcurve.h" 

namespace ito33
{

namespace finance
{


/// The simplest possible yield curve: constant flat rate.
class ITO33_DLLDECL YieldCurveFlat : public YieldCurve
{
public:
  /**
      Constructor initializes the new yield curve with the constant rate.
  
      @param dRate the constant zero coupon rate value
   */
  YieldCurveFlat(double dRate);

  // default copy ctor, assignment operator and dtor are ok

  /// Returns the constant rate of this yield curve.
  double GetRate() const { return m_dRate; }

  // implement the base class pure virtuals
  virtual void GetZeroRates(const double *pdDays,
                            double *pdValues,
                            size_t nPoints) const;

  virtual shared_ptr<YieldCurve> Perturb(double dShift) const;

    // implement base class pure virtuals
  virtual void Visit(YieldCurveVisitor& visitor) const;

  virtual XML::Tag Dump(const char *name, XML::Tag& tagParent) const;

private:
  double m_dRate;
};

} // namespace finance

} // namespace ito33

#endif // _ITO33_FINANCE_YIELDCURVE_FLAT_H_
