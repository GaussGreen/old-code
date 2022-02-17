/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/yieldcurve_zerocoupon.cpp
// Purpose:     YieldCurveZeroCoupon class implementation
// Created:     2004/06/17
// RCS-ID:      $Id: yieldcurve_zerocoupon.cpp,v 1.18 2006/08/23 10:23:05 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include <cmath>

#include "ito33/beforestd.h"
#include <algorithm>
#include "ito33/afterstd.h"

#include "ito33/date.h"
#include "ito33/useexception.h"
#include "ito33/sharedptr.h"
#include "ito33/constants.h"

#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/extrapolationmode.h"

#include "ito33/finance/error.h"
#include "ito33/finance/yieldcurve_zerocoupon.h"
#include "ito33/finance/yieldcurve_annuallycompounded.h"
#include "ito33/finance/yieldcurve_visitor.h"

#include "ito33/xml/write.h"

#include "ito33/xml/finance/yieldcurve.h"
#include "ito33/xml/finance/common.h"

extern const ito33::finance::Error
  ITO33_YIELDCUVE_ANNUALLYCCOMPOUNDED,
  ITO33_YIELDCURVE_INVALID_RATE,
  ITO33_YIELDCURVE_SETLEGS_BADSIZE;

namespace ito33
{

namespace finance
{

extern XML::Tag 
Dump(const YieldCurveLeg& leg, const char* legName, XML::Tag& tagParent);

void YieldCurveZeroCoupon::DoValidate()
{
  CheckReferenceDate();

  m_pImpl = shared_ptr<YieldCurveAnnuallyCompounded>
            (new YieldCurveAnnuallyCompounded(GetReferenceDate(),
                                              m_legs.size()));

  Date reference(GetReferenceDate());

  for ( Legs::const_iterator iter = m_legs.begin();
        iter != m_legs.end();
        ++iter)
  {
    // get the day diff
    Date dt(reference);
    switch ( iter->GetMaturityUnit() )
    {
      case TimeUnit_Day:
        dt.AddDays(iter->GetMaturityDuration());
        break;

      case TimeUnit_Month:
        dt.AddMonths(iter->GetMaturityDuration());
        break;

      case TimeUnit_Year:
        dt.AddYears(iter->GetMaturityDuration());
        break;

      case TimeUnit_Max:
        FAIL("invalid maturity value");
        break;
    }
    
    size_t nDays = Date::DaysDiff(reference, dt);

    double dDayDiffFrom = Date::DaysDiffWithDayCount(reference, dt, m_dcc);
    double dDayDiffTo = nDays * ONEDAY; // convention of yieldCurveAnnuallyCompounded

    // ( 1 + R )^T = ( 1 + r )^t in both day count convention
    double dRate = pow( 1. + iter->GetRate(), dDayDiffTo / dDayDiffFrom ) - 1.;

    //  convert rate to that expressed on ACT/365.
    m_pImpl->AddLeg(nDays, dRate);
  }

  m_pImpl->Validate();
}

void YieldCurveZeroCoupon::SetDayCountConvention(Date::DayCountConvention dcc)
{
  ito33::Validate(dcc);

  m_dcc = dcc;
}

void YieldCurveZeroCoupon::AddLeg(size_t nMaturityDuration,
                                  TimeUnit maturityUnit,
                                  double dRate)
{
  m_legs.push_back(YieldCurveLeg(dRate, nMaturityDuration, maturityUnit));
}

// ----------------------------------------------------------------------------
// YieldCurveZeroCoupon accessors
// ----------------------------------------------------------------------------

void
YieldCurveZeroCoupon::GetZeroRates(const double *pdDays,
                                   double *pdValues,
                                   size_t nPoints) const
{
  Validate();

  return m_pImpl->GetZeroRates(pdDays, pdValues, nPoints);
}

shared_ptr<YieldCurve>
YieldCurveZeroCoupon::Perturb(double dShift) const
{
  Validate();

  shared_ptr<YieldCurveZeroCoupon>
    pyc(new YieldCurveZeroCoupon(GetReferenceDate(), m_legs.size()));

  pyc->SetDayCountConvention(GetDayCountConvention());

  for ( Legs::const_iterator iter = m_legs.begin();
        iter != m_legs.end();
        ++iter)
    pyc->AddLeg(iter->GetMaturityDuration(),
                iter->GetMaturityUnit(),
                iter->GetRate() + dShift);

  return pyc;
}

XML::Tag
YieldCurveZeroCoupon::Dump(const char *name, XML::Tag& tagParent) const
{
  XML::Tag tagYC(name, tagParent);
  XML::Tag tagYCVisitor(XML_TAG_YIELDCURVEZEROCOUPON_ROOT,tagYC);

  tagYCVisitor.Element(XML_TAG_YIELDCURVE_REFERENCE_DATE)(GetReferenceDate());

  for ( size_t n = 0; n < m_legs.size(); ++n )
  {
    finance::Dump(m_legs[n], XML_TAG_ZEROCOUPONRATE_ROOT, tagYCVisitor); 
  }

  return tagYC;
}

void YieldCurveZeroCoupon::Visit(YieldCurveVisitor& visitor) const
{
  visitor.OnYieldCurveZeroCoupon(*this);
}

} // namespace finance

} // namespace ito33
