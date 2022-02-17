/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/yieldcurve_annuallycompounded.cpp
// Purpose:     YieldCurveAnnuallyCompounded class implementation
// Created:     2004/06/17
// RCS-ID:      $Id: yieldcurve_annuallycompounded.cpp,v 1.26 2006/08/21 14:28:03 zhang Exp $
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

#include "ito33/useexception.h"
#include "ito33/sharedptr.h"
#include "ito33/constants.h"

#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/extrapolationmode.h"

#include "ito33/finance/error.h"
#include "ito33/finance/yieldcurve_annuallycompounded.h"
#include "ito33/finance/yieldcurve_visitor.h"

#include "ito33/xml/write.h"

#include "ito33/xml/finance/yieldcurve.h"
#include "ito33/xml/finance/common.h"

extern const ito33::finance::Error
  ITO33_YIELDCUVE_ANNUALLYCCOMPOUNDED,
  ITO33_YIELDCURVE_INVALID_RATE,
  ITO33_YIELDCURVE_SETLEGS_BADSIZE;

using namespace ito33;
using namespace ito33::finance;

// ============================================================================
// YieldCurveAnnuallyCompounded implementation
// ============================================================================
// ============================================================================
// inline functions implementation
// ============================================================================


// ----------------------------------------------------------------------------
// private functions
// ----------------------------------------------------------------------------

// this function is used for equality of two dates expressed as doubles (we
// can't compare them using ==, of course...)
static inline bool
AreDaysEqual(const ZeroRate& leg1,
             const ZeroRate& leg2)
{
  return leg1.GetDuration() == leg2.GetDuration();
}

// ----------------------------------------------------------------------------
// YieldCurveAnnuallyCompounded operations
// ----------------------------------------------------------------------------

bool YieldCurveAnnuallyCompounded::CheckData()
{
  // we should have at least 2 legs
  if ( m_legs.size() < 2 )
    return false;

  // sort the legs
  const YieldCurveAnnuallyCompounded::Legs::iterator 
    pBegin = m_legs.begin(),
    pEnd = m_legs.end();

  std::sort(pBegin, pEnd);

  // now we can check for duplicates
  if ( std::adjacent_find(pBegin, pEnd, AreDaysEqual) != pEnd )
    return false;

  // note that we don't test any more for negavitve values here because
  // Leg does it in its ctor

  return true;
}

void YieldCurveAnnuallyCompounded::DoValidate()
{
  CHECK_COND( CheckData(), ITO33_YIELDCUVE_ANNUALLYCCOMPOUNDED );

  // prepare the data in a format such that GetZeroRates() can work efficiently
  // with them
  size_t nLegs = m_legs.size();

  m_pdDays.resize(nLegs);
  m_pdRates.resize(nLegs);

  for ( size_t n = 0; n < m_legs.size(); ++n )
  {
    m_pdDays[n] = m_legs[n].GetDuration() * ONEDAY;
    m_pdRates[n] = m_legs[n].GetRate();
  }
}

// ----------------------------------------------------------------------------
// YieldCurveAnnuallyCompounded accessors
// ----------------------------------------------------------------------------

void
YieldCurveAnnuallyCompounded::GetZeroRates(const double *pdDays,
                                           double *pdValues,
                                           size_t nPoints) const
{
  Validate();

  ASSERT_MSG
  (
    m_pdDays[0] >= 0., 
    "The maturity must be greater than the reference date of the yield curve"
  ); 

  numeric::Interpolate
  (
    &(m_pdDays[0]), &(m_pdRates[0]), m_pdDays.size(), 
    pdDays, pdValues, nPoints,
    numeric::ExtrapolationMode_Constant, numeric::ExtrapolationMode_Constant
  );
}

void YieldCurveAnnuallyCompounded::SetLegs
            (
              const std::vector<size_t>& pnDays,
              const std::vector<double>& pdRates
            )
{
  if ( pnDays.size() != pdRates.size() )
    ThrowInvalidSetLegsData();

  Clear();

  m_legs.resize(pnDays.size());
  for(size_t n = 0; n < pnDays.size(); n++)
  {
    if(pdRates[n] < 0)
      ThrowInvalidRate();

    m_legs[n] = ZeroRate(pdRates[n], pnDays[n]);
  }

  // invalidate first before Validate()
  Invalidate();

  Validate();
}

shared_ptr<YieldCurve>
YieldCurveAnnuallyCompounded::Perturb(double dShift) const
{
  Validate();

  shared_ptr<YieldCurveAnnuallyCompounded>
    pyc(new YieldCurveAnnuallyCompounded(GetReferenceDate(), m_pdDays.size()));

  for ( Legs::const_iterator iter = m_legs.begin();
        iter != m_legs.end();
        ++iter)
    pyc->AddLeg(iter->GetDuration(), iter->GetRate() + dShift);

  return pyc;
}

XML::Tag
YieldCurveAnnuallyCompounded::Dump(const char *name, XML::Tag& tagParent) const
{
  XML::Tag tagYC(name, tagParent);
  XML::Tag tagYCVisitor(XML_TAG_YIELDCURVEANNUALLYCOMPOUNDED_ROOT,tagYC);

  tagYCVisitor.Element(XML_TAG_YIELDCURVE_REFERENCE_DATE)(GetReferenceDate());

  for ( size_t n = 0; n < m_legs.size(); ++n )
  {
    XML::Tag tagLeg(XML_TAG_YIELDCURVE_LEG, tagYCVisitor);

    tagLeg.Element(XML_TAG_YIELDCURVE_LEG_DAY)
                  (m_legs[n].GetDuration());

    tagLeg.Element(XML_TAG_FINANCE_VALUE)(m_legs[n].GetRate());
 
  }

  return tagYC;
}

void YieldCurveAnnuallyCompounded::Visit(YieldCurveVisitor& visitor) const
{
  visitor.OnYieldCurveAnnuallyCompounded(*this);
}

////////////////////////// exception thrown /////////////////////////////////////
void YieldCurveAnnuallyCompounded::ThrowInvalidRate()
{
  throw EXCEPTION(ITO33_YIELDCURVE_INVALID_RATE);
}

void YieldCurveAnnuallyCompounded::ThrowInvalidSetLegsData()
{
  throw EXCEPTION(ITO33_YIELDCURVE_SETLEGS_BADSIZE);
}

