/////////////////////////////////////////////////////////////////////////////
// Name:        src/finance/yieldcurve_flat.cpp
// Purpose:     YieldCurveFlat class implementation
// Author:      Vadim Zeitlin, ZHANG Yunzhi
// Created:     12.02.04
// RCS-ID:      $Id: yieldcurve_flat.cpp,v 1.16 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#if ( defined(_MSC_VER) && _MSC_VER >= 1400 )
// still want to use fill_n but don't know how to avoid the compile warning
// warning C4996: 'std::_Fill_n' was declared deprecated
// so just disable the warning
#pragma warning(disable:4996)
#endif

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/debug.h"
#include "ito33/date.h"

#include "ito33/beforestd.h"
#include <algorithm> // for fill_n
#include "ito33/afterstd.h"

#include "ito33/useexception.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/error.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/yieldcurve_visitor.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/yieldcurve.h"
#include "ito33/xml/finance/common.h"

using namespace ito33;
using namespace ito33::finance;

extern const ito33::finance::Error ITO33_NEG_FLATRATE;

// ============================================================================
// YieldCurveFlat implementation
// ============================================================================

YieldCurveFlat::YieldCurveFlat(double dRate)
: YieldCurve(Date(2000, Date::Jan, 1)), m_dRate(dRate) 
  {
    if (m_dRate < 0.)
      throw EXCEPTION_MSG
            (
              ITO33_NEG_FLATRATE,
              TRANS("YieldCurveFlat: Rate must be non negatif.")
            );
  }

void
YieldCurveFlat::GetZeroRates(const double *pdDates,
                             double *pdValues,
                             size_t nPoints) const
{
  CHECK_VOID( pdDates && pdValues,
              "YieldCurveFlat::GetZeroRates(): NULL parameter" );

  std::fill_n(pdValues, nPoints, m_dRate);
}

shared_ptr<YieldCurve>
YieldCurveFlat::Perturb(double dShift) const
{
  return make_ptr( new YieldCurveFlat(m_dRate + dShift) );
}

XML::Tag
YieldCurveFlat::Dump(const char *name, XML::Tag& tagParent) const
{

  XML::Tag tagYC(name, tagParent);
  XML::Tag tagYCVisitor(XML_TAG_YIELDCURVEFLAT_ROOT,tagYC);
  tagYCVisitor.Element(XML_TAG_FINANCE_FLAT)(m_dRate);

  return tagYC;
}


void YieldCurveFlat::Visit(YieldCurveVisitor& visitor) const
{
  visitor.OnYieldCurveFlat(*this);
}

