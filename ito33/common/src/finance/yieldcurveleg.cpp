/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/yieldcurveleg.cpp
// Purpose:     YieldCurveLeg class implementation
// Author:      ZHANG Yunzhi
// Created:     2005/04/11
// RCS-ID:      $Id: yieldcurveleg.cpp,v 1.4 2006/08/23 10:23:05 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/timeunit.h"
#include "ito33/finance/yieldcurveleg.h"
#include "ito33/finance/cashrate.h"
#include "ito33/finance/swaprate.h"

#include "ito33/xml/write.h"

#include "ito33/xml/finance/timeunit.h"
#include "ito33/xml/finance/daycountconvention.h"
#include "ito33/xml/finance/frequency.h"
#include "ito33/xml/finance/yieldcurve.h"

extern const ito33::finance::Error
  ITO33_YIELDCURVELEG_RATE,
  ITO33_YIELDCURVELEG_MATURITY,
  ITO33_INVALID_TIMEUNIT;

namespace ito33
{

namespace finance
{

YieldCurveLeg::YieldCurveLeg(double dRate,
              size_t nMaturityDuration,
              TimeUnit maturityUnit
              )
              :m_dRate(dRate),
               m_nMaturityDuration(nMaturityDuration),
               m_maturityUnit(maturityUnit)
{
  CHECK_COND(nMaturityDuration > 0, ITO33_YIELDCURVELEG_MATURITY);

  CHECK_COND(dRate >= 0, ITO33_YIELDCURVELEG_RATE);

  CHECK_COND(maturityUnit < TimeUnit_Max, ITO33_INVALID_TIMEUNIT);
}

extern void DumpPart(const YieldCurveLeg& part, XML::Tag& tagObject)
{
  tagObject.Element(XML_TAG_SWAPCURVELEG_RATE)(part.GetRate());
  tagObject.Element(XML_TAG_SWAPCURVELEG_MATURITY_DURATION)
                   ( part.GetMaturityDuration() );
  tagObject.Element(XML_TAG_SWAPCURVELEG_MATURITY_UNIT)( 
                    XML::GetNameOfTimeUnit(part.GetMaturityUnit()));
}

/**
    Writes YieldCurveLeg object in tagParent and with tag name as
    legName. 
    
    @returns tag object.
 */
extern XML::Tag Dump(const YieldCurveLeg& leg, const char* legName,
                     XML::Tag& tagParent)
{
  ito33::XML::Tag tagCashRate(legName, tagParent);

  DumpPart(leg, tagCashRate);

  return tagCashRate;
}

} // namespace finance

} // namespace ito33
