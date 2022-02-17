/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/bondlike/conversionpricereset.cpp
// Purpose:     Characteristics of the conversion price reset date
// Author:      Yann and David
// Created:     2004/10/20
// RCS-ID:      $Id: conversionpricereset.cpp,v 1.7 2006/08/19 22:43:51 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"
#include "ito33/date.h"

#include "ito33/finance/bondlike/conversionpricereset.h"
#include "ito33/finance/bondlike/bonderror.h"
#include "ito33/finance/bondlike/resetflooredby.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/bondlike/conversionpricereset.h"
#include "ito33/xml/finance/bondlike/common.h"
#include "ito33/xml/finance/common.h"

extern const ito33::finance::BondError
    ITO33_BONDLIKE_INVALID_RESET_CAP_RATE,
    ITO33_BONDLIKE_INVALID_RESET_FLOOR_RATE,
    ITO33_BONDLIKE_INVALID_RESET_MULTIPLIER;

namespace ito33
{

namespace finance
{

ConversionPriceReset::ConversionPriceReset(Date resetDate, double dFloor)
   : m_resetDate(resetDate), m_dFloor(dFloor),
     m_dCap(1.0), m_dMultiplier(1.)
{
  CHECK_COND
  (
    m_dFloor <= 1.0 && m_dFloor >= 0.0,
    ITO33_BONDLIKE_INVALID_RESET_FLOOR_RATE
  );
} 

void ConversionPriceReset::SetMultiplier(double dMultiplier)
{
  CHECK_COND
  (
    dMultiplier >= 0.0 && dMultiplier <= 10.0,
    ITO33_BONDLIKE_INVALID_RESET_MULTIPLIER
  );

  m_dMultiplier = dMultiplier;   
}

void ConversionPriceReset::SetCap(double dCap)
{
  CHECK_COND
  (
    dCap >= 1.0 && dCap <= 10.0,
    ITO33_BONDLIKE_INVALID_RESET_CAP_RATE
  );

  m_dCap = dCap;
}


XML::Tag ConversionPriceReset::Dump(ito33::XML::Tag& tagParent) const
{
  XML::Tag tagMe(XML_TAG_CONVPRICERESET_ROOT, tagParent);
  
  tagMe.Element(XML_TAG_FINANCE_DATE)(m_resetDate);
  tagMe.Element(XML_TAG_CONVPRICERESET_CAPRATE)(m_dCap);
  tagMe.Element(XML_TAG_CONVPRICERESET_FLOORRATE)(m_dFloor);
  tagMe.Element(XML_TAG_CONVPRICERESET_MULTIPLIER)(m_dMultiplier);

  return tagMe;
}


} // namespace finance

} // namespace ito33

