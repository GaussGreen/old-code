/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/hazardratelinear.cpp
// Purpose:     linear (space only) hazard rate class
// Author:      David 
// Created:     03/12/19
// RCS-ID:      $Id: hazardratelinear.cpp,v 1.17 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <math.h>

#include "ito33/ihg/hazardratelinear.h"
#include "ito33/ihg/hazardrate_visitor.h"

#include "ito33/xml/write.h"

#include "ihg/xml/hazardrate.h"

namespace ito33
{

namespace ihg
{

void HazardRateLinear::GetHazardRates(double /* dTime */, 
                                      const double*  pdSpots,
                                      double *pdValues, 
                                      size_t nNumber) const
{
  size_t n;

  for( n = 0; n < nNumber; n++)
  {
    pdValues[n] = m_dSlope * pdSpots[n] + m_dB;

    if (pdValues[n] < 0.0) pdValues[n] = 0.0;
  }
}

void HazardRateLinear::Dump(ito33::XML::Tag& tagParent) const
{
  ito33::XML::Tag tag(XML_TAG_HAZARDRATELINEAR_ROOT, tagParent);
  /* ctor: HazardRateLinear(double dValue0, double dS0, double dValueS0) */
  /* m_dB = dValue0;                         */
  /*  m_dSlope = (dValueS0 - dValue0) / dS0; */
  tag.Element(XML_TAG_HAZARDRATE_B)(m_dB);
  tag.Element(XML_TAG_HAZARDRATE_SLOPE)(m_dSlope);
}

void HazardRateLinear::Visit(ito33::ihg::HazardRateVisitor &visitor) const
{
  visitor.OnHazardRateLinear(*this);
}

} // namespace ihg

} // namespace ito33
