/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/hazardratedecay.cpp
// Purpose:     exponentially decaying (space only) hazard rate class
// Author:      David (based on idea in calldef project)
// Created:     03/12/19
// RCS-ID:      $Id: hazardratedecay.cpp,v 1.16 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <math.h>

#include "ito33/ihg/hazardratedecay.h"
#include "ito33/ihg/hazardrate_visitor.h"

#include "ito33/xml/write.h"

#include "ihg/xml/hazardrate.h"

namespace ito33
{

namespace ihg
{

void HazardRateDecay::GetHazardRates(double /* dTime */, 
                                     const double*  pdSpots,
                                     double *pdValues, size_t nNumber) const
{
  size_t n;

  for( n = 0; n < nNumber; n++)
    pdValues[n] = m_dAlpha * exp( -pdSpots[n] / m_dS0);

}

void HazardRateDecay::Dump(ito33::XML::Tag& tagParent) const
{
  ito33::XML::Tag tag(XML_TAG_HAZARDRATEDECAY_ROOT, tagParent);
  /* The analytic formula is h(S) = alpha * exp(-S/S0) */
  /* ctor: HazardRateDecay(double dAlpha, double dS0) */
  tag.Element(XML_TAG_HAZARDRATE_ALPHA)(m_dAlpha);
  tag.Element(XML_TAG_HAZARDRATE_S0)(m_dS0);
}

void HazardRateDecay::Visit(ito33::ihg::HazardRateVisitor& visitor) const
{
  visitor.OnHazardRateDecay(*this);
}

} // namespace ihg


} // namespace ito33
