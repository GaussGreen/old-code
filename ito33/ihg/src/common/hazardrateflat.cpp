/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/hazardrateflat.cpp
// Purpose:     Flat (constant) hazard rate class
// Author:      David 
// Created:     2005/05/05
// RCS-ID:      $Id: hazardrateflat.cpp,v 1.19 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/ihg/hazardrateflat.h"
#include "ito33/ihg/hazardrate_visitor.h"

#include "ito33/xml/write.h"

#include "ihg/xml/hazardrate.h"

extern const ito33::Error ITO33_BAD_PARAM;

namespace ito33
{

namespace ihg
{


HazardRateFlat::HazardRateFlat(double dValue) 
{
  if ( dValue < 0 || dValue > 10 )
  {
    throw EXCEPTION_MSG
          (
            ITO33_BAD_PARAM,
            TRANS("Flat hazard rate cannot be negative or greater than 1000%!")
          );
  }

  m_dValue = dValue;
}

double HazardRateFlat::GetValueAtTime(double /* dTime */) const
{
  return m_dValue;
}

shared_ptr<HazardRate> HazardRateFlat::Perturb(double dShift)
{
  return shared_ptr<HazardRate>(new HazardRateFlat(m_dValue + dShift));
}

void HazardRateFlat::GetHazardRates(double /* dTime */, 
                                    const double* /* pdS */,
                                    double *pdValues, 
                                    size_t nNumber) const
{
  size_t nIdx;

  for ( nIdx = 0; nIdx < nNumber; nIdx++)
    pdValues[nIdx] = m_dValue;
}

void HazardRateFlat::Dump(ito33::XML::Tag& tagParent) const
{
  ito33::XML::Tag tag(XML_TAG_HAZARDRATEFLAT_ROOT, tagParent);
  /* ctor: HazardRateFlat(double dValue) */
  tag.Element(XML_TAG_HAZARDRATE_FLAT)(m_dValue);
}

void HazardRateFlat::Visit(ito33::ihg::HazardRateVisitor& visitor) const
{
  visitor.OnHazardRateFlat(*this);
}


} // namespace ihg

} // namespace ito33
