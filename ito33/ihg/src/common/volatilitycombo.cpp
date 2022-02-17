/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/volatiliycombo.cpp
// Purpose:     Implementations of volatilitycombo class
// Author:      David
// Created:     2004/06/02
// RCS-ID:      $Id: volatilitycombo.cpp,v 1.11 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/src/common/volatilitycombo.cpp
 */

#include "ito33/sharedptr.h"
#include "ito33/arraycheckers.h"

#include "ito33/numeric/piecewiseconstantfunctor.h"

#include "ito33/ihg/spotcomponent.h"
#include "ito33/ihg/volatilitycombo.h"
#include "ito33/ihg/volatility_visitor.h"

#include "ito33/xml/write.h"
#include "ito33/xml/write_vector.h"

#include "ihg/xml/spotcomponent.h"
#include "ihg/xml/volatility.h"

extern const ito33::Error ITO33_NULL_PARAM;

namespace ito33
{

namespace ihg
{


VolatilityCombo::VolatilityCombo(shared_ptr<SpotComponent> pSpotComponent)
  : VolatilityWithTimeComponent(),
    m_pSpotComponent(pSpotComponent)
{
}

VolatilityCombo::VolatilityCombo
                 (shared_ptr<SpotComponent> pSpotComponent,
                  const Date* pDates, const double* pdValues, 
                  size_t nNbTimes)
                : VolatilityWithTimeComponent(pDates, pdValues, nNbTimes),
                  m_pSpotComponent(pSpotComponent)
{
  CheckNonNegativity(pdValues, nNbTimes);

  if ( !m_pSpotComponent )
    throw EXCEPTION_MSG
          (
            ITO33_NULL_PARAM,
            TRANS("Setting invalid spot component for volatility")
          );
}

VolatilityCombo::VolatilityCombo
    (shared_ptr<SpotComponent> pSpotComponent,
     const std::vector<Date>& dates, const std::vector<double>& values)
   : VolatilityWithTimeComponent(dates, values),
     m_pSpotComponent(pSpotComponent)
{
  if ( !m_pSpotComponent )
    throw EXCEPTION_MSG
          (
            ITO33_NULL_PARAM,
            TRANS("Setting invalid spot component for volatility")
          );
}

void VolatilityCombo::GetVols
     (double dTime, const double *pdS, double *pdVols, size_t nNbS) const
{
  // Get the first spot component values
  m_pSpotComponent->GetValues(pdS, pdVols, nNbS);

  // Get the time compoent value
  double dTmp = (*m_pTimeComponent)(dTime);

  for (size_t nIdxS = 0; nIdxS < nNbS; nIdxS++)
    pdVols[nIdxS] *= dTmp;
}

void VolatilityCombo::Dump(ito33::XML::Tag& tagParent) const
{
  ito33::XML::Tag tag(XML_TAG_VOLATILITYCOMBO_ROOT, tagParent);

  // Split the remaining output into two separate tags 
  // Output the spot component first
  tag.Element(XML_TAG_SPOTCOMPONENT_ROOT, *m_pSpotComponent);

  // Output the time component second
  ito33::XML::Tag tagTime(XML_TAG_VOLATILITY_TIMECOMPONENT, tag); 

  DumpVector( tagTime, m_pDates,
              XML_TAG_FINANCE_DATES, XML_TAG_FINANCE_DATE);
  DumpVector( tagTime, m_pTimeComponent->GetY(),
              XML_TAG_FINANCE_VALUES, XML_TAG_FINANCE_VALUE);
}

void VolatilityCombo::Visit(VolatilityVisitor& visitor) const
{
  visitor.OnVolatilityCombo(*this);
}


} // namespace ihg

} // namespace ito33
