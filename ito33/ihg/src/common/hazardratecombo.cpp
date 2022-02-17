/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/hazardratecombo.cpp
// Purpose:     Combination of two hazard rate classes
// Author:      David 
// Created:     2005/06/02
// RCS-ID:      $Id: hazardratecombo.cpp,v 1.22 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/arraycheckers.h"

#include "ito33/numeric/piecewiseconstantfunctor.h"

#include "ito33/finance/modelparametersconsumer.h"

#include "ito33/ihg/spotcomponent.h"
#include "ito33/ihg/hazardratecombo.h"
#include "ito33/ihg/hazardrate_visitor.h"

#include "ito33/xml/write.h"
#include "ito33/xml/write_vector.h"
#include "ito33/xml/finance/common.h"

#include "ihg/xml/hazardrate.h"
#include "ihg/xml/spotcomponent.h"

namespace ito33
{

namespace ihg
{


HazardRateCombo::HazardRateCombo
                 (const shared_ptr<SpotComponent>& pSpotComponent,
                  const Date* pDates, const double* pdValues, 
                  size_t nNbTimes)
                : HazardRateWithTimeComponent(pDates, pdValues, nNbTimes),
                  m_pSpotComponent(pSpotComponent)
{
  CheckNonNegativity(pdValues, nNbTimes);
}

HazardRateCombo::HazardRateCombo
                 (const shared_ptr<SpotComponent>& pSpotComponent,
                  const std::vector<Date>& dates, 
                  const std::vector<double>& values)
                : HazardRateWithTimeComponent(dates, values),
                  m_pSpotComponent(pSpotComponent)
{
  CheckNonNegativity(values);
}

void HazardRateCombo::GetHazardRates(double dTime, const double* pdS,
                                     double *pdValues, size_t nNumber) const
{
  // Get the first spot component values
  m_pSpotComponent->GetValues(pdS, pdValues, nNumber);

  // Get the time compoent value
  double dTmp = (*m_pTimeComponent)(dTime);

  for (size_t nIdxS = 0; nIdxS < nNumber; nIdxS++)
    pdValues[nIdxS] *= dTmp;
}

void HazardRateCombo::Dump(ito33::XML::Tag& tagParent) const
{
  ito33::XML::Tag tag(XML_TAG_HAZARDRATECOMBO_ROOT, tagParent);

  // Split the remaining output into two separate tags 
  // Output the spot component first
  tag.Element(XML_TAG_SPOTCOMPONENT_ROOT, *m_pSpotComponent);

  // Output the time component second
  ito33::XML::Tag tagTime(XML_TAG_HAZARDRATE_TIMECOMPONENT, tag); 

  DumpVector( tagTime, m_pDates,
              XML_TAG_FINANCE_DATES, XML_TAG_FINANCE_DATE);
  DumpVector( tagTime, m_pTimeComponent->GetY(),
              XML_TAG_FINANCE_VALUES, XML_TAG_FINANCE_VALUE);

}

void HazardRateCombo::GetModelParameters
  (finance::ModelParametersConsumer& visitor) const
{
  // Actually, the only spot component we support is "power"
  m_pSpotComponent->GetModelParameters(visitor, MODEL_PARAM_NAME_HR_TIMESPOT);

  visitor.OnTimeComponentValues(MODEL_PARAM_NAME_HR_TIMESPOT, m_pDates,
    m_pTimeComponent->GetY() );
}


void HazardRateCombo::Visit(ito33::ihg::HazardRateVisitor &visitor) const
{
  visitor.OnHazardRateCombo(*this);
}

} // namespace ihg

} // namespace ito33
