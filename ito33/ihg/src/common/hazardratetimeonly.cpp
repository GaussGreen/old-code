/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/hazardratetimeonly.cpp
// Purpose:     implementation of HazardRateTimeOnly class
// Author:      (z)
// Created:     03/11/14
// RCS-ID:      $Id: hazardratetimeonly.cpp,v 1.31 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/arraycheckers.h"
#include "ito33/sharedptr.h"

#include "ito33/numeric/piecewiseconstantfunctor.h"

#include "ito33/finance/modelparametersconsumer.h"

#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardrate_visitor.h"

#include "ito33/xml/write.h"
#include "ito33/xml/write_vector.h"
#include "ito33/xml/finance/common.h"

#include "ihg/xml/hazardrate.h"

namespace ito33
{

namespace ihg
{

  
HazardRateTimeOnly::HazardRateTimeOnly(const Date *pDates,
                                       const double *pdValues,
                                       size_t nNbTimes)
  : HazardRateWithTimeComponent
    ( 
      pDates, 
      CheckNonNegativity(pdValues, nNbTimes),
      nNbTimes
    )
{
}

HazardRateTimeOnly::HazardRateTimeOnly(const std::vector<Date>& dates,
                                       const std::vector<double>& values)
  : HazardRateWithTimeComponent( dates, CheckNonNegativity(values) )
{ 
}


shared_ptr<HazardRate> HazardRateTimeOnly::Perturb(double dShift)
{
  std::vector<double> pdNewV( GetTimeComponentValues() );

  for(size_t n = 0; n < pdNewV.size(); n++)
    pdNewV[n] += dShift;

  return shared_ptr<HazardRate>(new HazardRateTimeOnly(m_pDates, pdNewV) );
}

double HazardRateTimeOnly::GetValueAtTime(double dTime) const
{
  return (*m_pTimeComponent)(dTime);
}

void HazardRateTimeOnly::GetHazardRates(double dTime, 
                                        const double *  /* pdSpots */,
                                        double *pdValues, 
                                        size_t nNbS) const
{
  double dAlpha = GetValueAtTime(dTime);
  
  for (size_t n = 0; n < nNbS; n++)
    pdValues[n] = dAlpha;
}



void HazardRateTimeOnly::Dump(ito33::XML::Tag& tagParent) const
{
  ito33::XML::Tag tag(XML_TAG_HAZARDRATETIMEONLY_ROOT, tagParent);
 
  DumpVector( tag, m_pDates,
              XML_TAG_FINANCE_DATES, XML_TAG_FINANCE_DATE);
  DumpVector( tag, m_pTimeComponent->GetY(),
              XML_TAG_FINANCE_VALUES, XML_TAG_FINANCE_VALUE);


} //end Dump

void HazardRateTimeOnly::GetModelParameters
  (finance::ModelParametersConsumer& visitor) const
{
  visitor.OnTimeComponentValues(MODEL_PARAM_NAME_HR_TIMEONLY , m_pDates,
    m_pTimeComponent->GetY() );
}

void HazardRateTimeOnly::Visit(ito33::ihg::HazardRateVisitor& visitor) const
{
  visitor.OnHazardRateTimeOnly(*this);
}



} // namespace ihg

} // namespace ito33
