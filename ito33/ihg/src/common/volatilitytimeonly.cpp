/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/volatilitytimeonly.cpp
// Purpose:     implementation of HazardRateTimeOnly class
// Author:      (z), David
// Created:     2004/05/19
// RCS-ID:      $Id: volatilitytimeonly.cpp,v 1.18 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 - 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <vector>
#include "ito33/afterstd.h"

#include "ito33/arraycheckers.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/modelparametersconsumer.h"

#include "ito33/numeric/piecewiseconstantfunctor.h"

#include "ito33/ihg/volatilitytimeonly.h"
#include "ito33/ihg/volatility_visitor.h"

#include "ito33/xml/write.h"
#include "ito33/xml/write_vector.h"

#include "ihg/xml/volatility.h"

namespace ito33
{

namespace ihg
{


VolatilityTimeOnly::VolatilityTimeOnly
   (const Date *pDates, const double *pdValues, size_t nNbTimes)
  : VolatilityWithTimeComponent(pDates, pdValues, nNbTimes)
{
}

VolatilityTimeOnly::VolatilityTimeOnly
   (const std::vector<Date>& dates, const std::vector<double>& values)
  : VolatilityWithTimeComponent( dates, values )
{ 
}

void VolatilityTimeOnly::GetVols(double dTime, const double* /* pdS */,
                                 double *pdVols, size_t nNbS) const
{
  double dVol = (*m_pTimeComponent)(dTime);
  
  for (size_t n = 0; n < nNbS; n++)
    pdVols[n] = dVol;
}

shared_ptr<Volatility> VolatilityTimeOnly::Perturb(double dShift) const
{
  std::vector<double> pdNewV( GetTimeComponentValues() );

  for(size_t n = 0; n < pdNewV.size(); n++)
    pdNewV[n] += dShift;

  return shared_ptr<Volatility>(
            new VolatilityTimeOnly(m_pDates, pdNewV));
}

void VolatilityTimeOnly::Dump(ito33::XML::Tag& tagParent) const
{
  ito33::XML::Tag tag(XML_TAG_VOLATILITYTIMEONLY_ROOT, tagParent);
 
  DumpVector( tag, m_pDates,
              XML_TAG_FINANCE_DATES, XML_TAG_FINANCE_DATE);
  DumpVector( tag, m_pTimeComponent->GetY(),
              XML_TAG_FINANCE_VALUES, XML_TAG_FINANCE_VALUE);

} //end Dump

void VolatilityTimeOnly::Visit(VolatilityVisitor& visitor) const
{
  visitor.OnVolatilityTimeOnly(*this);
}

void 
VolatilityTimeOnly::GetModelParameters(finance::ModelParametersConsumer& visitor) const
{
  visitor.OnTimeComponentValues(MODEL_PARAM_NAME_VOL_TIMEONLY,
    m_pDates, m_pTimeComponent->GetY() );
}

} // namespace ihg

} // namespace ito33
