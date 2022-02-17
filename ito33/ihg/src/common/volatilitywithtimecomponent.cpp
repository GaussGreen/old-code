/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/volatilitywithtimecomponent.cpp
// Purpose:     implementation for VolatilityWithTimeComponent class
// Created:     2005/03/04
// RCS-ID:      $Id: volatilitywithtimecomponent.cpp,v 1.8 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/dateutils.h"
#include "ito33/array.h"
#include "ito33/useexception.h"

#include "ito33/numeric/piecewiseconstantfunctor.h"
#include "ito33/numeric/mesh/specialtimes.h"

#include "ito33/ihg/volatilitywithtimecomponent.h"

extern const ito33::Error ITO33_BAD_PARAM, ITO33_BAD_DATA;

namespace ito33
{

namespace ihg
{


VolatilityWithTimeComponent::VolatilityWithTimeComponent()
{
  // Need to allocate, even if empty, in case other functions try to access
  // the dates or times
  m_pTimeComponent = make_ptr( new numeric::PiecewiseConstantFunctor
                                   ( NULL, NULL, 0 ) );
}

VolatilityWithTimeComponent::VolatilityWithTimeComponent
    (const Date* pDates, const double* pdValues, size_t nNbTimes) 
   : m_pDates(pDates, pDates + nNbTimes)
{
  Array<double> pdTimes(nNbTimes);
  size_t nIdx;

  for ( nIdx = 0; nIdx < nNbTimes; nIdx++)
  {
    CHECK_COND_MSG(pDates[nIdx].IsValid(), ITO33_BAD_DATA,
                   "invalid date array of the volatility.");
           
    pdTimes[nIdx] = GetDoubleFrom(pDates[nIdx]);

    if ( pdValues[nIdx] > 5 || pdValues[nIdx] < 0 )
    {
      throw EXCEPTION_MSG
            (
              ITO33_BAD_PARAM,
              TRANS("Time component volatility cannot be negative " 
                    "or greater than 500%!")
            );
    }
  }

  m_pTimeComponent = make_ptr( new numeric::PiecewiseConstantFunctor
                                   ( pdTimes.Get(), pdValues, nNbTimes ) );
}

VolatilityWithTimeComponent::VolatilityWithTimeComponent
(const std::vector<Date>& dates, const std::vector<double>& values):
  m_pDates(dates)
{
  Array<double> pdTimes(dates.size());
  size_t nIdx;

  for ( nIdx = 0; nIdx < dates.size(); nIdx++)
  {
        
    CHECK_COND_MSG(dates[nIdx].IsValid(), ITO33_BAD_DATA,
                   "invalid date array of the volatility.");

    pdTimes[nIdx] = GetDoubleFrom(dates[nIdx]);

    if ( values[nIdx] > 5 || values[nIdx] < 0 )
    {
      throw EXCEPTION_MSG
            (
              ITO33_BAD_PARAM,
              TRANS("Time component volatility cannot be negative " 
                    "or greater than 500%!")
            );
    }

  }

  m_pTimeComponent = make_ptr( new numeric::PiecewiseConstantFunctor
                               ( pdTimes.Get(), &values[0], values.size() ) );
}

const std::vector<double>&
VolatilityWithTimeComponent::GetTimeComponentValues() const
{
  return m_pTimeComponent->GetY();
}

void VolatilityWithTimeComponent::ResetTimeComponent 
     (const Date* pDates, const double *pdValues, size_t nNbTimes)
{
  m_pDates = std::vector<Date>(pDates, pDates + nNbTimes);

  Array<double> pdTimes(nNbTimes);
  size_t nIdx;

  for ( nIdx = 0; nIdx < nNbTimes; nIdx++)
    pdTimes[nIdx] = GetDoubleFrom(pDates[nIdx]);

  m_pTimeComponent = make_ptr( new numeric::PiecewiseConstantFunctor
                                   ( pdTimes.Get(), pdValues, nNbTimes ) );
}

void VolatilityWithTimeComponent::GetSpecialTimes
     (numeric::mesh::SpecialTimes& specialTimes) const
{
  specialTimes.clear();

  const std::vector<double>& pdTimes = m_pTimeComponent->GetX();

  // use standard refine level for these time points
  size_t nIdx;
  for ( nIdx = 0; nIdx < pdTimes.size(); nIdx++)
    specialTimes.push_back( numeric::mesh::SpecialTime(pdTimes[nIdx]) );                    
}
 
} // namespace ihg

} // namespace ito33
