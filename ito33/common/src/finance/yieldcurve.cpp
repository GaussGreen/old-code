/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/yieldcurve.cpp
// Purpose:     YieldCurve classes implementation
// Author:      Vadim Zeitlin
// Created:     12.02.04
// RCS-ID:      $Id: yieldcurve.cpp,v 1.25 2006/08/21 14:28:03 zhang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include <cmath>

#include "ito33/sharedptr.h"
#include "ito33/dateutils.h"
#include "ito33/array.h"
#include "ito33/useexception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/yieldcurve.h"

extern const ito33::finance::Error ITO33_YIELDCURVE_REFERENCEDATE_UNDEF;
extern const ito33::Error ITO33_BAD_DATE;

namespace ito33
{

namespace finance
{


YieldCurve::YieldCurve(const Date& referenceDate)
  : m_ReferenceDate(referenceDate), 
    m_dReferenceDate( GetDoubleFrom(m_ReferenceDate) ),
    m_bValidated(false)
{
  CHECK_COND(referenceDate.IsValid(), ITO33_BAD_DATE);
}


YieldCurve::YieldCurve()
  : m_dReferenceDate( 1.e99 ),
    m_bValidated(false)
{
}


const Date& YieldCurve::GetReferenceDate() const
{
  CheckReferenceDate();

  return m_ReferenceDate;
}

void YieldCurve::SetReferenceDate(const Date& referenceDate)
{
  CHECK_COND(referenceDate.IsValid(), ITO33_BAD_DATE);

  m_ReferenceDate = referenceDate;
  m_dReferenceDate = GetDoubleFrom(m_ReferenceDate);

  // invalidate the yield curve because of reference date change
  Invalidate();
}

void YieldCurve::Validate() const
{
  if ( !m_bValidated )
  {
    const_cast<YieldCurve *>(this)->DoValidate();
  
    const_cast<YieldCurve *>(this)->m_bValidated = true;
  }
  //else: nothing to do, we had been already validated
}

void YieldCurve::CheckReferenceDate() const
{
  CHECK_COND(m_ReferenceDate.IsValid(), ITO33_YIELDCURVE_REFERENCEDATE_UNDEF);
}

double YieldCurve::GetContinuousRate(double dMaturity) const
{
  CheckReferenceDate();

  double dRate = GetZeroRate(dMaturity - m_dReferenceDate);

  return log(1. + dRate);
}


double YieldCurve::GetContinuousRate(Date maturity) const
{
  return GetContinuousRate( GetDoubleFrom(maturity) );
}


void YieldCurve::GetDiscountFactor(const double *pdMaturities, 
                                   double *pdDF, 
                                   size_t nNb) const
{
  GetCompoundFactor(pdMaturities, pdDF, nNb);
    
  for (size_t nI = 0; nI < nNb; nI++)
    pdDF[nI] = 1. / pdDF[nI];
}

void YieldCurve::GetCompoundFactor(const double *pdMaturities, 
                                   double *pdCF, 
                                   size_t nNb) const
{
  CheckReferenceDate();

  size_t
    nI;

  Array<double> 
    pdMaturitiesTmp(nNb);

  for(nI = 0; nI < nNb; nI++)
    pdMaturitiesTmp[nI] = pdMaturities[nI] - m_dReferenceDate;

  GetZeroRates(pdMaturitiesTmp.Get(), pdCF, nNb);
    
  for (nI = 0; nI < nNb; nI++)
    pdCF[nI] = pow(1. + pdCF[nI], pdMaturitiesTmp[nI]); 
}

double YieldCurve::GetForwardDiscountFactor(double dMaturity1, 
                                            double dMaturity2) const
{
  ASSERT( dMaturity1 <= dMaturity2 );

  double
    pdMaturitiesTmp[2],
    dFDF;

  pdMaturitiesTmp[0] = dMaturity1;
  pdMaturitiesTmp[1] = dMaturity2;

  GetForwardDiscountFactor(pdMaturitiesTmp, &dFDF, 2);

  return dFDF;
}

double YieldCurve::GetForwardDiscountFactor(Date date1, Date date2) const
{
  return GetForwardDiscountFactor(GetDoubleFrom(date1), GetDoubleFrom(date2));
}

//*********************************************************************

void YieldCurve::GetForwardDiscountFactor(const double *pdMaturities, 
                                          double *pdFDF, 
                                          size_t nNb) const
{
  Array<double> 
    pdCFTmp(nNb);

  GetCompoundFactor(pdMaturities, pdCFTmp.Get(), nNb);
  
  for(size_t nI = 0; nI < nNb - 1; nI++)
    pdFDF[nI] = pdCFTmp[nI] / pdCFTmp[nI + 1];
}


} // namespace finance

} // namespace ito33
