/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/varianceswapterms.cpp
// Purpose:     Implementation of variance swap class
// Created:     2006/07/20
// RCS-ID:      $Id: varianceswapterms.cpp,v 1.3 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"
#include "ito33/enum_values_names.h"

#include "ito33/finance/error.h"
#include "ito33/finance/varianceswapterms.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/swaptype.h"
#include "ito33/xml/finance/returntype.h"
#include "ito33/xml/finance/swappayofftype.h"
#include "ito33/xml/finance/varianceswap.h"

extern const ito33::finance::Error 
  ITO33_VARIANCESWAP_INVALID_MATURITY,
  ITO33_VARIANCESWAP_INVALID_STARTOFSAMPLINGPERIOD,
  ITO33_VARIANCESWAP_INVALID_NBSAMPLINGRETURNS,
  ITO33_VARIANCESWAP_START_AFTER_MATURITY,
  ITO33_VARIANCESWAP_CAP_MULTIPLIER_TOO_SMALL,
  ITO33_VARIANCESWAP_CAP_MULTIPLIER_TOO_LARGE,
  ITO33_VARIANCESWAP_NON_STRICTLY_POSITIVE_ANNUAL_RETURN_FREQUENCY,
  ITO33_VARIANCESWAP_NON_STRICTLY_POSITIVE_CORRIDOR,
  ITO33_VARIANCESWAP_UP_CORRIDOR_INVALID,
  ITO33_VARIANCESWAP_DOWN_CORRIDOR_INVALID;

namespace ito33
{

namespace finance
{

VarianceSwapTerms::VarianceSwapTerms(Date maturityDate,
                                     SwapType swapType,
                                     Date startOfSamplingPeriod,
                                     size_t nNbSamplingReturns)
  : m_maturityDate(maturityDate),
    m_swapType(swapType),
    m_startOfSamplingPeriod(startOfSamplingPeriod),
    m_nNbSamplingReturns(nNbSamplingReturns),
    m_dCapMultiplier(-1.),
    m_nAnnualReturnFrequency(252),
    m_dUpCorridorBarrier(0.0),
    m_dDownCorridorBarrier(0.0)
{
  //default to log returns
  m_returnType = finance::Return_Log;

  // standard input checks
  CHECK_COND( m_maturityDate.IsValid(), ITO33_VARIANCESWAP_INVALID_MATURITY );

  CHECK_COND( m_startOfSamplingPeriod.IsValid(), 
              ITO33_VARIANCESWAP_INVALID_STARTOFSAMPLINGPERIOD );

  CHECK_COND( m_startOfSamplingPeriod < m_maturityDate, 
              ITO33_VARIANCESWAP_START_AFTER_MATURITY );

  CHECK_COND( m_nNbSamplingReturns > 0, 
              ITO33_VARIANCESWAP_INVALID_NBSAMPLINGRETURNS );
}

void VarianceSwapTerms::SetCapMultiplier(double dCapMultiplier)
{
  CHECK_COND( dCapMultiplier >= 1., 
    ITO33_VARIANCESWAP_CAP_MULTIPLIER_TOO_SMALL);

  CHECK_COND( dCapMultiplier <= 10.0, 
    ITO33_VARIANCESWAP_CAP_MULTIPLIER_TOO_LARGE);

  m_dCapMultiplier = dCapMultiplier;
}

void VarianceSwapTerms::SetReturnType(ReturnType returnType)
{
  m_returnType = returnType;
}

void VarianceSwapTerms::SetAnnualReturnFrequency(size_t nAnnualReturnFrequency)
{
  CHECK_COND( nAnnualReturnFrequency > 0,
    ITO33_VARIANCESWAP_NON_STRICTLY_POSITIVE_ANNUAL_RETURN_FREQUENCY);

  m_nAnnualReturnFrequency = nAnnualReturnFrequency;
}

void VarianceSwapTerms::SetUpCorridorBarrier(double dUpCorridorBarrier)
{
  CHECK_COND( dUpCorridorBarrier > 0.0, 
    ITO33_VARIANCESWAP_NON_STRICTLY_POSITIVE_CORRIDOR);

  if ( m_dDownCorridorBarrier > 0.0 )
    CHECK_COND( dUpCorridorBarrier < m_dDownCorridorBarrier, 
                ITO33_VARIANCESWAP_UP_CORRIDOR_INVALID);

  m_dUpCorridorBarrier = dUpCorridorBarrier;
}

void VarianceSwapTerms::SetDownCorridorBarrier(double dDownCorridorBarrier)
{
  CHECK_COND( dDownCorridorBarrier > 0.0, 
    ITO33_VARIANCESWAP_NON_STRICTLY_POSITIVE_CORRIDOR);

  if ( m_dUpCorridorBarrier > 0.0 )
    CHECK_COND( dDownCorridorBarrier > m_dUpCorridorBarrier, 
                ITO33_VARIANCESWAP_DOWN_CORRIDOR_INVALID);

  m_dDownCorridorBarrier = dDownCorridorBarrier;
}

void VarianceSwapTerms::DumpMe(XML::Tag& tagParent) const
{
  tagParent.Element(XML_TAG_FINANCE_MATURITY)(m_maturityDate);

  tagParent.Element(XML_TAG_VARIANCESWAP_STARTOFPERIOD)
                   (m_startOfSamplingPeriod);

  tagParent.Element(XML_TAG_VARIANCESWAP_NBSAMPLINGRETURNS)
                   (m_nNbSamplingReturns);

  tagParent.Element(XML_TAG_VARIANCESWAP_SWAP_TYPE)
                   (
                     GetNameFromEnumValue(
                       m_swapType,
                       SIZEOF(g_swapTypes),
                       g_swapTypes)
                   );

  tagParent.Element(XML_TAG_VARIANCESWAP_RETURN_TYPE)
                   (
                     GetNameFromEnumValue(
                       m_returnType,
                       SIZEOF(g_returnTypes),
                       g_returnTypes)
                   );  

  if ( m_nAnnualReturnFrequency != 252 )
    tagParent.Element(XML_TAG_VARIANCESWAP_ANNUALRETURNFREQUENCY)
                     (m_nAnnualReturnFrequency);

  if ( m_dCapMultiplier > 0 )
    tagParent.Element(XML_TAG_VARIANCESWAP_CAP_MULTIPLIER)
                     (m_dCapMultiplier);

  if ( m_dUpCorridorBarrier > 0 )
    tagParent.Element(XML_TAG_VARIANCESWAP_UP_CORRIDOR)
                     (m_dUpCorridorBarrier);

  if ( m_dDownCorridorBarrier > 0 )
    tagParent.Element(XML_TAG_VARIANCESWAP_DOWN_CORRIDOR)
                     (m_dDownCorridorBarrier);
}


} // namespace finance

} // namespace ito33
