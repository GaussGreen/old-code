/////////////////////////////////////////////////////////////////////////////
// Name:        finance/bondlike/pepsaveragingperiod.cpp
// Purpose:     class for defining peps averaging period
// Author:      ITO33 Canada
// Created:     May 10, 2005
// RCS-ID:      $Id: pepsaveragingperiod.cpp,v 1.7 2006/08/19 22:43:51 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file finance/bondlike/pepsaveragingperiod.cpp
 */


#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/bondlike/pepsaveragingperiod.h"
#include "ito33/finance/bondlike/bonderror.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/bondlike/pepsaveragingperiod.h"

extern const ito33::finance::BondError 
   ITO33_GENERALIZEDPEPSLIKE_AVG_PERIOD_TOO_MANY_DAYS,
   ITO33_GENERALIZEDPEPSLIKE_NOT_STOCK_AVERAGING,
   ITO33_GENERALIZEDPEPSLIKE_NOT_CONVERSION_RATIO_AVERAGING,
   ITO33_GENERALIZEDPEPSLIKE_NEGATIVE_AVERAGE,
   ITO33_GENERALIZEDPEPSLIKE_AVG_PERIOD_START_DATE_AFTER_END_DATE,
   ITO33_GENERALIZEDPEPSLIKE_INVALID_NBSAMPLINAVERAGES,
   ITO33_GENERALIZEDPEPSLIKE_ZERO_NB_SAMPLES_USED,
   ITO33_GENERALIZEDPEPSLIKE_NO_CURRENT_AVERAGE_SET,
   ITO33_GENERALIZEDPEPSLIKE_INVALID_NB_SAMPLES_USED;

namespace ito33
{

namespace finance
{

/* static */ 
shared_ptr<PEPSAveragingPeriod> 
PEPSAveragingPeriod::CreateWithStock(Date averageStartDate, 
  Date averageEndDate, size_t nNumberOfSamplingAverages)
{
  shared_ptr<PEPSAveragingPeriod> pAvgPeriod
    ( 
      new PEPSAveragingPeriod(averageStartDate, averageEndDate, 
                              nNumberOfSamplingAverages)
    );
  
  pAvgPeriod->m_bHasStockAveraging = true;

  return pAvgPeriod;
}

/* static */
shared_ptr<PEPSAveragingPeriod> 
PEPSAveragingPeriod::CreateWithConversionRatio(Date averageStartDate, 
  Date averageEndDate, size_t nNumberOfSamplingAverages)
{
  shared_ptr<PEPSAveragingPeriod> pAvgPeriod
    (
       new PEPSAveragingPeriod(averageStartDate, averageEndDate, 
                            nNumberOfSamplingAverages)
     ); 

  pAvgPeriod->m_bHasStockAveraging = false;

  return pAvgPeriod;
}

PEPSAveragingPeriod::PEPSAveragingPeriod(Date averageStartDate, 
  Date averageEndDate, size_t nNumberOfSamplingAverages)
  : m_averageStartDate(averageStartDate), 
    m_averageEndDate(averageEndDate),
    m_nNumberOfSamplingAverages(nNumberOfSamplingAverages),
    m_nNbSamplesUsed(0),
    m_dCurrentAverage(-1.),
    m_bHasStockAveraging(false) 
{

  CHECK_COND(  m_nNumberOfSamplingAverages > 0,              
    ITO33_GENERALIZEDPEPSLIKE_INVALID_NBSAMPLINAVERAGES );

  CHECK_COND( averageStartDate < averageEndDate, 
    ITO33_GENERALIZEDPEPSLIKE_AVG_PERIOD_START_DATE_AFTER_END_DATE);
  
  CHECK_COND( Date::DaysDiff(averageStartDate, averageEndDate) <= 30, 
    ITO33_GENERALIZEDPEPSLIKE_AVG_PERIOD_TOO_MANY_DAYS);
}


double PEPSAveragingPeriod::GetCurrentStockAverage() const
{
  CHECK_COND( m_bHasStockAveraging, 
    ITO33_GENERALIZEDPEPSLIKE_NOT_STOCK_AVERAGING);

  return m_dCurrentAverage;
}

double PEPSAveragingPeriod::GetCurrentConversionRatioAverage() const
{
   CHECK_COND( !m_bHasStockAveraging, 
    ITO33_GENERALIZEDPEPSLIKE_NOT_CONVERSION_RATIO_AVERAGING);

   return m_dCurrentAverage;
}

void PEPSAveragingPeriod::SetCurrentStockAverage(double dCurrentStockAverage, 
                                                 size_t nNbSamplesUsed)
{
  CHECK_COND( m_bHasStockAveraging, 
    ITO33_GENERALIZEDPEPSLIKE_NOT_STOCK_AVERAGING);

  CHECK_COND(dCurrentStockAverage > 0., 
    ITO33_GENERALIZEDPEPSLIKE_NEGATIVE_AVERAGE);

  CHECK_COND( nNbSamplesUsed > 0, 
    ITO33_GENERALIZEDPEPSLIKE_ZERO_NB_SAMPLES_USED);

  m_dCurrentAverage = dCurrentStockAverage;  

  m_nNbSamplesUsed = nNbSamplesUsed;
}

void PEPSAveragingPeriod::SetCurrentConversionRatioAverage(
  double dCurrentConversionRatioAverage, size_t nNbSamplesUsed)
{   
  CHECK_COND( !m_bHasStockAveraging, 
    ITO33_GENERALIZEDPEPSLIKE_NOT_CONVERSION_RATIO_AVERAGING);
  
  CHECK_COND( dCurrentConversionRatioAverage > 0., 
    ITO33_GENERALIZEDPEPSLIKE_NEGATIVE_AVERAGE);

  CHECK_COND( nNbSamplesUsed > 0,
    ITO33_GENERALIZEDPEPSLIKE_ZERO_NB_SAMPLES_USED);

  m_dCurrentAverage = dCurrentConversionRatioAverage;

  m_nNbSamplesUsed = nNbSamplesUsed;
}

void PEPSAveragingPeriod::ValidateWith(Date valuationDate) const
{
  // The current average has to be set if 
  // valuation date is after the begining of the averaging period
  if ( valuationDate > m_averageStartDate )
  {  
    double dCurrentAverage = -1;
    
    if ( HasStockAveraging() )
      dCurrentAverage = GetCurrentStockAverage();
    else
      dCurrentAverage = GetCurrentConversionRatioAverage();

    CHECK_COND( dCurrentAverage > 0,
                ITO33_GENERALIZEDPEPSLIKE_NO_CURRENT_AVERAGE_SET); 

    // check that if we the valuationDate is before
    // the end of the averaging period the user did 
    // not enter a number of sample used greater than
    // it needs to be
    if ( valuationDate < m_averageEndDate )        
      CHECK_COND( m_nNbSamplesUsed < m_nNumberOfSamplingAverages, 
                  ITO33_GENERALIZEDPEPSLIKE_INVALID_NB_SAMPLES_USED );
  }
}

XML::Tag PEPSAveragingPeriod::Dump(XML::Tag& tagParent) const
{
  XML::Tag tagMe(XML_TAG_FINANCE_BONDLIKE_PEPS_AVERAGING_PERIOD, tagParent);

  tagMe.Element(XML_TAG_FINANCE_BONDLIKE_PEPS_START_DATE)(m_averageStartDate);

  tagMe.Element(XML_TAG_FINANCE_BONDLIKE_PEPS_END_DATE)(m_averageEndDate);
  
  tagMe.Element(XML_TAG_FINANCE_BONDLIKE_PEPS_STOCK_AVERAGING)(m_bHasStockAveraging);

  tagMe.Element(XML_TAG_FINANCE_BONDLIKE_PEPS_NB_SAMPLING_AVERAGES)
    (m_nNumberOfSamplingAverages);

  if ( m_dCurrentAverage > 0 )
  {
    tagMe.Element( XML_TAG_FINANCE_BONDLIKE_PEPS_CURRENT_AVERAGE)
      (m_dCurrentAverage);
    tagMe.Element(XML_TAG_FINANCE_BONDLIKE_PEPS_NB_SAMPLES_USED)
      (m_nNbSamplesUsed);
  }
 
  return tagMe;
}


} // namespace finance

} // namespace ito33
