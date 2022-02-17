/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/exoticoption/asianoption.cpp
// Purpose:     Implementation of shared ptr for Option class
// Author:      ITO 33 Canada
// Created:     March 29, 2005
// RCS-ID:      $Id: asianoption.cpp,v 1.8 2006/08/19 22:40:30 wang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"
#include "ito33/enum_values_names.h"

#include "ito33/finance/exoticoption/asianoption.h"

#include "ito33/finance/derivative_visitor.h"
#include "ito33/finance/optionerror.h"
#include "ito33/finance/optionlike.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/frequency.h"
#include "ito33/xml/finance/exoticoption/asianoption.h"

extern const ito33::finance::OptionError 
      ITO33_OPTION_NEGATIVE_STRIKE, ITO33_ASIANOPTION_NEGATIVE_CURRENTAVERAGE,
      ITO33_ASIANOPTION_AVG_START_DATE_AFTER_MATURITY_DATE,
      ITO33_ASIANOPTION_AVG_END_DATE_BEFORE_AVG_START_DATE,
      ITO33_ASIANOPTION_AVG_END_DATE_AFTER_MATURITY_DATE,
      ITO33_ASIANOPTION_AVG_END_DATE_NOT_SET,
      ITO33_ASIANOPTION_INVALID_NBSAMPLINAVERAGES,
      ITO33_ASIANOPTION_INVALID_NB_SAMPLES_USED,
      ITO33_ASIANOPTION_ZERO_NB_SAMPLES_USED;

namespace ito33
{

namespace finance
{

AsianOption::AsianOption(double dFixedStrike, Date maturityDate,
         OptionType optionType, ExerciseType exerciseType, 
         Date averageStartDate, size_t nNumberOfSamplingAverages)
         :OptionLike(maturityDate, optionType, exerciseType),
         m_dFixedStrike(dFixedStrike),
         m_dCurrentAverage(0.),
         m_averageStartDate(averageStartDate),
         m_nNumberOfSamplingAverages(nNumberOfSamplingAverages),
         m_averageEndDate(maturityDate),
         m_nNbSamplesUsed(0)
{
  
  CHECK_COND( m_dFixedStrike > 0, ITO33_OPTION_NEGATIVE_STRIKE);

  CHECK_COND(  m_nNumberOfSamplingAverages > 0,              
    ITO33_ASIANOPTION_INVALID_NBSAMPLINAVERAGES );

  CHECK_COND( m_maturityDate > m_averageStartDate, 
    ITO33_ASIANOPTION_AVG_START_DATE_AFTER_MATURITY_DATE);
}

AsianOption::AsianOption(Date maturityDate, OptionType optionType, 
        ExerciseType exerciseType, Date averageStartDate, 
        size_t nNumberOfSamplingAverages)
        :OptionLike(maturityDate, optionType, exerciseType),
        m_dFixedStrike(-1.),
        m_dCurrentAverage(0.),
        m_averageStartDate(averageStartDate),
        m_nNumberOfSamplingAverages(nNumberOfSamplingAverages),
        m_averageEndDate(maturityDate),
        m_nNbSamplesUsed(0)
{
  
  CHECK_COND(  m_nNumberOfSamplingAverages > 0,              
    ITO33_ASIANOPTION_INVALID_NBSAMPLINAVERAGES );

  CHECK_COND( m_maturityDate >  m_averageStartDate, 
    ITO33_ASIANOPTION_AVG_START_DATE_AFTER_MATURITY_DATE);
}

void AsianOption::SetCurrentAverage(double dCurrentAverage, size_t 
                                    nNbSamplesUsed)
{
  CHECK_COND( dCurrentAverage > 0, 
    ITO33_ASIANOPTION_NEGATIVE_CURRENTAVERAGE);
  
  CHECK_COND( nNbSamplesUsed < m_nNumberOfSamplingAverages,
              ITO33_ASIANOPTION_INVALID_NB_SAMPLES_USED);

  CHECK_COND( nNbSamplesUsed > 0,
              ITO33_ASIANOPTION_ZERO_NB_SAMPLES_USED);

  m_dCurrentAverage = dCurrentAverage;

  m_nNbSamplesUsed = nNbSamplesUsed;
}

void AsianOption::SetAverageEndDate(Date averageEndDate)
{
  m_averageEndDate = averageEndDate;
  
  CHECK_COND( m_averageEndDate > m_averageStartDate,
    ITO33_ASIANOPTION_AVG_END_DATE_BEFORE_AVG_START_DATE);

  CHECK_COND( m_maturityDate >= m_averageEndDate,
    ITO33_ASIANOPTION_AVG_END_DATE_AFTER_MATURITY_DATE);
}

double AsianOption::GetFixedStrike() const
{
  return m_dFixedStrike;
}

void AsianOption::Visit(DerivativeVisitor& visitor) const
{
  visitor.OnAsianOption(*this);
}

XML::Tag AsianOption::Dump(XML::Tag& tagParent) const
{
  // this is pretty straightforward: just dump all our data under a single
  // option tag
  XML::Tag tagOption(XML_TAG_ASIAN_OPTION_ROOT, tagParent);
 
  DumpMe(tagOption);

  //strike output is optional in case of floating strike option
  if ( m_dFixedStrike > 0 )
    tagOption.Element(XML_TAG_FINANCE_STRIKE)(m_dFixedStrike);

  tagOption.Element(XML_TAG_ASIAN_OPTION_AVG_START_DATE)(m_averageStartDate);

  if ( m_averageEndDate !=  m_maturityDate)
    tagOption.Element(XML_TAG_ASIAN_OPTION_AVG_END_DATE)(m_averageEndDate);

  tagOption.Element(XML_TAG_ASIAN_OPTION_NB_SAMPLING_AVERAGES)
    (m_nNumberOfSamplingAverages);

  if ( m_dCurrentAverage > 0 )
  {
    tagOption.Element(XML_TAG_ASIAN_OPTION_CURRENT_AVERAGE)(m_dCurrentAverage);
    tagOption.Element(XML_TAG_ASIAN_OPTION_NB_SAMPLES_USED)(m_nNbSamplesUsed);
  }

  return tagOption;
}


} // namespace finance

} // namespace ito33
