/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/bondlike/generalizedpepslike.cpp
// Purpose:     financial class for PERCS-like instrument
// Created:     2004/08/17 
// RCS-ID:      $Id: generalizedpepslike.cpp,v 1.11 2006/08/19 22:43:51 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"
#include "ito33/error.h"

#include "ito33/finance/sessiondata.h"

#include "ito33/finance/derivative_visitor.h"
#include "ito33/finance/derivative_modifyingvisitor.h"

#include "ito33/finance/bondlike/bondliketerms.h"
#include "ito33/finance/bondlike/pepsaveragingperiod.h"
#include "ito33/finance/bondlike/bonderror.h"
#include "ito33/finance/bondlike/generalizedpepslike.h"
#include "ito33/finance/bondlike/callschedule.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/bondlike/generalizedpepslike.h"
#include "ito33/xml/finance/bondlike/common.h"

extern const ito33::finance::Error   
  ITO33_AVERAGINGPERIOD_ALREADY_SET,
  ITO33_AVERAGINGPERIOD_NULL;

extern const ito33::finance::BondError 
  ITO33_BONDLIKE_NULL_CALLSCHEDULE,
  ITO33_BONDLIKE_NULL_BONDLIKETERMS,
  ITO33_MANDATORY_GUARANTEED_YIELD_TO_CALL_INVALID,
  ITO33_GENERALIZEDPEPSLIKE_CALL_ALREADY_SET,
  ITO33_GENERALIZEDPEPSLIKE_DOWNSIDE_CONVERSION_RATIO,
  ITO33_GENERALIZEDPEPSLIKE_UPSIDE_CONVERSION_RATIO,
  ITO33_GENERALIZEDPEPSLIKE_LOWER_HIGHER_STRIKE,
  ITO33_NULL_GENERALIZEDPEPSLIKE_CALL,
  ITO33_GENERALIZEDPEPSLIKE_AVERAGINGPERIOD_AVG_END_DATE_INVALID,
  ITO33_GENERALIZEDPEPSLIKE_NO_CURRENT_AVERAGE_SET;


namespace ito33
{

namespace finance
{

GeneralizedPEPSLike::GeneralizedPEPSLike(
                      const shared_ptr<BondLikeTerms>& pBondLikeTerms,
                      double dDownsideConversionRatio,
                      double dLowerStrike,
                      double dUpsideBaseConversionRatio,
                      double dHigherStrike)
    : ConvertibleLike(pBondLikeTerms),
      m_dDownsideConversionRatio(dDownsideConversionRatio), 
      m_dLowerStrike(dLowerStrike),
      m_dUpsideBaseConversionRatio(dUpsideBaseConversionRatio),
      m_dHigherStrike(dHigherStrike),
      m_bHasOptionalConversion(false),
      m_bCallSet(false)
{
  CHECK_PTR
  (
    pBondLikeTerms,
    ITO33_BONDLIKE_NULL_BONDLIKETERMS
  );

  CHECK_COND
  (
    dDownsideConversionRatio > 0.,
    ITO33_GENERALIZEDPEPSLIKE_DOWNSIDE_CONVERSION_RATIO
  );

  CHECK_COND
  (
    dUpsideBaseConversionRatio > 0,
    ITO33_GENERALIZEDPEPSLIKE_UPSIDE_CONVERSION_RATIO
  );

  CHECK_COND
  (
    dLowerStrike > 0 && dHigherStrike > dLowerStrike,
    ITO33_GENERALIZEDPEPSLIKE_LOWER_HIGHER_STRIKE
  );
}

void GeneralizedPEPSLike::SetCallFixedCash
                              (const shared_ptr<CallSchedule>& pCallSchedule)
{
  CHECK_COND(!m_bCallSet, ITO33_GENERALIZEDPEPSLIKE_CALL_ALREADY_SET);

  m_pCallSchedule = CHECK_PTR(pCallSchedule, ITO33_BONDLIKE_NULL_CALLSCHEDULE);

  // Check if any of the call periods have a yield to call
  CHECK_COND( !m_pCallSchedule->HasYield(), 
              ITO33_MANDATORY_GUARANTEED_YIELD_TO_CALL_INVALID );

  m_bCallSet = true;
}

void GeneralizedPEPSLike::SetAveragingPeriod
    (const shared_ptr<PEPSAveragingPeriod>& pAveragingPeriod)
{

  CHECK_COND(!m_pAveragingPeriod, ITO33_AVERAGINGPERIOD_ALREADY_SET);

  m_pAveragingPeriod = CHECK_PTR(pAveragingPeriod, ITO33_AVERAGINGPERIOD_NULL);

  // Check that the maturity date and the end of averaging period
  // are no more than 5 days apart.
  CHECK_COND( Date::DaysDiff( m_pBondLikeTerms->GetMaturityDate(), 
    m_pAveragingPeriod->GetAverageEndDate() ) <= 5, 
    ITO33_GENERALIZEDPEPSLIKE_AVERAGINGPERIOD_AVG_END_DATE_INVALID);

}

void GeneralizedPEPSLike::SetGeneralizedPEPSLikeCall
        (const shared_ptr<GeneralizedPEPSLikeCall>& pGeneralizedPEPSLikeCall)
{
  CHECK_COND(!m_bCallSet, ITO33_GENERALIZEDPEPSLIKE_CALL_ALREADY_SET);

  CHECK_PTR
        (
          pGeneralizedPEPSLikeCall,
          ITO33_NULL_GENERALIZEDPEPSLIKE_CALL
        );

  m_pGeneralizedPEPSLikeCall = pGeneralizedPEPSLikeCall;

  m_bCallSet = true;
}

void GeneralizedPEPSLike::ValidateWith(const SessionData& sessionData) const
{
  ConvertibleLike::ValidateWith(sessionData);

  if ( m_pAveragingPeriod )
    m_pAveragingPeriod->ValidateWith( sessionData.GetValuationDate() );
}

void GeneralizedPEPSLike::Visit(DerivativeVisitor& visitor) const
{
   visitor.OnGeneralizedPEPSLike(*this);
}

void GeneralizedPEPSLike::Visit(DerivativeModifyingVisitor& visitor)
{
   visitor.OnGeneralizedPEPSLike(*this);
}

XML::Tag GeneralizedPEPSLike::Dump(XML::Tag& tagParent) const
{ 
  XML::Tag tagPEPS(XML_TAG_GENERALIZED_PEPSLIKE_ROOT, tagParent);

  ConvertibleLike::DumpMe(tagPEPS);

  tagPEPS.Element(XML_TAG_GENERALIZED_PEPSLIKE_LOWER_STRIKE)
                 (m_dLowerStrike);
  tagPEPS.Element(XML_TAG_GENERALIZED_PEPSLIKE_HIGHER_STRIKE)
                 (m_dHigherStrike);
  tagPEPS.Element(XML_TAG_GENERALIZED_PEPSLIKE_DOWNSIDE_CONVERSION_RATIO)
                 (m_dDownsideConversionRatio);
  tagPEPS.Element(XML_TAG_GENERALIZED_PEPSLIKE_UPSIDE_BASE_CONVERSION_RATIO)
                 (m_dUpsideBaseConversionRatio);

  // optional conversions
  tagPEPS.Element(XML_TAG_GENERALIZED_PEPSLIKE_HAS_OPTIONAL_CONVERSION)
                 (m_bHasOptionalConversion);

  // optional calls
  if ( m_pCallSchedule && !( m_pCallSchedule->GetAll().empty() ) )
  {
    tagPEPS.Element(*m_pCallSchedule);
  }
  
  if( m_pGeneralizedPEPSLikeCall)
  {
    tagPEPS.Element(*m_pGeneralizedPEPSLikeCall);
  }

  if ( m_pAveragingPeriod )
  {
    m_pAveragingPeriod->Dump(tagPEPS);
  }

  return tagPEPS;
}


} // namespace finance

} // namespace ito33
