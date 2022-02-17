/////////////////////////////////////////////////////////////////////////////
// Name:        finance/bondlike/callschedule.cpp
// Purpose:     class for standard call of bond
// Created:     2004/08/17 
// RCS-ID:      $Id: callschedule.cpp,v 1.43 2006/08/19 22:43:51 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <algorithm>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/bondlike/bonderror.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/callperiod.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/bondlike/callschedule.h"
#include "ito33/xml/finance/bondlike/common.h"
#include "ito33/xml/finance/common.h"

extern const ito33::Error ITO33_NULL_PARAM;

extern const ito33::finance::BondError
  ITO33_BONDLIKE_CALLS_SOFTAFTERHARD,
  ITO33_BONDLIKE_CALLS_WRONGSCHEDULE,
  ITO33_BONDLIKE_SOFTCALL_BADTRIGGER,
  ITO33_BONDLIKE_CALL_WRONG_START_END,
  ITO33_BONDLIKE_CALL_BADSTRIKE,
  ITO33_BONDLIKE_CALL_GUARANTEED_YIELD_TOO_SMALL,
  ITO33_BONDLIKE_CALL_GUARANTEED_YIELD_TOO_LARGE,
  ITO33_BONDLIKE_CALL_STRIKE_AND_YIELD,
  ITO33_BONDLIKE_CALLPERIOD_NOT_STRIKE,
  ITO33_BONDLIKE_CALLPERIOD_NOT_YIELD,
  ITO33_BONDLIKE_CALL_INCONSISTENT_YIELD;

namespace ito33
{

namespace finance
{

double CallPeriod::GetStrike() const
{
  CHECK_COND( !HasYield(), ITO33_BONDLIKE_CALLPERIOD_NOT_STRIKE);

  return m_dStrike;
}

double CallPeriod::GetGuaranteedYield() const
{
  CHECK_COND( HasYield(), ITO33_BONDLIKE_CALLPERIOD_NOT_YIELD);

  return m_dYield;
}

CallPeriod::CallPeriod(Date startDate, Date endDate)
                     : m_startDate(startDate), m_endDate(endDate),
                       m_dTriggerRate(0.), m_dStrike(0), m_dYield(0.)
{
  CHECK_COND(startDate <= endDate, ITO33_BONDLIKE_CALL_WRONG_START_END);
}

/* static */ shared_ptr<CallPeriod>
CallPeriod::CreateWithStrike
            (Date startDate, Date endDate, double dStrike)
{
  CHECK_COND(dStrike > 0, ITO33_BONDLIKE_CALL_BADSTRIKE);
  
  shared_ptr<CallPeriod> pCallPeriod( new CallPeriod(startDate, endDate) );

  pCallPeriod->m_dStrike = dStrike;
  pCallPeriod->m_bHasYield = false;


  return pCallPeriod;
}

/* static */ shared_ptr<CallPeriod>
CallPeriod::CreateWithYield
            (Date startDate, Date endDate, double dYield)
{
  CHECK_COND( dYield >= -0.1, 
              ITO33_BONDLIKE_CALL_GUARANTEED_YIELD_TOO_SMALL);

  CHECK_COND( dYield <= .2, 
              ITO33_BONDLIKE_CALL_GUARANTEED_YIELD_TOO_LARGE);

  shared_ptr<CallPeriod> pCallPeriod( new CallPeriod(startDate, endDate) );

  pCallPeriod->m_dYield = dYield;
  pCallPeriod->m_bHasYield = true;

  return pCallPeriod;
}

void CallPeriod::SetTrigger(double dTriggerRate)
{
  CHECK_COND(dTriggerRate > 0, ITO33_BONDLIKE_SOFTCALL_BADTRIGGER);

  m_dTriggerRate = dTriggerRate;
}

void CallSchedule::AddCallPeriod(const shared_ptr<CallPeriod>& pCallPeriod)
{
  CHECK_PTR(pCallPeriod, ITO33_NULL_PARAM);

  Elements::iterator iterBefore, iterNext;

  for ( iterBefore = m_calls.end(), iterNext = m_calls.begin();
        iterNext != m_calls.end() 
        && pCallPeriod->GetEndDate() > (*iterNext)->GetStartDate();
        iterBefore = iterNext, iterNext++)
    ;

  // at this moment call.endDate <= iterNext->startDate or iterNext == end()
  // and we compare call to iterBefore to check schedule error
  if(
      // wrong schedule with the next one : equal
      iterNext != m_calls.end() 
      && pCallPeriod->GetStartDate() == (*iterNext)->GetEndDate() 
    ||
      // wrong schedule with the last one
      (iterBefore != m_calls.end() 
      && pCallPeriod->GetStartDate() < (*iterBefore)->GetEndDate())
      // we have always endDate > iterBefore->Start, thus eqaul is not possible
    )
  {
    throw EXCEPTION(ITO33_BONDLIKE_CALLS_WRONGSCHEDULE);
  }

  if ( pCallPeriod->IsSoft() )
  {
    CHECK_COND
    (
      iterBefore == m_calls.end() || (*iterBefore)->IsSoft(),
      ITO33_BONDLIKE_CALLS_SOFTAFTERHARD
    );
  }
  else
  {
    CHECK_COND
    (
      iterNext == m_calls.end() || !( (*iterNext)->IsSoft() ),
      ITO33_BONDLIKE_CALLS_SOFTAFTERHARD
    );
  }

  m_calls.insert(iterNext, pCallPeriod);
}

void CallSchedule::Validate() const
{
  // If any call period has a yield, then KeepAccrued must be true
  // and ForfeitCoupon must be false.
  CHECK_COND( 
    !HasYield() || (m_bKeepAccrued == true && m_bForfeitCoupon == false),
    ITO33_BONDLIKE_CALL_INCONSISTENT_YIELD);
}

bool CallSchedule::HasYield() const
{
  Elements::const_iterator iter;
  for (iter = m_calls.begin(); iter != m_calls.end(); ++iter)
  {
    if ( (*iter)->HasYield() )
      return true;
  }

  return false;
}

bool CallSchedule::HasSoftCall() const
{
  Elements::const_iterator iter;
  for (iter = m_calls.begin(); iter != m_calls.end(); ++iter)
  {
    if ( (*iter)->IsSoft() )
      return true;
  }

  return false;
}

bool CallSchedule::HasMakeWhole() const 
{
  return m_makeWholeType != MakeWholeType_Max
    && !m_calls.empty() && (*m_calls.begin())->GetTrigger() > 0;
}

XML::Tag CallSchedule::Dump(XML::Tag& tagParent) const
{
  XML::Tag tagMe(XML_TAG_BONDLIKE_CALLSCHEDULE_ROOT, tagParent);

  Call::DumpMe(tagMe);

  XML::Tag tagPeriods(XML_TAG_BONDLIKE_CALLSCHEDULE_CALLPERIODS, tagMe);

  Elements::const_iterator iter;
  for(iter = m_calls.begin(); iter != m_calls.end(); ++iter)
    tagPeriods.Element(*(*iter));

  return tagMe;
}

XML::Tag CallPeriod::Dump(XML::Tag& tagParent) const
{
  XML::Tag tagMe(XML_TAG_BONDLIKE_CALLPERIOD_ROOT, tagParent);

  tagMe.Element(XML_TAG_BONDLIKE_STARTDATE)(m_startDate);

  tagMe.Element(XML_TAG_BONDLIKE_ENDDATE)(m_endDate);

  if ( m_bHasYield )
    tagMe.Element(XML_TAG_FINANCE_YIELD)(m_dYield);
  else
    tagMe.Element(XML_TAG_FINANCE_STRIKE)(m_dStrike);

  if (m_dTriggerRate > 0)
    tagMe.Element(XML_TAG_BONDLIKE_TRIGGERRATE)(m_dTriggerRate);

  return tagMe;
}


} // namespace finance

} // namespace ito33
