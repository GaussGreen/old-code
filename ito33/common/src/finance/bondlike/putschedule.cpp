/////////////////////////////////////////////////////////////////////////////
// Name:        finance/bondlike/putschedule.cpp
// Purpose:     class for standard put provision of a bond
// Author:      ZHANG Yunzhi
// Created:     2004 may 6
// RCS-ID:      $Id: putschedule.cpp,v 1.31 2006/08/19 22:43:51 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/bondlike/bonderror.h"
#include "ito33/finance/bondlike/putperiod.h"
#include "ito33/finance/bondlike/putschedule.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/bondlike/putschedule.h"
#include "ito33/xml/finance/bondlike/common.h"
#include "ito33/xml/finance/common.h"

extern const ito33::finance::BondError 
  ITO33_BONDLIKE_PUTPERIOD_NOT_STRIKE,
  ITO33_BONDLIKE_PUTPERIOD_NOT_YIELD,
  ITO33_BONDLIKE_PUT_SCHEDULE_DUPLICATEDDATE,
  ITO33_BONDLIKE_PUT_SCHEDULE_BADSTRIKE,
  ITO33_BONDLIKE_PUT_SCHEDULE_GUARANTEED_YIELD_TOO_SMALL,
  ITO33_BONDLIKE_PUT_SCHEDULE_GUARANTEED_YIELD_TOO_LARGE,
  ITO33_BONDLIKE_PUT_SCHEDULE_INCONSISTENT_YIELD;

namespace ito33
{

namespace finance
{

//)))))))))))))))))))))))))))))) PUTPERIOD PART )))))))))))))))))))))))))))
double PutPeriod::GetStrike() const
{
  CHECK_COND( !HasYield(), ITO33_BONDLIKE_PUTPERIOD_NOT_STRIKE);

  return m_dStrike;
}

double PutPeriod::GetGuaranteedYield() const
{
  CHECK_COND( HasYield(), ITO33_BONDLIKE_PUTPERIOD_NOT_YIELD);

  return m_dYield;
}

//==============================PutSchedule part============================

void PutSchedule::AddPutWithStrike(Date putDate, double dPutStrike)
{
  CHECK_COND(
      m_puts.find(putDate) == m_puts.end(),
      ITO33_BONDLIKE_PUT_SCHEDULE_DUPLICATEDDATE);

  CHECK_COND(
      dPutStrike > 0,
      ITO33_BONDLIKE_PUT_SCHEDULE_BADSTRIKE);

  m_puts[putDate] = PutData(dPutStrike, false);
  m_pPutPeriods.clear();
}


void PutSchedule::AddPutWithYield(Date putDate, double dPutYield)
{
  CHECK_COND(
      m_puts.find(putDate) == m_puts.end(),
      ITO33_BONDLIKE_PUT_SCHEDULE_DUPLICATEDDATE);

  CHECK_COND(
      dPutYield >= -0.1,
      ITO33_BONDLIKE_PUT_SCHEDULE_GUARANTEED_YIELD_TOO_SMALL);

  CHECK_COND( 
      dPutYield <= .2, 
      ITO33_BONDLIKE_PUT_SCHEDULE_GUARANTEED_YIELD_TOO_LARGE);

  m_puts[putDate] = PutData(dPutYield, true);
  m_pPutPeriods.clear();
}

bool PutSchedule::HasYield() const
{
  Elements::const_iterator iter;
  for(iter = m_puts.begin(); iter != m_puts.end(); ++iter)
  {
    if( iter->second.bYield )
      return true;
  }
  return false;
}

void PutSchedule::Validate() const
{
  // If any put period has a yield, then KeepAccrued must be true
  // and ForfeitCoupon must be false.
  CHECK_COND( 
    !HasYield() || (m_bKeepAccrued == true && m_bForfeitCoupon == false),
    ITO33_BONDLIKE_PUT_SCHEDULE_INCONSISTENT_YIELD);
}

XML::Tag PutSchedule::Dump(XML::Tag& tagParent) const
{
  XML::Tag tagMe(XML_TAG_BONDLIKE_PUTSCHEDULE_ROOT, tagParent);

  tagMe.Element(XML_TAG_BONDLIKE_KEEPACCRUED)(m_bKeepAccrued);
  tagMe.Element(XML_TAG_BONDLIKE_FORFEITCOUPON)(m_bForfeitCoupon);

  XML::Tag tagPeriods(XML_TAG_BONDLIKE_PUTCHEDULE_PUTS, tagMe);

  Elements::const_iterator iter;

  for (iter = m_puts.begin(); iter != m_puts.end(); ++iter)
  {
    XML::Tag tagPut(XML_TAG_BONDLIKE_PUTSCHEDULE_PUT_ROOT, tagPeriods);

    tagPut.Element(XML_TAG_FINANCE_DATE)(iter->first);

    if ( iter->second.bYield )
      tagPut.Element(XML_TAG_FINANCE_YIELD)(iter->second.dYield);
    else
      tagPut.Element(XML_TAG_FINANCE_STRIKE)(iter->second.dStrike);     
  }

  return tagMe;
}

const PutSchedule::CollectionElements& PutSchedule::GetPutPeriods() const
{
  if(m_pPutPeriods.empty())
    CreateCollectionElements();

  return m_pPutPeriods;
}

void PutSchedule::CreateCollectionElements() const
{
  Elements::const_iterator iter;
  PutSchedule *pCF = const_cast<PutSchedule *>(this);

  pCF->m_pPutPeriods.clear();
  for (iter = m_puts.begin(); iter != m_puts.end(); ++iter)
  {
    pCF->m_pPutPeriods.push_back
            (
              PutPeriod
              (
               iter->first, 
               iter->second.bYield? iter->second.dYield : iter->second.dStrike,
               iter->second.bYield
              )
            );
  }
}

} // namespace finance

} // namespace ito33
