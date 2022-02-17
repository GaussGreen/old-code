/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/bondlike/call.cpp
// Purpose:     base call class for bond-like
// Author:      Wang
// Created:     2004/09/05 
// RCS-ID:      $Id: call.cpp,v 1.21 2006/08/19 22:43:51 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file fcommon/src/finance/bondlike/call.cpp
 */

#include "ito33/useexception.h"

#include "ito33/finance/bondlike/bonderror.h"
#include "ito33/finance/bondlike/call.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/bondlike/call.h"
#include "ito33/xml/finance/bondlike/common.h"
#include "ito33/xml/finance/common.h"

#include "ito33/xml/finance/bondlike/trigger.h"

extern const ito33::finance::BondError
  ITO33_BONDLIKE_MAKEWHOLETYPE_UNDEFINED,
  ITO33_BONDLIKE_MAKEWHOLETYPE_NOTPREMIUM,
  ITO33_BONDLIKE_MAKEWHOLETYPE_NOTCOUPON,
  ITO33_BONDLIKE_MAKEWHOLETYPE_INITIAL_PREMIUM,
  ITO33_BONDLIKE_CALL_NOTICEPERIOD,
  ITO33_BONDLIKE_CALL_TRIGGERPERIOD,
  ITO33_BONDLIKE_TRIGGERASPERCENTAGEOF;

namespace ito33
{

namespace finance
{

void Call::SetNoticePeriod(size_t nNoticePeriod)
{
  CHECK_COND( nNoticePeriod >= 1 && nNoticePeriod <= 120,
              ITO33_BONDLIKE_CALL_NOTICEPERIOD );

  m_nNoticePeriod = nNoticePeriod;
}

void Call::SetPremiumMakeWhole(double dInitialPremium)
{
  CHECK_COND(dInitialPremium > 0, 
             ITO33_BONDLIKE_MAKEWHOLETYPE_INITIAL_PREMIUM);

  m_makeWholeType = MakeWholeType_Premium;
  m_dMakeWholePremium = dInitialPremium;
}

void Call::SetTriggerCheckPeriod(size_t nTriggerPeriod, size_t nTriggerHistory)
{
  CHECK_COND(nTriggerPeriod <= 30, ITO33_BONDLIKE_CALL_TRIGGERPERIOD);

  m_nTriggerPeriod   = nTriggerPeriod;
  m_nTriggerHistory  = nTriggerHistory;
}

void Call::SetTriggerAsPercentageOf(TriggerAsPercentageOf asPercentageOf)
{
  CHECK_COND
  ( 
    asPercentageOf < TriggerAsPercentageOf_Max,
    ITO33_BONDLIKE_TRIGGERASPERCENTAGEOF
  );

  m_triggerAsPercentageOf     = asPercentageOf;
}

MakeWholeType Call::GetMakeWholeType() const
{
  CHECK_COND(HasMakeWhole(), ITO33_BONDLIKE_MAKEWHOLETYPE_UNDEFINED);

  return m_makeWholeType;
}

double Call::GetMakeWholePremium() const
{
  CHECK_COND(m_makeWholeType == MakeWholeType_Premium,
              ITO33_BONDLIKE_MAKEWHOLETYPE_NOTPREMIUM);

  return m_dMakeWholePremium;
}

bool Call::IsPVCouponMakeWhole() const
{
  CHECK_COND(m_makeWholeType == MakeWholeType_Coupon,
              ITO33_BONDLIKE_MAKEWHOLETYPE_NOTCOUPON);

  return m_bPVCouponMakeWhole;
}

void Call::DumpMe(ito33::XML::Tag& tagRoot) const
{
  tagRoot.Element(XML_TAG_BONDLIKE_KEEPACCRUED)(m_bKeepAccrued);
  tagRoot.Element(XML_TAG_BONDLIKE_FORFEITCOUPON)(m_bForfeitCoupon);

  if ( HasNoticePeriod() )
    tagRoot.Element(XML_TAG_BONDLIKE_NOTICEPERIOD)(m_nNoticePeriod);

  if ( HasMakeWhole() )
  {
    XML::Tag tagMakeWhole(XML_TAG_BONDLIKE_CALL_MAKEWHOLE_ROOT, tagRoot);


    if(m_makeWholeType == MakeWholeType_Coupon)
    {
      tagMakeWhole.Element(XML_TAG_FINANCE_TYPE)
                          (XML_VALUE_BONDLIKE_CALL_MAKEWHOLE_TYPE_COUPON);
      tagMakeWhole.Element(XML_TAG_BONDLIKE_CALL_MAKEWHOLE_PVCOUPON)
                          (m_bPVCouponMakeWhole);
    }
    else
    {
      tagMakeWhole.Element(XML_TAG_FINANCE_TYPE)
                          (XML_VALUE_BONDLIKE_CALL_MAKEWHOLE_TYPE_PREMIUM);
      tagMakeWhole.Element(XML_TAG_BONDLIKE_CALL_MAKEWHOLE_PREMIUM)
                          (m_dMakeWholePremium);
    }
  }

  if (m_nTriggerPeriod)
  {
    tagRoot.Element(XML_TAG_BONDLIKE_CALL_TRIGGERPERIOD)(m_nTriggerPeriod);
    tagRoot.Element(XML_TAG_BONDLIKE_CALL_TRIGGERHISTORY)(m_nTriggerHistory);
  }

  tagRoot.Element (XML_TAG_BONDLIKE_TAPO)
                    (
                      GetNameFromEnumValue(
                        m_triggerAsPercentageOf,
                        SIZEOF(g_triggerAsPercentageOfs),
                        g_triggerAsPercentageOfs)
                    );
}


} // namespace finance

} // namespace ito33
