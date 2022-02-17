/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/bondlike/attachedwarrantconvertiblebond.cpp
// Purpose:     financial convertible bond with warrant class
// Author:      Ito 33
// Created:     2004/12/31
// RCS-ID:      $Id: attachedwarrantconvertiblebond.cpp,v 1.10 2006/08/19 22:43:51 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/sessiondata.h"

#include "ito33/finance/bondlike/bonderror.h"
#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/sharedependentconversion.h"
#include "ito33/finance/bondlike/attachedwarrantconvertiblebond.h"

#include "ito33/finance/derivative_visitor.h"
#include "ito33/finance/derivative_modifyingvisitor.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/bondlike/attachedwarrantconvertiblebond.h"

extern const ito33::finance::BondError 
  ITO33_AWCB_RESET_BEFORE_ISSUE,
  ITO33_AWCB_RESET_AFTER_MATURITY,
  ITO33_AWCB_CROSSCURRENCY_UNDEFINED_STRIKE,
  ITO33_AWCB_NULL_CONVERSION,
  ITO33_AWCB_INVALID_CALL_NOTICE,
  ITO33_AWCB_NO_CURRENT_CONVERSION_RATIO_SET;

namespace ito33
{

namespace finance
{


AttachedWarrantConvertibleBond::AttachedWarrantConvertibleBond(
    const shared_ptr<BondTerms>& pBondTerms,
    const shared_ptr<ShareDependentConversion>& pShareDepConversion)
  : CBBase(pBondTerms),
    m_pShareDepConversion
    (
      CHECK_PTR
      (
       pShareDepConversion,
       ITO33_AWCB_NULL_CONVERSION
       )
    )
{ 
  // Make sure the conversion specification is consistent with the issue and
  // maturity dates
  CHECK_COND
  (
    pShareDepConversion->GetStartDate() >= pBondTerms->GetIssueDate(),
    ITO33_AWCB_RESET_BEFORE_ISSUE
  );

  CHECK_COND
  (
    pShareDepConversion->GetEndDate() <= pBondTerms->GetMaturityDate(),
    ITO33_AWCB_RESET_AFTER_MATURITY
  );
}

void AttachedWarrantConvertibleBond::ValidateWith
     (const SessionData& sessionData) const
{
  CBBase::ValidateWith(sessionData);

  // If fixed strike is not given and in case of cross currency, 
  // the computed conversion price as seen in the code 
  // is not in the same unit of the underlying currency.
  // In this case we force the user to enter a fixed strike.
  if ( m_pShareDepConversion->GetFixedStrike() < 0. )
    CHECK_COND( !IsCrossCurrency(sessionData), 
                ITO33_AWCB_CROSSCURRENCY_UNDEFINED_STRIKE );

  // Do not allow call notice for contracts with a reset date
  // where the valuation date is before the reset date
  if ( m_pShareDepConversion->HasResetDate() )
  {
    CHECK_COND( !GetCallSchedule() 
                || ! GetCallSchedule()->HasNoticePeriod(),
                ITO33_AWCB_INVALID_CALL_NOTICE);

    // The current conversion ratio must be set if
    // the valuation date is after the reset date
    if (    GetSessionData()->GetValuationDate()
          > m_pShareDepConversion->GetResetDate() )
    {
      CHECK_COND( m_pShareDepConversion->GetCurrentRatio() > 0.,
                  ITO33_AWCB_NO_CURRENT_CONVERSION_RATIO_SET);
    }
  }
}

void AttachedWarrantConvertibleBond::Visit(DerivativeVisitor& visitor) const
{
  visitor.OnAttachedWarrantConvertibleBond(*this);
}

void AttachedWarrantConvertibleBond::Visit(DerivativeModifyingVisitor& visitor)
{
  visitor.OnAttachedWarrantConvertibleBond(*this);
}

XML::Tag AttachedWarrantConvertibleBond::Dump(XML::Tag& tagParent) const
{
  XML::Tag tagCB(XML_TAG_ATTACHEDWARRANTCONVERTIBLEBOND_ROOT, tagParent);

  // dump the base class parameters
  DumpMe(tagCB);

  m_pShareDepConversion->Dump(tagCB);

  return tagCB;
}


} // namespace finance

} // namespace ito33
