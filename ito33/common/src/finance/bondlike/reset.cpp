/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/bondlike/reset.cpp
// Purpose:     financial resettable convertible bond class
// Author:      Yann and David
// Created:     2004/10/05
// RCS-ID:      $Id: reset.cpp,v 1.27 2006/08/19 22:43:51 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/sessiondata.h"

#include "ito33/finance/bondlike/bonderror.h"
#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/resetconversionschedule.h"
#include "ito33/finance/bondlike/conversionpricereset.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/reset.h"

#include "ito33/finance/derivative_visitor.h"
#include "ito33/finance/derivative_modifyingvisitor.h"


#include "ito33/xml/write.h"
#include "ito33/xml/finance/bondlike/reset.h"

extern const ito33::finance::BondError 
  ITO33_BONDLIKE_RESET_AND_COCO,
  ITO33_BONDLIKE_RESET_AND_COCALL,
  ITO33_BONDLIKE_RESET_BEFORE_ISSUE,
  ITO33_BONDLIKE_RESET_AFTER_MATURITY,
  ITO33_BONDLIKE_RESET_DATE_AFTER_MATURITY,
  ITO33_BONDLIKE_EMPTY_RESET_DATES,
  ITO33_BONDLIKE_RESET_INCONSISTENCY_INITIAL_CURRENT_PRICE,
  ITO33_BONDLIKE_RESET_OID,
  ITO33_BONDLIKE_NULL_RESET_CONVERSION_SCHEDULE,
  ITO33_BONDLIKE_RESET_FIXEDRATE_NOT_SET;

namespace ito33
{

namespace finance
{

Reset::Reset(const shared_ptr<BondTerms>& pBondTerms,
             const shared_ptr<ResetConversionSchedule>& pResetConversionSchedule)
  : CBBase(pBondTerms),
    m_pResetConversionSchedule
    (
      CHECK_PTR
        (
          pResetConversionSchedule,
          ITO33_BONDLIKE_NULL_RESET_CONVERSION_SCHEDULE
        )
    )
{ 
  // Make sure the reset period is consistent with the issue and
  // maturity dates
  CHECK_COND
  (
    pResetConversionSchedule->GetStartDate() >= pBondTerms->GetIssueDate(),
    ITO33_BONDLIKE_RESET_BEFORE_ISSUE
  );

  CHECK_COND
  (
    pResetConversionSchedule->GetEndDate() <= pBondTerms->GetMaturityDate(),
    ITO33_BONDLIKE_RESET_AFTER_MATURITY
  );
}

void Reset::CheckTriggerPeriod() const
{
  CHECK_COND
  (
    !m_pCallSchedule || m_pCallSchedule->GetTriggerPeriod() == 0,
    ITO33_BONDLIKE_RESET_AND_COCALL
  );
}

void Reset::SetCallSchedule(const shared_ptr<CallSchedule>& pCallSchedule)
{
  CBBase::SetCallSchedule(pCallSchedule);

  CheckTriggerPeriod();
}

void Reset::Validate() const
{
  // Call function in base class
  CBBase::Validate();

  CheckTriggerPeriod();  
  
  // For now, make sure at least one reset date is specified. Otherwise,
  // a normal bond should be used for pricing
  CHECK_COND
  (
    !m_pResetConversionSchedule->GetAll().empty(),
    ITO33_BONDLIKE_EMPTY_RESET_DATES
  );

  // Cannot have a reset date on the maturity date.  Find the last
  // reset date and compare. A reset date on the valuation date is 
  // valid, but will be ignored during pricing.
  CHECK_COND
  (
    m_pResetConversionSchedule->GetAll().back()->GetDate() < GetMaturityDate(),
    ITO33_BONDLIKE_RESET_DATE_AFTER_MATURITY
  );
}

void Reset::ValidateWith(const SessionData& sessionData) const
{
  CBBase::ValidateWith(sessionData);
 
  // Fixed FX rate must be set for a cross currency reset
  if ( IsCrossCurrency(sessionData) )
    CHECK_COND( GetFixedFXRate() > 0, ITO33_BONDLIKE_RESET_FIXEDRATE_NOT_SET );

  const ResetConversionSchedule::Elements&
    resetSchedule = GetResetConversionSchedule()->GetAll();
  
  if ( sessionData.GetValuationDate() < resetSchedule.front()->GetDate() )
  {
    CHECK_COND
      (
        fabs(GetResetConversionSchedule()->GetInitialConversionPrice() -
        GetResetConversionSchedule()->GetCurrentConversionPrice()) < 1.e-4,
        ITO33_BONDLIKE_RESET_INCONSISTENCY_INITIAL_CURRENT_PRICE
      );
  }
}

void Reset::Visit(DerivativeVisitor& visitor) const
{
  visitor.OnReset(*this);
}

void Reset::Visit(DerivativeModifyingVisitor& visitor)
{
  visitor.OnReset(*this);
}

XML::Tag Reset::Dump(XML::Tag& tagParent) const
{
  XML::Tag tagReset(XML_TAG_RESET_ROOT, tagParent);

  DumpMe(tagReset);

  m_pResetConversionSchedule->Dump(tagReset);

  return tagReset;
}


} // namespace finance

} // namespace ito33
