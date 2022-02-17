/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/bondlike/pepslike.cpp
// Purpose:     financial class for PERCS-like instrument
// Author:      Wang
// Created:     2004/08/17 
// RCS-ID:      $Id: pepslike.cpp,v 1.17 2006/08/19 22:43:51 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/derivative_visitor.h"
#include "ito33/finance/bondlike/bonderror.h"

#include "ito33/finance/bondlike/pepslike.h"
#include "ito33/finance/bondlike/callfixedshare.h"
#include "ito33/finance/bondlike/callschedule.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/bondlike/pepslike.h"
#include "ito33/xml/finance/bondlike/common.h"

extern const ito33::Error ITO33_BAD_PARAM;

extern const ito33::finance::BondError 
  ITO33_PEPSLIKE_CALL_ALREADY_SET,
  ITO33_BONDLIKE_NULL_CALLSCHEDULE,
  ITO33_PEPSLIKE_CALL_FIXED_SHARE_RATIO,
  ITO33_PEPSLIKE_CALL_FIXED_SHARE_TRIGGER,
  ITO33_MANDATORY_GUARANTEED_YIELD_TO_CALL_INVALID;

namespace ito33
{

namespace finance
{

PEPSLike::PEPSLike(const shared_ptr<BondLikeTerms>& pBondLikeTerms,
          double dMaxConversionRatio,
          double dMinConversionRatio)
    : ConvertibleLike(pBondLikeTerms),
      m_dMaxConversionRatio(dMaxConversionRatio), 
      m_dMinConversionRatio(dMinConversionRatio),
      m_bHasOptionalConversion(false),
      m_bCallSet(false)
{
  CHECK_PTR_MSG
  (
    pBondLikeTerms,
    ITO33_BAD_PARAM,
    "PEPS definition: Invalid BondLikeTerms."
  );

  CHECK_COND_MSG
  (
    dMinConversionRatio > 0.,
    ITO33_BAD_PARAM,
    "PEPS definition: The minumum conversion ratio at maturity must be" 
    "a positive number."
  );

  CHECK_COND_MSG
  (
    dMaxConversionRatio >= dMinConversionRatio,
    ITO33_BAD_PARAM,
    "PEPS definition: The maximum conversion ratio at maturity must be "
    "greater or equal to the minimum conversion ratio."
  );

}

void PEPSLike::SetCallFixedCash(const shared_ptr<CallSchedule>& pCallSchedule)
{
  CHECK_COND(!m_bCallSet, ITO33_PEPSLIKE_CALL_ALREADY_SET);

  m_pCallSchedule = CHECK_PTR(pCallSchedule, ITO33_BONDLIKE_NULL_CALLSCHEDULE);

  // Check if any of the call periods have a yield to call
  CHECK_COND( !m_pCallSchedule->HasYield(), 
              ITO33_MANDATORY_GUARANTEED_YIELD_TO_CALL_INVALID );

  m_bCallSet = true;
}


void PEPSLike::SetCallFixedShare
                  (const shared_ptr<CallFixedShare>& pCallFixedShare)
{
  CHECK_COND(!m_bCallSet, ITO33_PEPSLIKE_CALL_ALREADY_SET);

  CHECK_PTR_MSG
        (
          pCallFixedShare,
          ITO33_BAD_PARAM,
          "PEPS-like definition: Invalid fixed share call."
        );

  CHECK_COND(pCallFixedShare->GetRatio() == m_dMinConversionRatio,
             ITO33_PEPSLIKE_CALL_FIXED_SHARE_RATIO);

  CHECK_COND(pCallFixedShare->GetTrigger() >= 1,
             ITO33_PEPSLIKE_CALL_FIXED_SHARE_TRIGGER);

  m_pCallFixedShare = pCallFixedShare;

  m_bCallSet = true;
}


void PEPSLike::Visit(DerivativeVisitor& visitor) const
{
  visitor.OnPEPSLike(*this);
}


XML::Tag PEPSLike::Dump(XML::Tag& tagParent) const
{
  XML::Tag tagPEPS(XML_TAG_PEPSLIKE_ROOT, tagParent);

  ConvertibleLike::DumpMe(tagPEPS);

  tagPEPS.Element(XML_TAG_PEPSLIKE_MIN_CONVERSION_RATIO)
                 (m_dMinConversionRatio);

  tagPEPS.Element(XML_TAG_PEPSLIKE_MAX_CONVERSION_RATIO)
                 (m_dMaxConversionRatio);

  // optional calls
  if ( m_pCallSchedule && !( m_pCallSchedule->GetAll().empty() ) )
  {
    tagPEPS.Element(XML_TAG_PEPSLIKE_CALL_TYPE)
                   (XML_VALUE_MANDATORY_FIXED_CASH);
    tagPEPS.Element(*m_pCallSchedule);
  }
  else if ( m_pCallFixedShare )
  {
    tagPEPS.Element(XML_TAG_PEPSLIKE_CALL_TYPE)
                   (XML_VALUE_MANDATORY_FIXED_SHARE);
    tagPEPS.Element(*m_pCallFixedShare);
  }

  // optional conversions
  tagPEPS.Element(XML_TAG_PEPSLIKE_HAS_OPTIONAL_CONVERSION)
                 (m_bHasOptionalConversion);

  return tagPEPS;
}


} // namespace finance

} // namespace ito33
