/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/bondlike/percslike.cpp
// Purpose:     financial class for PERCS-like instrument
// Author:      Wang
// Created:     2004/08/17 
// RCS-ID:      $Id: percslike.cpp,v 1.19 2006/08/19 22:43:51 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file common/src/finance/bondlike/percslike.cpp
 */

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/derivative_visitor.h"

#include "ito33/finance/bondlike/bonderror.h"
#include "ito33/finance/bondlike/percslike.h"
#include "ito33/finance/bondlike/callschedule.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/bondlike/percslike.h"

extern const ito33::finance::BondError 
  ITO33_BONDLIKE_NULL_CALLSCHEDULE,
  ITO33_BONDLIKE_NULL_BONDLIKETERMS,
  ITO33_PERCSLIKE_CAP_PRICE,
  ITO33_PERCSLIKE_MAX_CONVERSION_RATIO,
  ITO33_MANDATORY_GUARANTEED_YIELD_TO_CALL_INVALID;

namespace ito33
{

namespace finance
{


PERCSLike::PERCSLike(const shared_ptr<BondLikeTerms>& pBondLikeTerms,
                     double dCapPrice,
                     double dMaxConversionRatio)
                   : ConvertibleLike(pBondLikeTerms),
                     m_dCapPrice(dCapPrice),
                     m_dMaxConversionRatio(dMaxConversionRatio)
{
  CHECK_PTR
  (
    pBondLikeTerms,
    ITO33_BONDLIKE_NULL_BONDLIKETERMS
  );

  CHECK_COND
  (
    dCapPrice > 0.,
    ITO33_PERCSLIKE_CAP_PRICE
  );

  CHECK_COND
  (
    dMaxConversionRatio > 0.,
    ITO33_PERCSLIKE_MAX_CONVERSION_RATIO
  );
}

void PERCSLike::SetCallSchedule(const shared_ptr<CallSchedule>& pCallSchedule) 
{ 
  m_pCall = CHECK_PTR(pCallSchedule, ITO33_BONDLIKE_NULL_CALLSCHEDULE);
 
  // Check if any of the call periods have a yield to call
  CHECK_COND( !m_pCall->HasYield(), 
              ITO33_MANDATORY_GUARANTEED_YIELD_TO_CALL_INVALID );
}

XML::Tag PERCSLike::Dump(XML::Tag& tagParent) const
{
  XML::Tag tagPERCS(XML_TAG_PERCSLIKE_ROOT, tagParent);

  ConvertibleLike::DumpMe(tagPERCS);

  tagPERCS.Element(XML_TAG_PERCSLIKE_CAP_PRICE)
                 (m_dCapPrice);

  tagPERCS.Element(XML_TAG_PERCSLIKE_CONVERSION_RATIO_AT_MATURITY)
                 (m_dMaxConversionRatio);

  // optional calls
  if ( m_pCall && !( m_pCall->GetAll().empty() ) )
    tagPERCS.Element(*m_pCall);

  return tagPERCS;
}


void PERCSLike::Visit(DerivativeVisitor& visitor) const
{
  visitor.OnPERCSLike(*this);
}


} // namespace finance

} // namespace ito33
