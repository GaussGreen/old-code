/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/bondlike/convertiblebond.cpp
// Purpose:     financial convertible bond class
// Author:      ZHANG Yunzhi
// Created:     2004 may 6
// RCS-ID:      $Id: convertiblebond.cpp,v 1.24 2006/08/19 22:43:51 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/derivative_visitor.h"
#include "ito33/finance/derivative_modifyingvisitor.h"

#include "ito33/finance/bondlike/bonderror.h"
#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/finance/bondlike/conversionschedule.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/bondlike/convertiblebond.h"

extern const ito33::Error ITO33_BAD_PARAM;

extern const ito33::finance::BondError
  ITO33_BONDLIKE_NULL_CONVERSION,
  ITO33_BONDLIKE_NO_CONVERSION;

namespace ito33
{

namespace finance
{

ConvertibleBond::ConvertibleBond
                 (const shared_ptr<BondTerms>& pBondTerms,
                  const shared_ptr<ConversionSchedule>& pConversionSchedule)
                : CBBase(pBondTerms)
{   
  CHECK_PTR(pConversionSchedule, ITO33_BONDLIKE_NULL_CONVERSION);
 
  m_pConversionSchedule = pConversionSchedule;
}

void ConvertibleBond::Validate() const
{
  CBBase::Validate();

  CHECK_COND
  (
    ! ( m_pConversionSchedule->GetAll().empty() ),
    ITO33_BONDLIKE_NO_CONVERSION
  );
}

void ConvertibleBond::ValidateWith(const SessionData& sessionData) const
{
  CBBase::ValidateWith(sessionData);

  m_pConversionSchedule->ValidateWith(sessionData.GetValuationDate());
}

void ConvertibleBond::Visit(DerivativeVisitor& visitor) const
{
  visitor.OnConvertibleBond(*this);
}

void ConvertibleBond::Visit(DerivativeModifyingVisitor& visitor)
{
  visitor.OnConvertibleBond(*this);
}

XML::Tag ConvertibleBond::Dump(XML::Tag& tagParent) const
{
  XML::Tag tagCB(XML_TAG_CONVERTIBLEBOND_ROOT, tagParent);

  DumpMe(tagCB);

  tagCB.Element(*m_pConversionSchedule);

  return tagCB;
}


} // namespace fianance

} // namespace ito33

