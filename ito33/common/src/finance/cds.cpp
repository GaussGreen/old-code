/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/cds.cpp
// Purpose:     implement financial cds class
// Created:     2004/03/01
// RCS-ID:      $Id: cds.cpp,v 1.34 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2003 - 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/cds.h"
#include "ito33/finance/cashflowstream_uniform.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/derivative.h"
#include "ito33/xml/finance/cds.h"

#include "ito33/finance/derivative_visitor.h"
#include "ito33/finance/derivative_modifyingvisitor.h"


extern const ito33::finance::Error ITO33_CDS_INVALID_SPREADSTREAM;

namespace ito33
{

namespace finance
{

CDS::CDS(double dRecoveryRate,
         const shared_ptr<CashFlowStreamUniform>& pSpreadStream) 
       : CDSLike(dRecoveryRate)
{
  // Check and set the spread stream in the base class
  CHECK_PTR(pSpreadStream, ITO33_CDS_INVALID_SPREADSTREAM);

  m_pSpreadStream = pSpreadStream;
}


double CDS::GetSpread() const
{
  return m_pSpreadStream->GetAnnualPaymentAmount();
}


Date CDS::GetMaturityDate() const
{
  return m_pSpreadStream->GetLastPaymentDate();
}


void CDS::Visit(DerivativeVisitor& visitor) const
{
  visitor.OnCDS(*this);
}


void CDS::Visit(DerivativeModifyingVisitor& visitor)
{
  visitor.OnCDS(*this);
}


XML::Tag CDS::Dump(ito33::XML::Tag& tagParent) const 
{
  XML::Tag tagCDS(XML_TAG_CDS_ROOT, tagParent);
	
  DumpMe(tagCDS);

  // dump the cash stream of payments
  tagCDS.Element(XML_TAG_CDS_SPREADSTREAM, *m_pSpreadStream);

  // dump the recovery rate
  tagCDS.Element(XML_TAG_FINANCE_RECOVERYRATE)(m_dRecoveryRate); 

  // dump the market price, if set
  DumpMarketPrice(tagCDS);

  return tagCDS;
}


} // namespace finance

} // namespace ito33
