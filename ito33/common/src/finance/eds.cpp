/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/eds.cpp
// Purpose:     Implement financial EDS class
// Created:     2005/01/26
// RCS-ID:      $Id: eds.cpp,v 1.19 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/eds.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/barrier.h"
#include "ito33/xml/finance/eds.h"

#include "ito33/finance/derivative_visitor.h"
#include "ito33/finance/derivative_modifyingvisitor.h"

extern const ito33::Error ITO33_NULL_PARAM;
extern const ito33::finance::Error
  ITO33_NO_MARKETPRICE,
  ITO33_INVALID_RECOVERYRATE_1,
  ITO33_INVALID_BARRIER;

namespace ito33
{

namespace finance
{


EDS::EDS(double dRecoveryRate, 
         const shared_ptr<CashFlowStreamUniform>& pSpreadStream,
         double dBarrier)
       : Derivative(),
         m_dRecoveryRate(dRecoveryRate),
         m_pSpreadStream(pSpreadStream),
         m_dBarrier(dBarrier)
{
  CHECK_PTR_MSG(pSpreadStream,
                ITO33_NULL_PARAM,
                "EDS definition: Setting invalid spread stream.");

  CHECK_COND_1(dRecoveryRate >= 0. && dRecoveryRate <= 1., 
             ITO33_INVALID_RECOVERYRATE_1,
             dRecoveryRate);

  CHECK_COND(dBarrier > 0, ITO33_INVALID_BARRIER);
}

Date EDS::GetMaturityDate() const
{
  return m_pSpreadStream->GetLastPaymentDate();
}

void EDS::Visit(DerivativeVisitor& visitor) const
{
  visitor.OnEDS(*this);
}

void EDS::Visit(DerivativeModifyingVisitor& visitor)
{
  visitor.OnEDS(*this);
}


XML::Tag EDS::Dump(ito33::XML::Tag& tagParent) const 
{
  XML::Tag tagEDS(XML_TAG_EDS_ROOT,tagParent);
	
  DumpMe(tagEDS);

  // dump the cash stream of payments
  tagEDS.Element(XML_TAG_EDS_SPREADSTREAM, *m_pSpreadStream);

  // dump the recovery rate
  tagEDS.Element(XML_TAG_FINANCE_RECOVERYRATE)(m_dRecoveryRate); 

  tagEDS.Element(XML_TAG_BARRIER)(m_dBarrier);

  DumpMarketPrice(tagEDS);

  return tagEDS;
}


} // namespace finance

} // namespace ito33
