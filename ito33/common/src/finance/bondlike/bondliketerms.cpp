/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/bindkind/bondliketerms.cpp
// Purpose:     base Characteristics class for Bond and Peps-like
// Author:      Wang
// Created:     2004/08/13
// RCS-ID:      $Id: bondliketerms.cpp,v 1.28 2006/08/19 22:43:51 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"


#include "ito33/finance/error.h"
#include "ito33/finance/cashflowstream.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/bondlike/bondliketerms.h"

#include "ito33/finance/bondlike/bonderror.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/bondlike/bondliketerms_all.h"
#include "ito33/xml/finance/bondlike/common.h"
#include "ito33/xml/finance/common.h"

extern const ito33::Error 
  ITO33_BAD_PARAM,
  ITO33_NULL_PARAM;

extern const ito33::finance::Error 
  ITO33_ISSUEDATE_AFTER_MATURITYDATE,
  ITO33_INVALID_RECOVERYRATE_1;

extern const ito33::finance::BondError
  ITO33_BONDLIKE_TERMS_COUPON_ISSUEDATE,
  ITO33_BONDLIKE_TERMS_ISSUEPRICE,
  ITO33_BONDLIKE_TERMS_BAD_NOMINAL;

namespace ito33
{

namespace finance
{


BondLikeTerms::BondLikeTerms(Date issueDate, double dIssuePrice, 
                             Date maturityDate, double dNominal,
                             double dRecoveryRate)
                           : m_issueDate(issueDate),
                             m_dIssuePrice(dIssuePrice), 
                             m_maturityDate(maturityDate), 
                             m_dNominal(dNominal),
                             m_dRecoveryRate(dRecoveryRate)
{
  CHECK_COND
  (
    m_issueDate < m_maturityDate,
    ITO33_ISSUEDATE_AFTER_MATURITYDATE
  );

  CHECK_COND
  ( 
    m_dIssuePrice > 0 && m_dIssuePrice <= 2, 
    ITO33_BONDLIKE_TERMS_ISSUEPRICE
  );

  CHECK_COND_1
  (
    m_dRecoveryRate >= 0 && m_dRecoveryRate <= 1,
    ITO33_INVALID_RECOVERYRATE_1,
    m_dRecoveryRate
  );

  CHECK_COND
  (
    m_dNominal > 0,
    ITO33_BONDLIKE_TERMS_BAD_NOMINAL
  );
} 

void BondLikeTerms::SetCashDistribution
     (const shared_ptr<finance::CashFlowStream>& pCashDistribution)
{
  CHECK_COND
    (
      pCashDistribution->GetContractingDate() == m_issueDate,
      ITO33_BONDLIKE_TERMS_COUPON_ISSUEDATE
    );

  m_pCashDistribution 
      = CHECK_PTR_MSG
        (
          pCashDistribution,
          ITO33_NULL_PARAM,
          "Bond-like security terms definition: Invalid cash distribution."
        );

}

ito33::XML::Tag BondLikeTerms::Dump(ito33::XML::Tag& tagParent) const
{  
  XML::Tag tagMe(XML_TAG_BONDLIKETERMS_ROOT, tagParent);
  
  BondLikeTerms::DumpMe(tagMe);

  return tagMe;
}

// this function is called by subclasses of BondLikeTerms
void BondLikeTerms::DumpMe(ito33::XML::Tag& tagParent) const
{
	tagParent.Element(XML_TAG_BONDLIKE_ISSUEDATE)( GetIssueDate() );

	tagParent.Element(XML_TAG_BONDLIKE_ISSUEPRICE)( GetIssuePrice() );
	tagParent.Element(XML_TAG_FINANCE_MATURITY)( GetMaturityDate() );
  tagParent.Element(XML_TAG_BONDLIKE_NOMINAL)(m_dNominal);

  tagParent.Element(XML_TAG_FINANCE_RECOVERYRATE)( GetRecoveryRate() );

  if (GetCashDistribution())
	  tagParent.Element(XML_TAG_BONDLIKETERMS_CASHDISTRIBUTION,
                      *GetCashDistribution());
}


} // namespace finance

} // namespace ito33
