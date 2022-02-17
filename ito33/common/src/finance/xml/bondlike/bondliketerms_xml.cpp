/////////////////////////////////////////////////////////////////////////////
// Name:        bondliketerms_xml.cpp
// Purpose:     Restore bondliketerms objects from XML document
// Author:      ZHANG Yunzhi
// Created:     2004-09-03
// RCS-ID:      $Id: bondliketerms_xml.cpp,v 1.37 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/enum_values_names.h"

#include "ito33/finance/issuer.h"
#include "ito33/finance/bondlike/bondliketerms.h"
#include "ito33/finance/bondlike/bondterms.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/frequency.h"
#include "ito33/xml/finance/daycountconvention.h"
#include "ito33/xml/finance/cashflowstream_all.h"
#include "ito33/xml/finance/floatingrates.h"
#include "ito33/xml/finance/bondlike/bondliketerms_all.h"
#include "ito33/xml/finance/bondlike/common.h"

using namespace ito33::finance;

namespace ito33
{

namespace XML
{ 
  
/**
    Exception thrown when repayment information is not found in XML document.
 */
class MissingRepaymentsInBondTerms: public Exception
{
public:

  MissingRepaymentsInBondTerms(const char *filename,
                       size_t line,
                       const char *function)
    : Exception(ITO33_BAD_DATA,
                "Repayments/Redemption data is missed.",
                filename,
                line,
                function)
  {
  }
};

#define THROW_MISSING_REPAYMETNS_IN_BONDTERMS_EXCEPTION                       \
  throw MissingRepaymentsInBondTerms(__FILE__, __LINE__, __FUNCTION__)

/**
    Exception thrown when last repayment date (for bondTerms) is
    not maturity date.
 */
class InconsistentMaturityInBondTerms: public Exception
{
public:

  InconsistentMaturityInBondTerms(const char *filename,
                       size_t line,
                       const char *function)
    : Exception(ITO33_BAD_DATA,
                "Maturity date is not last redemption date.",
                filename,
                line,
                function)
  {
  }
};

#define THROW_INCONSISTENT_MATURITY_IN_BONDTERMS_EXCEPTION                    \
  throw InconsistentMaturityInBondTerms(__FILE__, __LINE__, __FUNCTION__)

void GetOptionalDataFromBondLikeTermsTo
          (
            const BondLikeTerms& bondLikeTerms, BondLikeTerms& terms
          )
{
  if ( bondLikeTerms.GetCashDistribution() )
    terms.SetCashDistribution(bondLikeTerms.GetCashDistribution());
}

BondLikeTerms* GetBondLikeTermsInNode(const xml::node& node)
{
  Date maturityDate
    = GetDateFromName(node, XML_TAG_FINANCE_MATURITY);

  double dRecoveryRate
    = GetDoubleFromName(node, XML_TAG_FINANCE_RECOVERYRATE);

  Date issueDate
    = GetDateFromName(node, XML_TAG_BONDLIKE_ISSUEDATE);
  
  double dIssuePrice
    = GetDoubleFromName(node, XML_TAG_BONDLIKE_ISSUEPRICE);

  double dNominal
    = GetDoubleFromName(node, XML_TAG_BONDLIKE_NOMINAL);

  BondLikeTerms 
    *pTerms = new BondLikeTerms(issueDate, dIssuePrice, maturityDate,
                                dNominal, dRecoveryRate);

  // get cash distribution, if any
  xml::node::const_iterator
    pNodeStream = node.find(XML_TAG_BONDLIKETERMS_CASHDISTRIBUTION);
  if ( pNodeStream != node.end() )
    pTerms->SetCashDistribution(
                GetCashFlowStreamInNode(*pNodeStream));

  return pTerms;
}

shared_ptr<BondLikeTerms> GetBondLikeTermsFromNode(const xml::node& node)
{
  return make_ptr( GetBondLikeTermsInNode
                  ( GetNodeByName(node, XML_TAG_BONDLIKETERMS_ROOT) ) );
}

shared_ptr<BondTerms> GetBondTermsFromNode(const xml::node& node)
{
  xml::node nodeRoot = GetNodeByName(node, XML_TAG_BONDTERMS_ROOT);

  BondLikeTerms* pBondLikeTerms = GetBondLikeTermsInNode(nodeRoot);

  double dRedemptionPrice
    = GetDoubleFromName(nodeRoot, XML_TAG_BONDTERMS_REDEMPTIONPRICE);
  
  shared_ptr<BondTerms> pTerms;

  pTerms = make_ptr(new BondTerms
                       (
                         pBondLikeTerms->GetIssueDate(),
                         pBondLikeTerms->GetIssuePrice(),
                         pBondLikeTerms->GetMaturityDate(),
                         pBondLikeTerms->GetNominal(),
                         dRedemptionPrice,
                         pBondLikeTerms->GetRecoveryRate()
                       )
                  ); 

  xml::node::const_iterator
    pnodeFound = nodeRoot.find(XML_TAG_FLOATINGRATES_ROOT);

  if ( pnodeFound != nodeRoot.end() )
    pTerms->SetFloatingRates
            ( GetFloatingRatesFromNode(*pnodeFound) );

  GetOptionalDataFromBondLikeTermsTo(*pBondLikeTerms, *pTerms);

  // Check for CashPayToZero yield
  pnodeFound = nodeRoot.find(XML_TAG_BONDTERMS_CASHPAYTOZEROACCRETIONRATE);

  if ( pnodeFound != nodeRoot.end() )
  {
    double dAccretionRateOfCashPayToZero = GetDoubleFromNode(*pnodeFound);

    pTerms->SetCashPayToZero(dAccretionRateOfCashPayToZero);
  }

  // get the compounding frequency if any
  pnodeFound = nodeRoot.find(XML_TAG_BONDTERMS_COMPOUNDING_FREQUENCY);

  if ( pnodeFound != nodeRoot.end() )
  {
    finance::Frequency
      freq = GetEnumFromNode(*pnodeFound,
                             SIZEOF(g_frequencys),
                             g_frequencys);

    pTerms->SetYieldCompoundingFrequency(freq);
  }

  // get the yield day count convention if any
  pnodeFound = nodeRoot.find(XML_TAG_BONDTERMS_YIELD_DCC);

  if ( pnodeFound != nodeRoot.end() )
  {
    Date::DayCountConvention
      dcc = GetEnumFromNode(*pnodeFound,
                            SIZEOF(g_dayCountConventions),
                            g_dayCountConventions);

    pTerms->SetYieldDayCountConvention(dcc);
  }

  // Check for yield to maturity
  pnodeFound = nodeRoot.find(XML_TAG_BONDTERMS_OIDYIELD);

  if ( pnodeFound != nodeRoot.end() )
  {
    double dYield = GetDoubleFromNode(*pnodeFound);

    pTerms->SetAccretingBond(dYield);
  }

  delete pBondLikeTerms;

  return pTerms;
}

} // namespace XML

} // namespace ito33
