/////////////////////////////////////////////////////////////////////////////
// Name:        testattachedwarrantconvertiblebond.cpp
// Purpose:     Acceptance test for attached warrant convertible bond
// Author:      ITO 33
// Created:     15/03/2005
// RCS-ID:      $Id: testattachedwarrantconvertiblebond.cpp,v 1.8 2006/08/19 23:22:40 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------
#include "ito33/date.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/bondlike/sharedependentconversion.h"
#include "ito33/finance/bondlike/attachedwarrantconvertiblebond.h"
#include "ito33/finance/bondlike/bondterms.h"

#include "ito33/tests/utilexml.h"

#include "ito33/tests/testattachedwarrantconvertiblebond.h"

#include "ito33/xml/write.h"

using namespace ito33;
using namespace ito33::finance;

// ----------------------------------------------------------------------------
// Attached warrant convertible bond
// ----------------------------------------------------------------------------

void 
 AttachedWarrantConvertibleBondTest::ConversionStartDateBeforeBondStartDate()
{
  
  Date issueDate    = Date(2005, Date::Mar, 15);
  Date maturityDate = Date(2006, Date::Mar, 15);

  double dIssuePrice     = 1.;
  double dParValue       = 110.;
  double dRedemptionRate = 1;
  double dRecoveryRate   = 0.;

  shared_ptr<BondTerms> pBondTerms( new BondTerms(issueDate, dIssuePrice,
    maturityDate, dParValue, dRedemptionRate, dRecoveryRate) );
 

  Date startConvDate = issueDate;
  startConvDate.AddDays(-1);
  Date endConvDate(2006, Date::Mar, 1);

  double dBaseRatio = 1.;
  double dIncrementalShareFactor = 1.;

  shared_ptr<ShareDependentConversion>
    pConv( new ShareDependentConversion(
      startConvDate, endConvDate, dBaseRatio, dIncrementalShareFactor) );

  shared_ptr<AttachedWarrantConvertibleBond> 
    pWarrant( new AttachedWarrantConvertibleBond(pBondTerms,pConv) );

}

void AttachedWarrantConvertibleBondTest::ConversionEndDateAfterBondEndDate()
{

  Date issueDate    = Date(2005, Date::Mar, 15);
  Date maturityDate = Date(2006, Date::Mar, 15);

  double dIssuePrice     = 1.;
  double dParValue       = 110.;
  double dRedemptionRate = 1;
  double dRecoveryRate   = 0.;

  shared_ptr<BondTerms> pBondTerms( new BondTerms(issueDate, dIssuePrice,
    maturityDate, dParValue, dRedemptionRate, dRecoveryRate) );
 

  Date startConvDate = issueDate;
  Date endConvDate   = maturityDate;
  endConvDate.AddDays(+1);

  double dBaseRatio = 1.;
  double dIncrementalShareFactor = 1.;

  shared_ptr<ShareDependentConversion>
    pConv( new ShareDependentConversion(
      startConvDate, endConvDate, dBaseRatio, dIncrementalShareFactor ) );

  shared_ptr<AttachedWarrantConvertibleBond> 
    pWarrant( new AttachedWarrantConvertibleBond(pBondTerms,pConv) );
}


void AttachedWarrantConvertibleBondTest::EmptySharedDependentConversion()
{

  Date issueDate    = Date(2005, Date::Mar, 15);
  Date maturityDate = Date(2006, Date::Mar, 15);

  double dIssuePrice     = 1.;
  double dParValue       = 110.;
  double dRedemptionRate = 1;
  double dRecoveryRate   = 0.;

  shared_ptr<BondTerms> pBondTerms( new BondTerms(issueDate, dIssuePrice,
    maturityDate, dParValue, dRedemptionRate, dRecoveryRate) );

  shared_ptr<AttachedWarrantConvertibleBond> 
    pWarrant( new AttachedWarrantConvertibleBond(pBondTerms,shared_ptr<ShareDependentConversion>() ) );
}

void AttachedWarrantConvertibleBondTest::Dump()
{
  Date issueDate    = Date(2005, Date::Mar, 15);
  Date maturityDate = Date(2006, Date::Mar, 15);

  double dIssuePrice     = 1.;
  double dParValue       = 110.;
  double dRedemptionRate = 1;
  double dRecoveryRate   = 0.;

  shared_ptr<BondTerms> pBondTerms( new BondTerms(issueDate, dIssuePrice,
    maturityDate, dParValue, dRedemptionRate, dRecoveryRate) );
 

  Date startConvDate = issueDate;
  Date endConvDate   = maturityDate;
  Date resetDate(2006, Date::Feb, 1);

  double dBaseRatio = 1.;
  double dIncrementalShareFactor = 1.;
  double dCapRatio  = 1.;

  shared_ptr<ShareDependentConversion>
    pConv( new ShareDependentConversion(
      startConvDate, endConvDate, dBaseRatio, dIncrementalShareFactor) );

  pConv->SetCapRatio(dCapRatio);

  pConv->SetResetDate( resetDate );

  shared_ptr<AttachedWarrantConvertibleBond> 
    pWarrant( new AttachedWarrantConvertibleBond(pBondTerms,pConv) );

  std::ostringstream oss;
 
  ExpectedXML expected(oss,
            "<?xml version=\"1.0\"?>"
            "<root>\n"
            "<attached_warrant_convertible_bond>\n"
            "<bond_terms>\n"
            "<issue_date>2005-03-15</issue_date>\n"
            "<issue_price>1</issue_price>\n"
            "<maturity>2006-03-15</maturity>\n"
            "<nominal>110</nominal>\n"
            "<recovery_rate>0</recovery_rate>\n"
            "<redemption_price>1</redemption_price>\n"
            "</bond_terms>\n"
            "<new_share>0</new_share>\n"
            "<share_dependent_conversion>\n"
            "<keep_accrued>0</keep_accrued>\n"
            "<forfeit_coupon>0</forfeit_coupon>\n"
            "<start_date>2005-03-15</start_date>\n"
            "<end_date>2006-03-15</end_date>\n"
            "<reset_date>2006-02-01</reset_date>\n"
            "<base_ratio>1</base_ratio>\n"
            "<cap_ratio>1</cap_ratio>\n"
            "<incremental_share_factor>1</incremental_share_factor>\n"
            "</share_dependent_conversion>\n"
            "</attached_warrant_convertible_bond>\n"
            "</root>"
         );

  ito33::XML::RootTag root("root",oss);

  pWarrant->Dump(root);

}

