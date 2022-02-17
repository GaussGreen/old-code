/////////////////////////////////////////////////////////////////////////////
// Name:        testsharedependentconversion.cpp
// Purpose:     Acceptance test for shared dependent conversion 
// Author:      ITO 33
// Created:     14/03/2005
// RCS-ID:      $Id: testsharedependentconversion.cpp,v 1.4 2006/08/19 23:22:41 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------
#include "ito33/date.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/bondlike/sharedependentconversion.h"
#include "ito33/finance/bondlike/cocotype.h"

#include "ito33/tests/utilexml.h"


#include "ito33/tests/testsharedependentconversion.h"

#include "ito33/xml/write.h"

using namespace ito33;
using namespace ito33::finance;

// ----------------------------------------------------------------------------
// Shared dependent conversion
// ----------------------------------------------------------------------------

void ShareDependentConversionTest::EndDateBeforeStartDate()
{
  Date startDate(2005, Date::Mar, 1);
  Date endDate(2005, Date::Feb, 1);

  double dBaseRatio = 1.;
  double dIncrementalShareFactor = 1.;

  shared_ptr<ShareDependentConversion>
    pConv( new finance::ShareDependentConversion
               (startDate, endDate, dBaseRatio, dIncrementalShareFactor) );

}

void ShareDependentConversionTest::BaseRatioNegative()
{

  Date startDate(2005, Date::Jan, 1);
  Date endDate(2006, Date::Dec, 1);

  double dBaseRatio = -1.;
  double dIncrementalShareFactor = 1.;

  shared_ptr<ShareDependentConversion>
    pConv( new finance::ShareDependentConversion
               ( startDate, endDate, dBaseRatio, dIncrementalShareFactor) );
}
 
void ShareDependentConversionTest::ResetDateBeforeStartDate()
{
  Date startDate(2005, Date::Jan, 1);
  Date endDate(2006, Date::Dec, 1);
  Date resetDate = startDate;
  resetDate.AddDays(-1);

  double dBaseRatio = 1.;
  double dIncrementalShareFactor = 1.;

  shared_ptr<ShareDependentConversion>
    pConv( new finance::ShareDependentConversion
               ( startDate, endDate,  dBaseRatio, dIncrementalShareFactor ) );

  pConv->SetResetDate(resetDate);
}
 
void ShareDependentConversionTest::ResetDateAfterEndDate()
{
  Date startDate(2005, Date::Jan, 1);
  Date endDate(2006, Date::Dec, 1);
  Date resetDate = endDate;
  resetDate.AddDays(1);

  double dBaseRatio = 1.;
  double dIncrementalShareFactor = 1.;
 
  shared_ptr<ShareDependentConversion>
    pConv( new finance::ShareDependentConversion
               ( startDate, endDate, dBaseRatio, dIncrementalShareFactor ) );

  pConv->SetResetDate(resetDate);

}
 
void ShareDependentConversionTest::IncrementalShareFactorNegative()
{
  Date startDate(2005, Date::Jan, 1);
  Date endDate(2006, Date::Dec, 1);

  double dBaseRatio = 1.;
  double dIncrementalShareFactor = -1.;

  shared_ptr<ShareDependentConversion>
    pConv( new finance::ShareDependentConversion
               ( startDate, endDate, dBaseRatio, dIncrementalShareFactor ) );

}
 
void ShareDependentConversionTest::CapRatioNegative()
{
  Date startDate(2005, Date::Jan, 1);
  Date endDate(2006, Date::Dec, 1);
 
  double dBaseRatio = 1.;
  double dIncrementalShareFactor = 1.;
  double dCapRatio  = -1.;

  shared_ptr<ShareDependentConversion>
    pConv( new finance::ShareDependentConversion
               ( startDate, endDate, dBaseRatio, dIncrementalShareFactor ) );

  pConv->SetCapRatio(dCapRatio);
}

void ShareDependentConversionTest::CapRatioLessThanBaseRatio()
{
  Date startDate(2005, Date::Jan, 1);
  Date endDate(2006, Date::Dec, 1);
  
  double dBaseRatio = 10.;
  double dIncrementalShareFactor = 1.;
  double dCapRatio  = 1.;

  shared_ptr<ShareDependentConversion>
    pConv( new finance::ShareDependentConversion
               ( startDate, endDate, dBaseRatio, dIncrementalShareFactor ) );

  pConv->SetCapRatio(dCapRatio);
}

void ShareDependentConversionTest::FixedStrikeNegative()
{
  Date startDate(2005, Date::Jan, 1);
  Date endDate(2006, Date::Dec, 1);

  double dBaseRatio = 1.;
  double dIncrementalShareFactor = 1.;

  shared_ptr<ShareDependentConversion>
    pConv( new finance::ShareDependentConversion
               ( startDate, endDate, dBaseRatio, dIncrementalShareFactor) );

  pConv->SetFixedStrike(-1.);
}
 
void ShareDependentConversionTest::TriggerRateNegative()
{
  Date startDate(2005, Date::Jan, 1);
  Date endDate(2006, Date::Dec, 1);
 
  double dBaseRatio = 1.;
  double dIncrementalShareFactor = 1.;

  shared_ptr<ShareDependentConversion>
    pConv( new finance::ShareDependentConversion
               ( startDate, endDate, dBaseRatio, dIncrementalShareFactor ) );

  double dTrigger = -1.;
  double dExtremeTrigger = 1.;
  double dChangeRate = 0.;
  CoCoType coCoType = CoCoType_CheckAnyTimeAndConvertAsOfCheckDate;

  pConv->SetCoCo( dTrigger, coCoType, dChangeRate, dExtremeTrigger);

}
 
void ShareDependentConversionTest::ExtremeTriggerRateNegative()
{
  Date startDate(2005, Date::Jan, 1);
  Date endDate(2006, Date::Dec, 1);

  double dBaseRatio = 1.;
  double dIncrementalShareFactor = 1.;

  shared_ptr<ShareDependentConversion>
    pConv( new finance::ShareDependentConversion
               ( startDate, endDate, dBaseRatio, dIncrementalShareFactor ) );

  double dTrigger = 1.;
  double dExtremeTrigger = -1.;
  double dChangeRate = 0.;
  CoCoType coCoType = CoCoType_CheckAnyTimeAndConvertAsOfCheckDate;

  pConv->SetCoCo( dTrigger, coCoType, dChangeRate, dExtremeTrigger);
}
 
void ShareDependentConversionTest::ChangeRatePositive()
{
  /* 
     Change positive, i.e. dTriggerRate has to be less
     than extreme trigger rate
  */
  Date startDate(2005, Date::Jan, 1);
  Date endDate(2006, Date::Dec, 1);
 
  double dBaseRatio = 1.;
  double dIncrementalShareFactor = 1.;
 
  shared_ptr<ShareDependentConversion>
    pConv( new finance::ShareDependentConversion
               ( startDate, endDate, dBaseRatio, dIncrementalShareFactor ) );

  double dTrigger = 10.;
  double dExtremeTrigger = 1.;
  double dChangeRate = 1.;
  CoCoType coCoType = CoCoType_CheckAnyTimeAndConvertAsOfCheckDate;

  pConv->SetCoCo( dTrigger, coCoType, dChangeRate, dExtremeTrigger);
}
 
void ShareDependentConversionTest::ChangeRateNegative()
{
  /*
    Change rate negative. trigger rate has to be greater
    than extreme trigger rate
  */
  Date startDate(2005, Date::Jan, 1);
  Date endDate(2006, Date::Dec, 1);

  double dBaseRatio = 1.;
  double dIncrementalShareFactor = 1.;

  shared_ptr<ShareDependentConversion>
    pConv( new finance::ShareDependentConversion
               ( startDate, endDate, dBaseRatio, dIncrementalShareFactor ) );

  double dTrigger = 1.;
  double dExtremeTrigger = 10.;
  double dChangeRate = -1.;
  CoCoType coCoType = CoCoType_CheckAnyTimeAndConvertAsOfCheckDate;

  pConv->SetCoCo( dTrigger, coCoType, dChangeRate, dExtremeTrigger);
}
 
void ShareDependentConversionTest::Dump()
{

  Date startDate(2005, Date::Jan, 1);
  Date endDate(2006, Date::Dec, 1);
  Date resetDate(2005, Date::May, 12);

  double dBaseRatio = 1.;
  double dIncrementalShareFactor = 1.;
  double dCapRatio  = 1.;

  shared_ptr<ShareDependentConversion>
    pConv( new finance::ShareDependentConversion
               ( startDate, endDate, dBaseRatio, dIncrementalShareFactor ) );

  pConv->SetCapRatio(dCapRatio);

  pConv->SetResetDate(resetDate);

  double dTrigger = 1.0;
  double dExtremeTrigger = 1.0;
  double dChangeRate = 0.0;
  CoCoType coCoType = CoCoType_CheckAnyTimeAndConvertAsOfCheckDate;

  pConv->SetCoCo( dTrigger, coCoType, dChangeRate, dExtremeTrigger);

  std::ostringstream oss;

  ExpectedXML expected(oss,
    "<?xml version=\"1.0\"?>"
    "<root>\n"
    "<share_dependent_conversion>\n"
    "<keep_accrued>0</keep_accrued>\n"
    "<forfeit_coupon>0</forfeit_coupon>\n"
    "<start_date>2005-01-01</start_date>\n"
    "<end_date>2006-12-01</end_date>\n"
    "<reset_date>2005-05-12</reset_date>\n"
    "<base_ratio>1</base_ratio>\n"
    "<cap_ratio>1</cap_ratio>\n"
    "<incremental_share_factor>1</incremental_share_factor>\n"
    "<trigger_rate>1</trigger_rate>\n"
    "<change_rate_in_trigger>0</change_rate_in_trigger>\n"
    "<extreme_trigger_rate>1</extreme_trigger_rate>\n"
    "<coco_type>check_any_time_and_convert_as_of_check_date</coco_type>\n"
    "</share_dependent_conversion>\n"
    "</root>\n");


  ito33::XML::RootTag root("root",oss);

  pConv->Dump(root);
}
