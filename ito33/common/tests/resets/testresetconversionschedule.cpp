/////////////////////////////////////////////////////////////////////////////
// Name:        testresetconversionpriceschedule.cpp
// Purpose:     Acceptance test for reset conversion price schedule. 
// Author:      Yann d'Halluin
// Created:     17/11/2004
// RCS-ID:      $Id: testresetconversionschedule.cpp,v 1.7 2006/08/19 23:22:40 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------
#include "ito33/date.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/bondlike/conversionpricereset.h"
#include "ito33/finance/bondlike/resetconversionschedule.h"
#include "ito33/finance/bondlike/resetflooredby.h"

#include "ito33/tests/utilexml.h"

#include "ito33/tests/testresetconversionschedule.h"

#include "ito33/xml/write.h"

using namespace ito33;
using namespace ito33::finance;

// ----------------------------------------------------------------------------
// Conversion schedule reset
// ----------------------------------------------------------------------------

void  ConversionScheduleResetTest::StartDateBeforeEndDate()
{

  Date startDate(2003,Date::Jan,1);
  Date endDate(2002,Date::Jan,1);

  ResetFlooredBy flooredBy = ResetFlooredBy_PrevailingConversionPrice;

  double dInitialConvPrice = 1.0;
  double dCurrentConvPrice = 1.0;

  ResetConversionSchedule ResetConvSchedule(startDate, endDate, dInitialConvPrice, 
                          dCurrentConvPrice, flooredBy);

 }//ConversionScheduleResetTest::StartDateBeforeEndDate()
 

void  ConversionScheduleResetTest::CurrentRatioNegative()
{
  Date startDate(2002,Date::Jan,1);
  Date endDate(2003,Date::Jan,1);

  ResetFlooredBy flooredBy = ResetFlooredBy_PrevailingConversionPrice;

  double dInitialConvPrice = 1.0;
  double dCurrentConvPrice = -1.0;

  ResetConversionSchedule ResetConvSchedule(startDate, endDate, dInitialConvPrice, 
                          dCurrentConvPrice, flooredBy);
}//ConversionScheduleResetTest::CurrentRatioNegative()


 void ConversionScheduleResetTest::InitialRatioNegative()
 {

  Date startDate(2002,Date::Jan,1);
  Date endDate(2003,Date::Jan,1);

  ResetFlooredBy flooredBy = ResetFlooredBy_PrevailingConversionPrice;

  double dInitialConvPrice = -1.0;
  double dCurrentConvPrice = 1.0;

  ResetConversionSchedule ResetConvSchedule(startDate, endDate, dInitialConvPrice, 
                          dCurrentConvPrice, flooredBy);

 }//ConversionScheduleResetTest::InitialRatioNegative()


void ConversionScheduleResetTest::ResetDateOutsideSchedule() 
{
  double dCapRate    = 1.0;
  double dFloorRate  = 1.0;
  double dMultiplier = 1.0;
  Date resetDate(2004, Date::Mar, 1);

  shared_ptr<ConversionPriceReset> 
    pConvPriceReset( new ConversionPriceReset(resetDate,dFloorRate) );

  pConvPriceReset->SetCap(dCapRate);
  pConvPriceReset->SetMultiplier(dMultiplier);


  Date startDate(2002,Date::Jan,1);
  Date endDate(2003,Date::Jan,1);

  ResetFlooredBy flooredBy = ResetFlooredBy_PrevailingConversionPrice;

  double dInitialConvPrice = 1.0;
  double dCurrentConvPrice = 1.0;

  ResetConversionSchedule ResetConvSch(startDate, endDate, dInitialConvPrice, 
                          dCurrentConvPrice, flooredBy);

  ResetConvSch.AddConversionPriceReset(pConvPriceReset);

 } //ConversionScheduleResetTest::ResetDateOutsideSchedule()
 
 void ConversionScheduleResetTest::Dump()
 {

  
  Date startDate(2002,Date::Jan,1);
  Date endDate(2003,Date::Jan,1);

  ResetFlooredBy flooredBy = ResetFlooredBy_PrevailingConversionPrice;

  double dInitialConvPrice = 1.0;
  double dCurrentConvPrice = 1.0;

  ResetConversionSchedule ResetConvSch(startDate, 
    endDate, dInitialConvPrice, dCurrentConvPrice, flooredBy);


  std::ostringstream oss;

    ExpectedXML expected(oss,
    "<?xml version=\"1.0\"?>"
    "<root>\n"
    "<reset_conversion_schedule>\n"  
    "<keep_accrued>0</keep_accrued>\n"
    "<forfeit_coupon>0</forfeit_coupon>\n"
    "<conversion_start>2002-01-01</conversion_start>\n"
    "<conversion_end>2003-01-01</conversion_end>\n"
    "<initial_conversion_price>1</initial_conversion_price>\n"
    "<current_conversion_price>1</current_conversion_price>\n"
    "<cash_value>0</cash_value>\n"
    "<floored_by>prevailing_conversion_price</floored_by>\n"
    "<conversion_price_resets/>\n"
    "</reset_conversion_schedule>\n"
   "</root>\n");


  ito33::XML::RootTag root("root",oss);

  ResetConvSch.Dump(root);

 }//ConversionScheduleResetTest::Dump()
