/////////////////////////////////////////////////////////////////////////////
// Name:        testconversionpricereset.cpp
// Purpose:     Acceptance test conversion price resets. 
// Author:      Yann d'Halluin
// Created:     17/11/2004
// RCS-ID:      $Id: testconversionpricereset.cpp,v 1.6 2006/08/19 23:22:40 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------
#include "ito33/date.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/bondlike/conversionpricereset.h"

#include "ito33/tests/utilexml.h"


#include "ito33/tests/testconversionpricereset.h"

#include "ito33/xml/write.h"

using namespace ito33;
using namespace ito33::finance;

// ----------------------------------------------------------------------------
// Conversion Price reset
// ----------------------------------------------------------------------------

void  ConversionPriceResetTest::FloorRateTooLarge()
{

  double dFloorRate = 2.0;
  Date resetDate(2003, Date::Mar, 1);

  shared_ptr<ConversionPriceReset> 
    pConversionPriceReset( new ConversionPriceReset(resetDate, dFloorRate) );

 }//ConversionPriceResetTest::FloorRateTooLarge()
 
void  ConversionPriceResetTest::FloorRateTooSmall()
{
  double dFloorRate = -0.1;
  Date resetDate(2003, Date::Mar, 1);

  shared_ptr<ConversionPriceReset> 
    pConversionPriceReset( new ConversionPriceReset(resetDate, dFloorRate) );

 }//ConversionPriceResetTest::FloorRateTooSmall()
 

void ConversionPriceResetTest::CapRateTooLarge()
{
  double dCapRate   = 20000;
  double dFloorRate = 1.0;
  Date resetDate(2003, Date::Mar, 1);

  shared_ptr<ConversionPriceReset> 
    pConversionPriceReset( new ConversionPriceReset(resetDate, dFloorRate) );

  pConversionPriceReset->SetCap(dCapRate);

}//ConversionPriceResetTest::CapRateTooLarge()
 
void  ConversionPriceResetTest::CapRateTooSmall()
{
  double dCapRate   = 0.2;
  double dFloorRate = 1.0;
  Date resetDate(2003, Date::Mar, 1);

  shared_ptr<ConversionPriceReset> 
    pConversionPriceReset( new ConversionPriceReset(resetDate, dFloorRate) );

  pConversionPriceReset->SetCap(dCapRate);

}//ConversionPriceResetTest::CapRateTooSmall()


 void ConversionPriceResetTest::MultiplierTooLarge()
 {
  double dFloorRate  = 1.0;
  double dMultiplier = 20.0;

  Date resetDate(2003, Date::Mar, 1);

  shared_ptr<ConversionPriceReset> 
    pConversionPriceReset( new ConversionPriceReset(resetDate, dFloorRate) );

  pConversionPriceReset->SetMultiplier(dMultiplier);
 }//ConversionPriceResetTest::MultiplierTooLarge()


 void ConversionPriceResetTest::MultiplierTooSmall()
 {
  double dFloorRate  = 1.0;
  double dMultiplier = -.1;

  Date resetDate(2003, Date::Mar, 1);

  shared_ptr<ConversionPriceReset> 
    pConversionPriceReset( new ConversionPriceReset(resetDate, dFloorRate) );

  pConversionPriceReset->SetMultiplier(dMultiplier);

 }//RConversionPriceResetTest::MultiplierTooSmall()

 void ConversionPriceResetTest::Dump()
 {
  double dCapRate    = 1.0;
  double dFloorRate  = 1.0;
  double dMultiplier = 1.0;
  Date resetDate(2003, Date::Mar, 1);

  shared_ptr<ConversionPriceReset> 
    pConversionPriceReset( new ConversionPriceReset(resetDate, dFloorRate) );

  pConversionPriceReset->SetCap(dCapRate);
  pConversionPriceReset->SetMultiplier(dMultiplier);

   std::ostringstream oss;

    ExpectedXML expected(oss,
    "<?xml version=\"1.0\"?>"
    "<root>\n"
    "<conversion_price_reset>\n"
    "<date>2003-03-01</date>\n"
    "<cap_rate>1</cap_rate>\n"
    "<floor_rate>1</floor_rate>\n"
    "<multiplier>1</multiplier>\n"
    "</conversion_price_reset>\n"
   "</root>\n");

  ito33::XML::RootTag root("root",oss);

  pConversionPriceReset->Dump(root);

 }//ConversionPriceResetTest::Dump()

void ConversionPriceResetTest::CapRateDefaultToOne()
{

  double dFloorRate  = 1.0;
  Date resetDate(2003, Date::Mar, 1);

  shared_ptr<ConversionPriceReset> 
    pConversionPriceReset( new ConversionPriceReset(resetDate, dFloorRate) );

 ITO33_ASSERT_DOUBLES_EQUAL( pConversionPriceReset->GetCap(),1.0);

}//ConversionPriceResetTest::ResetCapRateDefaultToOne()

void ConversionPriceResetTest::MultiplierDefaultToOne()
{

  double dFloorRate  = 1.0;
  Date resetDate(2003, Date::Mar, 1);

  shared_ptr<ConversionPriceReset> 
    pConversionPriceReset( new ConversionPriceReset(resetDate, dFloorRate) );

  ITO33_ASSERT_DOUBLES_EQUAL( pConversionPriceReset->GetMultiplier(), 1.0);
}//ConversionPriceResetTest::ResetMultiplierDefaultToOne()

