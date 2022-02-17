/////////////////////////////////////////////////////////////////////////////
// Name:        tests/bsutil_cppunit/testBSut.cpp
// Purpose:     test file for black-scholes utility function
// Author:      Laurence
// Created:     22/09/2003
// RCS-ID:      $Id: testBSut.cpp,v 1.14 2004/10/05 09:13:49 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/numeric/schemetype.h"
#include "ito33/common.h"
#include "ito33/exception.h"

#include "ito33/finance/blackscholesutils.h"
#include "ito33/finance/blackscholes.h"
#include "ito33/cppunit.h"

#include "ito33/tests/testBSut.h"

#include <cmath>

using namespace ito33;
using namespace ito33::finance;

// Results turn out to fit maples with 1e-4 precision because of
// the approximation of the cumulated normal. This is due to the 
// values choice on Spot and Strike that are aroud 1e^2

//--------------------------------------------------------------------------
//         tests for function CalculateCallPrice
//--------------------------------------------------------------------------            

//call price for variable spot
void BSutileTest::CallPriceFS()
{
  // Verifies results of function CalculateCallPrice on some values computed with
  // Mapple

  double dSpot = 100.;
  double dStrike = 101.;
  double dRate = 0.04;
  double dFRate = 0.05;
  double dVol = 0.2;
  double dTminust = 1.;
  double dPriceSpot[12];
  
  dPriceSpot[0] =   6.2667538588129;
  dPriceSpot[1] =   6.7329458516061;
  dPriceSpot[2] =   7.2181084999532;
  dPriceSpot[3] =   7.7220306772342;
  dPriceSpot[4] =   8.2444605318868;
  dPriceSpot[5] =   8.7851083567527;
  dPriceSpot[6] =   9.3436497193865;
  dPriceSpot[7] =   9.9197286056206;
  dPriceSpot[8] =  10.5129606099219;
  dPriceSpot[9] =  11.1229361291230;
  dPriceSpot[10] = 11.7492235290028;
  dPriceSpot[11] = 12.3913722562394;

  
    // Different values of price for variable Spot
  for (int i = 0; i < 12; i++)
  CPPUNIT_ASSERT_DOUBLES_EQUAL( CalculateCallPrice(dSpot-1+i, dStrike, dVol,
    dTminust, dRate, dFRate) , dPriceSpot[i] , 1e-4 );

} // CallPriceFS


//call price for variable strike
void BSutileTest::CallPriceFE()
{
  // Verifies results of function CalculateCallPrice on some values computed with
  // Mapple

  double dSpot = 100.;
  double dStrike = 101.;
  double dRate = 0.04;
  double dFRate = 0.05;
  double dVol = 0.2;
  double dTminust = 1.;
  double dPriceStrike[12];

  dPriceStrike[0] =  7.5792856779446; 
  dPriceStrike[1] =  7.1466420693023;  
  dPriceStrike[2] =  6.7329458516061;
  dPriceStrike[3] =  6.3378471589460; 
  dPriceStrike[4] =  5.9609614762926;
  dPriceStrike[5] =  5.6018732733383;  
  dPriceStrike[6] =  5.2601396385689;  
  dPriceStrike[7] =  4.9352938770853;  
  dPriceStrike[8] =  4.6268490395374;  
  dPriceStrike[9] =  4.3343013534503;  
  dPriceStrike[10] = 4.0571335321544;  
  dPriceStrike[11] = 3.7948179403859;

    // Different values of price for variable Strike
  for (int j = 0; j < 12; j++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL( CalculateCallPrice(dSpot, dStrike-2+j, dVol,
    dTminust, dRate, dFRate) , dPriceStrike[j] , 1e-4 );

}//CallPriceFE 

//call price for variable rate
void BSutileTest::CallPriceFR()
{
  // Verifies results of function CalculateCallPrice on some values computed with
  // Mapple

  double dSpot = 100.;
  double dStrike = 101.;
  double dFRate = 0.05;
  double dVol = 0.2;
  double dTminust = 1.;
  double dPriceRate[21];

  dPriceRate[0] =  6.3320345363910; 
  dPriceRate[1] =  6.3714518134498;
  dPriceRate[2] =  6.4110191507084;
  dPriceRate[3] =  6.4507364382129;
  dPriceRate[4] =  6.4906035613861;
  dPriceRate[5] =  6.5306204010320;
  dPriceRate[6] =  6.5707868333428;
  dPriceRate[7] =  6.6111027298918;
  dPriceRate[8] =  6.6515679576461;
  dPriceRate[9] =  6.6921823789653;
  dPriceRate[10] = 6.7329458516061;
  dPriceRate[11] = 6.7738582287263;
  dPriceRate[12] = 6.8149191935889;
  dPriceRate[13] = 6.8561290860695;
  dPriceRate[14] = 6.8974872496559;
  dPriceRate[15] = 6.9389936844590;
  dPriceRate[16] = 6.9806482207164;
  dPriceRate[17] = 7.0224506840985;
  dPriceRate[18] = 7.0644008957153;
  dPriceRate[19] = 7.1064986721234;
  dPriceRate[20] = 7.1487438253332;

    // Different values of price for variable Rate
  for (int k = 0; k < 21; k++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL( CalculateCallPrice(dSpot, dStrike, dVol,
    dTminust, 0.03+k/1000., dFRate) , dPriceRate[k] , 1e-4 );

}//CallPriceFR 

void BSutileTest::CallPriceSpot0()
{
  double dSpot = 1.;
  double dStrike = 101.;
  double dRate = 0.04;
  double dFRate = 0.05;
  double dVol = 0.2;
  double dTminust = 1.;

  CPPUNIT_ASSERT_DOUBLES_EQUAL( CalculateCallPrice(dSpot, dStrike, dVol,
    dTminust, dRate, dFRate) , 0. , 1e-7 );
}//CallPriceSpot0


//Check the behaviour of CalculateCallPrice for rate = 0
void BSutileTest::CallPriceRate0()
{
  double dSpot = 101.;
  double dStrike = 101.;
  double dRate = 0.0;
  double dFRate = 0.0;
  double dVol = 0.2;
  double dTminust = 1.;

  CPPUNIT_ASSERT_DOUBLES_EQUAL( CalculateCallPrice(dSpot, dStrike, dVol,
    dTminust, dRate, dFRate) , 8.045223129959854256 , 1e-4 );

}//CallPriceRate0

//Check the behaviour of CalculateCallPrice for spot groing towards infinity
void BSutileTest::CallPriceInfty()
{
  double dSpot = 400.;
  double dStrike = 100.;
  double dRate = 0.0;
  double dFRate = 0.0;
  double dVol = 0.1;
  double dTminust = 1.;
  double dSminusE = dSpot - dStrike;

  CPPUNIT_ASSERT_DOUBLES_EQUAL( CalculateCallPrice(dSpot, dStrike, dVol,
    dTminust, dRate, dFRate) , dSminusE , 1e-7 );

}//CallPriceInfty

//--------------------------------------------------------------------------
//         tests for function CalculatePutPrice
//--------------------------------------------------------------------------  


//put price for variable spot
void BSutileTest::PutPriceFS()
{
  // Verifies results of function CalculatePutPrice on some values computed with
  // Octave

  double dSpot = 100.;
  double dStrike = 101.;
  double dRate = 0.04;
  double dFRate = 0.05;
  double dVol = 0.2;
  double dTminust = 1.;
  double dPriceSpot[12];

  dPriceSpot[0] =  9.13477418762685;
  dPriceSpot[1] =  8.64973675591936;
  dPriceSpot[2] =  8.18366996980784;
  dPriceSpot[3] =  7.73636273254613;
  dPriceSpot[4] =  7.30756316269786;
  dPriceSpot[5] =  6.89698156306305;
  dPriceSpot[6] =  6.50429350119618;
  dPriceSpot[7] =  6.12914296292954;
  dPriceSpot[8] =  5.77114554273016;
  dPriceSpot[9] =  5.42989163743151;
  dPriceSpot[10] = 5.1049496128097;
  dPriceSpot[11] = 4.7958689155455;

    // Different values of price for variable Spot
  for (int i = 0; i < 12; i++)
  CPPUNIT_ASSERT_DOUBLES_EQUAL( CalculatePutPrice(dSpot-1+i, dStrike, dVol,
    dTminust, dRate, dFRate) , dPriceSpot[i] , 1e-4 );

} // PutPriceFS


//put price for variable strike
void BSutileTest::PutPriceFE()
{
// Verifies results of function CalculatePutPrice on some values computed with
// Mapple

  double dSpot = 100.;
  double dStrike = 101.;
  double dRate = 0.04;
  double dFRate = 0.05;
  double dVol = 0.2;
  double dTminust = 1.;
  double dPriceStrike[12];

  dPriceStrike[0] = 7.57449770395316;
  dPriceStrike[1] = 8.10264353446321;
  dPriceStrike[2] = 8.64973675591936;
  dPriceStrike[3] = 9.21542750241159;
  dPriceStrike[4] = 9.7993312589105;
  dPriceStrike[5] = 10.4010324951085;
  dPriceStrike[6] = 11.0200882994914;
  dPriceStrike[7] = 11.6560319771602;
  dPriceStrike[8] = 12.3083765787646;
  dPriceStrike[9] = 12.9766183318298;
  dPriceStrike[10] = 13.6602399496862;
  dPriceStrike[11] = 14.35871379707;



// Different values of price for variable Strike
  for (int j = 0; j < 12; j++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL( CalculatePutPrice(dSpot, dStrike-2+j, dVol,
    dTminust, dRate, dFRate) , dPriceStrike[j] , 1e-4 );

}//PutPriceFE 

//put price for variable rate
void BSutileTest::PutPriceFR()
{
  // Verifies results of function CalculatePutPrice on some values computed with
  // Mapple

  double dSpot = 100.;
  double dStrike = 101.;
  double dFRate = 0.05;
  double dVol = 0.2;
  double dTminust = 1.;
  double dPriceRate[21];

  dPriceRate[0]  = 9.22409097472183;
  dPriceRate[1]  = 9.165542244057;
  dPriceRate[2]  = 9.10724149063599;
  dPriceRate[3]  = 9.04918850663371;
  dPriceRate[4]  = 8.99138307970284;
  dPriceRate[5]  = 8.9338249929758;
  dPriceRate[6]  = 8.87651402506685;
  dPriceRate[7]  = 8.8194499500748;
  dPriceRate[8]  = 8.76263253758587;
  dPriceRate[9]  = 8.70606155267696;
  dPriceRate[10] = 8.64973675591936;
  dPriceRate[11] = 8.59365790338263;
  dPriceRate[12] = 8.53782474663905;
  dPriceRate[13] = 8.48223703276823;
  dPriceRate[14] = 8.42689450436217;
  dPriceRate[15] = 8.37179689953071;
  dPriceRate[16] = 8.31694395190718;
  dPriceRate[17] = 8.26233539065455;
  dPriceRate[18] = 8.20797094047187;
  dPriceRate[19] = 8.15385032160102;
  dPriceRate[20] = 8.09997324983389;


    // Different values of price for variable Rate
  for (int k = 0; k < 21; k++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL( CalculatePutPrice(dSpot, dStrike, dVol,
    dTminust, 0.03+k/1000., dFRate) , dPriceRate[k] , 1e-4 );

}//PutPriceFR 


//Chek the behaviour of CalculatePutPrice for values of the spot in 
//the vicinity of 0 it should be E-S for 0 rate
void BSutileTest::PutPriceSpot0()
{
  double dSpot = 1.;
  double dStrike = 101.;
  double dRate = 0.0;
  double dFRate = 0.0;
  double dVol = 0.2;
  double dTminust = 1.;

  CPPUNIT_ASSERT_DOUBLES_EQUAL( CalculatePutPrice(dSpot, dStrike, dVol,
    dTminust, dRate, dFRate) , 100. , 1e-7 );
}



void BSutileTest::PutPriceRate0()
{
  double dSpot = 101.;
  double dStrike = 101.;
  double dRate = 0.0;
  double dFRate = 0.0;
  double dVol = 0.2;
  double dTminust = 1.;

  CPPUNIT_ASSERT_DOUBLES_EQUAL( CalculatePutPrice(dSpot, dStrike, dVol,
    dTminust, dRate, dFRate) , 8.04522312995985 , 1e-4 );

}

//Check the behaviour of CalculatePutPrice for spot groing towards infinity
void BSutileTest::PutPriceInfty()
{
  double dSpot = 400.;
  double dStrike = 100.;
  double dRate = 0.0;
  double dFRate = 0.0;
  double dVol = 0.1;
  double dTminust = 1.;
//  double dSminusE = dSpot - dStrike;

  CPPUNIT_ASSERT_DOUBLES_EQUAL( CalculatePutPrice(dSpot, dStrike, dVol,
    dTminust, dRate, dFRate) , 0. , 1e-7 );

}//PutPriceInfty

//--------------------------------------------------------------------------
//         tests for function CalculateVegaCall
//--------------------------------------------------------------------------            

void BSutileTest::VerifVega()
{
//this function verifies that the implemented Vega function is indeed the
//derivative of the calculatecallprice function with respect to the 
//volatility
  double dSpot = 99.;
  double dStrike = 101.;
  double dRate = 0.04;
  double dFRate = 0.05;
  double dVol = 1.;
  double dTminust = 1.;

  double
  dPriceP,
  dPriceM,
  dDeriv,
  dH=1e-3;

  dPriceP=CalculateCallPrice(dSpot, dStrike, dVol+dH, dTminust, dRate, dFRate);
  dPriceM=CalculateCallPrice(dSpot, dStrike, dVol-dH, dTminust, dRate, dFRate);
  dDeriv=(dPriceP-dPriceM)/(2*dH);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( CalculateVegaCall(dSpot, dStrike, 
     dTminust, dRate, dFRate, dVol), dDeriv , 1e-3 );

//no better than 1e-3 have we been able to obtain precision for these values.
//This is most probably due to the lack of precision in the computation of the
//"cumulatednormal"
}//VerifVega

void BSutileTest::VegaCall()
{
  double dSpot = 100.;
  double dStrike = 101.;
  double dRate = 0.04;
  double dFRate = 0.05;
  double dVol = 0.2;
  double dTminust = 1.;
  double dVega[12];

  dVega[0] =  37.5221396806431;
  dVega[1] =  37.9485624092755;
  dVega[2] =  38.2801690851176;
  dVega[3] =  38.5173141797042;
  dVega[4] =  38.6610326811679;
  dVega[5] =  38.7129949061797;
  dVega[6] =  38.6754565229442;
  dVega[7] =  38.5512049072491;
  dVega[8] =  38.3435029101557;
  dVega[9] =  38.0560310560297;
  dVega[10] = 37.6928291165853;
  dVega[11] = 37.2582379236555;

// Different values of vega for variable Spot
  for (int i = 0; i < 12; i++)
  CPPUNIT_ASSERT_DOUBLES_EQUAL( CalculateVegaCall(dSpot-1+i, dStrike, 
    dTminust, dRate, dFRate, dVol) , dVega[i] , 1e-8 );

//Here we can have greater precision thanks to the exact calculation of
//the normal

} // VegaCall

//--------------------------------------------------------------------------
//         tests for function CalculateImpliedVolatility
//--------------------------------------------------------------------------            

void BSutileTest::ImpliedVol()
{
  double dSpot = 100.;
  double dStrike = 101.;
  double dRate = 0.04;
  double dFRate = 0.05;
  double dVol = 0.2;
  double dTminust = 1.;
  double dPriceSpot[12];
  double dPriceVol[12];
  double dTVol[12] = {.005, .01, .1, .101, .3, .4, .45, .6, 1., 2., 4.33, 9.};
  
  dPriceSpot[0] =   6.2667538588129;
  dPriceSpot[1] =   6.7329458516061;
  dPriceSpot[2] =   7.2181084999532;
  dPriceSpot[3] =   7.7220306772342;
  dPriceSpot[4] =   8.2444605318868;
  dPriceSpot[5] =   8.7851083567527;
  dPriceSpot[6] =   9.3436497193865;
  dPriceSpot[7] =   9.9197286056206;
  dPriceSpot[8] =  10.5129606099219;
  dPriceSpot[9] =  11.1229361291230;
  dPriceSpot[10] = 11.7492235290028;
  dPriceSpot[11] = 12.3913722562394;

    
    // Different values of price for variable Spot
  double dRes;
  for (int i = 0; i < 12; i++)
  {
    dRes = CalculateImpliedVolatility(dSpot-1+i, dStrike, 
    dPriceSpot[i], dTminust, dRate, dFRate);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( dRes , dVol , 1e-4 );
  }

  
  for (int j = 0; j < 12; j++)
  {
    dPriceVol[j] = CalculateCallPrice(dSpot, dStrike, dTVol[j], dTminust, 
    dRate, dFRate);
    dRes = CalculateImpliedVolatility(dSpot, dStrike, 
    dPriceVol[j], dTminust, dRate, dFRate);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( dRes, dTVol[j] , 1e-4 );
  }
  

}//impliedvol


//--------------------------------------------------------------------------
//         tests for function CalculateDelta
//--------------------------------------------------------------------------            

void BSutileTest::Delta()
{
  double dSpot = 100.;
  double dStrike = 101.;
  double dRate = 0.04;
  double dFRate = 0.05;
  double dVol = 0.2;
  double dTminust = 1.;
  double dDeltaSpot[12];
  double dDeltaSpotPut[12];

  dDeltaSpot[0] = 0.456647069999871;  
  dDeltaSpot[1] = 0.475708955888213;  
  dDeltaSpot[2] = 0.494581091053224; 
  dDeltaSpot[3] = 0.513221253606381;  
  dDeltaSpot[4] = 0.531590121397856; 
  dDeltaSpot[5] = 0.549651403614733;
  dDeltaSpot[6] = 0.567371932705965;  
  dDeltaSpot[7] = 0.584721718871737;  
  dDeltaSpot[8] = 0.601673969694852;
  dDeltaSpot[9] = 0.618205077748027;
  dDeltaSpot[10] = 0.634294579189002;
  dDeltaSpot[11] = 0.649925086461712;

  dDeltaSpotPut[0] = -0.494582354500843;
  dDeltaSpotPut[1] = -0.475520468612501;
  dDeltaSpotPut[2] = -0.45664833344749;
  dDeltaSpotPut[3] = -0.438008170894333;
  dDeltaSpotPut[4] = -0.419639303102858;
  dDeltaSpotPut[5] = -0.401578020885981;
  dDeltaSpotPut[6] = -0.383857491794749;
  dDeltaSpotPut[7] = -0.366507705628978;
  dDeltaSpotPut[8] = -0.349555454805862;
  dDeltaSpotPut[9] = -0.333024346752687;
  dDeltaSpotPut[10] = -0.316934845311712;
  dDeltaSpotPut[11] = -0.301304338039002;

     
    // Different values of delta for variable Spot with a call
  for (int i = 0; i < 12; i++)
  CPPUNIT_ASSERT_DOUBLES_EQUAL( CalculateDelta(1, dSpot-1+i, dStrike, 
   dRate, dFRate, dVol, dTminust) , dDeltaSpot[i], 1e-4 );

      // Different values of delta for variable Spot with a put
  for (int j = 0; j < 12; j++)
  CPPUNIT_ASSERT_DOUBLES_EQUAL( CalculateDelta(0, dSpot-1+j, dStrike, 
   dRate, dFRate, dVol, dTminust) , dDeltaSpotPut[j], 1e-4 );


}//CalculateDelta

void BSutileTest::DeltaTime()
{
  double dSpot = 100.;
  double dStrike = 101.;
  double dRate = 0.04;
  double dFRate = 0.05;
  double dVol = 0.2;
  double dTminust = 7.;
  double dDeltaSpot = 0.384179218863;
 
  CPPUNIT_ASSERT_DOUBLES_EQUAL( CalculateDelta(1, 
  dSpot, dStrike, dRate, dFRate, dVol, dTminust) , dDeltaSpot, 1e-4 );

}



//--------------------------------------------------------------------------
//         tests for function GetStrikeFromDelta
//--------------------------------------------------------------------------            

void BSutileTest::GetStrike()
{
  double dSpot = 100.;
  double dStrike = 101.;
  double dRate = 0.04;
  double dFRate = 0.05;
  double dVol = 0.2;
  double dTminust = 1.;
  double dDeltaSpot[12];
 
  dDeltaSpot[0] = 0.456647069999871;  
  dDeltaSpot[1] = 0.475708955888213;  
  dDeltaSpot[2] = 0.494581091053224; 
  dDeltaSpot[3] = 0.513221253606381;  
  dDeltaSpot[4] = 0.531590121397856; 
  dDeltaSpot[5] = 0.549651403614733;
  dDeltaSpot[6] = 0.567371932705965;  
  dDeltaSpot[7] = 0.584721718871737;  
  dDeltaSpot[8] = 0.601673969694852;
  dDeltaSpot[9] = 0.618205077748027;
  dDeltaSpot[10] = 0.634294579189002;
  dDeltaSpot[11] = 0.649925086461712;
     
    // Different values of delta for variable Spot with a call
  for (int i = 0; i < 12; i++)
  CPPUNIT_ASSERT_DOUBLES_EQUAL( GetStrikeFromDelta(dTminust, 
  dSpot-1+i, dDeltaSpot[i], dVol, dRate, dFRate) , dStrike, 1e-4 );

}

/***********************************************************************/

void BSutileTest::GetStrikeTime()
{
  double dSpot = 100.;
  double dStrike = 101.;
  double dRate = 0.04;
  double dFRate = 0.05;
  double dVol = 0.2;
  double dTminust = 7.;
  double dDeltaSpot = 0.384179218863;
 
  CPPUNIT_ASSERT_DOUBLES_EQUAL( GetStrikeFromDelta(dTminust, 
  dSpot, dDeltaSpot, dVol, dRate, dFRate) , dStrike, 1e-4 );

}
