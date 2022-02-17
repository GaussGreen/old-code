/////////////////////////////////////////////////////////////////////////////
// Name:        tests/ihg/testvolatility.cpp
// Purpose:     testing volatility
// Author:      Yann d'Halluin
// Created:     11/07/2004
// RCS-ID:      $Id: testvolatility.cpp,v 1.10 2006/06/15 09:27:37 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/beforestd.h"
#include <iostream>
#include <vector>
#include <stdlib.h>
#include "ito33/afterstd.h"

#include "math.h"

#include "ito33/cppunit.h"
#include "ito33/exception.h"
#include "ito33/array.h"
#include "ito33/date.h"
#include "ito33/dateutils.h"

#include "ito33/xml/write.h"

#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/volatilitypower.h"
#include "ito33/ihg/volatilitytanh.h"
#include "ito33/ihg/volatilitytimeonly.h"
#include "ito33/ihg/parametrization.h"

#include "ihg/tests/testvolatility.h"
#include "ihg/tests/testparametrization.h"

#include "ito33/tests/utilexml.h"

using namespace ito33;
using namespace ito33::ihg;

// ----------------------------------------------------------------------------
// VolatilityFlat tests
// ----------------------------------------------------------------------------

void VolatilityFlatTest::Value()
{
  double dVol = 0.4;
  VolatilityFlat vol(dVol);

 double dValue = vol.GetValue();

 ITO33_ASSERT_DOUBLES_EQUAL( dValue, dVol);

}

void VolatilityFlatTest::GetVolsSquared()
{

 double dVol = 0.3;
 double dVolSquared = dVol *dVol;
 VolatilityFlat vol(dVol);

 //Check the GetVolsSquared function for
 //different spot values
 double dValue;
 double dSpot = 100.0;
 vol.GetVolsSquared(0.0,&dSpot,&dValue,1);

 ITO33_ASSERT_DOUBLES_EQUAL(dValue, dVolSquared);

 double pdSpots[2] = {300.0, 500.0};
 double pdValues[2];

 vol.GetVolsSquared(3.0,pdSpots,pdValues,2);
 ITO33_ASSERT_DOUBLES_EQUAL(pdValues[0], dVolSquared);
 ITO33_ASSERT_DOUBLES_EQUAL(pdValues[1], dVolSquared);

}

void VolatilityFlatTest::GetVols()
{
  double dVol = 0.3;
  VolatilityFlat vol(dVol);
  
  // Check the GetVols function for various times
  // and array sizes
  double dValue;
  double dSpot = 100.0;
  vol.GetVols(0.0, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(dValue, dVol);

  vol.GetVols(1.0, &dSpot, &dValue, 1);
  ITO33_ASSERT_DOUBLES_EQUAL(dValue, dVol);

  double pdSpots[2] = {100.0, 200.0};
  double pdValues[2];
  
  vol.GetVols(0.0, pdSpots, pdValues, 2);
  ITO33_ASSERT_DOUBLES_EQUAL(pdValues[0], dVol);
  ITO33_ASSERT_DOUBLES_EQUAL(pdValues[1], dVol);

  vol.GetVols(100.1, pdSpots, pdValues, 2);
  ITO33_ASSERT_DOUBLES_EQUAL(pdValues[0], dVol);
  ITO33_ASSERT_DOUBLES_EQUAL(pdValues[1], dVol);

}


void VolatilityFlatTest::VolatilityNegative()
{
  // Cannot have negative volatility
 VolatilityFlat vol(-0.1);
}

void VolatilityFlatTest::VolatilityTooLarge()
{
  // Cannot be too large
  VolatilityFlat vol(10.0);
}

void VolatilityFlatTest::SerializationOutput()
{
  
  // Cannot be too large
  VolatilityFlat vol(.2);

  finance::ModelParametersConsumerTest visitor;

  vol.GetModelParameters(visitor);

  CPPUNIT_ASSERT( visitor.m_pParameters.size() == 1 );

  ITO33_ASSERT_DOUBLES_EQUAL( visitor.m_pParameters[0].value , .2 );

  CPPUNIT_ASSERT( visitor.m_sCategoryName == MODEL_PARAM_NAME_VOL_FLAT);

  CPPUNIT_ASSERT( visitor.m_pParameters[0].name == SCALAR_MODEL_PARAM_NAME_VALUE);

}


void VolatilityFlatTest::Dump()
{
  ExpectedXML expected(m_oss,
                "<?xml version=\"1.0\"?>"
                "<root>\n"
                "<volatility_flat>\n"
                "<flat>0.3</flat>\n"
                "</volatility_flat>\n"
                "</root>"
              );

  ito33::XML::RootTag root("root",m_oss);

  VolatilityFlat vol(.3);
  vol.Dump(root);
}


void VolatilityPowerTest::AlphaNegative()
{
  VolatilityPower vol(-1., -1., 100.);
}

void VolatilityPowerTest::SpotNegative()
{
  VolatilityPower vol(1.0, -1., -1.);
}
  
void VolatilityPowerTest::GetVolsSquared()
{
  double dAlpha        = .2;
  double dBeta         = 1.5;
  double dS0           = 52.3;
  double dAlphaSquared = dAlpha * dAlpha;

  VolatilityPower vol(dAlpha, dBeta, dS0);

  size_t nSize = 5;
  std::vector<double> pdS(nSize);

  pdS[0] = 10.0;
  pdS[1] = 20.0;
  pdS[2] = 40.0;
  pdS[3] = 80.0;
  pdS[4] = 160.0;

  std::vector<double> pdVolsSquared(nSize);
  vol.GetVolsSquared(0.0,&pdS[0],&pdVolsSquared[0],nSize);

  for ( size_t nIdx = 0 ; nIdx < nSize ; nIdx++)
  {
    double dTmp = 1./pdS[nIdx];
    double dVolSquared = dAlphaSquared * pow(dS0 * dTmp, 2.*dBeta);

    ITO33_ASSERT_DOUBLES_EQUAL(dVolSquared, pdVolsSquared[nIdx]);
  }
}


void VolatilityPowerTest::GetVols()
{
   double dAlpha        = .2;
  double dBeta         = 1.5;
  double dS0           = 52.3;

  VolatilityPower vol(dAlpha, dBeta, dS0);

  size_t nSize = 5;
  std::vector<double> pdS(nSize);

  pdS[0] = 10.0;
  pdS[1] = 20.0;
  pdS[2] = 40.0;
  pdS[3] = 80.0;
  pdS[4] = 160.0;

  std::vector<double> pdVols(nSize);
  vol.GetVols(0.0,&pdS[0],&pdVols[0],nSize);

  for ( size_t nIdx = 0 ; nIdx < nSize ; nIdx++)
  {
    double dTmp = 1./pdS[nIdx];
    double dVol = dAlpha * pow(dS0 * dTmp, dBeta);

    ITO33_ASSERT_DOUBLES_EQUAL(dVol, pdVols[nIdx]);
  }

}

void VolatilityPowerTest::SerializationOutput()
{
  double dAlpha        = .2;
  double dBeta         = 1.5;
  double dS0           = 52.3;

  VolatilityPower vol(dAlpha, dBeta, dS0);

  finance::ModelParametersConsumerTest visitor;

  vol.GetModelParameters(visitor);

  CPPUNIT_ASSERT( visitor.m_pParameters.size() == 3);

  CPPUNIT_ASSERT( visitor.m_sCategoryName == MODEL_PARAM_NAME_VOL_POWER);

  CPPUNIT_ASSERT( visitor.m_pParameters[0].name == SCALAR_MODEL_PARAM_NAME_ALPHA);
  CPPUNIT_ASSERT( visitor.m_pParameters[1].name == SCALAR_MODEL_PARAM_NAME_BETA); 
  CPPUNIT_ASSERT( visitor.m_pParameters[2].name == SCALAR_MODEL_PARAM_NAME_S0);

  CPPUNIT_ASSERT( visitor.m_pParameters[0].value == .2 );
  CPPUNIT_ASSERT( visitor.m_pParameters[1].value == 1.5 );
  CPPUNIT_ASSERT( visitor.m_pParameters[2].value == 52.3 );


}

void VolatilityPowerTest::Dump()
{ 
   VolatilityPower vol(1., 1., 1.);

   ExpectedXML expected(m_oss,
                "<?xml version=\"1.0\"?>"
                "<root>\n"
                "<volatility_power>\n"
                "<alpha>1</alpha>\n"
                "<beta>1</beta>\n"
                "<S0>1</S0>\n"
                "</volatility_power>\n"
                "</root>"
              );

  ito33::XML::RootTag root("root",m_oss);

  vol.Dump(root);
}

void VolatilityTanhTest::LeftNegative()
{
  double dLeft  = 1.0;
  double dRight = 2.0;
  double dScale = .2;
  double dS0    = 53.2;

  VolatilityTanh vol(dLeft, -dRight, dScale, dS0);
}

void VolatilityTanhTest::RightNegative()
{
  double dLeft  = 1.0;
  double dRight = 2.0;
  double dScale = .2;
  double dS0    = 53.2;

  VolatilityTanh vol(dLeft, -dRight, dScale, dS0);
}

void VolatilityTanhTest::ScaleNegative()
{
  double dLeft  = 1.0;
  double dRight = 2.0;
  double dScale = .2;
  double dS0    = 53.2;

  VolatilityTanh vol(dLeft, dRight, -dScale, dS0);
}

void VolatilityTanhTest::SpotNegative()
{
  double dLeft  = 1.0;
  double dRight = 2.0;
  double dScale = .2;
  double dS0    = 53.2;

  VolatilityTanh vol(dLeft, dRight, dScale, -dS0);
}

void VolatilityTanhTest::Dump()
{
  double dLeft  = 1.0;
  double dRight = 2.0;
  double dScale = .2;
  double dS0    = 53.2;

  VolatilityTanh vol(dLeft, dRight, dScale, dS0);
 
  ExpectedXML expected(m_oss,
                "<?xml version=\"1.0\"?>"
                "<root>\n"
                "<volatility_tanh>\n"
                "<left>1</left>\n"
                "<right>2</right>\n"
                "<scale>0.2</scale>\n"
                "<S0>53.2</S0>\n"
                "</volatility_tanh>\n"
                "</root>"
              );

  ito33::XML::RootTag root("root",m_oss);

  vol.Dump(root);

}

void VolatilityTanhTest::GetVolsSquared()
{
  double dLeft  = 1.0;
  double dRight = 2.0;
  double dScale = .2;
  double dS0    = 53.2;

  VolatilityTanh vol(dLeft, dRight, dScale, dS0);

  size_t nSize = 5;
  std::vector<double> pdS(nSize);

  pdS[0] = 10.0;
  pdS[1] = 20.0;
  pdS[2] = 40.0;
  pdS[3] = 80.0;
  pdS[4] = 160.0;

  std::vector<double> pdVolsSquared(nSize);
  vol.GetVolsSquared(0.0, &pdS[0], &pdVolsSquared[0], nSize);
  
  for (size_t nIdx = 0; nIdx < nSize; nIdx++)
  {
    double dVolSquared = 0.5 * (dRight - dLeft) 
                              * ( 1. + tanh( dScale * (pdS[nIdx] - dS0) ) ) 
                          + dLeft; 

    dVolSquared *= dVolSquared;

    ITO33_ASSERT_DOUBLES_EQUAL(dVolSquared, pdVolsSquared[nIdx]);
  }

}


void VolatilityTanhTest::GetVols()
{
  double dLeft  = 1.0;
  double dRight = 2.0;
  double dScale = .2;
  double dS0    = 53.2;

  VolatilityTanh vol(dLeft, dRight, dScale, dS0);

  size_t nSize = 5;
  std::vector<double> pdS(nSize);

  pdS[0] = 10.0;
  pdS[1] = 20.0;
  pdS[2] = 40.0;
  pdS[3] = 80.0;
  pdS[4] = 160.0;

  std::vector<double> pdVols(nSize);
  vol.GetVols(0.0, &pdS[0], &pdVols[0], nSize);

    for (size_t nIdx = 0; nIdx < nSize; nIdx++)
  {
    double dVol  = 0.5 * (dRight - dLeft) 
                       * ( 1. + tanh( dScale * (pdS[nIdx] - dS0) ) )
                 + dLeft; 

    ITO33_ASSERT_DOUBLES_EQUAL(dVol, pdVols[nIdx]);
  }
}


void VolatilityTanhTest::SerializationOutput()
{
  double dLeft  = 1.0;
  double dRight = 2.0;
  double dScale = .2;
  double dS0    = 53.2;

  VolatilityTanh vol(dLeft, dRight, dScale, dS0);

  finance::ModelParametersConsumerTest visitor;

  vol.GetModelParameters( visitor );

  CPPUNIT_ASSERT( visitor.m_pParameters.size() == 4 );

  CPPUNIT_ASSERT( visitor.m_sCategoryName == MODEL_PARAM_NAME_VOL_TANH);

  CPPUNIT_ASSERT( visitor.m_pParameters[0].name == SCALAR_MODEL_PARAM_NAME_LEFT_LIMIT);
  CPPUNIT_ASSERT( visitor.m_pParameters[1].name == SCALAR_MODEL_PARAM_NAME_RIGHT_LIMIT);
  CPPUNIT_ASSERT( visitor.m_pParameters[2].name == SCALAR_MODEL_PARAM_NAME_SCALE );
  CPPUNIT_ASSERT( visitor.m_pParameters[3].name == SCALAR_MODEL_PARAM_NAME_S0 );

  CPPUNIT_ASSERT( visitor.m_pParameters[0].value == 1.0);
  CPPUNIT_ASSERT( visitor.m_pParameters[1].value == 2.0);
  CPPUNIT_ASSERT( visitor.m_pParameters[2].value == .2);
  CPPUNIT_ASSERT( visitor.m_pParameters[3].value == 53.2);


}

void VolatilityTimeOnlyTest::VolatilityTooLarge()
{
  std::vector<Date> pDates;
  std::vector<double> pdValues;

  pDates.push_back( Date(2005,Date::Jan,1) );
  pdValues.push_back( 6. );

  VolatilityTimeOnly vol(pDates, pdValues);

}

void VolatilityTimeOnlyTest::SerializationOutput()
{  
  std::vector<Date> pDates;
  std::vector<double> pdValues;

  pDates.push_back( Date(2005, Date::Jan, 1) );
  pDates.push_back( Date(2005, Date::Jan, 2) );
  pDates.push_back( Date(2005, Date::Jan, 3) );
  pDates.push_back( Date(2005, Date::Jan, 4) );
  pDates.push_back( Date(2005, Date::Jan, 5) );
  pDates.push_back( Date(2005, Date::Jan, 6) );
  pDates.push_back( Date(2005, Date::Jan, 7) );
  pDates.push_back( Date(2005, Date::Jan, 8) );
  pDates.push_back( Date(2005, Date::Jan, 9) );
  pDates.push_back( Date(2005, Date::Jan, 10) );
  pdValues.push_back( .1 );
  pdValues.push_back( .2 );
  pdValues.push_back( .3 );
  pdValues.push_back( .4 );
  pdValues.push_back( .5 );
  pdValues.push_back( .6 );
  pdValues.push_back( .7 );
  pdValues.push_back( .8 );
  pdValues.push_back( .9 );
  pdValues.push_back( 1. );


  ihg::VolatilityTimeOnly vol(pDates, pdValues);

  finance::ModelParametersConsumerTest visitor;

  vol.GetModelParameters(visitor);

  CPPUNIT_ASSERT( visitor.m_pdVal.size() == 10 );

  CPPUNIT_ASSERT( visitor.m_pDates.size() == 10 );

  CPPUNIT_ASSERT( visitor.m_sCategoryName == MODEL_PARAM_NAME_VOL_TIMEONLY );

  size_t nIdx;

  for ( nIdx = 0 ; nIdx < visitor.m_pDates.size(); nIdx++ )
  {
    CPPUNIT_ASSERT( visitor.m_pDates[nIdx] == Date(2005, Date::Jan, nIdx+1 ) );
    ITO33_ASSERT_DOUBLES_EQUAL( visitor.m_pdVal[nIdx], (nIdx+1)/10. );
  }
}

void VolatilityTimeOnlyTest::Dump()
{
  std::vector<Date> pDates;
  std::vector<double> pdValues;

  pDates.push_back( Date(2005, Date::Jan, 1) );
  pDates.push_back( Date(2005, Date::Jan, 2) );
  pDates.push_back( Date(2005, Date::Jan, 3) );
  pDates.push_back( Date(2005, Date::Jan, 4) );
  pDates.push_back( Date(2005, Date::Jan, 5) );
  pdValues.push_back( .1 );
  pdValues.push_back( .2 );
  pdValues.push_back( .3 );
  pdValues.push_back( .4 );
  pdValues.push_back( .5 );



  ihg::VolatilityTimeOnly vol(pDates, pdValues);


  ExpectedXML expected(m_oss,
                "<?xml version=\"1.0\"?>"
                "<root>\n"
                "<volatility_time_only>\n"
                "<dates>\n"
                "<date>2005-01-01</date>\n"
                "<date>2005-01-02</date>\n"
                "<date>2005-01-03</date>\n"
                "<date>2005-01-04</date>\n"
                "<date>2005-01-05</date>\n"
                "</dates>\n"
                "<values>\n"
                "<value>0.1</value>\n"
                "<value>0.2</value>\n"
                "<value>0.3</value>\n"
                "<value>0.4</value>\n"
                "<value>0.5</value>\n"
                "</values>\n"
                "</volatility_time_only>\n"
                "</root>"
              );

  ito33::XML::RootTag root("root",m_oss);

  vol.Dump(root);  

}