
/////////////////////////////////////////////////////////////////////////////
// Name:       ihg/tests/testsuite/option/src/cppunit_forwardoption_test.cpp
// Purpose:     Base class for testing IHG projects
// Author:      Yann d'Halluin
// Created:     2005/06/13
// RCS-ID:      $Id: cppunit_forwardoption_test.cpp,v 1.3 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <string>
#include <iostream>
#include "ito33/afterstd.h"

#include <math.h>

#include "ito33/exception.h"
#include "forwardoptiontester.h"

#include "ito33/ihg/hazardratecallback.h"
#include "ito33/ihg/volatilitycallback.h"

#include "cppunit_forwardoption_test.h"

using namespace std; 


std::ostream* outputStream;
ito33::ihg::ForwardOptionTester forwardOptionTester;

void __stdcall parametricHR(double dTime, const double *pdS, double *pdValues, 
                            size_t nNbS, int i);

void __stdcall parametricVol(double dTime, const double *pdS, double *pdValues, 
                             size_t nNbS, int i);

namespace ito33
{

namespace ihg
{

namespace test
{

void CppUnitForwardOptionTest::Test1()
{
  string sYieldCurveFile    = "..\\txtfiles\\constyieldmedium.txt";
  string sForeignCurveFile  = "..\\txtfiles\\constyieldzero.txt";
  string sDividendFile      = "..\\txtfiles\\zerodividends.txt";
  string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  string sHazardRateFile    = "..\\txtfiles\\consthrzero.txt";
  string sMeshParamFile     = "..\\txtfiles\\nonuniformmeshes.txt";
  string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  string sOptionFile        = "txtfiles\\forwardoption1.txt";

  forwardOptionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  forwardOptionTester.SetSpotSharePrice(100.0);
  
  double dPrice = 13.7163503;

  CPPUNIT_ASSERT( forwardOptionTester.RunPriceConvergenceTests(dPrice) );
}


void CppUnitForwardOptionTest::Test2()
{
  string sYieldCurveFile    = "..\\txtfiles\\constyieldmedium.txt";
  string sForeignCurveFile  = "..\\txtfiles\\constyieldsmall.txt";
  string sDividendFile      = "..\\txtfiles\\zerodividends.txt";
  string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  string sHazardRateFile    = "..\\txtfiles\\consthrzero.txt";
  string sMeshParamFile     = "..\\txtfiles\\nonuniformmeshes.txt";
  string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  string sOptionFile        = "txtfiles\\forwardoption1.txt";

  forwardOptionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  forwardOptionTester.SetSpotSharePrice(100.0);
  
  double dPrice = 12.544088388;
  
  CPPUNIT_ASSERT( forwardOptionTester.RunPriceConvergenceTests(dPrice) );

}


void CppUnitForwardOptionTest::Test3()
{
  string sYieldCurveFile    = "..\\txtfiles\\constyieldmedium.txt";
  string sForeignCurveFile  = "..\\txtfiles\\constyieldsmall.txt";
  string sDividendFile      = "..\\txtfiles\\mixeddividends.txt";
  string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  string sHazardRateFile    = "..\\txtfiles\\consthrzero.txt";
  string sMeshParamFile     = "..\\txtfiles\\nonuniformmeshes.txt";
  string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  string sOptionFile        = "txtfiles\\forwardoption1.txt";

  forwardOptionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  forwardOptionTester.SetSpotSharePrice(100.0);
  
  double dPrice = 9.4709471;
  //double dDelta = 0.47360002;
  //double dGamma = 0.012719321;
  //double dVega = 37.123866;
  CPPUNIT_ASSERT( forwardOptionTester.RunPriceConvergenceTests(dPrice) );

}


void CppUnitForwardOptionTest::Test4()
{
  string sYieldCurveFile    = "..\\txtfiles\\constyieldmedium.txt";
  string sForeignCurveFile  = "..\\txtfiles\\constyieldzero.txt";
  string sDividendFile      = "..\\txtfiles\\zerodividends.txt";
  string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  string sHazardRateFile    = "..\\txtfiles\\consthrmedium.txt";
  string sMeshParamFile     = "..\\txtfiles\\nonuniformmeshes.txt";
  string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  string sOptionFile        = "txtfiles\\forwardoption1.txt";

  forwardOptionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  forwardOptionTester.SetSpotSharePrice(100.0);
  
  double dPrice = 18.815157;
  //double dDelta = 0.73041516;
  //double dGamma = 0.011013006;
  //double dVega = 33.039085;
  CPPUNIT_ASSERT( forwardOptionTester.RunPriceConvergenceTests(dPrice) );
}

void CppUnitForwardOptionTest::Test5()
{
  string sYieldCurveFile    = "..\\txtfiles\\constyieldmedium.txt";
  string sForeignCurveFile  = "..\\txtfiles\\constyieldzero.txt";
  string sDividendFile      = "..\\txtfiles\\zerodividends.txt";
  string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  string sHazardRateFile    = "..\\txtfiles\\consthrzero.txt";
  string sMeshParamFile     = "..\\txtfiles\\uniformmeshes.txt";
  string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  string sOptionFile        = "txtfiles\\forwardoption1.txt";

  forwardOptionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  forwardOptionTester.SetSpotSharePrice(100.0);

  forwardOptionTester.SetHazardRate(ito33::shared_ptr<ito33::ihg::HazardRate>
    (new ito33::ihg::HazardRateCallBack(parametricHR,0)));

  double dPrice = 23.625661;
  forwardOptionTester.RunPriceConvergenceTests(dPrice);

  sMeshParamFile = "..\\txtfiles\\nonuniformmeshes.txt";
  forwardOptionTester.ReadMeshParams(sMeshParamFile);
  
  CPPUNIT_ASSERT( forwardOptionTester.RunPriceTests(dPrice, 1.e-3) );

}

void CppUnitForwardOptionTest::Test6()
{
  string sYieldCurveFile    = "..\\txtfiles\\constyieldmedium.txt";
  string sForeignCurveFile  = "..\\txtfiles\\constyieldzero.txt";
  string sDividendFile      = "..\\txtfiles\\zerodividends.txt";
  string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  string sHazardRateFile    = "..\\txtfiles\\consthrzero.txt";
  string sMeshParamFile     = "..\\txtfiles\\uniformmeshes.txt";
  string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  string sOptionFile        = "txtfiles\\forwardoption1.txt";

  forwardOptionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  forwardOptionTester.SetSpotSharePrice(100.0);
  
  forwardOptionTester.SetVolatility(ito33::shared_ptr<ito33::ihg::Volatility>
    (new ito33::ihg::VolatilityCallBack(parametricVol,0)));

  // need 6 tests to get the implicit convergence down to 2
  double dPrice = 8.117522;
  forwardOptionTester.RunPriceConvergenceTests(dPrice, 6);

  sMeshParamFile = "..\\txtfiles\\nonuniformmeshes.txt";
  forwardOptionTester.ReadMeshParams(sMeshParamFile);
  
  CPPUNIT_ASSERT( forwardOptionTester.RunPriceTests(dPrice, 1.3e-3) );

}

void CppUnitForwardOptionTest::Test7()
{
  string sYieldCurveFile    = "..\\txtfiles\\constyieldmedium.txt";
  string sForeignCurveFile  = "..\\txtfiles\\constyieldsmall.txt";
  string sDividendFile      = "..\\txtfiles\\mixeddividends.txt";
  string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  string sHazardRateFile    = "..\\txtfiles\\consthrzero.txt";
  string sMeshParamFile     = "..\\txtfiles\\uniformmeshes.txt";
  string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  string sOptionFile        = "txtfiles\\forwardoption1.txt";

  forwardOptionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  forwardOptionTester.SetSpotSharePrice(100.0);
  
  forwardOptionTester.SetVolatility(ito33::shared_ptr<ito33::ihg::Volatility>
    (new ito33::ihg::VolatilityCallBack(parametricVol,0)));

  forwardOptionTester.SetHazardRate(ito33::shared_ptr<ito33::ihg::HazardRate>
    (new ito33::ihg::HazardRateCallBack(parametricHR,0)));

  double dPrice = 14.075379;
  forwardOptionTester.RunPriceConvergenceTests(dPrice, 6);

  sMeshParamFile = "..\\txtfiles\\nonuniformmeshes.txt";
  forwardOptionTester.ReadMeshParams(sMeshParamFile);
  
  CPPUNIT_ASSERT( forwardOptionTester.RunPriceTests(dPrice, 1.1e-3) );

}


void CppUnitForwardOptionTest::Test8()
{
  string sYieldCurveFile    = "..\\txtfiles\\constyieldmedium.txt";
  string sForeignCurveFile  = "..\\txtfiles\\constyieldzero.txt";
  string sDividendFile      = "..\\txtfiles\\zerodividends.txt";
  string sVolatilityFile    = "..\\txtfiles\\constvolmedium.txt";
  string sHazardRateFile    = "..\\txtfiles\\polynomhrlarge.txt";
  string sMeshParamFile     = "..\\txtfiles\\uniformmeshes.txt";
  string sNumParamFile      = "..\\txtfiles\\implicitnumparams.txt";
  string sValuationDateFile = "..\\txtfiles\\valuationdate.txt";
  string sOptionFile        = "txtfiles\\forwardoption1.txt";

  forwardOptionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
    sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
    sNumParamFile, sValuationDateFile, sOptionFile);

  forwardOptionTester.SetSpotSharePrice(100.0);

  double dPrice = 22.835853;
  forwardOptionTester.RunPriceConvergenceTests(dPrice);

  sMeshParamFile = "..\\txtfiles\\nonuniformmeshes.txt";
  forwardOptionTester.ReadMeshParams(sMeshParamFile);
  
  CPPUNIT_ASSERT( forwardOptionTester.RunPriceTests(dPrice, 1.e-3) );

}


} //end test
} // end ihg
}//end ito33