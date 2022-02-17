/////////////////////////////////////////////////////////////////////////////
// Name:        tests/simulatedannealing/testsa.cpp
// Purpose:     unit tests for simulated annealing
// Author:      ITO 33
// Created:     April 19, 2005
// RCS-ID:      $Id: testsa.cpp,v 1.7 2006/08/19 23:22:41 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "ito33/common.h"
#include "ito33/exception.h"
#include "ito33/cppunit.h"
#include "ito33/useexception.h"
#include "ito33/date.h"
#include "ito33/sharedptr.h"
#include "ito33/link.h"

#include "ito33/numeric/sa.h"
#include "ito33/numeric/asa.h"
#include "ito33/numeric/numericerror.h"

#include "ito33/tests/testsa.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/option.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/computationalflags.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardrateflat.h"

ITO33_FORCE_LINK_MODULE(IHGPriceOption);

using namespace ito33;
using namespace ito33::numeric;
using namespace ito33::finance;
using namespace ito33::ihg;

extern const ito33::Error ITO33_BAD_PARAM;



//global test functions
void linear(const std::vector<double> &pdParam, double &dF)
{
  double dFunc1 = pdParam[0] - pdParam[1];
  double dFunc2 = 2.0*pdParam[0] + 3.0*pdParam[1] - 2.0;

  dF = 0.5 * (dFunc1 * dFunc1 + dFunc2 * dFunc2);
}

void rosenbrock(const std::vector<double> &pdParam, double &dF)
{
  size_t nIdx;
  size_t nSize = pdParam.size();
  
  dF = 0.0;
  //Rosenbrock test function

  for ( nIdx = 1 ; nIdx < nSize; nIdx++)
  {
    dF += 100.*(pdParam[nIdx]-pdParam[nIdx-1])*(pdParam[nIdx]-pdParam[nIdx-1])
      + (1-pdParam[nIdx-1])*(1-pdParam[nIdx-1]);
  }
 
}

void sphere(const std::vector<double> &pdParam, double &dF)
{
  dF = 0.0;

  size_t nIdx;
  size_t nSize = pdParam.size();

  for (nIdx = 0 ; nIdx < nSize; nIdx++)
    dF = dF + pdParam[nIdx]*pdParam[nIdx];

}

void plateau(const std::vector<double> &pdParam, double &dF)
{
  size_t nIdx;
  size_t nSize = pdParam.size();
  dF = 0.;

  for (nIdx = 0 ; nIdx < nSize; nIdx++)
    dF = dF + floor( pdParam[nIdx] );

  dF = 6.* nSize + dF;
}

void calibration(const std::vector<double> &pdParam, double &dF)
{
 
  Date valuationDate(2003, Date::Feb, 1);

  shared_ptr<Numeraire> pCurrency( new Numeraire("EUR") );

  shared_ptr<Equity> pEquity(new Equity(50., pCurrency));
  shared_ptr<YieldCurve> pyf(new YieldCurveFlat(0.01));
    
  pEquity->SetBorrowCurve(pyf);
  
  shared_ptr<YieldCurve> pyc( new YieldCurveFlat(0.01) );
  shared_ptr<RateData> pRateData(new RateData() );
  pRateData->SetYieldCurve(pCurrency, pyc) ;
     
  shared_ptr<SessionData> 
    pSessionData(new SessionData(pRateData, pEquity, valuationDate));

  shared_ptr<Option> opt1(new Option(45, 
                                    Date(2005, Date::Feb, 1),
                                    Option_Call,
                                    ExerciseType_European
                                    )
                          );

  shared_ptr<Option> opt2(new Option(55, 
                                    Date(2006, Date::Feb, 1),
                                    Option_Call,
                                    ExerciseType_European
                                    )
                         );

   
  opt1->SetSessionData(pSessionData);
  opt2->SetSessionData(pSessionData);

  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);
    
  pModel->SetVolatility( shared_ptr<ihg::Volatility>(new VolatilityFlat(0.23456)) );
  pModel->SetHazardRate( shared_ptr<ihg::HazardRate>(new HazardRateFlat(0.123456)) );

  shared_ptr<ModelOutput> output1 = pModel->Compute(*opt1);
  shared_ptr<ModelOutput> output2 = pModel->Compute(*opt2);

  double dPrice1 = output1->GetPrice();
  double dPrice2 = output2->GetPrice();


  //Price using guesses from sa

  pModel->SetVolatility( shared_ptr<ihg::Volatility>(new VolatilityFlat( pdParam[0] )) );
  pModel->SetHazardRate( shared_ptr<ihg::HazardRate>(new HazardRateFlat( pdParam[1] )) );

  shared_ptr<ModelOutput> outputParam1 = pModel->Compute(*opt1);
  shared_ptr<ModelOutput> outputParam2 = pModel->Compute(*opt2);

  double dParamPrice1 = outputParam1->GetPrice();
  double dParamPrice2 = outputParam2->GetPrice();

  dF = (dPrice1 - dParamPrice1)/dParamPrice1*(dPrice1 - dParamPrice1)/dParamPrice1
       + (dPrice2 - dParamPrice2)/dParamPrice2*(dPrice2 - dParamPrice2)/dParamPrice2;

}

void SATest::operator () (const std::vector<double> &pdParam, double &dF)
{
 
  switch(m_test)  
  {
  case LINEAR:
    linear(pdParam,dF);
    break;
  case ROSENBROCK:
    rosenbrock(pdParam,dF);
    break;
  case SPHERE:
    sphere(pdParam,dF);
    break;
  case PLATEAU:
    plateau(pdParam,dF);
    break;
  case CALIBRATION:
    calibration(pdParam,dF);
    break;

  default:
    CPPUNIT_ASSERT(false);
    
  }

}


void SATest::TestLinear()
{
  m_test = LINEAR;

  size_t nSize = 2;
  size_t nIdx;

  std::vector<double> pdParamLowerBound(nSize);
  std::vector<double> pdParamUpperBound(nSize);
  std::vector<double> pdParam(nSize);

  for ( nIdx = 0 ; nIdx < nSize ; nIdx++ )
  {
    pdParamLowerBound[nIdx] = -10.;
    pdParamUpperBound[nIdx] = 10.;
  }

  // linear problem only needs 1 iteration
  numeric::SA solver(pdParamLowerBound, pdParamUpperBound);

  pdParam[0] = 3.;
  pdParam[1] = -2.1;

  numeric::NumericError err = solver(*this, pdParam);

  CPPUNIT_ASSERT(err == ITO33_NO_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.4, pdParam[0], 1e-4);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.4, pdParam[1], 1e-4);
  CheckObjectiveFunction(pdParam);

}


void SATest::TestRosenbrock()
{

  m_test = ROSENBROCK;

  size_t nSize = 4;
  size_t nIdx;

  std::vector<double> pdParamLowerBound(nSize);
  std::vector<double> pdParamUpperBound(nSize);
  std::vector<double> pdParam(nSize);

  for ( nIdx = 0 ; nIdx < nSize ; nIdx++ )
  {
    pdParamLowerBound[nIdx] = -10.;
    pdParamUpperBound[nIdx] = 10.;
    pdParam[nIdx] = 5.2;
  }

  numeric::SA solver(pdParamLowerBound, pdParamUpperBound, 1.e-6);

  numeric::NumericError err = solver(*this, pdParam);

  CPPUNIT_ASSERT(err == ITO33_NO_ERROR);

  for ( nIdx = 0 ; nIdx < nSize ; nIdx++ )
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, pdParam[nIdx], 1e-3);

  CheckObjectiveFunction(pdParam);
 
}

void SATest::TestSphere()
{
  m_test = SPHERE;

  size_t nSize = 2;
  size_t nIdx;

  std::vector<double> pdParamLowerBound(nSize);
  std::vector<double> pdParamUpperBound(nSize);
  std::vector<double> pdParam(nSize);

  for ( nIdx = 0 ; nIdx < nSize ; nIdx++ )
  {
    pdParamLowerBound[nIdx] = -10.;
    pdParamUpperBound[nIdx] = 10.;
    pdParam[nIdx] = -2.0;
  }

  numeric::SA solver(pdParamLowerBound, pdParamUpperBound);

  numeric::NumericError err = solver(*this, pdParam);

  CPPUNIT_ASSERT(err == ITO33_NO_ERROR);

  for ( nIdx = 0 ; nIdx < nSize ; nIdx++ )
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, pdParam[nIdx], 1e-4);

  CheckObjectiveFunction(pdParam);

}

void SATest::TestPlateau()
{
  m_test = PLATEAU;

  size_t nSize = 4;
  size_t nIdx;

  std::vector<double> pdParamLowerBound(nSize);
  std::vector<double> pdParamUpperBound(nSize);
  std::vector<double> pdParam(nSize);

  for ( nIdx = 0 ; nIdx < nSize ; nIdx++ )
  {
    pdParamLowerBound[nIdx] = -5.12;
    pdParamUpperBound[nIdx] = 5.12;
    pdParam[nIdx] = 4.0;
  }

  // linear problem only needs 1 iteration
  numeric::SA solver(pdParamLowerBound, pdParamUpperBound);

  numeric::NumericError err = solver(*this, pdParam);

  CPPUNIT_ASSERT(err == ITO33_NO_ERROR);

  CheckObjectiveFunction(pdParam);

}


void SATest::TestOptionCalibration()
{
 //try to find hazard rate and volatility

  m_test = CALIBRATION;

  size_t nSize = 2;
  size_t nIdx;

  std::vector<double> pdParamLowerBound(nSize);
  std::vector<double> pdParamUpperBound(nSize);
  std::vector<double> pdParam(nSize);

  for ( nIdx = 0 ; nIdx < nSize ; nIdx++ )
  {
    pdParamLowerBound[nIdx] = 0.;
    pdParamUpperBound[nIdx] = 1.;
    pdParam[nIdx] = .5;
  }

  // linear problem only needs 1 iteration
  numeric::SA solver(pdParamLowerBound, pdParamUpperBound, 1.e-6);

 // solver.SetOutputToTrue();
  solver.SetNumberOfStepsBeforeAlteringDirection(10);
  solver.SetNEpsilon(4);

  numeric::NumericError err = solver(*this, pdParam);

  CPPUNIT_ASSERT(err == ITO33_NO_ERROR);

  CheckObjectiveFunction(pdParam);
}

void SATest::CheckObjectiveFunction(const std::vector<double> &pdParam)
{
  double dF = 0.0;
  SATest::operator ()(pdParam,dF);

  CPPUNIT_ASSERT( dF <= 1.e-6 );
}

//___________________________ASA testing___________________________________//

void ASATest::operator () (const std::vector<double> &pdParam, double &dF)
{
 
  switch(m_test)  
  {
  case LINEAR:
    linear(pdParam,dF);
    break;
  case ROSENBROCK:
    rosenbrock(pdParam,dF);
    break;
  case SPHERE:
    sphere(pdParam,dF);
    break;
  case PLATEAU:
    plateau(pdParam,dF);
    break;
  case CALIBRATION:
    calibration(pdParam,dF);
    break;

  default:
    CPPUNIT_ASSERT(false);
    
  }

}





void ASATest::TestRosenbrock()
{

  m_test = ROSENBROCK;

  size_t nSize = 2;
  size_t nIdx;

  std::vector<double> pdParamLowerBound(nSize);
  std::vector<double> pdParamUpperBound(nSize);
  std::vector<double> pdParam(nSize);

  for ( nIdx = 0 ; nIdx < nSize ; nIdx++ )
  {
    pdParamLowerBound[nIdx] = -10.;
    pdParamUpperBound[nIdx] = 10.;
    pdParam[nIdx] = 5.;
  }

  numeric::ASA solver(pdParamLowerBound, pdParamUpperBound, 1.e-7);
  solver.SetLimitAcceptance(50000);
  solver.SetLimitGenerated(50000);
  solver.SetSamplingNumber(5);
  solver.SetCostParameterScaleRatio(.08);
  solver.SetReannealingFrequency(10);



  numeric::NumericError err = solver(*this, pdParam);

  CheckObjectiveFunction(pdParam);

  CPPUNIT_ASSERT( err == ITO33_NO_ERROR );

  for ( nIdx = 0 ; nIdx < nSize ; nIdx++ )
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, pdParam[nIdx], 1e-3);

 
 
}

void ASATest::TestSphere()
{
  m_test = SPHERE;

  size_t nSize = 2;
  size_t nIdx;

  std::vector<double> pdParamLowerBound(nSize);
  std::vector<double> pdParamUpperBound(nSize);
  std::vector<double> pdParam(nSize);

  for ( nIdx = 0 ; nIdx < nSize ; nIdx++ )
  {
    pdParamLowerBound[nIdx] = -5.12;
    pdParamUpperBound[nIdx] = 5.12;
    pdParam[nIdx] = 1.0;
  }

  numeric::ASA solver(pdParamLowerBound, pdParamUpperBound);

  numeric::NumericError err = solver(*this, pdParam);

  CheckObjectiveFunction(pdParam);

  CPPUNIT_ASSERT( err == ITO33_NO_ERROR );

  for ( nIdx = 0 ; nIdx < nSize ; nIdx++ )
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, fabs(pdParam[nIdx]), 1e-3);

  

}

void ASATest::TestPlateau()
{
  m_test = PLATEAU;

  size_t nSize = 4;
  size_t nIdx;

  std::vector<double> pdParamLowerBound(nSize);
  std::vector<double> pdParamUpperBound(nSize);
  std::vector<double> pdParam(nSize);

  for ( nIdx = 0 ; nIdx < nSize ; nIdx++ )
  {
    pdParamLowerBound[nIdx] = -5.12;
    pdParamUpperBound[nIdx] = 5.12;
    pdParam[nIdx] = 4.0;
  }

  // linear problem only needs 1 iteration
  numeric::ASA solver(pdParamLowerBound, pdParamUpperBound);

  numeric::NumericError err = solver(*this, pdParam);

  CheckObjectiveFunction(pdParam);

  CPPUNIT_ASSERT(err == ITO33_NO_ERROR);

  

}


void ASATest::TestOptionCalibration()
{
 //try to find hazard rate and volatility

  m_test = CALIBRATION;

  size_t nSize = 2;
  size_t nIdx;

  std::vector<double> pdParamLowerBound(nSize);
  std::vector<double> pdParamUpperBound(nSize);
  std::vector<double> pdParam(nSize);

  for ( nIdx = 0 ; nIdx < nSize ; nIdx++ )
  {
    pdParamLowerBound[nIdx] = 0.;
    pdParamUpperBound[nIdx] = 1.;
    pdParam[nIdx] = .5;
  }

  // linear problem only needs 1 iteration
  numeric::ASA solver(pdParamLowerBound, pdParamUpperBound);
  solver.SetOutput("asa_out");
  solver.SetLimitAcceptance(50);
  solver.SetLimitGenerated(300);
  //solver.SetSamplingNumber(5);
  //solver.SetCostParameterScaleRatio(.1);
  //solver.SetReannealingFrequency(20);

  numeric::NumericError err = solver(*this, pdParam);

  CheckObjectiveFunction(pdParam);

  CPPUNIT_ASSERT( err == ITO33_NO_ERROR );

}

void ASATest::CheckObjectiveFunction(const std::vector<double> &pdParam)
{
  double dF = 0.0;
  ASATest::operator ()(pdParam,dF);

  //for (size_t nIdx = 0; nIdx < pdParam.size();nIdx++)
  // std::cout << "["<<nIdx<<"] ="<< pdParam[nIdx] << std::endl;

  //std::cout <<"Objective: " << dF << std::endl;
  //CPPUNIT_ASSERT( dF <= 1.e-6 );
}

