/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/option/testsuite/optioninterface.cpp
// Purpose:     Base class for testing option
// Author:      Ito33 Canada
// Created:     2005/06/14
// RCS-ID:      $Id: optioninterface.cpp,v 1.11 2006/08/22 10:10:43 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/date.h"
#include "ito33/sharedptr.h"
#include "ito33/dateutils.h"
#include "ito33/useexception.h"

#include "ito33/ihg/theoreticalmodel.h"

#include "ito33/ihg/volatility.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardrate.h"
#include "ito33/ihg/hazardrateflat.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/issuer.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/optiontype.h"
#include "ito33/finance/moneymarket.h"
#include "ito33/finance/numeraire_list.h"
#include "ito33/finance/option.h"

#include "optioninterface.h"
#include "testparam.h"
#include "utiltest.h"

extern const ito33::Error ITO33_UNEXPECTED;

namespace ito33 
{

namespace ihg
{

namespace test
{

  
void OptionInterface::SetMaturityDate(Date maturityDate)
{
 
  m_pOption = make_ptr( new finance::Option
                            (m_pOption->GetStrike(),
                             maturityDate,
                             m_pOption->GetOptionType(),
                             m_pOption->GetExerciseType() ) );

}

void OptionInterface::SetOptionType(finance::OptionType optionType)
{

  m_pOption =  make_ptr( new finance::Option
                             (m_pOption->GetStrike(),
                              m_pOption->GetMaturityDate(),
                              optionType,
                              m_pOption->GetExerciseType() ) );

}

void OptionInterface::SetExerciseType(finance::ExerciseType exerciseType)
{

  m_pOption =  make_ptr( new finance::Option
                             (m_pOption->GetStrike(),
                              m_pOption->GetMaturityDate(),
                              m_pOption->GetOptionType(),
                              exerciseType ) );

}

void OptionInterface::SetValuationDate(Date valuationDate)
{
  m_pSessionData =  make_ptr( new finance::SessionData
                                  (m_pSessionData->GetRateData(),
                                   m_pSessionData->GetEquity(),
                                   valuationDate) );
}

void OptionInterface::SetDividends(shared_ptr<finance::Dividends> pDividends)
{
  m_pSessionData->GetEquity()->SetDividends(pDividends);
}

void OptionInterface::SetYieldCurve(shared_ptr<finance::YieldCurve> pYC)
{
  m_pSessionData->SetYieldCurve(pYC);
}

void OptionInterface::SetForeignCurve(shared_ptr<finance::YieldCurve> pYC)
{
 m_pSessionData->GetEquity()->SetBorrowCurve(pYC);
}

void OptionInterface::SetVolatility(shared_ptr<Volatility> pVol)
{
  m_pModel->SetVolatility(pVol);
}

void OptionInterface::SetHazardRate(shared_ptr<HazardRate> pHR)
{
  m_pModel->SetHazardRate(pHR);
}

void OptionInterface::SetSpotSharePrice(double dSpot)
{
  m_pSessionData->GetEquity()->SetSpotSharePrice(dSpot);
}

void OptionInterface::SetStrike(double dStrike)
{
  m_pOption =  make_ptr( new finance::Option
                             (dStrike,
                              m_pOption->GetMaturityDate(),
                              m_pOption->GetOptionType(),
                              m_pOption->GetExerciseType()) );
}

void OptionInterface::SetDebugOutputFile(const char *debugOutputFile)
{
   m_pModel->SetDebugOutputFile(debugOutputFile);
}

void OptionInterface::Solve()
{
  m_pOption->SetSessionData(m_pSessionData);

  shared_ptr<finance::ModelOutput> pOutput =  m_pModel->Compute( *m_pOption );

  m_dPrice = pOutput->GetPrice();
  m_dGamma = pOutput->GetGamma();
  m_dDelta = pOutput->GetDelta();
  m_dTheta = pOutput->GetTheta();
}

void OptionInterface::CreateDebugOutputFile(const char *fileName, 
                                            std::string sTestTitle)
{
  char outputFileName[1024];
  OutPutFileName(fileName, sTestTitle.c_str(), outputFileName);
  SetDebugOutputFile( outputFileName );
  Solve();
}


void OptionInterface::GetTestParameters(double &dStep, double &dMax,
                 double &dMin, TestParam testParam, TestType testType)
{
  switch (testType)
  {
  case SPOT: 
    dStep = testParam.m_dStep * GetSpotSharePrice();
    dMax  = testParam.m_dSMax * GetSpotSharePrice();
    dMin  = testParam.m_dSMin * GetSpotSharePrice();
    break;

  case VOL:
    dStep = testParam.m_dVolStep;
    dMax  = testParam.m_dVolMax;
    dMin  = testParam.m_dVolMin;
    break;

  case HAZARDRATE:
    dStep = testParam.m_dHazardRateStep;
    dMax  = testParam.m_dHazardRateMax;
    dMin  = testParam.m_dHazardRateMin;
    break;

  case YIELDRATE:
    dStep = testParam.m_dYieldRateStep;
    dMax  = testParam.m_dYieldRateMax;
    dMin  = testParam.m_dYieldRateMin;
    break;

  case FOREIGNRATE:
    dStep = testParam.m_dForeignRateStep;
    dMax  = testParam.m_dForeignRateMax;
    dMin  = testParam.m_dForeignRateMin;
    break;

  case STRIKE:
    dStep = testParam.m_dStrikeStep * GetStrike();
    dMax  = testParam.m_dStrikeMax  * GetStrike();
    dMin  = testParam.m_dStrikeMin  * GetStrike(); 
    break;

  case MATURITY: 
    dStep = testParam.m_dMaturityStep;
    dMax  = testParam.m_dMaturityMax;
    dMin  = testParam.m_dMaturityMin;
    break;


  default:
    throw EXCEPTION_MSG
       (
         ITO33_UNEXPECTED,
         "Comparison test: this test does not exists."
       );
  
  }
}

void OptionInterface::SetTestParameter(double dCurrentParam, TestType testType)
{
 
  switch (testType)
    {
    case SPOT:
      SetSpotSharePrice( dCurrentParam );
      break;
      
    case VOL:
      SetVolatility( shared_ptr<Volatility>(new VolatilityFlat(dCurrentParam)) );
      break;

    case HAZARDRATE:
      SetHazardRate( shared_ptr<HazardRate>(new HazardRateFlat(dCurrentParam)) );
      break;
    
    case YIELDRATE:
      SetYieldCurve( shared_ptr<finance::YieldCurve>(new finance::YieldCurveFlat(dCurrentParam)) );
      break;

    case FOREIGNRATE:
      SetForeignCurve( shared_ptr<finance::YieldCurve>(new finance::YieldCurveFlat(dCurrentParam)) );
      break;

    case MATURITY:
      {
        SetValuationDate( Date(2003, Date::Feb, 1) );
        Date maturityDate = GetValuationDate();
        maturityDate.AddDays((int)dCurrentParam);
        // std::cout << MatDate << std::endl;
        SetMaturityDate( maturityDate );
        break;
      }

    case STRIKE:    
      SetStrike(dCurrentParam);
      break;

    default:
      throw EXCEPTION_MSG
        (
          ITO33_UNEXPECTED,
          "Comparison test: this test does not exists."
        );
    
    }//end switch statement
}

} //namespace test

} // namespace ihg

} // namespace ito33

