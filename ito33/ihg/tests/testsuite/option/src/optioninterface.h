/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/testsuite/option/src/optioninterface.h
// Purpose:     Base class for testing IHG projects
// Author:      Ito33Canada 
// Created:     2005/06/09
// RCS-ID:      $Id: optioninterface.h,v 1.9 2006/08/22 10:10:43 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/tests/testsuite/option/src/optioninterface.h
    @brief Base class for testing IHG projects

 
*/

#ifndef _IHG_TESTS_TESTSUITE_OPTION_SRC_OPTIONINTERFACE_H_
#define _IHG_TESTS_TESTSUITE_OPTION_SRC_OPTIONINTERFACE_H_

#include "ito33/string.h"
#include "ito33/date.h"
#include "ito33/sharedptr.h"
#include "ito33/common.h"

#include "ito33/finance/option.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/optiontype.h"
#include "ito33/finance/exercisetype.h"
#include "ito33/finance/yieldcurve.h"
#include "ito33/finance/issuer.h"
#include "ito33/finance/moneymarket.h"
#include "ito33/finance/derivative.h"


#include "ito33/ihg/theoreticalmodel.h"

#include "testparam.h"

namespace ito33 
{

namespace finance
{
  class Dividends;
}

namespace ihg
{
  class Volatility;
  class HazardRate;

namespace test
{


class OptionInterface
{

protected:
  shared_ptr<finance::SessionData> m_pSessionData;
  shared_ptr<finance::Option>  m_pOption;
  shared_ptr<TheoreticalModel> m_pModel;

  double m_dPrice;
  double m_dDelta;
  double m_dGamma;
  double m_dTheta;

public:
  OptionInterface(shared_ptr<finance::SessionData> pSessionData,
    shared_ptr<finance::Option> pOption,
    shared_ptr<TheoreticalModel> pModel)
    : m_pSessionData(pSessionData),
      m_pOption(pOption),
      m_pModel(pModel),
      m_dPrice(0.0), m_dDelta(0.),
      m_dGamma(0.0), m_dTheta(0.)
  {
  }

  OptionInterface(shared_ptr<finance::SessionData> pSessionData,
    shared_ptr<TheoreticalModel> pModel)
    :m_pSessionData(pSessionData),
      m_pModel(pModel),
      m_dPrice(0.0), m_dDelta(0.),
      m_dGamma(0.0), m_dTheta(0.)
  {

  }

  shared_ptr<finance::Option> GetOption() const
  {
    return m_pOption;
  }

  shared_ptr<TheoreticalModel> GetModel() const
  {
    return m_pModel;
  }

  virtual ~OptionInterface() {}

  double GetPrice() const
  {
    return m_dPrice;
  }

  double GetDelta() const 
  {
    return m_dDelta; 
  }

  double GetGamma() const
  {
    return m_dGamma;
  }

  double GetTheta() const
  {
    return m_dTheta;
  }

  double GetSpotSharePrice() const
  {
    return m_pSessionData->GetEquity()->GetSpotSharePrice();
  }

  virtual double GetStrike() const
  {
    return m_pOption->GetStrike();
  }


  virtual shared_ptr<finance::Derivative> GetDerivative() const
  {
    return m_pOption;
  }

  Date GetValuationDate() const
  {
    return m_pSessionData->GetValuationDate();
  }

  shared_ptr< finance::YieldCurve > GetYieldCurve() const
  {
    return m_pSessionData->GetYieldCurve();
  }

  shared_ptr< finance::YieldCurve > GetForeignCurve() const
  {
    return m_pSessionData->GetEquity()->GetBorrowCurve();
  }

  virtual Date GetMaturityDate() const
  {
    return m_pOption->GetMaturityDate();
  }

  finance::OptionType GetOptionType() const
  {
    return m_pOption->GetOptionType();
  }

  void SetMaturityDate(Date maturityDate);

  void SetValuationDate(Date valuationDate);

  void SetOptionType(finance::OptionType optionType);

  void SetExerciseType(finance::ExerciseType exerciseType);

  void SetDividends(shared_ptr<finance::Dividends> pDividends);

  void SetYieldCurve(shared_ptr<finance::YieldCurve> pYC);

  void SetForeignCurve(shared_ptr<finance::YieldCurve> pYC);

  void SetVolatility(shared_ptr<Volatility> pVol);

  void SetSpotSharePrice(double dSpot);

  void SetHazardRate(shared_ptr<HazardRate> pHR);

  void SetStrike(double dStrike);

  void SetDebugOutputFile(const char *debugOutputFile);

  void CreateDebugOutputFile(const char* fileName, std::string sTestTitle);

  void GetTestParameters(double &dStep, double &dMax, double &dMin,
                         TestParam testParam, TestType testType);

  void SetTestParameter(double dCurrentParam, TestType testType);

  virtual void Solve();

};



} //namespace test

} // namespace ihg

} // namespace ito33



#endif
