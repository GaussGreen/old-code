/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/option/testsuite/utiltest.h
// Purpose:     util function for tests
// Author:      Yann d'Halluin David Pooley
// Created:     2004/05/02
// RCS-ID:      $Id: utiltest.h,v 1.3 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/tests/option/testsuite/utiltest.h
    @brief set of functions to help when writting tests results
**/

#ifndef _IHG_TESTS_TESTSUITE_UTILTEST_H_
#define _IHG_TESTS_TESTSUITE_UTILTEST_H_

#include "ito33/beforestd.h"
#include <vector>
#include "ito33/afterstd.h"

#include "ito33/string.h"
#include "ito33/sharedptr.h"
#include "ito33/date.h"

#include "ito33/useexception.h"
#include "ito33/xml/write.h"

namespace ito33
{
  namespace XML 
  {
    /// allow to display spot and a value next to each other
    struct ParameterOneValue
    {
      ParameterOneValue() {}

      void init(std::string sParameterTag, std::string sValueTag)
      {
        m_sParameterTag = sParameterTag;
        m_sValueTag     = sValueTag;
      }

      std::string m_sParameterTag;
      std::string m_sValueTag;
      std::vector<double> m_pdParameter;
      std::vector<double> m_pdValue;
    };
  
  template <> Tag MakeNamedTag(const char *name, 
    const ParameterOneValue &parameterOneValue, Tag& parent);

    struct ParameterTwoValues
    {
      ParameterTwoValues() {}

      void init(std::string parameterTag, 
        std::string value1Tag, std::string value2Tag)
      {
        m_parameterTag = parameterTag;
        m_value1Tag    = value1Tag;
        m_value2Tag    = value2Tag; 
      }

      std::string m_parameterTag;
      std::string m_value1Tag;
      std::string m_value2Tag;

      std::vector<double> m_pdParameter;
      std::vector<double> m_pdValue1;
      std::vector<double> m_pdValue2;
    };


   template <> Tag MakeNamedTag(const char *name, 
     const ParameterTwoValues &parameterTwoValues,  Tag& parent);

  } //end XML

namespace finance
{
  class Dividends;
  class YieldCurve;
}

  
namespace ihg
{

  class Volatility;
  class HazardRate;

namespace test
{

  //class TestParam;

  enum HazardRateType { DECAYHR,FLATHR, POWERHR, TIMEONLYHR};
  enum YieldType      { FLAT, ZEROCOUPON, NOYIELD};
  enum DividendType   { NODIVIDEND, CASH, YIELD};
  enum VolatilityType { FLATVOL};
  enum TestType       { VOL, SPOT ,STRIKE, MATURITY, HAZARDRATE, YIELDRATE, 
                        DIVIDENDRATE, FOREIGNRATE };
  enum TestDirection  { INCREASE, DECREASE };
  enum SpotAt         { NORMAL, ATK, LESSK, GREATERK};
  enum RATIO          { LINEAR, QUADRATIC};

/// Select a volatility
ito33::shared_ptr<ito33::ihg::Volatility> GetVolatility(
      VolatilityType volatilityType);

/// Select an hazard rate
ito33::shared_ptr<ito33::ihg::HazardRate> GetHazardRate(
       HazardRateType hazardRateType,double dSpot,ito33::Date ValuationDate);

///  Select a dividend type
ito33::shared_ptr<ito33::finance::Dividends> GetDividend(
       DividendType dividendType,ito33::Date date,double Spot);


/// Select a yield curve
ito33::shared_ptr<ito33::finance::YieldCurve> GetYieldCurve(
     YieldType yieldType,ito33::Date date);

    
/// Select a foreign curve
ito33::shared_ptr<ito33::finance::YieldCurve> GetForeignCurve(
      YieldType foreignType,ito33::Date date);

///Generate xml files
void GenerateXMLTestingFiles(size_t &nFileCounter);


inline void OutPutFileName(const char* filename, const char *testname,
                           char *outputFileName)
{
  std::string strFileName(filename);
  std::string::size_type nBegin = strFileName.find_last_of('/');

  if(nBegin == std::string::npos)
    nBegin = 1;
  else
    nBegin++;

  std::string::size_type nEnd = strFileName.find_last_of('.');

  std::string t = strFileName.substr(nBegin, nEnd - nBegin);

  sprintf(outputFileName,"failedcases/%s_%s.xml", t.c_str(), testname);

  
}

} //end namespace test
} //end namespace ihg
} //endnamespace ito33


#endif
