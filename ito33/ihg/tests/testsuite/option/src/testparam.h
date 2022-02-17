/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/option/testsuite/testparam.h
// Purpose:     define the different parameters for the tests
// Author:      Yann d'Halluin
// Created:     2004/05/02
// RCS-ID:      $Id: testparam.h,v 1.2 2005/06/27 16:20:57 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/tests/option/test/testsuite/testparam.h
    @brief set of control parameters to be used for different tests
**/
#ifndef _IHG_TESTS_TESTSUITE_TESTPARAM_H_
#define _IHG_TESTS_TESTSUITE_TESTPARAM_H_

#include "utiltest.h"
#include "ito33/xml/write.h"
#include <fstream>
#include "ito33/constants.h"

namespace ito33
{

namespace ihg
{

namespace test
{

  /**
  @brief define different control parameters to
         be used for the differents tests
*/
  class TestParam {

  public:
    ///Numerical check
    double m_dEPSCheck;

     /// [S,S+dStep] stock
    double m_dStep;

    /// [,dSMax*Spot]; as a percentage of the current spot
    double m_dSMax;

    // [dSMin*Spot,]; as a percentage of the current spot
    double m_dSMin;
  
     /// [dVol,dVol+dStep] volatility
    double m_dVolStep;

    /// [,dVolMax]; as a maximum value for the volatility
    double m_dVolMax;

    // [dVolMin,]; as a minimum value for the  volatility
    double m_dVolMin;

    ///[dYield,dYield+Step]
    double m_dYieldRateStep;

    /// [,dMax]; as a maximum value for the yield rate
    double m_dYieldRateMax;

    // [dMin,]; as a minimum value for yield rate
    double m_dYieldRateMin;

    ///[dForeignRate,dForeignRate+Step]
    double m_dForeignRateStep;

    /// [,dMax]; as a maximum value for the foreign rate
    double  m_dForeignRateMax;

    // [dMin,]; as a minimum value for foreign rate
    double m_dForeignRateMin;

    ///[dHazard,dHazardRate+Step]
    double m_dHazardRateStep;

    /// [,dMax]; as a maximum value for the hazard rate
    double  m_dHazardRateMax;

    /// [dMin,]; as a minimum value for hazard rate
    double m_dHazardRateMin;

    ///number of iteration for convergence testing
    size_t m_nIterConvergence;

    ///number of grid point to start with
    size_t m_nGridSize;

    ///number of timestep to start with
    size_t m_nTimeSize;

    /// maturity step in years
    double m_dMaturityStep;

    /// max maturity
    double m_dMaturityMax;

    /// min maturity
    double m_dMaturityMin;

    /// strike steps
    double m_dStrikeStep;

    /// Maximum strike
    double m_dStrikeMax;

    /// Minimum strike
    double m_dStrikeMin;

    ///File currently being tested
    char buffer[1024];

    public:
    /// Constructor
    TestParam() 
    {
      m_dEPSCheck = 1.e-2;

      m_dStep = .1;
      m_dSMax = 2.;
      m_dSMin = .1;

      m_dVolStep = .1;
      m_dVolMax = 2.;
      m_dVolMin = .01;
    
      m_dYieldRateStep = .02;
      m_dYieldRateMin  = 0;
      m_dYieldRateMax  = .2;

      m_dForeignRateStep = .02;
      m_dForeignRateMin  = 0;
      m_dForeignRateMax  = .2;

      m_dHazardRateStep = .1;
      m_dHazardRateMin  = 0;
      m_dHazardRateMax  = 1.0;

      m_nIterConvergence = 5;
      m_nGridSize = 101;
      m_nTimeSize = 51;

      m_dMaturityStep = 30*4; //in terms of days
      m_dMaturityMax  = 4*INVERSEONEDAY; 
      m_dMaturityMin  = 1;

      m_dStrikeStep = .1;
      m_dStrikeMax = 2.;
      m_dStrikeMin = .1;

      buffer[0]=0;
    } //end Constructor

    void SetFileName(const char *fileName)
    {
      buffer[0]=0;
      sprintf(buffer,fileName);
    }

    char *GetFileName()
    {
      return buffer;
    }

   
  };
} //end namespace test
} //end namespace ihg
} //end namespace ito33


#endif
