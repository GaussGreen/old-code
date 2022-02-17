/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/analysisdate_surface.h
// Purpose:     Tests on Analysis date and ComputeSurface flags
// Created:     2004/11/22
// RCS-ID:      $Id: analysisdate_surface.h,v 1.1 2004/11/29 09:27:02 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _IHG_TESTS_ANALYSISDATE_SURFACE_H_
#define _IHG_TESTS_ANALYSISDATE_SURFACE_H_

#include <cstddef>

#include "ito33/sharedptr.h"
#include "ito33/string.h"

namespace ito33
{
  class Date;

namespace finance
{
  class Derivative;
}

namespace ihg
{

  class TheoreticalModel;

  const std::string 
    CanCompute[] = {"Price", "Delta", "Gamma", "Theta", "Rho", "Vega", "Fugit"};  
  
  // Price, Delta, Gamma, Theta, Rho, Vega, Fugit
  const size_t NUMBEROFGREEK = SIZEOF(CanCompute);

// ----------------------------------------------------------------------------
// Analysis date tests
// ----------------------------------------------------------------------------

class AnalysisDateTest 
{

public:

  AnalysisDateTest(const TheoreticalModel& model,
                   const finance::Derivative* pDerivative)
                 : m_model(model), m_pDerivative(pDerivative)
  {
    // for option and cds, only fugit is not computed for now
    for (size_t nIdx = 0; nIdx < NUMBEROFGREEK - 1; nIdx++)
      m_pbShouldHaveDatas[nIdx] = true; 

    m_pbShouldHaveDatas[NUMBEROFGREEK - 1] = false;

    // by default, nothing is checked yet
    for (size_t nIdx = 0; nIdx < NUMBEROFGREEK; nIdx++)
      m_pbHasDatas[nIdx] = false;
  }

  void Report();


protected:

  bool NoDataForAnalysisDateBeforeValuationDate();
  bool NoDataForAnalysisDateAtMaturityDate();

  void DataAvailability();
  bool ConsistencyAtValuationDate();
  
  bool m_pbShouldHaveDatas[NUMBEROFGREEK];
  bool m_pbHasDatas[NUMBEROFGREEK];
  

private:

  TheoreticalModel m_model;
  const finance::Derivative* m_pDerivative;

  NO_COPY_CLASS(AnalysisDateTest);
};

class SurfaceTest 
{

public:

  SurfaceTest(const TheoreticalModel& model,
              const finance::Derivative* pDerivative)
            : m_model(model), m_pDerivative(pDerivative)
  {
    // for option and cds, only fugit is not computed for now
    for (size_t nIdx = 0; nIdx < NUMBEROFGREEK - 1; nIdx++)
      m_pbShouldHaveDatas[nIdx] = true; 

    m_pbShouldHaveDatas[NUMBEROFGREEK - 1] = false;

    // by default, nothing is checked yet
    for (size_t nIdx = 0; nIdx < NUMBEROFGREEK; nIdx++)
      m_pbHasDatas[nIdx] = false;
  }

  void Report();


protected:

  bool ConsistencyAtValuationDate();
  void DataAvailability();
  
  bool m_pbShouldHaveDatas[NUMBEROFGREEK];
  bool m_pbHasDatas[NUMBEROFGREEK];


private:

  TheoreticalModel m_model;
  const finance::Derivative* m_pDerivative;
  
  NO_COPY_CLASS(SurfaceTest);
};


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_TESTS_ANALYSISDATE_SURFACE_H_
