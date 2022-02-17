///////////////////////////////////////////////////////////////////////////////
// File:             ito33/numeric/calibrator_nag.h                          //
// Purpose:          wrapper class for NAG calibrator                        // 
// Created:          2004/12/8                                               //
// RCS-ID:           $Id: bisecnewton.h,v 1.9 2004/10/05 09:13:38 pedro E    //
// Copyright (C) Trilemma LLP  2004 -                                        //
///////////////////////////////////////////////////////////////////////////////

//#define STEPDEBUGNAG 

#ifndef _ITO33_NUMERIC_CALIBRATOR_NAG_H_
#define _ITO33_NUMERIC_CALIBRATOR_NAG_H_

#include <nag.h>
#include <nage04.h>

// include the nag dll
#ifdef _MSC_VER
  #pragma comment(lib, "nagc.lib")
#endif

#include "ito33/array.h"

namespace ito33
{

namespace numeric
{

#ifndef TRUE
  #define TRUE Nag_TRUE
#endif

#ifndef FALSE
  #define FALSE Nag_FALSE
#endif

template<class Target>
static void __stdcall 
ObjectifFunctionNAG(long /*nNbUnknown*/, 
                    double *pdX, double *pdObjectif, double* pdGradients, 
                    Nag_Comm *pComm)
{
  Target *pTarget = (Target *)pComm->p;

  double dObjectif;

  if ( pComm->flag == 2)
    pTarget->ComputeObjectif(pdX, dObjectif, pdGradients, true);
    //pTarget->ComputeObjectifWithDerivs(pdX, dObjectif, pdGradients);
  else
    pTarget->ComputeObjectif(pdX, dObjectif, pdGradients, false);

  *pdObjectif = dObjectif;
}

class CalibratorNAG
{
public:

  CalibratorNAG(size_t nNbParams, double *pdLowerBounds, double *pdUpperBounds,
                double dTolerance = 1.e-4, size_t nNbMaxIters = 40)
              : m_nNbParams(nNbParams), 
                m_dTolerance(dTolerance), m_nNbMaxIters(nNbMaxIters)
  {
    m_pdGradients = new double [m_nNbParams];
    
    m_pdLowerBounds = new double [m_nNbParams];
    m_pdUpperBounds = new double [m_nNbParams];
    for (size_t nIdx = 0; nIdx < m_nNbParams; nIdx++)
    {
      m_pdLowerBounds[nIdx] = pdLowerBounds[nIdx];
      m_pdUpperBounds[nIdx] = pdUpperBounds[nIdx];
    }
  }  

  void ResetBounds(double* pdLowerBounds, double* pdUpperBounds)
  {
    for (size_t nIdx = 0; nIdx < m_nNbParams; nIdx++)
    {
      m_pdLowerBounds[nIdx] = pdLowerBounds[nIdx];
      m_pdUpperBounds[nIdx] = pdUpperBounds[nIdx];
    }
  }

  template <class Target>
  void Calibrate(Target& target, double* pdX, bool bUserDerivs = false)
  {
    m_iErrorCode = 0;

    double dObjectif;
    
    Nag_E04_Opt nagOptions;   
    
    Nag_Comm nagComm;
    
    static NagError nagFail, nagFailOption;
  
    // Nag options
    e04xxc(&nagOptions);

    nagOptions.local_search = TRUE;

    nagOptions.optim_tol = m_dTolerance;
  
    nagOptions.max_iter = m_nNbMaxIters; 

    //nagOptions.linesearch_tol = 0.9;

#ifdef STEPDEBUGNAG
    //nagOptions.deriv_check = TRUE;
    nagOptions.deriv_check = FALSE;
     
    nagOptions.print_level = Nag_Soln_Iter_Full;
   
    nagOptions.list = TRUE;      
#else
    nagOptions.deriv_check = FALSE;  
      
    nagOptions.print_level = Nag_NoPrint;
     
    nagOptions.list = FALSE; 
#endif  
  
    nagComm.p = (void *)(&target);

    if ( bUserDerivs )
    {
/*
      nag_opt_nlp
      ( long(m_nNbParams), 0, 0, NULL, 0,
        m_pdLowerBounds.Get(), m_pdUpperBounds.Get(),
        ObjectifFunctionNAG<Target>,
        NULLFN,
        pdX, &dObjectif, m_pdGradients.Get(),
        &nagOptions, &nagComm, &nagFail
      );
*/

      nag_opt_bounds_deriv
      (
        long(m_nNbParams), ObjectifFunctionNAG<Target>, 
        Nag_Bounds, m_pdLowerBounds.Get(), m_pdUpperBounds.Get(),
        pdX, &dObjectif, m_pdGradients.Get(), 
        &nagOptions, &nagComm, &nagFail
      );

    }
    else
    {
      nag_opt_bounds_no_deriv
      (
        long(m_nNbParams), ObjectifFunctionNAG<Target>, 
        Nag_Bounds, m_pdLowerBounds.Get(), m_pdUpperBounds.Get(),
        pdX, &dObjectif, m_pdGradients.Get(), 
        &nagOptions, &nagComm, &nagFail
      );
    }

    e04xzc(&nagOptions, "all", &nagFailOption);
    
    m_iErrorCode = nagFail.code;

    m_errorMessage = nagFail.message;

    m_dObjectif = dObjectif;
  }
  
  int GetErrorCode() const { return m_iErrorCode; }

  char* GetErrorMessage() const { return m_errorMessage; }

  double GetObjectif() const { return m_dObjectif; }


private:

  size_t m_nNbParams;

  double m_dTolerance;

  size_t m_nNbMaxIters;

  Array<double> m_pdLowerBounds;
  Array<double> m_pdUpperBounds;

  Array<double> m_pdGradients;

  int m_iErrorCode;

  char* m_errorMessage;

  double m_dObjectif;

}; // class CalibratorNAG


} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_CALIBRATOR_NAG_H_
