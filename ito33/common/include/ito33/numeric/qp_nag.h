///////////////////////////////////////////////////////////////////////////////
// File:             ito33/numeric/qp_nag.h                                 
// Purpose:          wrapper class for NAG QP calibrator                
// Created:          2005/05/10                                            
// RCS-ID:           $Id: bisecnewton.h,v 1.9 2004/10/05 09:13:38 pedro
// Copyright         (C) 2005 -  Trilemma LLP                                  
///////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_NUMERIC_QP_NAG_H_
#define _ITO33_NUMERIC_QP_NAG_H_

#include "ito33/array.h"

#include <nag.h>
#include <nage04.h>

#ifdef _MSC_VER
  #pragma comment(lib, "nagc.lib")
#endif

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

/**
   Class to wrapp the QP method for minimization problem.

   The objectif function is 1 / 2 X^T A X + C^T X
 */
class QPMinimizerNAG
{
public:

  /**
     Ctor takes the dimension and the bounds

     @param nNbX Number of parameters
     @param pdLowerBounds The lower bound, default to none
     @param pdUpperBounds The upper bound, default to none
   */
  QPMinimizerNAG
  (size_t nNbX, double* pdLowerBounds = 0, double* pdUpperBounds = 0)
  {
    m_nNbX = nNbX;
    m_pdLowerBounds = Array<double>(m_nNbX);
    m_pdUpperBounds = Array<double>(m_nNbX);

    if ( pdLowerBounds )
      for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
        m_pdLowerBounds[nIdx] = pdLowerBounds[nIdx];
    else
      for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
        m_pdLowerBounds[nIdx] = -1.e20;

    if ( pdUpperBounds )
      for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
        m_pdUpperBounds[nIdx] = pdUpperBounds[nIdx];
    else
      for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
        m_pdUpperBounds[nIdx] = 1.e20;

    m_piKX = Array<long>(m_nNbX);
  }

  /**
     Minimizes the problem defined by the given matrix and vector.

     Note that we can't force the input parameters to be const since
     the NAG function doesn't allow.

     @param ppdMatrix The matrix
     @param pdC The vector
     @param pdX The solution
   */
  double operator()(double** ppdMatrix, double *pdC, double *pdX)
  {
    Nag_E04_Opt nagOptions;  
    static NagError nagFail, nagFailOption;
    double dObjectif;

    nag_opt_init(&nagOptions);
    
    nagOptions.prob = Nag_QP2; 

    nagOptions.print_level = Nag_NoPrint;
   
    nagOptions.list = FALSE; 

    nagOptions.ftol = 1.e-12;

    nagOptions.rank_tol = 1.e-12;

    nag_opt_lin_lsq(m_nNbX, m_nNbX, 0, (double *)0, 0, 
                    m_pdLowerBounds.Get(), m_pdUpperBounds.Get(), 
                    pdC, (double *)0, ppdMatrix[0], m_nNbX, m_piKX.Get(),
                    pdX, &dObjectif, &nagOptions, NAGCOMM_NULL, &nagFail);

    nag_opt_free(&nagOptions, "all", &nagFailOption);
   
    return dObjectif;
  }


private:

  size_t m_nNbX;

  Array<double> m_pdLowerBounds;

  Array<double> m_pdUpperBounds;

  Array<long> m_piKX;

}; // class QPMinimizerNAG


} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_QP_NAG_H_
