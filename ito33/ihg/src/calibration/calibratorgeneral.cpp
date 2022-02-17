/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/calibratorgeneral.cpp
// Purpose:     implementation of general calibration for the IHG model
// Created:     2005/07/04
// RCS-ID:      $Id: calibratorgeneral.cpp,v 1.14 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// Set to use NAG for the minimization. Otherwise use simulated annealing.
//#define USE_NAG_IHG

#include "ito33/sharedptr.h"
#include "ito33/exception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/derivative.h"
#include "ito33/finance/derivatives.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/numeric/exception.h"

#include "ito33/ihg/theoreticalmodel.h"

#include "ihg/translator.h"
#include "ihg/calibratorgeneral.h"

#ifdef USE_NAG_IHG
  #include "ito33/numeric/calibrator_nag.h"
#else
  #include "ito33/numeric/sa.h"
  #include "ito33/numeric/asa.h"
#endif

extern const ito33::finance::Error ITO33_CALIBRATION_FAIL;


namespace ito33
{

namespace ihg
{


void 
CalibratorGeneral::DoComputeObjectif
                   ( const double *pdX, double& dObjectif, 
                     double* pdGradients, bool bComputeGradient )
{

  // TODO: Add support for forward pricer
  ASSERT_MSG( !m_pForwardOption, 
              "Forward pricing not yet supported in IHG general calibration");

  // Get a pricing model from the current guess
  shared_ptr<TheoreticalModel> pModel( (*m_pTranslator)(pdX) );
  pModel->SetExternalFlagsToDefaults();

  size_t nIdx;
  
  dObjectif = 0.0;

  const finance::Derivatives::Elements& elements( m_pDerivatives->GetAll() );
  finance::Derivatives::Elements::const_iterator iter;

  for (iter = elements.begin(); iter != elements.end(); ++iter)
  {
    shared_ptr<finance::ModelOutput> output = pModel->Compute(*iter->first);

    double dComputedPrice = output->GetPrice();

    double dError = dComputedPrice - iter->first->GetMarketPrice();

    dObjectif += iter->second * dError * dError;
  }

  if (bComputeGradient)
  {
    // Compute the gradient by finite differences. Assume a positive
    // shift is always valid (i.e. assume upper bounds for calibration 
    // are less than than upper bounds in underlying classes)
    double dShift = 1.e-7;

    size_t nNbParams = m_pTranslator->GetNbParameters();
    std::vector<double> pdXShift(nNbParams);
    size_t nIdx = 0;
    for (nIdx = 0; nIdx < nNbParams; nIdx++)
      pdXShift[nIdx] = pdX[nIdx];

    for (nIdx = 0; nIdx < nNbParams; nIdx++)
    {
      pdXShift[nIdx] += dShift;
      shared_ptr<ihg::TheoreticalModel> pModel( (*m_pTranslator)(&pdXShift[0]));
      pdXShift[nIdx] -= dShift;

      double dObjectifShift = 0.0;

      for (iter = elements.begin(); iter != elements.end(); ++iter)
      {
        shared_ptr<finance::ModelOutput> output = pModel->Compute(*iter->first);
  
        double dComputedPrice = output->GetPrice();

        double dError = dComputedPrice - iter->first->GetMarketPrice();

        dObjectifShift += iter->second * dError * dError;
      }

      pdGradients[nIdx] = (dObjectifShift - dObjectif) / dShift;
    } // loop over parameters
  } // if sensitivity needed

  // Save the data if this is the best guess so far
  if (dObjectif < m_dObjectif)
  {
    m_dObjectif = dObjectif;

    for (nIdx = 0; nIdx < m_nNbUnknowns; nIdx++)
      m_pdFinalX[nIdx] = pdX[nIdx];
  }

} // ComputeObjectif


shared_ptr<ihg::TheoreticalModel>
CalibratorGeneral::Calibrate
(Translator* pTranslator, const finance::Derivatives& derivatives)
{
  // setup the translator
  m_pTranslator = pTranslator;

  // Get initial guess, bounds, etc
  InitializeArrays(*m_pTranslator);

  // Make the lists needed by the objective function. Do noy use forward 
  // option pricer
  ConstructDerivativeLists(derivatives, false);

  // Do the calibration. Objective function computes the gradients.
  RunCalibrator(false);

  // Return the calibrated model
  return GetLastCalibratedProcess();

} // CalibratorGeneralIHG::Calibrate()


void CalibratorGeneral::RunCalibrator(bool bUserGradients)
{

  #ifdef USE_NAG_IHG
    // Construct the nag calibrator
    double dTolerance = 1.e-6;
    size_t nMaxIterations = 200;
    numeric::CalibratorNAG calibrator(m_nNbUnknowns, &m_pdLowerBounds[0], 
                                      &m_pdUpperBounds[0], dTolerance,
                                      nMaxIterations);

    // Set to large value so the objective routine saves the new data
    m_dObjectif = 1.e99;

    // Actually do the calibration
    try
    {    
      calibrator.Calibrate(*this, &m_pdX[0], bUserGradients);
    }
    catch(ito33::numeric::Exception)
    {      
      throw EXCEPTION(ITO33_CALIBRATION_FAIL);
    }
  
  #else
    // Construct the simulated annealing calibrator
    //double dTolerance = 1.e-5;
    
    // just to shut up complier warning
    bUserGradients = false;

    numeric::ASA solver(m_pdLowerBounds, m_pdUpperBounds);
    //solver.SetOutput("asa_out");
    //solver.SetOutputDetailsTrue();
    //solver.SetSamplingNumber(1);
    //solver.SetCostParameterScaleRatio(1.);
    //solver.SetReannealingFrequency(10);
    //solver.SetLimitGenerated(5000);
    solver.SetLimitAcceptance(200); 

    // Set to large value so the objective routine saves the new data
    m_dObjectif = 1.e99;

    // Actually do the calibration
    try
    {    
      solver(*this, m_pdX);
    }
    catch(ito33::numeric::Exception)
    {
      throw EXCEPTION(ITO33_CALIBRATION_FAIL);
    }

    //std::cout <<"Func calls: " << solver.GetNbEvaluationFunction() <<std::endl;
  #endif
}


shared_ptr<ihg::TheoreticalModel>
CalibratorGeneral::GetLastCalibratedProcess() const
{
  return (*m_pTranslator)(m_pdFinalX);
}


void CalibratorGeneral::operator () 
(const std::vector<double>& pdParams, double& dF)
{
  // This function is called by the simulated annealing code. Just pass it
  // on to the objective function
  ComputeObjectif(&pdParams[0], dF, NULL, 0);
}


} // namespace ihg

} // namespace ito33
