/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/calibration/calibratorgeneral.cpp
// Purpose:     implementation of general calibration for the HG model
// Created:     2005/07/04
// RCS-ID:      $Id: calibratorgeneral.cpp,v 1.54 2006/08/19 23:44:26 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// 0 = simulated annealing
// 1 = nag
// 2 = mixed nag and simulated annealing
#define HG_CALIBRATION_METHOD 1

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/derivative.h"
#include "ito33/finance/option.h"
#include "ito33/finance/forwardoption.h"
#include "ito33/finance/derivatives.h"

#include "ito33/numeric/exception.h"

#include "ito33/hg/multioutput.h"

#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/underlyingprocess.h"

#include "hg/priceforwardoption.h"
#include "hg/translator.h"
#include "hg/calibratorgeneral.h"
#include "hg/numoutput.h"

#if HG_CALIBRATION_METHOD == 0
  #include "ito33/numeric/sa.h"
  #include "ito33/numeric/asa.h"
#elif HG_CALIBRATION_METHOD == 1
  #include "ito33/numeric/calibrator_nag.h"
#elif HG_CALIBRATION_METHOD == 2
  #include "ito33/numeric/sa.h"
  #include "ito33/numeric/asa.h"
  #include "ito33/numeric/calibrator_nag.h"
#endif

extern const ito33::finance::Error ITO33_CALIBRATION_FAIL;

namespace ito33
{

  using namespace finance;

namespace hg
{

void 
CalibratorGeneral::DoComputeObjectif(const double *pdX, double& dObjectif, 
                                   double* pdGradients, bool bComputeGradient)
{
  size_t nIdx;

  // Construct a pricing model from the current guess
  shared_ptr<UnderlyingProcess> pUnderlyingProcess( (*m_pTranslator)(pdX) );

  shared_ptr<TheoreticalModel> 
    pModel( new TheoreticalModel(pUnderlyingProcess) );
  
  shared_ptr<finance::ComputationalFlags> 
    pFlags(new finance::ComputationalFlags);

  pFlags->SetDiscretizationMethod(m_iDiscretizationMethod);

  // Make sure all problems use the same process params for mesh construction.
  // The params come from the initial guess.
  pModel->SetUnderlyingProcessForMesh(m_pUnderlyingProcessForMesh);

  // Need sensitivities to get the gradient
  if (bComputeGradient)
  {
    pFlags->SetSensitivityFlags( m_pTranslator->GetFlags() );

    for (nIdx = 0; nIdx < m_nNbUnknowns; nIdx++)
      pdGradients[nIdx] = 0.0;
  }

  pModel->SetExternalFlags(pFlags);

  // Actually price
  dObjectif = 0.0;

  if ( m_pDerivatives )
  {
    const Derivatives::Elements& elements( m_pDerivatives->GetAll() );
 
    Derivatives::Elements::const_iterator iter;

    for (iter = elements.begin(); iter != elements.end(); ++iter)
    {
      shared_ptr<finance::ModelOutput> output = pModel->Compute(*iter->first);

      double dComputedPrice = output->GetPrice();

      double dError = dComputedPrice - iter->first->GetMarketPrice();

      dObjectif += iter->second * dError * dError;

      shared_ptr<hg::NumOutput> 
        pNumOutput( static_pointer_cast<NumOutput>( output->GetNumOutput() ) );

      if (bComputeGradient)
      {
        std::vector<double> pdSensitivities = pNumOutput->GetSensitivities();

        for (nIdx = 0; nIdx < m_nNbUnknowns; nIdx++)
        {
          pdGradients[nIdx] += 2.0 * dError * pdSensitivities[nIdx]
                            * iter->second;
        }
      } // if computing gradients
    }
  }

  // Price option list using forward pricer
  if ( m_pForwardOption )
  {
    shared_ptr<MultiOutput> 
      forwardOutput( PriceForwardOption(*pModel, *m_pForwardOption) );
    
    if ( forwardOutput->IsSensitivityOnObjectif() )
    {
      dObjectif += 2. * forwardOutput->GetObjectif();

      if (bComputeGradient)
      {
        shared_ptr<NumOutput>
          pNumOutput( static_pointer_cast<NumOutput>
                      ( forwardOutput->GetNumOutput() ) );

        std::vector<double> pdSensitivities(pNumOutput->GetSensitivities());

        for (nIdx = 0; nIdx < m_nNbUnknowns; nIdx++)
        {
          pdGradients[nIdx] += 2.0 * pdSensitivities[nIdx];
        }
      } // if gradient requested
    }
    else
    { 
      std::vector<double> pdPrices;
      pdPrices = forwardOutput->GetPrices();

      std::vector< std::vector<double> > ppdSensitivities;

      if (bComputeGradient)
        ppdSensitivities = forwardOutput->GetMultiSensitivities();

      // loop over the options and update the objective function
      const ForwardOption::Elements& elements( m_pForwardOption->GetAll() );
      ForwardOption::Elements::const_iterator iter;
      size_t nIdxOp = 0;
      for (iter = elements.begin(); iter != elements.end(); ++iter, ++nIdxOp)
      {
        double dComputedPrice = pdPrices[nIdxOp];

        double dError = dComputedPrice - iter->m_pOption->GetMarketPrice();
        
        dObjectif += iter->m_dWeight * dError * dError;
        
        if (bComputeGradient)
        {
          for (nIdx = 0; nIdx < m_nNbUnknowns; nIdx++)
          {
            pdGradients[nIdx] += 2.0 * dError 
                               * ppdSensitivities[nIdxOp][nIdx]
                               * iter->m_dWeight;
          } // if gradient requested
        }
      }
    }
  } // if using forward pricer

  // Save the data if this is the best guess so far
  if (dObjectif < m_dObjectif)
  {
    m_dObjectif = dObjectif;

    for (nIdx = 0; nIdx < m_nNbUnknowns; nIdx++)
      m_pdFinalX[nIdx] = pdX[nIdx];
  }

} // ComputeObjectifWithDerivs



shared_ptr<UnderlyingProcess>
CalibratorGeneral::Calibrate(Translator* pTranslator, 
                             const finance::Derivatives& derivatives)
{
  // setup the translator
  m_pTranslator = pTranslator;

  // Get initial guess, bounds, etc
  InitializeArrays(*m_pTranslator);

  // Make the lists needed by the objective function
  ConstructDerivativeLists(derivatives, true);

  // Create an underlying process to be used by all meshes.
  // m_pdX should be the initial guess as set in InitializeArrays.
  m_pUnderlyingProcessForMesh = (*m_pTranslator)(m_pdX);

  // Do the calibration. Objective function does not compute gradients.
  RunCalibrator(true);

  // Return the calibrated model
  return GetLastCalibratedProcess();

} // CalibratorGeneralHG::Calibrate()


void CalibratorGeneral::RunCalibrator(bool bUserGradients)
{

  #if HG_CALIBRATION_METHOD == 1
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
    catch(const numeric::Exception&)
    {      
      throw EXCEPTION(ITO33_CALIBRATION_FAIL);
    }
 
  #elif HG_CALIBRATION_METHOD == 0
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
    catch(const numeric::Exception&)
    {
      throw EXCEPTION(ITO33_CALIBRATION_FAIL);
    }

    //std::cout <<"Func calls: " << solver.GetNbEvaluationFunction() <<std::endl;
  #elif HG_CALIBRATION_METHOD == 2
    // Use nag to get good initial guess for ASA. Then use ASA to get out
    // of local minima.  Finish with nag

    size_t nIdx;

    // Set to large value so the objective routine saves the new data
    m_dObjectif = 1.e99;

    double dTolerance = 1.e-6;
    size_t nMaxIterations = 10;
    numeric::CalibratorNAG calibrator(m_nNbUnknowns, &m_pdLowerBounds[0], 
                                      &m_pdUpperBounds[0], dTolerance,
                                      nMaxIterations);
    
    // Actually do the calibration
    try
    {    
      calibrator.Calibrate(*this, &m_pdX[0], bUserGradients);
    }
    catch(const numeric::Exception&)
    {      
      throw EXCEPTION(ITO33_CALIBRATION_FAIL);
    }

    // now use ASA
    numeric::ASA solver(m_pdLowerBounds, m_pdUpperBounds);
    solver.SetOutput("asa_out");
    solver.SetOutputDetailsTrue();
    solver.SetLimitGenerated(1000);

    for (nIdx = 0; nIdx < m_nNbUnknowns; nIdx++)
      m_pdX[nIdx] = m_pdFinalX[nIdx];

    try
    {    
      solver(*this, m_pdX);
    }
    catch(const numeric::Exception&)
    {
      throw EXCEPTION(ITO33_CALIBRATION_FAIL);
    }

    // finish with nag
    nMaxIterations = 100;
    numeric::CalibratorNAG calibrator2(m_nNbUnknowns, &m_pdLowerBounds[0], 
                                       &m_pdUpperBounds[0], dTolerance,
                                       nMaxIterations);
    
    for (nIdx = 0; nIdx < m_nNbUnknowns; nIdx++)
      m_pdX[nIdx] = m_pdFinalX[nIdx];

    try
    {    
      calibrator2.Calibrate(*this, &m_pdX[0], bUserGradients);
    }
    catch(const numeric::Exception&)
    {      
      throw EXCEPTION(ITO33_CALIBRATION_FAIL);
    }

  #endif
}


shared_ptr<UnderlyingProcess>
CalibratorGeneral::GetLastCalibratedProcess() const
{
  return (*m_pTranslator)(m_pdFinalX);
}


void CalibratorGeneral::operator () 
(const std::vector<double>& pdParams, double& dF)
{
  // This function is called by the simulated annealing code. Just pass it
  // on to the objective function
  ComputeObjectif(&pdParams[0], dF, NULL, false);
}


} // namespace hg

} // namespace ito33
