/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/calibration/calibrator_cdsrecovery.cpp
// Purpose:     general calibration for the HG model with CDS recovery
// Created:     2005/07/24
// RCS-ID:      $Id: calibrator_cdsrecovery.cpp,v 1.15 2006/08/19 23:44:26 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file tests/tmp/market/calibrator_cdsrecovery.cpp
   @brief General calibration for the HG model with CDS recovery
 */

// Set to use NAG for the minimization. Otherwise use simulated annealing.
#define USE_NAG_HG_RECOVERY

#include <cmath>

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/derivative.h"
#include "ito33/finance/option.h"
#include "ito33/finance/cds.h"
#include "ito33/finance/forwardoption.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/eds.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/impliedcdsspreads.h"
#include "ito33/finance/impliededsspreads.h"
#include "ito33/finance/derivatives.h"

#include "ito33/numeric/exception.h"

#include "ito33/hg/modeloutput.h"
#include "ito33/hg/multioutput.h"
#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/underlyingprocess.h"

#include "hg/numoutput.h"
#include "hg/priceforwardoption.h"
#include "hg/translator.h"
#include "hg/calibrator_cdsrecovery.h"

#ifdef USE_NAG_HG_RECOVERY
  #include "ito33/numeric/calibrator_nag.h"
#else
  #include "ito33/numeric/sa.h"
  #include "ito33/numeric/asa.h"
#endif

extern const ito33::finance::Error ITO33_CALIBRATION_FAIL;

namespace ito33
{

namespace hg
{

void 
CalibratorCDSRecovery::DoComputeObjectif
                       (const double *pdX, double& dObjectif, 
                        double* pdGradients, bool bComputeGradient)
{

  size_t nIdx;

  // Construct a pricing model from the current guess
  shared_ptr<UnderlyingProcess> pUnderlyingProcess( (*m_pTranslator)(pdX) );

  shared_ptr<TheoreticalModel> 
    pModel( new TheoreticalModel(pUnderlyingProcess) );

  // Make sure all problems use the same process params for mesh construction.
  // The params come from the initial guess.
  pModel->SetUnderlyingProcessForMesh(m_pUnderlyingProcessForMesh);

  // Need sensitivities to get the gradient
  if (bComputeGradient)
  {
    shared_ptr<finance::ComputationalFlags> 
      flags(new finance::ComputationalFlags);
    flags->SetSensitivityFlags( m_pTranslator->GetFlags() );
    pModel->SetExternalFlags(flags);

    for (nIdx = 0; nIdx < m_nNbUnknowns; nIdx++)
      pdGradients[nIdx] = 0.0;
  }

  // Actually price
  dObjectif = 0.0;
  
  const finance::Derivatives::Elements& elements( m_pDerivatives->GetAll() );
  finance::Derivatives::Elements::const_iterator iter;

  for (iter = elements.begin(); iter != elements.end(); ++iter)
  {

    shared_ptr<finance::Derivative> pDeriv = iter->first;
  
    if (dynamic_cast<finance::CDS *>( pDeriv.get() ) != 0)
    {
      // cds: objective can be based on implied spreads
      ObjectiveForCDS(pdX, dObjectif, pdGradients, bComputeGradient, 
                      pDeriv, iter->second, pModel);

    } 
    else if (dynamic_cast<finance::EDS *>( pDeriv.get() ) != 0)
    {      
      // eds: objective can be based on implied spreads
      ObjectiveForEDS(pdX, dObjectif, pdGradients, bComputeGradient, 
                      pDeriv, iter->second, pModel);

    }
    else
    {
      // Usual case: objective function based on price
      shared_ptr<finance::ModelOutput> output = pModel->Compute(*pDeriv);

      double dComputedPrice = output->GetPrice();

      const double dMarketPrice = pDeriv->GetMarketPrice();

      double dError = dComputedPrice - dMarketPrice;

      dObjectif += iter->second * dError * dError;

      if (bComputeGradient)
      {
        shared_ptr<NumOutput> 
          pNO( static_pointer_cast<NumOutput>( output->GetNumOutput() ) );

        std::vector<double> pdSensitivities = pNO->GetSensitivities();

        for (nIdx = 0; nIdx < m_nNbUnknowns - 1; nIdx++)
        {
          pdGradients[nIdx] += 2.0 * dError * pdSensitivities[nIdx]
                             * iter->second;
        }
      } // if computing gradients

    } // usual deriv (not cds or eds)
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
          pNO( static_pointer_cast<NumOutput>(forwardOutput->GetNumOutput()) );

        std::vector<double> pdSensitivities(pNO->GetSensitivities());

        for (nIdx = 0; nIdx < m_nNbUnknowns - 1; nIdx++)
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
      const finance::ForwardOption::Elements& 
        elements( m_pForwardOption->GetAll() );

      finance::ForwardOption::Elements::const_iterator iter;

      size_t nIdxOp = 0;
      for (iter = elements.begin(); iter != elements.end(); ++iter, ++nIdxOp)
      {
        const double dMarketPrice = iter->m_pOption->GetMarketPrice();

        double dComputedPrice = pdPrices[nIdxOp];

        double dError = dComputedPrice - dMarketPrice;
        
        dObjectif += iter->m_dWeight * dError * dError;
        
        if (bComputeGradient)
        {
          for (nIdx = 0; nIdx < m_nNbUnknowns - 1; nIdx++)
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


void CalibratorCDSRecovery::ObjectiveForCDS(
                                 const double *pdX, 
                                 double& dObjectif,
                                 double* pdGradients, 
                                 bool bComputeGradient,
                                 shared_ptr<finance::Derivative>& pDerivative,
                                 double dWeight,
                                 shared_ptr<TheoreticalModel>& pModel)
{
  // If CDS, create a new CDS with the specified recovery.  Also,
  // compare spreads instead of prices

  // Convert back to original CDS
  shared_ptr<finance::CDS> 
    pOriginalCDS( dynamic_pointer_cast<finance::CDS>(pDerivative) );

  // Create a new CDS
  double dRecoveryRate = pdX[m_nNbUnknowns - 1];

  shared_ptr<finance::CDS> 
    pNewCDS(new finance::CDS(dRecoveryRate, pOriginalCDS->GetSpreadStream() ));
    
  pNewCDS->SetSessionData(pOriginalCDS->GetSessionData());
  pNewCDS->SetMarketPrice(pOriginalCDS->GetMarketPrice());

  if ( m_bUseSpreads )
  {
    // Extract info needed to compute implied spreads, then compute 
    // the implied spread
    Date valuationDate = pNewCDS->GetSessionData()->GetValuationDate();
    Date firstDate     = pNewCDS->GetSpreadStream()->GetFirstPaymentDate();
    Date::DayCountConvention dayCount = 
      pNewCDS->GetSpreadStream()->GetDayCountConvention();
    finance::Frequency freq = pNewCDS->GetSpreadStream()->GetPaymentFrequency();

    finance::ImpliedCDSSpreads impliedSpreads(valuationDate,
                                              firstDate, 
                                              dayCount,
                                              freq, 
                                              dRecoveryRate);

    std::vector<Date> pDates(1);
    pDates[0] = pNewCDS->GetMaturityDate();

    std::vector<double> pdSpreads = 
      impliedSpreads.Compute(pModel, pNewCDS->GetSessionData(), pDates);

    // Compute error based on market spread and implied spread
    double dMarketSpread = pNewCDS->GetSpreadStream()->GetAnnualPaymentAmount();

    double dComputedSpread = pdSpreads[0];

    double dError = dComputedSpread - dMarketSpread;

    dObjectif += dWeight * dError * dError;

    // Compute gradients, if requested
    if (bComputeGradient)
    {
      // Use finite differences for the spread gradients
      double dShift = 1.e-6;

      // Copy the param vector (since it is const)
      std::vector<double> pdXCopy(m_nNbUnknowns, 0.0);
      for (size_t nIdx = 0; nIdx < m_nNbUnknowns; nIdx++)
        pdXCopy[nIdx] = pdX[nIdx];

      for (size_t nIdx = 0; nIdx < m_nNbUnknowns - 1; nIdx++)
      {
        pdXCopy[nIdx] = pdX[nIdx] + dShift;

        shared_ptr<UnderlyingProcess> 
          pUnderlyingProcessNew( (*m_pTranslator)(&pdXCopy[0]) );

        shared_ptr<TheoreticalModel> 
          pModelNew( new TheoreticalModel(pUnderlyingProcessNew) );

        pdSpreads = 
          impliedSpreads.Compute(pModelNew, pNewCDS->GetSessionData(), pDates);

        double dComputedSpreadShift = pdSpreads[0];

        double dGradient = (dComputedSpreadShift - dComputedSpread) / dShift;

        pdGradients[nIdx] += 2.0 * dError * dGradient * dWeight;

        pdXCopy[nIdx] = pdX[nIdx];

      } // loop over gradients

      // Also need gradient w.r.t. cds recovery
      finance::ImpliedCDSSpreads impliedSpreads2(valuationDate,
                                                 firstDate, 
                                                 dayCount,
                                                 freq, 
                                                 dRecoveryRate + dShift);

      pdSpreads = 
        impliedSpreads2.Compute(pModel, pNewCDS->GetSessionData(), pDates);

      double dComputedSpreadShift = pdSpreads[0];

      double dGradient = (dComputedSpreadShift - dComputedSpread) / dShift;

      pdGradients[m_nNbUnknowns - 1] += 2.0 * dError * dGradient * dWeight;

    } // if gradient needed
  }
  else
  {
    // Objective based on price
    shared_ptr<finance::ModelOutput> output = pModel->Compute(*pNewCDS);

    double dComputedPrice = output->GetPrice();

    const double dMarketPrice = pNewCDS->GetMarketPrice();

    double dError = dComputedPrice - dMarketPrice;

    dObjectif += dWeight * dError * dError;

    if (bComputeGradient)
    {
      shared_ptr<NumOutput> 
        pNO( static_pointer_cast<NumOutput>(output->GetNumOutput()) );

      std::vector<double> pdSensitivities = pNO->GetSensitivities();

      for (size_t nIdx = 0; nIdx < m_nNbUnknowns - 1; nIdx++)
      {
        pdGradients[nIdx] += 2.0 * dError * pdSensitivities[nIdx]
                           * dWeight;
      }

      // Just use finite differences for cds recovery sensitivity
      double dShift = 1.e-7;

      shared_ptr<finance::CDS> pCDSShift
        (new finance::CDS( pdX[m_nNbUnknowns - 1] + dShift, 
                           pOriginalCDS->GetSpreadStream() ) );

      pCDSShift->SetSessionData( pOriginalCDS->GetSessionData() );

      output = pModel->Compute(*pCDSShift);

      double dShiftedPrice = output->GetPrice();

      double dSensitivity = (dShiftedPrice - dComputedPrice) / dShift;

      pdGradients[m_nNbUnknowns - 1] += 2.0 * dError * dSensitivity
                                      * dWeight;

    } // if computing gradients

  } // if objective on spreads or price

}

void CalibratorCDSRecovery::ObjectiveForEDS(const double *pdX, 
                                 double& dObjectif,
                                 double* pdGradients, 
                                 bool bComputeGradient,
                                 shared_ptr<finance::Derivative>& pDerivative,
                                 double dWeight,
                                 shared_ptr<TheoreticalModel>& pModel)
{
  // Convert back to original EDS
  shared_ptr<finance::EDS> 
    pEDS( dynamic_pointer_cast<finance::EDS>(pDerivative) );

  if ( m_bUseSpreads )
  {
    // Extract info needed to compute implied spreads, then compute 
    // the implied spread
    Date valuationDate = pEDS->GetSessionData()->GetValuationDate();
    Date firstDate     = pEDS->GetSpreadStream()->GetFirstPaymentDate();
    Date::DayCountConvention dayCount = 
      pEDS->GetSpreadStream()->GetDayCountConvention();
    finance::Frequency freq = pEDS->GetSpreadStream()->GetPaymentFrequency();

    double dRecoveryRate = pEDS->GetRecoveryRate();
    double dBarrier      = pEDS->GetBarrier();

    finance::ImpliedEDSSpreads impliedSpreads(valuationDate, firstDate, 
                          dayCount, freq, dBarrier, dRecoveryRate);

    std::vector<Date> pDates(1);
    pDates[0] = pEDS->GetMaturityDate();

    std::vector<double> pdSpreads = 
      impliedSpreads.Compute(pModel, pEDS->GetSessionData(), pDates);

    // Compute error based on market spread and implied spread
    double dMarketSpread = pEDS->GetSpreadStream()->GetAnnualPaymentAmount();
    
    double dComputedSpread = pdSpreads[0] ;

    double dError = dComputedSpread - dMarketSpread;

    dObjectif += dWeight * dError * dError;

    // Compute gradients, if requested
    if (bComputeGradient)
    {
      // Use finite differences for the spread gradients
      double dShift = 1.e-6;

      // Copy the param vector (since it is const)
      std::vector<double> pdXCopy(m_nNbUnknowns, 0.0);
      for (size_t nIdx = 0; nIdx < m_nNbUnknowns; nIdx++)
        pdXCopy[nIdx] = pdX[nIdx];

      for (size_t nIdx = 0; nIdx < m_nNbUnknowns - 1; nIdx++)
      {
        pdXCopy[nIdx] = pdX[nIdx] + dShift;

        shared_ptr<UnderlyingProcess> 
          pUnderlyingProcessNew( (*m_pTranslator)(&pdXCopy[0]) );

        shared_ptr<TheoreticalModel> 
          pModelNew( new TheoreticalModel(pUnderlyingProcessNew) );

        pdSpreads = impliedSpreads.Compute(pModelNew, pEDS->GetSessionData(),pDates);

        double dComputedSpreadShift = pdSpreads[0];

        double dGradient = (dComputedSpreadShift - dComputedSpread) / dShift;

        pdGradients[nIdx] += 2.0 * dError * dGradient * dWeight;

        pdXCopy[nIdx] = pdX[nIdx];

      } // loop over gradients

    } // if gradient needed
  }
  else
  {
    // Objective based on price
    shared_ptr<finance::ModelOutput> output = pModel->Compute(*pEDS);

    double dComputedPrice = output->GetPrice();

    const double dMarketPrice = pEDS->GetMarketPrice();

    double dError = dComputedPrice - dMarketPrice;

    dObjectif += dWeight * dError * dError;

    if (bComputeGradient)
    {
      shared_ptr<NumOutput> 
        pNO( static_pointer_cast<NumOutput>( output->GetNumOutput() ) );

      std::vector<double> pdSensitivities = pNO->GetSensitivities();

      for (size_t nIdx = 0; nIdx < m_nNbUnknowns - 1; nIdx++)
      {
        pdGradients[nIdx] += 2.0 * dError * pdSensitivities[nIdx]
                           * dWeight;
      }
    } // if computing gradients

  }  // if objective based on spread or price

}

void CalibratorCDSRecovery::RunCalibrator(bool bUserGradients)
{

  #ifdef USE_NAG_HG_RECOVERY
    // Construct the nag calibrator
    double dTolerance = 1.e-6;
    size_t nMaxIterations = 100;
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


shared_ptr<UnderlyingProcess>
CalibratorCDSRecovery::Calibrate(Translator* pTranslator, 
                             const finance::Derivatives& derivatives)
{
  // setup the translator
  m_pTranslator = pTranslator;

  // Get initial guess, bounds, etc
  InitializeArrays(m_pTranslator);

  // Make the lists needed by the objective function
  ConstructDerivativeLists(derivatives, true);
 
  // Compute scale factors for relative error calculations
  AddRelTolToWeights();

  // Create an underlying process to be used by all meshes.
  // m_pdX should be the initial guess as set in InitializeArrays.
  m_pUnderlyingProcessForMesh = (*m_pTranslator)(m_pdX);

  // Do the calibration. Objective function does compute gradients.
  RunCalibrator(true);

  // Return the calibrated model
  return GetLastCalibratedProcess();

} 


shared_ptr<UnderlyingProcess>
CalibratorCDSRecovery::GetLastCalibratedProcess() const
{
  return (*m_pTranslator)(m_pdFinalX);
}


double CalibratorCDSRecovery::GetCDSRecovery() const
{
  return m_pdFinalX[m_nNbUnknowns-1];
}


void CalibratorCDSRecovery::operator () 
(const std::vector<double>& pdParams, double& dF)
{
  // This function is called by the simulated annealing code. Just pass it
  // on to the objective function
  ComputeObjectif(&pdParams[0], dF, NULL, false);
}


void CalibratorCDSRecovery::InitializeArrays(const Translator* pTranslator)
{
  // Save the number of unknowns.  Needed by the objective function
  m_nNbUnknowns = pTranslator->GetNbParameters() + 1;

  // Start setting up the data needed by the optimization routine.
  m_pdX.resize(m_nNbUnknowns);
  m_pdFinalX.resize(m_nNbUnknowns);
  m_pdLowerBounds.resize(m_nNbUnknowns);
  m_pdUpperBounds.resize(m_nNbUnknowns);

  // Get the process bounds in temp vectors.  Copy below.
  std::vector<double> pdLowerBoundsTmp = pTranslator->GetLowerBounds();
  std::vector<double> pdUpperBoundsTmp = pTranslator->GetUpperBounds();

  // Initial guess for the calibration.
  // The order of this constuction must be 'undone' by the objective function.
  // This is made easier by the translator class
  pTranslator->GetParameters(&m_pdX[0]);

  // copy into final arrays.  
  for (size_t nIdx = 0; nIdx < m_nNbUnknowns - 1; nIdx++)
  {
    m_pdLowerBounds[nIdx] = pdLowerBoundsTmp[nIdx];
    m_pdUpperBounds[nIdx] = pdUpperBoundsTmp[nIdx];
  }

  // add on the extra cds recovery param
  m_pdLowerBounds[m_nNbUnknowns - 1] = 1.e-8;
  m_pdUpperBounds[m_nNbUnknowns - 1] = 0.999;
  m_pdX[m_nNbUnknowns - 1] = 0.2;
  
}

void CalibratorCDSRecovery::AddRelTolToWeights()
{

  // Pre-compute scale factors for relative tolerances
  if ( m_pDerivatives )
  {
    finance::Derivatives::Elements& elements( m_pDerivatives->GetAll() );
    finance::Derivatives::Elements::iterator iter;

    for (iter = elements.begin(); iter != elements.end(); ++iter)
    {
      // if cds or eds, check for calibration against spreads. Otherwise just 
      // use relative tolerance based on price.
      double dMarketValue;
      if (dynamic_cast<finance::CDS *>( iter->first.get() ) != 0)
      {

        shared_ptr<finance::CDS> 
          pCDS = dynamic_pointer_cast<finance::CDS>(iter->first);

        if ( m_bUseSpreads )
          dMarketValue = pCDS->GetSpreadStream()->GetAnnualPaymentAmount();
        else
          dMarketValue = iter->first->GetMarketPrice();
      }
      else if (dynamic_cast<finance::EDS *>( iter->first.get() ) != 0)
      {

        shared_ptr<finance::EDS> 
          pEDS = dynamic_pointer_cast<finance::EDS>(iter->first);

        if ( m_bUseSpreads )
          dMarketValue = pEDS->GetSpreadStream()->GetAnnualPaymentAmount();
        else
          dMarketValue = iter->first->GetMarketPrice();
      }
      else
      {
        // normal contract
        dMarketValue = iter->first->GetMarketPrice();
      }

      double dScale = fabs(dMarketValue);
      if ( dScale > 1.e-5 )
        dScale = 1.0 / dScale;
      else
        dScale = 1.0;    

      iter->second *= dScale * dScale;

    } // loop over the derivative list
  }

  // Now handle the options to be priced by forward equation
  if ( m_pForwardOption )
  {
    finance::ForwardOption::Elements& elements( m_pForwardOption->GetAll() );

    finance::ForwardOption::Elements::iterator iter;

    for (iter = elements.begin(); iter != elements.end(); ++iter)
    {
      double dScale = fabs( iter->m_pOption->GetMarketPrice() );
      if ( dScale > 1.e-5 )
        dScale = 1.0 / dScale;
      else
        dScale = 1.0;    

      iter->m_dWeight *= dScale * dScale;
    } 
  }
}


} // namespace hg

} // namespace ito33
