#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"

#include "ito33/finance/cds.h"
#include "ito33/finance/option.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/hrtimeonlycalibrator_usecds.h"
#include "ito33/ihg/voltimeonlycalibrator_useoptions.h"

#include "ihg/calibrator.h"

namespace ito33
{

namespace ihg
{

void Calibrator::Calibrate(OptionElements pOptions, CDSElements pCDS,
                           shared_ptr<Volatility> pVol, shared_ptr<HazardRate> pHR)
{

  // Check data
  CheckCDS(pCDS);
  CheckOptions(pOptions);

  if ( pOptions.size() == 0 && !pVol )
  {
    throw EXCEPTION_MSG
          (
            ITO33_BAD_DATA,
            TRANS("Option list is empty and no volatility was specified in calibrator.")
          );
  }

  if ( pCDS.size() == 0 && !pHR )
  {
    throw EXCEPTION_MSG
          (
            ITO33_BAD_DATA,
            TRANS("CDS list is empty and no hazard rate was specified in calibrator.")
          );
  }

  if ( pCDS.size() == 0 && pCDS.size() == 0 )
  {
    throw EXCEPTION_MSG
          (
            ITO33_BAD_DATA,
            TRANS("CDS and options lists are empty in calibrator.")
          );
  }


  // Determine what type of calibration needs to be done, and do it
  if ( pOptions.size() == 0 )
    CalibrateHR(pCDS, pVol, pHR);
  else if ( pCDS.size() == 0 )
    CalibrateVol(pOptions, pHR, pVol);
  else
    CalibrateVolAndHR(pOptions, pCDS, pVol, pHR);

}


void Calibrator::CalibrateHR(CDSElements pCDSs, shared_ptr<Volatility> pVol, 
                   shared_ptr<HazardRate> pHR)
{
  m_pTheoreticalModel->SetVolatility(pVol);
  m_pTheoreticalModel->SetHazardRate(pHR);  

  // Do the calibration
  HazardRateTimeOnlyCalibratorUseCDS hrCalibrator(pCDSs, pVol);
  if ( pHR != NULL)
    hrCalibrator.SetExternalHazardRate(pHR);

  m_pHRTimeOnly = hrCalibrator.Calibrate();
}


void Calibrator::CalibrateVol(OptionElements pOptions, shared_ptr<HazardRate> pHR,
                    shared_ptr<Volatility> pVol)
{
  m_pTheoreticalModel->SetVolatility(pVol);
  m_pTheoreticalModel->SetHazardRate(pHR);  

  // Do the calibration
  VolTimeOnlyCalibratorUseOptions volCalibrator(pOptions, pHR);
  if ( pVol != NULL)
    volCalibrator.SetExternalVolatility(pVol);

  m_pVolTimeOnly = volCalibrator.Calibrate();
}


void Calibrator::CalibrateVolAndHR(OptionElements pOptions, CDSElements pCDSs, 
                 shared_ptr<Volatility> pVol, shared_ptr<HazardRate> pHR)
{

  // Setup a parametrization for the calibration
  shared_ptr<Volatility> pVolTmp = m_pTheoreticalModel->GetVolatility();
  shared_ptr<HazardRate> pHRTmp = m_pTheoreticalModel->GetHazardRate();

  // Setup the calibrarors
  VolTimeOnlyCalibratorUseOptions 
    volCalibrator(pOptions, pHRTmp);
  if ( pVol != NULL)
    volCalibrator.SetExternalVolatility(pVol);

  HazardRateTimeOnlyCalibratorUseCDS hrCalibrator(pCDSs, pVolTmp);  
  if ( pHR != NULL)
    hrCalibrator.SetExternalHazardRate(pHR);

  // Now do the calibration
  size_t nIteration = 0;
  double dError = 1.e100;
  while (nIteration < 40 )
  {
    nIteration++;
  
    // Fit the options
    m_pVolTimeOnly = volCalibrator.Calibrate();

    // set the fitted vol
    if ( pVol != NULL)
    {
      shared_ptr<Volatility> 
        pTmpVol( new VolatilityCombo(m_pVolTimeOnly, pVol) );

      m_pTheoreticalModel->SetVolatility(pTmpVol);
    }
    else
      m_pTheoreticalModel->SetVolatility(m_pVolTimeOnly);

    // Now fit the hazard rate to the cds
    m_pHRTimeOnly = hrCalibrator.Calibrate();

    // set the fitted hr
    if ( pHR != NULL)
    {
      shared_ptr<HazardRate> 
        pTmpHR( new HazardRateCombo(m_pHRTimeOnly, pHR) );

      m_pTheoreticalModel->SetHazardRate(pTmpHR);
    }
    else
      m_pTheoreticalModel->SetHazardRate(m_pHRTimeOnly);

    // now compute the error
    size_t nIdx;
    double dCurrentError = 0.0;

    for (nIdx = 0; nIdx < pOptions.size(); nIdx++)
    {
      shared_ptr<ihg::ModelOutput>  
        output = m_pTheoreticalModel->Compute(*pOptions[nIdx]);

      double dAbsError = fabs( output->GetPrice() 
                               - (*pOptions[nIdx]).GetPrice() );

      double dScale = (*pOptions[nIdx]).GetPrice();
      if (dScale < 1.0)
        dScale = 1.0;

      dCurrentError += (dAbsError / dScale) * (dAbsError / dScale);
    }
      
    for (nIdx = 0; nIdx < pCDSs.size(); nIdx++)
    {
      shared_ptr<ihg::ModelOutput> 
        output = m_pTheoreticalModel->Compute(*pCDSs[nIdx]);

      double dAbsError = fabs( output->GetPrice()
                              - (*pCDSs[nIdx]).GetPrice() );

      double dScale = (*pCDSs[nIdx]).GetPrice();
      if (dScale < 1.0)
        dScale = 1.0;

      dCurrentError += (dAbsError / dScale) * (dAbsError / dScale);
    }

    if (dCurrentError < 1.e-8)
      break;

    if (dCurrentError > dError && nIteration > 10)
      break;

    if ( fabs(dCurrentError - dError) < 1.e-10 && nIteration > 4)
      break;

    dError = dCurrentError;
    //std::cout << "Error = " << dError << std::endl;
  }

}



void Calibrator::CheckCDS(CDSElements& pCDS)
{

  if (pCDS.size() == 0)
    return;

  // Make sure the maturity dates are increasing, and the session datas
  // are the same
  CDSElements::iterator iterCDS = pCDS.begin();

  shared_ptr<finance::SessionData> 
    sessionData = (*iterCDS)->GetSessionData();

  Date lastMaturity = (*iterCDS)->GetMaturityDate();
  ++iterCDS;

  for (; iterCDS != pCDS.end(); ++iterCDS)
  {
    if ( (*iterCDS)->GetMaturityDate() < lastMaturity )
    {
      throw EXCEPTION_MSG
            (
              ITO33_BAD_DATA,
              TRANS("Maturities do not increase in list of CDS contracts in calibrator")
            );
    }

    if ( (*iterCDS)->GetSessionData() != sessionData )
    {
      throw EXCEPTION_MSG
            (
              ITO33_BAD_DATA,
              TRANS("Session datas are not the same in list of CDS contracts in calibrator")
            );
    }

    lastMaturity = (*iterCDS)->GetMaturityDate();

  }

}

void Calibrator::CheckOptions(OptionElements& pOptions)
{
  if (pOptions.size() == 0)
    return;

  // Make sure the maturity dates are increasing, and the session datas 
  // are the same
  OptionElements::iterator iterOption = pOptions.begin();

  shared_ptr<finance::SessionData> 
    sessionData = (*iterOption)->GetSessionData();

  Date lastMaturity = (*iterOption)->GetMaturityDate();
  ++iterOption;

  for (; iterOption != pOptions.end(); ++iterOption)
  {
    if ( (*iterOption)->GetMaturityDate() < lastMaturity )
    {
      throw EXCEPTION_MSG
            (
              ITO33_BAD_DATA,
              TRANS("Maturities do not increase in list of options contracts in calibrator")
            );
    }
  
    if ( (*iterOption)->GetSessionData() != sessionData )
    {
      throw EXCEPTION_MSG
            (
              ITO33_BAD_DATA,
              TRANS("Sessions are not the same in list of options contracts in calibrator")
            );
    }

    lastMaturity = (*iterOption)->GetMaturityDate();

  }

}


} // namespace ihg

} // namespace ito33
