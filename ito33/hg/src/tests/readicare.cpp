/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/tests/readicare.cpp
// Purpose:     Read ICARE calibration output file
// Created:     2005/06/20
// RCS-ID:      $Id: readicare.cpp,v 1.4 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <fstream>
#include "ito33/afterstd.h"

#include "ito33/vector.h"
#include "ito33/array.h"
#include "ito33/constants.h"

#include "ito33/finance/dividends.h"
#include "ito33/finance/yieldcurve_annuallycompounded.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/option.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/cds.h"
#include "ito33/finance/derivatives.h"

#include "ito33/hg/underlyingprocess.h"
#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/parametrization.h"

#include "hg/tests/readicare.h"

namespace ito33
{

namespace hg
{


void 
ReadICARE
(std::ifstream& sIn, shared_ptr<Parametrization>& pParamtrization,
 shared_ptr<finance::Derivatives>& pDerivatives)
{
  int iOutPutData;

  // Read and discard the output option, let the caller decide what to do
  sIn >> iOutPutData; 

  // Read the yield curve
  size_t nNbLegs;

  sIn >> nNbLegs;

  Array<double> pdYears(nNbLegs);
  Array<double> pdRates(nNbLegs);

  for (size_t nIdxLeg = 0; nIdxLeg < nNbLegs; nIdxLeg++)
    sIn >> pdYears[nIdxLeg] >> pdRates[nIdxLeg];
  
  // Read the borrow curve
  size_t nNbBorrowLegs;

  sIn >> nNbBorrowLegs;
  
  Array<double> pdBorrowYears(nNbBorrowLegs);
  Array<double> pdBorrowRates(nNbBorrowLegs);  

  for (size_t nIdxLeg = 0; nIdxLeg < nNbBorrowLegs; nIdxLeg++)
    sIn >>  pdBorrowYears[nIdxLeg] >> pdBorrowRates[nIdxLeg];  
   
  // Read the dividends
  size_t nNbDividends;
  int iYieldCash;
  int iDate;
  double dValue;
  double dPseudoYield;

  shared_ptr<finance::Dividends> pDividends(new finance::Dividends);

  sIn >> nNbDividends; 

  for (size_t nIdxDiv = 0; nIdxDiv < nNbDividends; nIdxDiv++)
  {
    sIn >> iYieldCash >> iDate >> dValue;
  
    if ( iYieldCash )
      pDividends->AddYield(iDate, dValue);
    else
    {
      sIn >> dPseudoYield;
        
      if ( dPseudoYield == 1)
        pDividends->Add(finance::Dividend::Cash, iDate, dValue);
      else
        pDividends->Add(finance::Dividend::PseudoCash, iDate, dValue, dPseudoYield);
    }
  }

  // Read the valuation date
  int iValuationDate;

  sIn >> iValuationDate;

  Date valuationDate(iValuationDate);

  // Read the spot share price
  double dSpot;

  sIn >> dSpot;
  
  // Create the borrow curve
  shared_ptr<finance::YieldCurveAnnuallyCompounded>
    pBorrowCurve( new finance::YieldCurveAnnuallyCompounded(valuationDate) );

  for (size_t nIdxLeg = 0; nIdxLeg < nNbBorrowLegs; nIdxLeg++)
    pBorrowCurve->AddLeg(int(pdBorrowYears[nIdxLeg] * INVERSEONEDAY + 0.04), 
                         pdBorrowRates[nIdxLeg]);
  //Create currency
  shared_ptr<finance::Numeraire> pCurrency(new finance::Numeraire("EUR") );

  // Create the equity
  shared_ptr<finance::Equity>
    pEquity( new finance::Equity(dSpot, pCurrency) );

  pEquity->SetBorrowCurve(pBorrowCurve);

  pEquity->SetDividends(pDividends);

  // The yield curve
  shared_ptr<finance::YieldCurveAnnuallyCompounded>
    pYieldCurve( new finance::YieldCurveAnnuallyCompounded(valuationDate) );

  for (size_t nIdxLeg = 0; nIdxLeg < nNbLegs; nIdxLeg++)
    pYieldCurve->AddLeg(int(pdYears[nIdxLeg] * INVERSEONEDAY + 0.04),
                        pdRates[nIdxLeg]);

  shared_ptr<finance::RateData> pRateData( new finance::RateData);
  pRateData->SetYieldCurve(pCurrency, pYieldCurve);

  // Create the session
  shared_ptr<finance::SessionData>
    pSessionData( new finance::SessionData(pRateData, pEquity, valuationDate) );

  // Read the Vanilla Smile
  std::vector< shared_ptr<finance::Option> > options;

  size_t nNbOptions;
  int iMaturityDate;
  double dStrike;
  double dBSVol;

  // Read the number of vanilla options
  sIn >> nNbOptions;
  
  // A temporary model to price the option
  std::vector<double> pdVolsTmp(1, 0.2), pdDefaultIntensitiesTmp(1);

  shared_ptr<UnderlyingProcess> 
    processTmp( new UnderlyingProcess(1, pdVolsTmp, pdDefaultIntensitiesTmp) );
  
  TheoreticalModel modelTmp(processTmp);

  for (size_t nIdxOp = 0; nIdxOp < nNbOptions; nIdxOp++)
  {
    sIn >> iMaturityDate >> dStrike >> dBSVol;

    shared_ptr<finance::Option>
      pOption( new finance::Option(dStrike, iMaturityDate, 
                                   finance::Option_Call,
                                   finance::ExerciseType_European) );

    pOption->SetSessionData(pSessionData);

    // use the temporary model to convert bs vol to market price
    processTmp->SetVolatilities( std::vector<double>(1, dBSVol) );

    pOption->SetMarketPrice( modelTmp.Compute(*pOption)->GetPrice() );

    options.push_back(pOption);
  }

  // Read the CDS curve
  std::vector< shared_ptr<finance::CDS> > cdss;

  size_t nNbCDS;

  int iFrequency;
  double dMaturity, dSpread, dRecoveryRate;

  // Read the number of CDS
  sIn >> nNbCDS;

  for (size_t nIdxCDS = 0; nIdxCDS < nNbCDS; nIdxCDS++)
  {
    sIn >> iFrequency >> dMaturity >> dSpread >> dRecoveryRate;
    
    // Note that there is no way to get the exact transformation
    // in some cases because of different conventions in hg and ICARE
    // but it shouldn't matter much
    double dAnnualAmount = dSpread * iFrequency;
    Date maturityDate = valuationDate;
    maturityDate.AddYears( int(dMaturity + 1.e-6) );
    
    Date firstCoupon = valuationDate;
    firstCoupon.AddMonths(12 / iFrequency);
    Date::DayCountConvention dcc = Date::DayCountConvention_Act365;

    shared_ptr<finance::CashFlowStreamUniform>
      pSpreads( new finance::CashFlowStreamUniform
                    ( valuationDate, firstCoupon, maturityDate, dAnnualAmount, 
                      dcc, static_cast<finance::Frequency>(iFrequency) ) );

    shared_ptr<finance::CDS>
      pCDS( new finance::CDS(dRecoveryRate, pSpreads) );

    pCDS->SetSessionData(pSessionData);

    // in ICARE, the cds price is just zero by default
    pCDS->SetMarketPrice(0.0);

    cdss.push_back(pCDS);
  }

  // Read the parametrization and the initial guess
  size_t nNbRegimes;
  size_t nNbJumps;

  std::vector<double> pdVols, pdDefaultIntensities;

  double dVol, dDefaultIntensity, dAmplitude, dIntensity;
  
  sIn >> nNbRegimes;

  Jumps jumpMatrix[NMAXREGIMES][NMAXREGIMES];

  for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
  {
    for (size_t nIdxR2 = 0; nIdxR2 < nNbRegimes; nIdxR2++)
    {
      Jumps& jumps = jumpMatrix[nIdxR1][nIdxR2];

      sIn >> nNbJumps;
      for (size_t nIdxJump = 0; nIdxJump < nNbJumps; nIdxJump++)
      {
        sIn >> dAmplitude >> dIntensity;

        jumps.Add(dIntensity, dAmplitude);
      }
    }  

    sIn >> dVol;
    pdVols.push_back(dVol);
  }

  int iHasDefault;

  // discard the flag indicating if the model has a default regime. It was used
  // by old files
  sIn >> iHasDefault;

  for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
  {
    sIn >> dDefaultIntensity;

    pdDefaultIntensities.push_back(dDefaultIntensity);
  }
  
  shared_ptr<UnderlyingProcess>
    pUP( new UnderlyingProcess(nNbRegimes, pdVols, pdDefaultIntensities) );

  for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
    for (size_t nIdxR2 = 0; nIdxR2 < nNbRegimes; nIdxR2++)
      pUP->SetJumps(nIdxR1, nIdxR2, jumpMatrix[nIdxR1][nIdxR2]);

  pParamtrization = shared_ptr<Parametrization>( new Parametrization(pUP) );

  // Read the parameters for the meshes and discard them
  size_t nNbT, nZoomT, nNbZoomT;
  size_t nNbS, nZoomS, nNbZoomS;
  double dRhoT, dRhoS;
  double dHorizon;

  sIn >> nNbT >> nZoomT >> nNbZoomT >> dRhoT;
  sIn >> nNbS >> nZoomS >> nNbZoomS >> dRhoS;

  sIn >> dHorizon; 
  
  // use or not the vega weighting for the vanilla smile,
  // not yet implemented in HG
  int iOptionWeighting;

  sIn  >> iOptionWeighting; 

  // Read and discard another user option on pricing or calibration
  // let the caller decide it
  int iUserTest;

  sIn >> iUserTest;
  
  // Read and discard the option on using the gradient in the calibration
  // We always use gradient 
  int iUseGradient;

  sIn >> iUseGradient;
  
  // Read the tolerance for the nag calibrator, not yet exported in HG
  double dToleranceForNag;  
  
  sIn >> dToleranceForNag;

  // Read the relative weight given for the CDS Smile   
  double dWeightCDS;

  sIn >> dWeightCDS;
  
  pDerivatives = shared_ptr<finance::Derivatives>(new finance::Derivatives);

  for (size_t nIdxOp = 0; nIdxOp < options.size(); nIdxOp++)
    pDerivatives->Add(options[nIdxOp]);

  for (size_t nIdxCDS = 0; nIdxCDS < cdss.size(); nIdxCDS++)
    pDerivatives->AddWithWeight(cdss[nIdxCDS], dWeightCDS);

  sIn.close();
}


} // namespace hg

} // namespace ito33
