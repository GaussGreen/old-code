/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/src/test_greeks.cpp
// Purpose:     Compare computed greeks in output wih financial level
//                greeks
// Created:     March 31, 2005
// RCS-ID:      $Id: test_greeks.cpp,v 1.14 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <cmath>
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/array.h"
#include "ito33/sharedptr.h"
#include "ito33/autoptr.h"
#include "ito33/constants.h"

#include "ito33/constants.h"
#include "ito33/numeric/predicatedouble.h"

#include "ito33/finance/derivative.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/option.h"
#include "ito33/finance/eds.h"
#include "ito33/finance/yieldcurve.h"
#include "ito33/finance/issuer.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/spotfxrates.h"
#include "ito33/finance/qualitycontrol.h"
#include "ito33/finance/computationalflags.h"

#include "ito33/finance/bondlike/reset.h"
#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/finance/bondlike/bond.h"
#include "ito33/finance/bondlike/pepslike.h"
#include "ito33/finance/bondlike/percslike.h"
#include "ito33/finance/bondlike/attachedwarrantconvertiblebond.h"

#include "ito33/finance/derivativevisitors/derivative_visitor_goodtype.h"

#include "ito33/ihg/theoreticalmodel.h"

#include "ito33/ihg/volatility.h"
#include "ito33/ihg/hazardrate.h"

#include "ito33/numeric/interpolation.h"

#include "ihg/xml/pricingreaderrecursive.h"

#include "ihg/tests/test_greeks.h"

using namespace ito33::finance;

namespace ito33
{

namespace ihg
{

TestGreeks::TestGreeks()
           : m_OutputFile("./greek_report.xml"),
             m_RootTag("root",m_OutputFile)
{
  m_pQualityControl = shared_ptr<finance::QualityControl>(new QualityControl);
}

double TestGreeks::GetSpot(const std::string &sFileName)
{
  shared_ptr<finance::Derivative> pDerivative;
  finance::ModelOutput output;
  shared_ptr<SessionData> pSessionData;
  shared_ptr<TheoreticalModel> pModel (new TheoreticalModel());

  ReadXMLFile( sFileName, pSessionData, pDerivative, pModel);

  return pSessionData->GetSpotSharePrice();

} //GetSpot

void TestGreeks::Run(const std::list<std::string> &fileList)
{
  std::list<std::string>::const_iterator iter;

  /*
  /// To refine the mesh
  shared_ptr<finance::QualityControl> pQualityControl(new QualityControl);  
  pQualityControl->SetComputationQuality(ComputationQuality_High);  
  SetQualityControl(pQualityControl);
  */

  for ( iter = fileList.begin(); iter != fileList.end(); ++iter )
  {
    std::cout << "Testing Greeks for: " << *iter << std::endl;

    //Get the spot so that the greeks are tested
    //for different values
    //double dSpot = GetSpot(*iter);

    //Step 1: Compare Vega
    CheckGreekVega(*iter);

    //Step 2: Compare Rho
    CheckGreekRho(*iter);

    //Step 3: Compare FXDelta
    CheckGreekFXDelta(*iter);

    //Step 3: Checks the convergence
    CheckGreekFXDeltaConv(*iter);
  }
   
}


void TestGreeks::CheckGreekVega(const std::string &sFileName)
{    
  const double dShiftVol = SHIFT;
  const double dInverseShiftVol = 1. / dShiftVol;

  shared_ptr<finance::Derivative> pDerivative;
  finance::ModelOutput output;
  shared_ptr<SessionData> pSessionData;
  shared_ptr<TheoreticalModel> pModel (new TheoreticalModel);

  ReadXMLFile( sFileName, pSessionData, pDerivative, pModel);
  
  // Set analysis date at the pricing date
  /*pModel->GetComputationalFlags()
        ->SetAnalysisDate(pSessionData->GetValuationDate());*/
  
  shared_ptr<ComputationalFlags> pFlags (new ComputationalFlags());
  pFlags->SetComputeSurface(true);
  pFlags->SetComputeVega(true);

  pDerivative->SetComputationalFlags(pFlags);

  PriceContract(*pDerivative, pModel, output);

  double dPrice = output.GetPrice();
  double dVega = 0.;
  
  Values
    pdOriginalSpots,
    pdOriginalPrices,
    pdOriginalGreeks;
    
  if( output.HasSpotAtAnalysisDate() )
    pdOriginalSpots = output.GetSpotsAtAnalysisDate();

  if( output.HasPriceAtAnalysisDate() )
    pdOriginalPrices = output.GetPricesAtAnalysisDate();

  if( output.HasVegaAtAnalysisDate() )
    pdOriginalGreeks = output.GetVegasAtAnalysisDate();

  if( output.HasVega() )
    dVega  = output.GetVega();
  else
  {
    std::cout << "Vega testing does not apply. " << std::endl;
    return;
  }

  ito33::XML::Tag test("test", m_RootTag);
  test.Element("file")(sFileName);
  test.Element("greek")("Vega");

  //shift the volatility
  shared_ptr<Volatility> pVol = pModel->GetVolatility()->Perturb(dShiftVol);
  pModel->SetVolatility(pVol);  

  // Deactivate the computation of the vega
  pDerivative->GetComputationalFlags()->SetComputeVega(false);

  PriceContract(*pDerivative, pModel, output);

  double dVegaShifted   = ( output.GetPrice() - dPrice )* dInverseShiftVol;
  double dError = fabs(dVegaShifted-dVega)/std::max(fabs(dVega),1.);

  //store the results
  test.Element("vega_computed_code")(dVega);
  test.Element("vega_manual_shifting")(dVegaShifted);
  test.Element("error")(dError);

  CheckResult(dError, dVega, dVegaShifted, 
    std::string("Vega"), sFileName);  
  
  if( output.HasVegaAtAnalysisDate() )
  {
    Values 
      pdShiftedSpots,
      pdShiftedPrices;   
  
    if( output.HasSpotAtAnalysisDate() )
      pdShiftedSpots = output.GetSpotsAtAnalysisDate();

    if( output.HasPriceAtAnalysisDate() )
      pdShiftedPrices = output.GetPricesAtAnalysisDate();
    
    CheckGreekSurfaceAtAnalysisDateResults( pdOriginalGreeks, pdOriginalPrices,
      pdOriginalSpots, pdShiftedPrices, pdShiftedSpots, dInverseShiftVol, 
      std::string("Vega"), sFileName );
  }

} //CheckGreekVega

void TestGreeks::CheckResult(double dError, double dOriginalGreek, 
        double dGreek, std::string sTestName, std::string sFileName)
{
  
  if( dError  > 5./100. )
  {
    std::cout << std::endl;
    std::cout.precision(12);
    std::cout << "Failed test for greek  : " << sTestName << std::endl;
    std::cout << "Filename               : " << sFileName << std::endl;
    std::cout << "Original Greek         : " << dOriginalGreek << std::endl;
    std::cout << "Computed shifted Greek : " << dGreek << std::endl;
    std::cout << "Error         : " << dError << std::endl;
  }

} // CheckResult

void TestGreeks::CheckGreekSurfaceAtAnalysisDateResults
                 ( const Values& pdOriginalGreeks, 
                   const Values& pdOriginalPrices, 
                   const Values& pdOriginalSpots,  
                   const Values& pdShiftedPrices, 
                   const Values& pdShiftedSpots,
                   double dInverseShift,
                   std::string sTestName, std::string sFileName )
{  
  double
    dShiftedGreek,
    dError,
    dMaxError = 0.;
  
  size_t
    nIdx,
    nNbOriginalSpots = pdOriginalSpots.size(), 
    nNbShiftedSpots = pdShiftedSpots.size();
 
  ASSERT_MSG(   (pdOriginalGreeks.size() == pdOriginalPrices.size() ) 
             && (pdOriginalGreeks.size() == pdOriginalSpots.size() )
             && (pdOriginalSpots.size() > 0)
             && (pdShiftedPrices.size() == pdShiftedSpots.size() )
             && (pdShiftedSpots.size() > 0), 
             "Invalid data : can't check greek surface at analysis date "
             "results"
            );

  Array<double> 
    pdInitialSpots(nNbOriginalSpots),
    pdPerturbedPricesOnInitialSpots(nNbOriginalSpots),
    pdPerturbedSpots(nNbShiftedSpots), 
    pdPerturbedPrices(nNbShiftedSpots);

  for(nIdx = 0; nIdx < nNbOriginalSpots; ++nIdx)
    pdInitialSpots[nIdx] = pdOriginalSpots[nIdx];

  for(nIdx = 0; nIdx < nNbShiftedSpots; ++nIdx)
  {
    pdPerturbedSpots[nIdx] = pdShiftedSpots[nIdx];
    pdPerturbedPrices[nIdx] = pdShiftedPrices[nIdx];
  }

  Interpolate(pdPerturbedSpots.Get(), pdPerturbedPrices.Get(), nNbShiftedSpots,
              pdInitialSpots.Get(), pdPerturbedPricesOnInitialSpots.Get(), 
              nNbOriginalSpots, numeric::ExtrapolationMode_Linear, 
              numeric::ExtrapolationMode_Linear);

  for(nIdx = 0; nIdx < nNbOriginalSpots; ++nIdx)
  {
    dShiftedGreek = 
      (pdPerturbedPricesOnInitialSpots[nIdx] - pdOriginalPrices[nIdx]) 
      * dInverseShift;
    
    dError = fabs(dShiftedGreek - pdOriginalGreeks[nIdx]) 
                   / std::max(fabs(pdOriginalGreeks[nIdx]),1.);
    
    if( dMaxError < dError )
      dMaxError = dError;
  }

  if( dMaxError  > 5./100. )
  {
    std::cout << std::endl;
    std::cout.precision(12);
    std::cout 
      << "Failed test for the surface at analysis date for the greek : " 
      << sTestName << std::endl;
    std::cout << "Filename               : " << sFileName << std::endl;
    std::cout << "Maximum error : " << dMaxError << std::endl;
  }

} // CheckGreekSurfaceAtAnalysisDateResults

void TestGreeks::CheckGreekRho(const std::string &sFileName)
{ 
  const double dShiftYC = SHIFT;
  const double dInverseShiftYC = 1 / dShiftYC;

  shared_ptr<finance::Derivative> pDerivative;
  finance::ModelOutput output;
  shared_ptr<SessionData> pSessionData;
  shared_ptr<TheoreticalModel> pModel (new TheoreticalModel());
  
  ReadXMLFile(sFileName, pSessionData, pDerivative, pModel);
  
  shared_ptr<ComputationalFlags> pFlags (new ComputationalFlags());
  pFlags->SetComputeSurface(true);
  pFlags->SetComputeRho(true);
  // Set analysis date at the pricing date
  pFlags->SetAnalysisDate(pSessionData->GetValuationDate());

  pDerivative->SetComputationalFlags(pFlags);

  finance::DerivativeVisitorGoodType visitor;
  pDerivative->Visit(visitor);

  PriceContract(*pDerivative, pModel, output );

  double dPrice = output.GetPrice();
  double 
    dRho = 0.,
    dUnderlyingRho = 0.;
   
  Values
    pdOriginalSpots, 
    pdOriginalPrices,
    pdOriginalRhos,
    pdOriginalUnderlyingRhos;
  
  if( output.HasPriceAtAnalysisDate() )
    pdOriginalPrices = output.GetPricesAtAnalysisDate();

  if( output.HasSpotAtAnalysisDate() )
    pdOriginalSpots = output.GetSpotsAtAnalysisDate();

  if( output.HasRhoAtAnalysisDate() )
    pdOriginalRhos = output.GetRhosAtAnalysisDate();

  if( output.HasUnderlyingRhoAtAnalysisDate() )
    pdOriginalUnderlyingRhos = output.GetUnderlyingRhosAtAnalysisDate();

  if( output.HasRho() )
  {
    dRho = output.GetRho(); 
    dUnderlyingRho = output.GetUnderlyingRho();
  }
  else
  {
    std::cout << "Can not computed Rho for this instrument." << std::endl;
    return;
  }

  ito33::XML::Tag test("test",m_RootTag);
  test.Element("file")(sFileName);  
  test.Element("greek")("Rho");
  
  if ( visitor.GetConvertibleBond() )
  {
    if ( visitor.GetConvertibleBond()->IsExchangeable() )
    {    
      std::cout << std::endl;
      std::cout << "Instrument is exchangeable. Rho is not supported." 
                <<  std::endl << std::endl;  
      return;
    }
  }

  if( visitor.GetGeneralizedPEPSLike() )
  {
    if ( visitor.GetGeneralizedPEPSLike()->IsExchangeable() )
    {          
      std::cout << std::endl 
                << "Instrument is exchangeable. Rho is not supported." 
                <<  std::endl << std::endl;      
      return;
    }  
  }

  if( visitor.GetGeneralizedPEPSLike() )
  {    
    if ( visitor.GetGeneralizedPEPSLike()->IsExchangeable() )
    {        
      
      std::cout << std::endl 
        << "Instrument is exchangeable. Rho is not supported." 
        <<  std::endl << std::endl;      
      return;
    }    
  }

  if( visitor.GetPEPSLike() )
  {
    if ( visitor.GetPEPSLike()->IsExchangeable() )
    {           
      std::cout << std::endl 
                << "Instrument is exchangeable. Rho is not supported." 
                <<  std::endl << std::endl;      
      return;
    }
  }

  if ( visitor.GetPERCSLike() ) 
  { 
    if ( visitor.GetPERCSLike()->IsExchangeable() )
    {    
      std::cout << std::endl 
                << "Instrument is exchangeable. Rho is not supported." 
                <<  std::endl << std::endl;  
      return;
    }
  }
      
  if ( visitor.GetReset() ) 
  {
    if ( visitor.GetReset()->IsExchangeable() )
    {
      std::cout << std::endl;
      std::cout << "Instrument is exchangeable. Rho is not supported." 
                <<  std::endl << std::endl;  
      return;
    }
  }

  if ( visitor.GetAttachedWarrantCB() ) 
  {
    if ( visitor.GetAttachedWarrantCB()->IsExchangeable() )
    {
      std::cout << "Instrument is exchangeable. Rho is not supported." 
                <<  std::endl;  
      return;
    }
  }
  
  if ( visitor.GetCBOption() )
  {
    if ( visitor.GetCBOption()->GetConvertibleBond()->IsExchangeable() )
    {    
      std::cout << std::endl;
      std::cout << "Instrument is exchangeable. Rho is not supported." 
                <<  std::endl << std::endl;  
      return;
    }
  }
    
  // deactivate the computation of the rho
  pDerivative->GetComputationalFlags()->SetComputeRho(false);
  
  double 
    dRhoShifted,
    dError;

  // "Underlying" rho **********************************************

  // Save the origin "underlying" YC
  shared_ptr<finance::YieldCurve> 
    pYC = pSessionData->GetYieldCurve();
  
  // Perturb the yield curve.
  shared_ptr<finance::YieldCurve> pYCNew( pYC->Perturb(dShiftYC) );

  pSessionData->SetYieldCurve(pYCNew);

  PriceContract(*pDerivative, pModel, output);

  dRhoShifted = ( output.GetPrice() - dPrice ) * dInverseShiftYC;
  dError = fabs(dRhoShifted - dUnderlyingRho)
         / std::max(fabs(dUnderlyingRho),1.);
  
  test.Element("underlying_rho_computed_code")(dUnderlyingRho);
  test.Element("underlying_rho_manual_shifting")(dRhoShifted);
  test.Element("error")(dError);

  CheckResult(dError, dUnderlyingRho, dRhoShifted, 
    std::string("UnderlyingRho"), sFileName);
  
  if( output.HasUnderlyingRhoAtAnalysisDate() )
  {    
    Values 
      pdShiftedSpots,
      pdShiftedPrices;   
  
    if( output.HasSpotAtAnalysisDate() )
      pdShiftedSpots = output.GetSpotsAtAnalysisDate();

    if( output.HasPriceAtAnalysisDate() )
      pdShiftedPrices = output.GetPricesAtAnalysisDate();

    CheckGreekSurfaceAtAnalysisDateResults( pdOriginalUnderlyingRhos, 
      pdOriginalPrices, pdOriginalSpots, pdShiftedPrices, pdShiftedSpots, 
      dInverseShiftYC, std::string("UnderlyingRho"), sFileName );
  }
  
  // Restore the origin "underlying" YC
  pSessionData->SetYieldCurve(pYC);
  
  // "Derivative" rho *********************************************

  if ( pDerivative->IsCrossCurrency() )
  {
    if( !pDerivative->GetNumeraire() )
    {
      std::cout << std::endl;
      std::cout << "Currency for a cross-currency instrument not defined!" 
                <<  std::endl << std::endl;  
      return;
    }

    // Save the origin "derivative" YC
    shared_ptr<Numeraire> pDerivCurrency = pDerivative->GetNumeraire();

    pYC = pSessionData->GetRateData()->GetYieldCurve( pDerivCurrency );
    
    pYCNew = pYC->Perturb(dShiftYC);

    pSessionData->GetRateData()->SetYieldCurve(pDerivCurrency, pYCNew);

    PriceContract(*pDerivative, pModel, output);
    
    // Restore the origin "derivative" YC
    pSessionData->GetRateData()->SetYieldCurve(pDerivCurrency, pYC);
  }  

  dRhoShifted = ( output.GetPrice() - dPrice ) * dInverseShiftYC;
  
  dError = fabs(dRhoShifted - dRho) / std::max(fabs(dRho),1.);
  
  test.Element("derivative_rho_computed_code")(dRho);
  test.Element("derivative_rho_manual_shifting")(dRhoShifted);
  test.Element("error")(dError);

  CheckResult(dError, dRho, dRhoShifted, std::string("DerivativeRho"), 
    sFileName);
  
  if( output.HasRhoAtAnalysisDate() )
  {    
    Values 
      pdShiftedSpots,
      pdShiftedPrices;   
  
    if( output.HasSpotAtAnalysisDate() )
      pdShiftedSpots = output.GetSpotsAtAnalysisDate();

    if( output.HasPriceAtAnalysisDate() )
      pdShiftedPrices = output.GetPricesAtAnalysisDate();
    CheckGreekSurfaceAtAnalysisDateResults( pdOriginalRhos, 
      pdOriginalPrices, pdOriginalSpots, pdShiftedPrices, pdShiftedSpots, 
      dInverseShiftYC, std::string("DerivativeRho"), sFileName );
  }
  
} //CheckGreekRho

void TestGreeks::CheckGreekFXDelta(const std::string &sFileName)
{
  shared_ptr<finance::Derivative> pDerivative;
  finance::ModelOutput output;
  shared_ptr<SessionData> pSessionData;
  shared_ptr<TheoreticalModel> pModel (new TheoreticalModel);
  
  ReadXMLFile( sFileName, pSessionData, pDerivative, pModel);
  pModel->SetDebugOutputFile("fxdelta.xml");
  
  if ( pDerivative->IsCrossCurrency() )
  {
    PriceContract(*pDerivative, pModel, output);

    double dPrice = output.GetPrice();
    double dFXDelta;

    if( output.HasFXDelta() )
      dFXDelta  = output.GetFXDelta();
    else
    {
      std::cout << "FX delta testing does not apply. " << std::endl;
      return;
    }

    std::cout << "***FXDelta***" << std::endl;
    // mannual fx delta
    ito33::XML::Tag test("test", m_RootTag);
    test.Element("file")(sFileName);
    test.Element("greek")("FXDelta");
  
    const double dShiftFXRate = SHIFT;
    const double dInverseShiftFXRate = 1. / dShiftFXRate;

    //shift the FX rate
    pDerivative->PerturbFXRate( dShiftFXRate );

    PriceContract(*pDerivative, pModel, output);

    //Restore FX rate
    pDerivative->PerturbFXRate( -dShiftFXRate );

    double 
      dFXDeltaShifted   = ( output.GetPrice() - dPrice ) * dInverseShiftFXRate;
    double 
      dError = fabs(dFXDeltaShifted-dFXDelta)/std::max(fabs(dFXDelta),1.);

    //store the results
    test.Element("fxdelta_computed_code")(dFXDelta);
    test.Element("fxdelta_manual_shifting")(dFXDeltaShifted);
    test.Element("error")(dError);

    CheckResult(dError, dFXDelta, dFXDeltaShifted, 
      std::string("FXDelta"), sFileName);
  }

} //CheckGreekFXDelta

void TestGreeks::CheckGreekFXDeltaConv(const std::string &sFileName)
{
  shared_ptr<finance::Derivative> pDerivative;
  finance::ModelOutput output;
  shared_ptr<SessionData> pSessionData;
  shared_ptr<TheoreticalModel> pModel (new TheoreticalModel);

  ReadXMLFile( sFileName, pSessionData, pDerivative, pModel);

  if ( pDerivative->IsCrossCurrency() )
  {  
    /// To refine the mesh
    shared_ptr<finance::QualityControl> pQualityControl(new QualityControl); 

    std::vector<double> FXDeltas;

    std::vector<ComputationQuality> CompQual;
    for (int i=0; i < ComputationQuality_Max; i++)
      CompQual.push_back((ComputationQuality) i); 

    for (std::vector<ComputationQuality>::const_iterator 
         iter = CompQual.begin();
         iter != CompQual.end();
         ++iter)
    { 
      pQualityControl->SetComputationQuality(*iter);
      pModel->SetQualityControl(pQualityControl);

      PriceContract(*pDerivative, pModel, output);

      if( output.HasFXDelta() )
        FXDeltas.push_back( output.GetFXDelta() );
      else
      {
        std::cout << "FX delta convergence testing does not apply." 
                  << std::endl;
        return;
      }
    }

    size_t 
      nIdx,
      nTests = FXDeltas.size();

    double
      dError,
      dPreviousFXDelta = FXDeltas[0];

    for ( nIdx = 0; nIdx < nTests; ++nIdx )
    {
      dError = fabs(FXDeltas[nIdx] - dPreviousFXDelta)
             / std::max(fabs(FXDeltas[nIdx]),1.);

      std::cout << FXDeltas[nIdx] << " " << dPreviousFXDelta 
        << " " << FXDeltas[nIdx] - dPreviousFXDelta << std::endl;
      
      if( dError  > 5./100. )
        std::cout << "Failed test for convergence FXDelta." << std::endl;

      dPreviousFXDelta = FXDeltas[nIdx];
    }

  }
}

/**
  Read xml file

  @param filename absolute path of file name
  @param sessionData
  @param derivative
  @param model
*/
void ReadXMLFile(std::string filename,
                 shared_ptr<finance::SessionData>& pSessionData,
                 shared_ptr<Derivative>& pDerivative,
                 shared_ptr<ito33::ihg::TheoreticalModel>& pModel)
{

 ihg::XML::PricingReaderRecursive reader(filename.c_str());

 pSessionData = reader.ReadSessionData();

 reader.ReadDerivative(pDerivative);
 
 pDerivative->SetSessionData(pSessionData);
 
 finance::DerivativeVisitorGoodType visitor;
 pDerivative->Visit(visitor);
 if( visitor.GetCBOption() )
   visitor.GetCBOption()->GetConvertibleBond()->SetSessionData(pSessionData);
 
 reader.ReadTheoreticalModel(pModel);
}

void PriceContract(finance::Derivative& pDerivative,
                   const shared_ptr<ito33::ihg::TheoreticalModel>& pModel,
                   finance::ModelOutput& output)
{
  try 
  {
    output = *pModel->Compute( pDerivative );
  }
  catch ( ito33::Exception& e )
  {
    printf("ITO33 exception:\n%s\n", e.GetFullMessage().c_str());
  }
  catch ( std::exception& e )
  {
    printf("std exception: %s\n", e.what());
  }
  catch ( ... )
  {
    puts("unknown exception!");
  }
  
}

} // namespace ihg

} // namespace ito33
