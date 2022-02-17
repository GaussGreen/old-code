///////////////////////////////////////////////////////////////////////////////
// File:             common/src/finance/impliedvolpriceconverter.cpp
// Purpose:          functions to get price from implied vol and vice versa
// Author:           ITO 33                                           
// Created:          2005/07/18
// RCS-ID:           $Id: impliedvolpriceconverter.cpp,v 1.6 2006/08/19 23:06:55 wang Exp $                                                  
// Copyright         (C) 2005 - 2006 Trilemma LLP
///////////////////////////////////////////////////////////////////////////////

/**

  @file ito33/finance/impliedvolpriceconverter.cpp

*/
#include "ito33/sharedptr.h"
#include "ito33/link.h"

#include "ito33/finance/option.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/impliedvolpriceconverter.h"

#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardrateflat.h"
#include "ito33/ihg/theoreticalmodel.h"

#include "ihg/impliedvol.h"

ITO33_FORCE_LINK_MODULE(IHGPriceOption);

namespace ito33
{

namespace finance
{

double GetPriceFromImpliedVol(const Option& option, double dImpliedVol)
{
  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

  shared_ptr<ihg::VolatilityFlat> pVol( new ihg::VolatilityFlat(dImpliedVol) );
  pModel->SetVolatility( pVol );

  shared_ptr<ihg::HazardRateFlat> pHR( new ihg::HazardRateFlat(0.0) );
  pModel->SetHazardRate( pHR );

  return pModel->Compute(option)->GetPrice();
}

double GetVegaFromImpliedVol(const Option& option, double dImpliedVol)
{
  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

  shared_ptr<ihg::VolatilityFlat> pVol( new ihg::VolatilityFlat(dImpliedVol) );
  pModel->SetVolatility( pVol );

  shared_ptr<ihg::HazardRateFlat> pHR( new ihg::HazardRateFlat(0.0) );
  pModel->SetHazardRate( pHR );

  shared_ptr<ComputationalFlags> pFlags(new ComputationalFlags);
  pFlags->SetComputeVega(true);

  pModel->SetExternalFlags(pFlags);

  return pModel->Compute(option)->GetVega();
}

double GetImpliedVolFromPrice(const Option& option, double dPrice)
{
  //make a local copy of the option
  //since implied vol works on the 
  //market price
  Option optionTmp(option);

  //set the market price to be dPrice
  optionTmp.SetMarketPrice(dPrice);

  //create an hazard rate of zero
  shared_ptr<ihg::HazardRateFlat> pHR( new ihg::HazardRateFlat( 0.0 ) );

  ihg::ImpliedVol impliedVol(optionTmp, pHR);

  return impliedVol.Compute();
}

} // namespace finance

} // namespace ito33
