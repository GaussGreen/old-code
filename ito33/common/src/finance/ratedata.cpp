/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/ratedata.cpp
// Purpose:     Implemenation of the RateData class
// Created:     2006/03/17
// RCS-ID:      $Id: ratedata.cpp,v 1.5 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/yieldcurve.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/spotfxrates.h"
#include "ito33/finance/moneymarket.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/ratedata.h"

extern const ito33::finance::Error ITO33_RATEDATA_INVALID_YIELD_CURVE,
                          ITO33_RATEDATA_MISSING_YIELDCURVE,
                          ITO33_NUMERAIRE_INVALID,
                          ITO33_RATEDATA_INVALID_SPOTFXRATES;


namespace ito33
{

namespace finance
{

// Construct with an empty SpotFXRates
RateData::RateData() : m_pSpotFXRates(new SpotFXRates)
{
}


void RateData::SetYieldCurve(const std::string& currency, 
                             const shared_ptr<YieldCurve>& pYieldCurve)
{
  SetYieldCurve(shared_ptr<Numeraire>(new Numeraire(currency)),
                pYieldCurve);
}


void RateData::SetYieldCurve(const shared_ptr<Numeraire>& pNumeraire, 
                             const shared_ptr<YieldCurve>& pYieldCurve)
{ 

  // standard checks for null data
  CHECK_COND( pYieldCurve, ITO33_RATEDATA_INVALID_YIELD_CURVE );
  CHECK_COND( pNumeraire, ITO33_NUMERAIRE_INVALID );

  // Check if the curve (money market) has already been set
  YCElements::iterator iter = m_MoneyMarkets.find( pNumeraire->GetCode() );

  if ( iter != m_MoneyMarkets.end() )
  {
    // update
    (iter->second)->SetYieldCurve(pYieldCurve);

    return;
  }

  // Not set, so create a new moneymarket and add to the map
  shared_ptr<MoneyMarket> pMarket(new MoneyMarket(pNumeraire, pYieldCurve));

  m_MoneyMarkets[pNumeraire->GetCode()] = pMarket;

}


void RateData::SetSpotFXRates(const shared_ptr<SpotFXRates>& pSpotFXRates)
{
  // standard check for null data
  CHECK_COND( pSpotFXRates, ITO33_RATEDATA_INVALID_SPOTFXRATES );

  m_pSpotFXRates = pSpotFXRates;
}


const shared_ptr<YieldCurve>& 
RateData::GetYieldCurve(const shared_ptr<Numeraire>& pNumeraire)
{

  // standard check for null data
  CHECK_COND( pNumeraire, ITO33_NUMERAIRE_INVALID );

  // find the money market
  YCElements::iterator iter = m_MoneyMarkets.find( pNumeraire->GetCode() );

  // check if the curve was found
  CHECK_COND(iter != m_MoneyMarkets.end(), ITO33_RATEDATA_MISSING_YIELDCURVE);

  return (iter->second)->GetYieldCurve();

}


void RateData::SetYieldCurveForMesh(const shared_ptr<Numeraire>& pNumeraire, 
                                    const shared_ptr<YieldCurve>& pYieldCurve)
{ 
  // Internal function, but still check.  Assume yield curve is valid
  ASSERT_MSG( pNumeraire, "Invalid numeraire in SetYieldCurveForMesh.");
  ASSERT_MSG( pYieldCurve, "Invalid yield curve in SetYieldCurveForMesh.");

  // Check if the curve (money market) has already been set
  YCElements::iterator iter = m_MoneyMarkets.find( pNumeraire->GetCode() );

  // if not set, something is wrong
  ASSERT_MSG( iter != m_MoneyMarkets.end(), 
              "Cannot find money market in SetYieldCurveForMesh" );

  (iter->second)->SetYieldCurveForMesh(pYieldCurve);

  return;
}


const shared_ptr<YieldCurve>& 
RateData::GetYieldCurveForMesh(const shared_ptr<Numeraire>& pNumeraire)
{

  // internal function, but still check
  ASSERT_MSG( pNumeraire, "Missing numeraire in SetYieldCurveForMesh.");

  // find the curve
  YCElements::iterator iter = m_MoneyMarkets.find( pNumeraire->GetCode() );

  // check if found
  ASSERT_MSG( iter != m_MoneyMarkets.end(), 
              "Cannot find money market in GetYieldCurveForMesh.");

  return (iter->second)->GetYieldCurveForMesh();

}


void RateData::Dump(XML::Tag& tagParent) const
{

  XML::Tag tagRateData(XML_TAG_RATEDATA_ROOT, tagParent);

  // output the exchange rates
  m_pSpotFXRates->Dump(tagRateData);

  // output the yield curves
  YCElements::const_iterator iterYC;

  for (iterYC = m_MoneyMarkets.begin(); iterYC != m_MoneyMarkets.end(); iterYC++)  
    (iterYC->second)->Dump(tagRateData);
  
}

} // namespace finance

} // namespace ito33
