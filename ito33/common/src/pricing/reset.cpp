/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/reset.cpp
// Author:      Yann and David
// Created:     2004/11/03
// RCS-ID:      $Id: reset.cpp,v 1.13 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/vector.h"
#include "ito33/autoptr.h"
#include "ito33/dateutils.h"

#include "ito33/finance/sessiondata.h"

#include "ito33/finance/bondlike/reset.h"
#include "ito33/finance/bondlike/conversionschedule.h"
#include "ito33/finance/bondlike/conversionperiod.h"
#include "ito33/finance/bondlike/conversionpricereset.h"
#include "ito33/finance/bondlike/resetconversionschedule.h"

#include "ito33/pricing/reset.h"
#include "ito33/pricing/cbconversions.h"

using namespace ito33::finance;

namespace ito33
{

  ITO33_IMPLEMENT_AUTOPTR(pricing::Reset);

namespace pricing
{

Reset::Reset(const finance::Reset& reset) : CB()
{
  //-------------------------------------------------------------------------
  // Construct a conversion schedule in order to construct
  // the cbconversions
  //-------------------------------------------------------------------------

  const shared_ptr<finance::ResetConversionSchedule>& pResetConvSched = 
      reset.GetResetConversionSchedule();

  shared_ptr<finance::ConversionSchedule> pConvSched(new
      finance::ConversionSchedule());
  
  m_dCurrentConvPrice = pResetConvSched->GetCurrentConversionPrice();
  m_dInitialConvPrice = pResetConvSched->GetInitialConversionPrice();

  // Check if cross currency.  Note that this function will be repeated
  // in GetCBBaseData below.  However, exchange rates are needed
  // to construct the conversion period, which must be done before
  // GetCBBaseData is called.
  GetCrossCurrencyData(*reset.GetSessionData(), reset.GetNumeraire() );

  //convert to the currency of the bond
  if ( m_bIsCrossCurrency )
  {
    double dFixedFXRate = reset.GetFixedFXRate();

    m_dCurrentConvPrice *= dFixedFXRate;
    m_dInitialConvPrice *= dFixedFXRate;
  }

  double dCurrentRatio = 
    reset.GetBondTerms()->GetNominal() / m_dCurrentConvPrice;

  shared_ptr<finance::ConversionPeriod> pConvPeriod(new
      finance::ConversionPeriod(
                                pResetConvSched->GetStartDate(),
                                pResetConvSched->GetEndDate(),
                                dCurrentRatio
                                ));

  pConvPeriod->SetCash(pResetConvSched->GetCash() );

  pConvSched->AddConversionPeriod(pConvPeriod);

  pConvSched->SetKeepAccrued(   pResetConvSched->GetKeepAccrued() );
  pConvSched->SetForfeitCoupon( pResetConvSched->GetForfeitCoupon() );
  
  //create the reset conversion schedule
  m_conversions = CBConversions(pConvSched, 
                                reset.GetSessionData()->GetValuationDate());

  // note: must construct m_conversions before calling
  GetCBBaseData(reset);
 
  //-------------------------------------------------------------------------
  //
  // Now store the reset information
  //
  //--------------------------------------------------------------------------
   
  // Just save the trivial paramaters. They should have already been
  // verified at the finance level.
  m_FlooredBy = pResetConvSched->GetResetFlooredBy();
  
  // Convert the reset dates into times. Assume the dates are already sorted.
  // Also save all the rate and multipliers
  finance::ResetConversionSchedule::Elements 
    conversionPriceResets = pResetConvSched->GetAll();

  size_t nNbTerms = conversionPriceResets.size();

  m_pdResetTimes.resize(nNbTerms);
  m_pdCapRates.resize(nNbTerms);
  m_pdFloorRates.resize(nNbTerms);
  m_pdMultipliers.resize(nNbTerms);

  finance::ResetConversionSchedule::Elements::const_iterator iterTerms;
  size_t nCounter = 0;

  for (iterTerms = conversionPriceResets.begin();
       iterTerms != conversionPriceResets.end();
       ++iterTerms)
  {
    double dTime              = GetDoubleFrom( (*iterTerms)->GetDate() );
    m_pdResetTimes[nCounter]  = dTime;
    m_pdCapRates[nCounter]    = (*iterTerms)->GetCap();
    m_pdFloorRates[nCounter]  = (*iterTerms)->GetFloor();
    m_pdMultipliers[nCounter] = (*iterTerms)->GetMultiplier();
    nCounter++;

  } //end loop over the conversion schedule

  // Must have at least one date if this function is called.
  // Again, this should have already been verified
  ASSERT_MSG( m_pdResetTimes.size() >= 1,
              "No reset dates specified in call to SetResetData"); 

} //Reset::Reset(const finance::Reset& reset)


void Reset::FlipFloorAndCap()
{

  // For 1D pricing, the variable is a function of the conversion ratio k.
  // For normal path dep pricing, we are capping and flooring the
  // conversion price N/k.  Thus, for 1D, need to flip the caps and floors.
  size_t nNbTerms = m_pdCapRates.size();
  size_t nIdx;
  for (nIdx = 0; nIdx < nNbTerms; nIdx++)
  {
    double dTmp = m_pdCapRates[nIdx];
    m_pdCapRates[nIdx] = 1.0/m_pdFloorRates[nIdx];
    m_pdFloorRates[nIdx] = 1.0/dTmp;
  }

} // Reset::FlipFloorAndCap()

} //namespace pricing

} //namespace ito33
