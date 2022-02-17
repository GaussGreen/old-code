///////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/sharedependentconversion.cpp
// Purpose:     Share dependent conversion
// Author:      Ito33
// Created:     2005/01/13
// RCS-ID:      $Id: sharedependentconversion.cpp,v 1.7 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
///////////////////////////////////////////////////////////////////////////

/**
   @todo We should probably need to use AddMonthsAdjustedForEndOfMonth 
         instead of AddMonth.
 */
#include "ito33/sharedptr.h"
#include "ito33/autoptr.h"
#include "ito33/dateutils.h"

#include "ito33/numeric/predicatetime.h"

#include "ito33/finance/bondlike/cocotype.h"
#include "ito33/finance/bondlike/bondliketerms.h"
#include "ito33/finance/bondlike/sharedependentconversion.h"
#include "ito33/finance/bondlike/bonderror.h"

#include "ito33/pricing/conversionprovisions.h"
#include "ito33/pricing/sharedependentconversion.h"
#include "ito33/pricing/cblikeparams.h"
#include "ito33/pricing/cbcalls.h"
#include "ito33/pricing/cbconversions.h"

extern const double DEFAULT_BIG_RATIO = 1.e99;

using namespace ito33::finance;

namespace ito33
{
  //Implementation of the AutoPtrDeleter for ShareDependentConversion class.
  ITO33_IMPLEMENT_AUTOPTR(pricing::ShareDependentConversion);
  
namespace pricing
{
   
ShareDependentConversion::ShareDependentConversion(
  const shared_ptr<finance::ShareDependentConversion> &pShareDeConv,
                  Date valuationDate):   
  CBConversions()
  {  
    if ( !pShareDeConv)
      return;

    // Get the data from the finance level class
    m_bKeepAccrued = pShareDeConv->GetKeepAccrued();
    m_bForfeitCoupon = pShareDeConv->GetForfeitCoupon();
    
    // share dependent info
    Date startDate             = pShareDeConv->GetStartDate();
    m_dBaseRatio               = pShareDeConv->GetBaseRatio();
    m_dIncrementalShareFactor  = pShareDeConv->GetIncrementalShareFactor();
    m_dFixedStrike             = pShareDeConv->GetFixedStrike();
    m_dCapRatio                = pShareDeConv->GetCapRatio();
    if ( m_dCapRatio < 0.0 )
      m_dCapRatio = DEFAULT_BIG_RATIO;

    m_bHasResetTime = pShareDeConv->HasResetDate();

    if ( m_bHasResetTime )
      m_dResetTime = GetDoubleFrom( pShareDeConv->GetResetDate() );

    // coco related info
    double dTriggerRate        = 0.0;
    double dExtremeTriggerRate = 0.0;
    double dRateChange         = 0.0;
    CoCoType cocoType          = CoCoType_Max;
    bool bIsLastTriggerMet     = false;    
    
    if ( pShareDeConv->HasCoCo() )
    {
      dTriggerRate        = pShareDeConv->GetTrigger();
      dExtremeTriggerRate = pShareDeConv->GetExtremeTrigger();
      dRateChange         = pShareDeConv->GetChangeRate();
      cocoType            = pShareDeConv->GetCoCoType();
      bIsLastTriggerMet   = pShareDeConv->GetIsLastTriggerConditionMet();
    }
         
    // Coco is only supported until the reset time. Construct 
    // conversion periods from start to reset time. If reset time
    // is undefined, go to the end time
    Date finalEndDate = pShareDeConv->GetEndDate();

    if ( m_bHasResetTime && valuationDate < pShareDeConv->GetResetDate() )
      finalEndDate = pShareDeConv->GetResetDate();    

    // if the trigger was met, and we convert as of check date, then
    // the coco feature is no longer needed 
    if (  bIsLastTriggerMet &&
         (cocoType == CoCoType_CheckQuarterlyAndConvertAsOfCheckDate ||
          cocoType == CoCoType_CheckAnyTimeAndConvertAsOfCheckDate) )
    {
      cocoType     = CoCoType_Max;
      dTriggerRate = 0.0; 
    }

    // Breakup the window if quarterly type coco was not removed due to 
    // lastTriggerMet above
    if ( dTriggerRate > 0.0 && 
         IsWindowTooLong(startDate, finalEndDate, cocoType,
                         dTriggerRate, dExtremeTriggerRate, dRateChange, 
                         bIsLastTriggerMet) ) 
    {
      BreakupWindow(startDate, finalEndDate, m_dBaseRatio, 0.0,
                    cocoType, dTriggerRate, dExtremeTriggerRate, dRateChange,
                    bIsLastTriggerMet, valuationDate);
    }
    else //regular conversion i.e no quarterly coco
    {

      if ( bIsLastTriggerMet == true &&
           cocoType == CoCoType_CheckQuarterlyAndConvertDuringNextQuarter &&
           valuationDate > pShareDeConv->GetStartDate() &&
           valuationDate < finalEndDate )
      {
        dTriggerRate = 0.0;
        cocoType = CoCoType_Max;
      }

      // start time ---> reset time (or end time)
      m_pdStartTimes.push_back( GetDoubleFrom( startDate ) );
      m_pdEndTimes.push_back( GetDoubleFrom( finalEndDate ) );
      m_pdCashs.push_back(0);
      m_pdRatios.push_back( m_dBaseRatio );
      m_pdTriggerRates.push_back( dTriggerRate );
      m_pCoCoTypes.push_back( cocoType );
      m_pbIsActive.push_back( true );

    }
   
    // only add this part is a reset time has been specified
    if ( m_bHasResetTime && valuationDate < pShareDeConv->GetResetDate())
    {

      //this coco type is accepted since it is all the time and no path
      //dependet feature are needed in this case for CoCo
      if ( cocoType != finance::CoCoType_CheckAnyTimeAndConvertOnCheckDate )
      {
        //coco after the reset date is not allowed
        dTriggerRate = 0;
        cocoType = finance::CoCoType_Max;
      }

      //add the remaining part from reset date to maturity date
      m_pdStartTimes.push_back( m_dResetTime );
      m_pdEndTimes.push_back( GetDoubleFrom( pShareDeConv->GetEndDate() ) );
      m_pdCashs.push_back(0);
      m_pdRatios.push_back(  m_dBaseRatio );
      m_pdTriggerRates.push_back( dTriggerRate );
      m_pCoCoTypes.push_back( cocoType );
      m_pbIsActive.push_back( true );
    }

    //adjust the number of conversions appropriately
    m_nNbConversions = m_pdStartTimes.size();

}//ShareDependentConversion::init()


bool ShareDependentConversion::GetGrossParities
     (const double* pdS, size_t nNbS, const double* pdNewShares, 
      double* pdValues) const
{
  
  size_t nIdxConversion = m_pParams->GetIndexConversion();
  
  if ( nIdxConversion == INVALIDINDEX )
    return false;

  //we are between reset time and maturity time
  //in this period we are solving for a regular cb conversion
  if ( HasResetTime() && 
      numeric::IsEqualOrAfter( m_pdStartTimes[nIdxConversion], m_dResetTime) )
    return CBConversions::GetGrossParities(pdS, nNbS, pdNewShares, pdValues);
 
  double dStrike = GetStrike();

  size_t nIdxS;
  for (nIdxS = 0; nIdxS < nNbS ; nIdxS++)
  {
    double dS = pdS[nIdxS];
    double dCurrentRatio;

    if ( dS > dStrike )
    {
      dCurrentRatio = m_dBaseRatio 
        + (dS - dStrike) / dS * m_dIncrementalShareFactor;  
      //conversion ratio can never be larger than cap ratio
      dCurrentRatio = std::min(dCurrentRatio, m_dCapRatio);
    }
    else //conversion ratio can never be less than base ratio
      dCurrentRatio = m_dBaseRatio;

    if ( m_pParams->GetCBLike().IsCrossCurrency() )    
      dCurrentRatio *= m_pParams->GetFXRate();    

     pdValues[nIdxS] = dCurrentRatio * pdNewShares[nIdxS];

  } //end loop

  return true;
}

double ShareDependentConversion::GetStrike() const
{ 
  double dStrike = m_dFixedStrike;

  if ( dStrike < 0. )
    dStrike = (m_pParams->GetClaim() - m_pParams->GetAccruedInterest()) 
                     / m_dBaseRatio;

  return dStrike;
}


} // namespace pricing

} // namespace ito33
