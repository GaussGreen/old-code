/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/mandatory.cpp
// Author:      Wang
// Created:     2004/08/16
// RCS-ID:      $Id: mandatory.cpp,v 1.39 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file common/src/pricing/mandatory.cpp
   @brief implementation of the contract for mandatory
 */

#include "ito33/sharedptr.h"
#include "ito33/autoptr.h"
#include "ito33/vector.h"
#include "ito33/dateutils.h"
#include "ito33/debug.h"

#include "ito33/finance/payoffconstant.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/finance/bondlike/bondliketerms.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/callfixedshare.h"
#include "ito33/finance/bondlike/convertiblelike.h"
#include "ito33/finance/bondlike/percslike.h"
#include "ito33/finance/bondlike/pepslike.h"
#include "ito33/finance/bondlike/generalizedpepslike.h"
#include "ito33/finance/bondlike/pepsaveragingperiod.h"

#include "ito33/pricing/cashflows.h"
#include "ito33/pricing/cbconversions.h"
#include "ito33/pricing/mandatory.h"
#include "ito33/pricing/callvariableshare.h"
#include "ito33/pricing/callfixedshare.h"

namespace ito33
{

  ITO33_IMPLEMENT_AUTOPTR(pricing::Mandatory);

  using namespace numeric;
  using namespace finance;

namespace pricing
{

Mandatory::Mandatory(const finance::PERCSLike& percs)
    : m_callType(MandatoryCallType_FixedCash), // only fixed cash call is allowed
      m_bHasStockAveraging(false),
      m_bHasAveragingPeriod(false)
{
  shared_ptr<BondLikeTerms> pblk( percs.GetBondLikeTerms() );
  
  GetBondLikeTermsData( *pblk, *percs.GetSessionData(), percs.GetNumeraire() );

  m_pCashFlows = make_ptr( new CashFlows(pblk->GetCashDistribution(),
                                         pblk->GetNominal()) );  

  m_calls = CBCalls( percs.GetCallSchedule() );

  m_conversions = MandatoryConversion(percs);

  GetConvertibleLikeData(percs);

  m_pPayoff = make_ptr( new PayoffConstant(0) );
}


Mandatory::Mandatory(const finance::PEPSLike& peps)
    : m_callType(MandatoryCallType_FixedCash),// it works for no call
    m_bHasStockAveraging(false),
    m_bHasAveragingPeriod(false)
{
  shared_ptr<BondLikeTerms> pblk( peps.GetBondLikeTerms() );
  
  GetBondLikeTermsData(*pblk, *peps.GetSessionData(), peps.GetNumeraire());

  m_pCashFlows = make_ptr( new CashFlows(pblk->GetCashDistribution(),
                                         pblk->GetNominal()) );  

  if(peps.GetCallFixedCash())
  {
    m_callType = MandatoryCallType_FixedCash;
    m_calls = CBCalls( peps.GetCallFixedCash() );
  }
  else if(peps.GetCallFixedShare())
  {
    m_callType = MandatoryCallType_FixedShare;
    m_callFixedShare = CallFixedShare
              ( peps.GetCallFixedShare(), peps.GetMaturityDate() );
  }

  m_conversions= MandatoryConversion(peps);

  GetConvertibleLikeData(peps);

  m_pPayoff = make_ptr( new PayoffConstant(0) );
}


Mandatory::Mandatory(const finance::GeneralizedPEPSLike& peps)
    : m_callType(MandatoryCallType_FixedCash), // it works for no call
      m_bHasStockAveraging(false),
      m_bHasAveragingPeriod(false)
{
  shared_ptr<BondLikeTerms> pblk( peps.GetBondLikeTerms() );
  
  GetBondLikeTermsData(*pblk, *peps.GetSessionData(), peps.GetNumeraire());
  
  std::vector<Date> repaymentDates(1, peps.GetMaturityDate());

  m_pCashFlows = make_ptr( new CashFlows(pblk->GetCashDistribution(),
                                         pblk->GetNominal()) );  

  if (peps.GetCallFixedCash())
  {
    m_callType = MandatoryCallType_FixedCash;
    m_calls = CBCalls( peps.GetCallFixedCash() );
  }
  else if(peps.GetGeneralizedPEPSLikeCall())
  {
    if(peps.GetGeneralizedPEPSLikeCall()->GetType()
        == finance::GeneralizedPEPSLikeCallType_FixedShare)
    {
      m_callType = MandatoryCallType_FixedShare;
      m_callFixedShare = CallFixedShare( peps );
    }
    else
    {
      m_callType = MandatoryCallType_VariableShare;
      m_callVariableShare = CallVariableShare( peps );
    }
  }

  m_conversions = MandatoryConversion(peps);

  GetConvertibleLikeData(peps);

  m_pPayoff = make_ptr( new PayoffConstant(0) );

  if ( peps.HasAveragingPeriod() )
  {
    shared_ptr<finance::PEPSAveragingPeriod> pAvgPeriod = peps.GetAveragingPeriod();

    m_dAverageStartTime = GetDoubleFrom( pAvgPeriod->GetAverageStartDate() );
    m_dAverageEndTime   = GetDoubleFrom( pAvgPeriod->GetAverageEndDate() );
    m_nNumberOfSamplingAverages = pAvgPeriod->GetNbSamplingAverages();

    m_bHasStockAveraging    = pAvgPeriod->HasStockAveraging();
    m_bHasAveragingPeriod   = true;

    m_nNbSamplesUsed = pAvgPeriod->GetNbSamplesUsed();

    if ( m_bHasStockAveraging )
      m_dCurrentAverage = pAvgPeriod->GetCurrentStockAverage();
    else
      m_dCurrentAverage = pAvgPeriod->GetCurrentConversionRatioAverage();

    // Note that if the current average has been set but the
    // valuation date is before the begining of the averaging period
    // the current average does not affect the results. The current
    // average set by the user if any should simply be ignored
    if ( peps.GetSessionData()->GetValuationDate() 
            <= pAvgPeriod->GetAverageStartDate() )
    {
      if ( m_bHasStockAveraging)
        m_dCurrentAverage = peps.GetSessionData()->GetSpotSharePrice();
      else
        m_dCurrentAverage = 1.0;

      m_nNbSamplesUsed = 0;
    }

  }//end averaging period

}


double Mandatory::GetClaim(double dTime, bool bPlus) const
{
  return GetNominal() + m_pCashFlows->GetAccruedInterest(dTime, bPlus);
}


} //namespace pricing

} //namespace ito33
