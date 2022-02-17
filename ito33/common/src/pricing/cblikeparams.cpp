/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/cblikeparams.cpp
// Purpose:     params class for CB-like
// Author:      Wang
// Created:     2004/08/19
// RCS-ID:      $Id: cblikeparams.cpp,v 1.91 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <algorithm>
#include "ito33/afterstd.h"

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"
#include "ito33/useexception.h"
#include "ito33/dateutils.h"

#include "ito33/finance/dividends.h"
#include "ito33/finance/qualitycontrol.h"
#include "ito33/finance/yieldcurve.h"
#include "ito33/finance/payoffdiscrete.h"

#include "ito33/numeric/numparams.h"

#include "ito33/pricing/dividendevents.h"
#include "ito33/pricing/paymentevent.h"
#include "ito33/pricing/eventmanager.h"
#include "ito33/pricing/cbevent.h"

#include "ito33/pricing/cashflows.h"

#include "ito33/pricing/cb.h"
#include "ito33/pricing/cbconversions.h"
#include "ito33/pricing/cbparams.h"
#include "ito33/pricing/cblikeparams.h"

extern const ito33::Error ITO33_UNEXPECTED;

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(pricing::CBLikeParams);

namespace pricing
{

using namespace finance;
using namespace numeric;
using namespace numeric::mesh;

typedef shared_ptr<CBEvent> CBEventPtr;

void CBLikeParams::InitMakeWholeData()
{  
  // TOFIX : have problem here when soft call period and cash distribution
  // end after maturity

  CallProvisions& calls = *(GetCalls());
  CashFlows& cashFlows = *(GetCashFlows());

  // calculate finalPermium for make-whole
  if ( calls.HasMakeWhole() )
  {
    size_t nIdx;
    double dTime = 0; // initialize to avoid warning, the value makes no sense
    for (nIdx = 0;
         nIdx < calls.GetNbCalls() && calls.GetTriggerRate(nIdx) > 0;
         nIdx++)
      dTime = calls.GetEndTime(nIdx);
    // after the loop, dTime will be the final time of all soft call periods
    // without call notice

    // so add call notice period if any, note that we have shifted 
    // original call period to the left by dCallNoticePeriod.
    double dCallNoticePeriod = calls.GetNoticePeriod();
    if ( dCallNoticePeriod > 0 ) 
      dTime += dCallNoticePeriod;

    size_t &nIdxAfter = m_nIndexAfterLastCouponMakeWhole; // shorter name
    for (nIdxAfter = 0;
            nIdxAfter < cashFlows.GetNbCashFlows() 
         && cashFlows.GetTime(nIdxAfter) <= dTime + TIMETOLERANCE;
         nIdxAfter++)
      ;

    if (calls.GetMakeWholeType() == MakeWholeType_Premium)
    {
      m_dMakeWholeFinalPremium = calls.GetMakeWholePremium();
      for (nIdx = 0; nIdx < nIdxAfter; nIdx++)
        m_dMakeWholeFinalPremium -= cashFlows.GetAmount(nIdx);

      // m_dMakeWholeFinalPremium can be negative as well as
      // for m_dMakeWholeCurrentPremium.

      // however GetCurrentMakeWholePositivePremium() will return
      // no negative value
    }
    else
    {
      m_dMakeWholeFinalPremium = 0;
    }
  }
}

void CBLikeParams::Init()
{
  Params::Init();

  m_CBEventManager.Init();

  // setup the Params pointer for call, conversion and put
  GetPuts()->SetParams(this);
  GetCalls()->SetParams(this);
  GetConversions()->SetParams(this);

  InitMakeWholeData();

  // Construct events  ===============================================

  // Dividend events
  ConstructDividendEvents();
  
  shared_ptr<CBEvent> pcbevent;
  
  // New share events
  if ( HasNewShare() )
  {
    Date TmpDate = m_cbLike.GetFiscalYearStartDate();

    if ( IsEqualOrBefore(GetDoubleFrom(TmpDate), m_dValuationTime) )
    {
      while( IsEqualOrBefore(GetDoubleFrom(TmpDate), m_dValuationTime) )
        TmpDate.AddYears(1);
    }
    else
    {
      while ( IsAfter(GetDoubleFrom(TmpDate), m_dValuationTime) )
        TmpDate.AddYears(-1);

      TmpDate.AddYears(1);
    }
    //Then, TmpDate is the first fiscal year after the valuation time

    while ( IsEqualOrBefore(GetDoubleFrom(TmpDate), m_dStoppingTime) )
    {
      pcbevent = make_ptr( new CBEvent( GetDoubleFrom(TmpDate),
                                        CBET_StartOfYear ) );
      m_CBEventManager.AddEvent(pcbevent);
      TmpDate.AddYears(1);
    }

    //Dividends added to m_CBEventManager     
    if ( m_pDividends )
    {     
      Dividends::Elements::iterator pDividend;
      
      for (pDividend = m_pDividends->GetAll().begin(); 
           pDividend != m_pDividends->GetAll().end(); 
           ++pDividend)
      {
        double dDividendTime = GetDoubleFrom(pDividend->date);

        // The dividend at valuation time will be taken into acount, but not the one
        // at the maturity
        if (   !IsBefore(dDividendTime, m_dValuationTime) 
            && IsBefore(dDividendTime, m_dStoppingTime) )
        {
          switch (pDividend->type)
          {
            case Dividend::PseudoCash :
            case Dividend::Cash :          
            case Dividend::Yield :          
              m_CBEventManager.AddEvent 
                  ( make_ptr( new CBEvent(dDividendTime, CBET_PayDividend) ) );
              break;
            default:
              throw EXCEPTION_MSG
                  (
                    ITO33_UNEXPECTED,
                    TRANS("Sorry, unknown dividend type.")
                  );
          }
        }
      }
    }
  }  

  size_t nIdx;

  // Cashflow events
  CashFlows& cashFlows = *(GetCashFlows());
  size_t nNbCashFlows = cashFlows.GetNbCashFlows();

  for (nIdx = 0; nIdx < nNbCashFlows; nIdx++)
  {
    double dTime = cashFlows.GetTime(nIdx);

    // Coupon at maturity will be taken into account
    if (   !IsBefore(dTime, m_dValuationTime) 
        && !IsAfter(dTime, m_dStoppingTime) )
    {               
      // add to base event manager (will be applied to price)
      m_eventManager.AddEvent
                     ( shared_ptr<Event>(
                        new PaymentEvent(dTime, cashFlows.GetAmount(nIdx))) );

      // add to cb event manager 
      // coupon event is continuous for computation of accrued interest 
      // and discret for computation of coupon amount
      m_CBEventManager.AddEvent
                       ( make_ptr( new CBEvent(dTime, CBET_PayCoupon, nIdx) ) );
    }
  }

  // Put events
  CBPuts& puts = *(GetPuts());
  size_t nNbPuts = puts.GetNbPuts();

  for (nIdx = 0; nIdx < nNbPuts; nIdx++)
  {
    double dTime = puts.GetTime(nIdx);
    if (    !IsBefore(dTime, m_dValuationTime)
          && !IsAfter(dTime, m_dStoppingTime) )
    {
      pcbevent = make_ptr( new CBEvent(dTime, CBET_MonoDatePut, nIdx) );
      m_CBEventManager.AddEvent(pcbevent);
    }
  }

  // Call events
  CallProvisions& calls = *(GetCalls());
  size_t nNbCalls = calls.GetNbCalls();

  for (nIdx = 0; nIdx < nNbCalls; nIdx++)
  {
    double dStartTime = calls.GetStartTime(nIdx);
    double dEndTime = calls.GetEndTime(nIdx);

    if ( !calls.IsContinuous(nIdx) ) // mono date event
    {
      if (   !IsBefore(dStartTime, m_dValuationTime)
          && IsBefore(dStartTime, m_dStoppingTime) )
      {
        pcbevent = make_ptr( new CBEvent(dStartTime, CBET_MonoDateCall, nIdx) );
        m_CBEventManager.AddEvent(pcbevent);
      }
    }
    else
    {
      if ( !IsAfter(dStartTime, m_dValuationTime) )
      {
        if ( IsAfter(dEndTime, m_dValuationTime) )
        {
          pcbevent = make_ptr(new CBEvent(m_dValuationTime, CBET_StartCall,
                                          nIdx) );
          m_CBEventManager.AddEvent(pcbevent);
        }
      }
      else if ( IsBefore(dStartTime, m_dStoppingTime) )
      {
        pcbevent = make_ptr( new CBEvent(dStartTime, CBET_StartCall, nIdx) );
        m_CBEventManager.AddEvent(pcbevent);
      }

      if ( !IsBefore(dEndTime, m_dStoppingTime) )
      {
        if ( IsBefore(dStartTime, m_dStoppingTime) )
        {
          pcbevent = make_ptr(new CBEvent(m_dStoppingTime, CBET_EndCall, nIdx));
          m_CBEventManager.AddEvent(pcbevent);
        }
      }
      else if ( IsAfter(dEndTime, m_dValuationTime) )
      {
        pcbevent = make_ptr( new CBEvent(dEndTime, CBET_EndCall, nIdx) );
        m_CBEventManager.AddEvent(pcbevent);      
      }
    }
  }

  // Conversion events
  ConversionProvisions& conv = *(GetConversions());
  size_t nNbConversions = conv.GetNbConversions();

  for (nIdx = 0; nIdx < nNbConversions; nIdx++)
  {
    double dStartTime = conv.GetStartTime(nIdx);
    double dEndTime = conv.GetEndTime(nIdx);

    if ( !conv.IsContinuous(nIdx) )
    {
      if (   !IsBefore(dStartTime, m_dValuationTime)
           && IsBefore(dStartTime, m_dStoppingTime) )
      {
        pcbevent = make_ptr( new CBEvent(dStartTime, CBET_MonoDateConversion,
                                         nIdx) );
        m_CBEventManager.AddEvent(pcbevent);
      }
    }
    else
    {
      if ( IsBefore(dStartTime, m_dValuationTime) )
      {
        if ( IsAfter(dEndTime, m_dValuationTime) )
        {
          pcbevent = make_ptr( new CBEvent(m_dValuationTime,
                                           CBET_StartConversion, nIdx) );
          m_CBEventManager.AddEvent(pcbevent);
        }
      }
      else if ( IsBefore(dStartTime, m_dStoppingTime) )
      {
        pcbevent = make_ptr( new CBEvent(dStartTime, CBET_StartConversion,
                                         nIdx) );
        m_CBEventManager.AddEvent(pcbevent);
      }

      if ( !IsBefore(dEndTime,  m_dStoppingTime) )
      {
        if ( IsBefore(dStartTime, m_dStoppingTime) )
        {
          pcbevent = make_ptr( new CBEvent(m_dStoppingTime,
                                           CBET_EndConversion, nIdx) );
          m_CBEventManager.AddEvent(pcbevent);
        }
      }
      else if ( IsAfter(dEndTime, m_dValuationTime) )
      {
        pcbevent = make_ptr( new CBEvent(dEndTime, CBET_EndConversion, nIdx) );
        m_CBEventManager.AddEvent(pcbevent);
      }
    }
  }

  // Add maturity date in the cb event manager.
  pcbevent = make_ptr( new CBEvent(m_dStoppingTime, CBET_Maturity) );
  m_CBEventManager.AddEvent(pcbevent);

  // initialize mini cb params that are going to be used for solving the
  // small cb for call notice. 
  if ( HasNoticePeriod() ) 
    SetupCallNoticeParams();
}

void CBLikeParams::UpdateDiscreteCBEventData()
{
  // Put event is discrete, the index is activated only in DoEvents
  m_nIdxPut = INVALIDINDEX;

  // the current coupon amount is 0. The value is changed only during DoEvents
  m_dCouponAmount = 0;
}

void CBLikeParams::SetInitialState(double dTime)
{
  StateValuesMustBeReset();
  
  Params::SetInitialState(dTime);

  StateValuesMustBeReset();

  m_CBEventManager.SetInitialState(dTime);

  UpdateDiscreteCBEventData();

  size_t nIdx, nNb;

  // Call index
  CallProvisions& calls = *(GetCalls());
  nNb = calls.GetNbCalls();
  for (nIdx = nNb - 1; nIdx < nNb; nIdx--)
  {
    if (    !IsBefore(dTime, calls.GetStartTime(nIdx)) 
         && !IsAfter(dTime, calls.GetEndTime(nIdx)) )
      break;
  }
  m_nIdxCall = nIdx;
  // if nothing found, at the end of loop, we have just
  // m_nIdx = INVALIDINDEX;

  // Conversion index
  ConversionProvisions& conv = *(GetConversions());
  nNb = conv.GetNbConversions();
  for (nIdx = nNb - 1; nIdx < nNb; nIdx--)
  {
    // Assume backward pricing. If dTime is at the boundary
    // of two windows, pick the left window
    if (    IsAfter(dTime, conv.GetStartTime(nIdx)) 
         && IsEqualOrBefore(dTime, conv.GetEndTime(nIdx)) )
      break;

    // For mandatories
    if ( AreTimesEqual(dTime, conv.GetStartTime(nIdx))
        && AreTimesEqual(dTime, conv.GetEndTime(nIdx)) )
        break;
  }
  m_nIdxConversion = nIdx;
  // if nothing found, at the end of loop, we have just
  // m_nIdx = INVALIDINDEX;

  // set the next coupon index, coupon index is used to update m_dSlopeCoupon
  // and for make whole premium
  CashFlows& cashFlows = *(GetCashFlows());
  nNb = cashFlows.GetNbCashFlows();

  for (m_nIdxNextCoupon = 0;
          m_nIdxNextCoupon < nNb 
       && !IsAfter(cashFlows.GetTime(m_nIdxNextCoupon), dTime);
       m_nIdxNextCoupon++)
    ;

  if ( m_nIdxNextCoupon == 0 || m_nIdxNextCoupon >= nNb ) // no next coupon
    m_dSlopeCoupon = 0;
  else
    m_dSlopeCoupon = cashFlows.GetAmount(m_nIdxNextCoupon)
                   / (  cashFlows.GetTime(m_nIdxNextCoupon) 
                      - cashFlows.GetTime(m_nIdxNextCoupon - 1) );
   
  // init current make-whole premium
  if ( calls.HasMakeWhole() )
  {
    m_dMakeWholeCurrentPremium = m_dMakeWholeFinalPremium;
    for (nIdx = m_nIdxNextCoupon; 
         nIdx < m_nIndexAfterLastCouponMakeWhole; 
         nIdx++)
      m_dMakeWholeCurrentPremium 
          += cashFlows.GetAmount(nIdx) *
             (
                calls.ShouldDiscountCouponForMakeWhole() 
              ? m_cbLike.GetDerivativeCurve()->GetForwardDiscountFactor
                      (m_dCurrentTime, cashFlows.GetTime(nIdx))
              : 1
             );
  } 

  // initialize the claim.
  UpdateClaim();

  if ( m_cbLike.IsCrossCurrency() )
    UpdateFXRate();

  UpdateAccruedInterest();

  // Initialize the conversion price
  UpdateConversionPrice();

  m_dFormerTime = dTime;
}

void CBLikeParams::GetInitialValues
     (const double* pdS, size_t nNbS, double *pdValues)
{
  ASSERT_MSG(m_cbLike.GetPayoff(), "payoff has not been set for contract!");

  m_cbLike.GetPayoff()->Get(pdS, pdValues, nNbS);
}

void CBLikeParams::Update(double dTime)
{
  StateValuesMustBeReset();

  Params::Update(dTime);

  m_CBEventManager.SetCurrentTime(dTime);
  
  UpdateDiscreteCBEventData();

  // according to the way the mesh is build
  // the index are only updated at the end of grid

  // Case where there is a change of grid and thus Update is called twice.
  // The second call will be at the end of a subgrid.   
  if ( AreTimesEqual(dTime, m_dFormerTime) )
  { 
    m_dDerivativeYieldCuveFDF = 1;

    // We update the continuous index for the new subgrid.
    size_t nIdx, nNb;

    // Update call index m_nIdxCall
    // the index is changed only in the following two cases:
    // - end of a call period, therefore m_nIdxCall is set to be its index
    // - end of a period without call, therefore m_nIdxCall
    //                                 is set to be invalid
    // otherwise, we leave m_nIdxCall unchanged.
    CallProvisions& calls = *(GetCalls());
    nNb = calls.GetNbCalls();
    for (nIdx = 0; nIdx < nNb; nIdx++)
    {
      if (    AreTimesEqual(calls.GetEndTime(nIdx), dTime)
           && calls.IsContinuous(nIdx) )
      {
        m_nIdxCall = nIdx;
        break;
      }
      else if ( AreTimesEqual(calls.GetStartTime(nIdx), dTime) )
      { 
        m_nIdxCall = INVALIDINDEX;
        break;
      }
    }

    // Update conversion index m_nIdxConversion
    // see the explanation for updating index of call
    ConversionProvisions& conv = *(GetConversions());
    nNb = conv.GetNbConversions();
    for (nIdx = 0; nIdx < nNb; nIdx++)
    {
      if (    AreTimesEqual(conv.GetEndTime(nIdx), dTime)
           && conv.IsContinuous(nIdx) )
      {
        m_nIdxConversion = nIdx;
        break;
      }
      else if ( AreTimesEqual(conv.GetStartTime(nIdx), dTime) )
      {
        m_nIdxConversion = INVALIDINDEX;
        break;
      }
    }

    // update m_dSlopeCoupon    
    CashFlows& cashFlows = *( GetCashFlows() );
    if (    m_nIdxNextCoupon == 0 
          || m_nIdxNextCoupon >= cashFlows.GetNbCashFlows() )
      m_dSlopeCoupon = 0;
    else
      m_dSlopeCoupon = cashFlows.GetAmount(m_nIdxNextCoupon)
                     / (  cashFlows.GetTime(m_nIdxNextCoupon)
                        - cashFlows.GetTime(m_nIdxNextCoupon - 1) );

    // We update the some initial values for the new subgrid.
    if (    GetCalls()->HasMakeWhole()
        && m_nIdxNextCoupon < m_nIndexAfterLastCouponMakeWhole )
    {     
      // this is not done in UpdateIndex() as that function works 
      // at beginning of each grid. the constraints updated after
      // UpdateIndex() depend on MakeWholePremium at t+.
      // 
      // here, we calculate MakeWholePremium at t-.
      if ( AreTimesEqual(dTime, GetCashFlows()->GetTime(m_nIdxNextCoupon)) )
        m_dMakeWholeCurrentPremium
          += GetCashFlows()->GetAmount(m_nIdxNextCoupon);
    }

    // update the claim at t-
    UpdateClaim(false);
  }   
  else // update inside a subgrid
  {
    m_dDerivativeYieldCuveFDF
       = m_cbLike.GetDerivativeCurve()->GetForwardDiscountFactor
                                                    (dTime, m_dFormerTime);
    if (    GetCalls()->HasMakeWhole()
        && m_nIdxNextCoupon < m_nIndexAfterLastCouponMakeWhole )
    {
      if ( GetCalls()->ShouldDiscountCouponForMakeWhole() )
        m_dMakeWholeCurrentPremium *= m_dDerivativeYieldCuveFDF;
    }

    UpdateClaim(true); 
  }

  UpdateAccruedInterest();

  if ( m_cbLike.IsCrossCurrency() )
    UpdateFXRate();

  // Common updates
  UpdateConversionPrice();

  m_dFormerTime = dTime;
}

bool CBLikeParams::UpdateMonoDateEventIndex()
{
  ASSERT_MSG(m_nIdxPut == INVALIDINDEX,
             "Index of put should be invalid in the PDE solving processus");

  bool bConstraintsNeedUpdate = false;

  const CBEvent* pCBEvent;

  while ( (pCBEvent = GetCBEventAtCurrentTime()) != 0 )
  {
    MonoDateEventIndexTreatment(pCBEvent, bConstraintsNeedUpdate);
  }

  return bConstraintsNeedUpdate;
}

void CBLikeParams::MonoDateEventIndexTreatment(const CBEvent* pCBEvent, 
                                               bool& bConstraintsNeedUpdate)
{
  switch (pCBEvent->m_cbeventType)
  {
    case CBET_MonoDatePut:
      m_nIdxPut = pCBEvent->m_nIndex;
      bConstraintsNeedUpdate = true;
      break;
    case CBET_MonoDateCall:
      m_nIdxCall = pCBEvent->m_nIndex;
      bConstraintsNeedUpdate = true;
      break;
    case CBET_MonoDateConversion:
      m_nIdxConversion = pCBEvent->m_nIndex;

      UpdateConversionPrice();

      bConstraintsNeedUpdate = true;
      break;        
    case CBET_PayCoupon:           
    {
      m_nIdxNextCoupon = pCBEvent->m_nIndex;

      CashFlows& cashFlows = *(GetCashFlows());

      if (    m_nIdxNextCoupon > 0 
            && m_nIdxNextCoupon < cashFlows.GetNbCashFlows() )
        m_dCouponAmount = cashFlows.GetAmount(m_nIdxNextCoupon);

      bConstraintsNeedUpdate = true;
      break;
    }
    case CBET_PayDividend:
    {
      if ( HasNewShare() )
        bConstraintsNeedUpdate = true;
      break;
    }
    case CBET_StartOfYear:
    {
      if ( HasNewShare() )
        m_bIsStartOfYear = true;
      break;
    }
    // This is necessary to treat correctly the maturity case.
    // See the comment at the begining of the function 
    // CBInstData::SetInitialValue() (in the cbinstdata.cpp file)
    case CBET_Maturity:
      bConstraintsNeedUpdate = true;
      break;
  }
}

void CBLikeParams::GetSpecialTimes(SpecialTimes& specialTimes) const
{
  SpecialTimes::iterator Result;

  // call the function in the base class
  Params::GetSpecialTimes(specialTimes);

  SpecialTimes specialTimesCB;

  // Add CBEvent times to the special times
  m_CBEventManager.GetSpecialTimes(specialTimesCB);
   
  SpecialTimes::const_iterator iter;
  for (iter = specialTimesCB.begin(); iter != specialTimesCB.end(); ++iter)
    specialTimes.push_back(*iter);

  // In case of call notice, add all (special_time - notice period)
  // times into the mesh.
  CallProvisions* pCalls = GetCalls();
  SpecialTimes callNoticePoints;
  if ( pCalls->HasNoticePeriod() )
  {
    double dNoticePeriod = pCalls->GetNoticePeriod();
    size_t nNbCalls = pCalls->GetNbCalls();
    
    for (size_t nIdxCall = 0; nIdxCall < nNbCalls; nIdxCall++)
    {
      double dStartTime = pCalls->GetStartTime(nIdxCall);
      double dEndTime = pCalls->GetEndTime(nIdxCall);

      // TODO: improve the efficiency
      SpecialTimes::const_iterator iter;
      for (iter = specialTimes.begin();
           iter != specialTimes.end();
           ++iter)
      {
        double dTime = iter->GetTime();
        double dTimeWithNotice = dTime - dNoticePeriod;
        if ( IsAfter(dTimeWithNotice, dStartTime)
          && IsBefore(dTimeWithNotice, dEndTime) )
        {
          callNoticePoints.push_back
              ( SpecialTime(dTimeWithNotice, iter->GetRefineLevel()) );
        }
      }  // loop through the existing special times
      
    } // loop over call windows

    // Add the new points
    SpecialTimes::const_iterator iterNew;
    for (iterNew = callNoticePoints.begin();
         iterNew != callNoticePoints.end();
         ++iterNew)
      specialTimes.push_back( *iterNew );
  } // if call notice


  // We suppose here that the obtained special times are always 
  // not bigger than the maturity of the contract  
  SpecialTime valuationTime(m_dValuationTime);
  Result = std::remove_if
           (specialTimes.begin(), specialTimes.end(),
            std::bind2nd(std::less<SpecialTime>(), valuationTime));

  //Use erase to remove no needed elements.  
  specialTimes.erase( Result, specialTimes.end() );

  // Add the valuation time 
  specialTimes.push_back(valuationTime);
}

void CBLikeParams::GetRecoveryValues
                   (const double *pdS, size_t nNbS, double *pdRecoveryValues)
{
  double dRecoveryRate = m_cbLike.GetRecoveryRate();
  double dValue = m_dClaim * dRecoveryRate;

  size_t n;
  for ( n = 0; n < nNbS; n++)
    pdRecoveryValues[n] = dValue;

  bool bConversionOnDefault;
  if ( m_cbLike.IsExchangeable(bConversionOnDefault) )
  {
    std::vector<double> pdConv(nNbS);

    // Remark: We pass, in the third argument of the function GetGrossParities,
    // the spot array because in this case (ie: Exchangeable), it is not 
    // possible to have new share feature.
    GetConversions()->GetGrossParities(pdS, nNbS, pdS, &pdConv[0]);

    if (bConversionOnDefault)
    {
      // max(Re * claim, rs)
      for (n = 0; n < nNbS; n++)
        pdRecoveryValues[n] = std::max(pdRecoveryValues[n], pdConv[n]);
    }
    else
    {
      // Re * max(claim, rs) =  max(Re * claim, Re * rs)
      for (n = 0; n < nNbS; n++)
        pdRecoveryValues[n] =
          std::max(pdRecoveryValues[n], dRecoveryRate * pdConv[n]);
    }
  }
}

bool CBLikeParams::GetCallConstraintValues
      (const double* pdS, size_t nNbS, double* pdValues,
       size_t& nIdxStartConversion, 
       const double* pdNewSharePrices, 
       bool bHasNoticePeriod)
{
//  if ( m_nIdxCall == INVALIDINDEX )
//    return false;
  
  return GetCalls()->GetCallConstraintValues
                     (pdS, nNbS, pdValues, nIdxStartConversion, 
                      pdNewSharePrices, bHasNoticePeriod);

//  return true;
}

bool CBLikeParams::GetConversionConstraintValues
                   (const double* pdS, size_t nNbS, 
                    const double* pdNewSharePrices, 
                    double* pdValues)
{
//  if ( m_nIdxConversion == INVALIDINDEX )
//    return false;

  return GetConversions()->GetConversionConstraintValues
                           ( pdS, nNbS, pdNewSharePrices, 
                             pdValues );
//  return true;
}

bool CBLikeParams::GetPutConstraintValues
                   (const double* pdS, size_t nNbS, double* pdValues)
{
  return GetPuts()->GetPutConstraintValues(pdS, nNbS, pdValues);
}

CBLikeParams* CBLikeParams::GetCallNoticeParams()
{
  // need to solve from t+NoticePeriod back to t
  double dNoticePeriod       =  GetCalls()->GetNoticePeriod();
  double valuationTimeNotice = m_dCurrentTime;
  double maturityTimeNotice  = valuationTimeNotice + dNoticePeriod;
  
  // ensure that we are not going beyond the maturity date
  ASSERT_MSG( maturityTimeNotice <= m_cbLike.GetMaturityTime(),
    "Call date with notice period is beyond the maturity date.");
  
  // Set the maturity time of the call notice CB contract
  m_pCallNoticeParams->GetCBLike().SetMaturityTime(maturityTimeNotice);
  m_pCallNoticeParams->SetStoppingTime(maturityTimeNotice);

  // Set the valuation time of the call notice params
  m_pCallNoticeParams->SetValuationTime(valuationTimeNotice);

  // Set the spot fx rate for the call notice CB contract
  m_pCallNoticeParams->GetCBLike().SetSpotFXRate(GetFXRate(m_dCurrentTime));

  return m_pCallNoticeParams.get();
}

void CBLikeParams::SetupCallNoticeParams()
{
  // Use the clone function to generate the call notice params.
  m_pCallNoticeParams = AutoPtr<CBLikeParams>( Clone() );

  /*
    When the issuer calls back, the holder of the security can convert during 
    the notice period, regardless the conversion trigger. To be confirmed.
  */
  m_pCallNoticeParams->GetConversions()->RemoveTriggers();

  /*
    By convention, when the issuer calls for early redemption, the holder of
    the security can't put the security.
  */
  m_pCallNoticeParams->GetPuts()->Clear();

  // Remove the call constraints of the mini cb
  m_pCallNoticeParams->GetCalls()->DeactivateAllCalls();


  // Prepare the NumParams for the call notice mini CB problem 
  shared_ptr<NumParams> pNumParams( new NumParams(*m_pNumParams) );
              
  /*
    Use FEW time points for the mini cb probleme, as the problem
    convergences very quickly on time. 
    
    CAN'T use fewer grid points for the mini cb problem
    as that decreases directly the global precison. For example, if
    we use half of number of grid points here, the precision of the 
    mini cb problem is at least half of that of the "big" problem
    without call notice. Thus the precision of "big" cb is destroied.
  */
  
  /*
    As the time mesh is refined at maturity of mini cb problem and
    that the call notice period is not long (at most 60), even we set
    number of time points as 3, we would get more about 10 points. That is
    enough to have very good precision.

    The precision of the whole problem with respected to time 
    is limited by some value as function of MaturityOfBigCB/InitialNbTimes.
    So setting 1000 time points for mini cb problem will multiply the
    run time by 100 but won't gain any precision.

    That works even when the maturity of the big cb is short. In fact,
    we prefer diminuating run time by 10
    than gaining one more precision after 4th decimal digit.
  */
  const size_t nNbTimeStepsMinForCallNotice = 3;

  /* 
    The following method is more precise, but much more time consuming
    size_t nNbTimeStepsMinForCallNotice = 
    m_pNumParams->GetNbTimeSteps() / (GetStoppingTime() - GetValuationTime())
                                   * 60. / 365;
    if (nNbTimeStepsMinForCallNotice < 3)
      nNbTimeStepsMinForCallNotice = 3;
  */

  pNumParams->SetNbTimeSteps(nNbTimeStepsMinForCallNotice);

  m_pCallNoticeParams->SetNumParams(pNumParams); 
}


void CBLikeParams::UpdateClaim(bool bPlus)
{
  // Update the claim value
  m_dClaim = m_cbLike.GetClaim(m_dCurrentTime, bPlus);

  #ifndef NDEBUG
    m_bClaimUpdated = true;
  #endif
}

double CBLikeParams::GetClaimFromYield(double dYield, bool bPlus) const
{
  return GetClaimFromYield(m_dCurrentTime,dYield,bPlus);
}//double CBLikeParams::GetClaimFromYield(double dYield, bool bPlus)

double 
CBLikeParams::GetClaimFromYield(double dTime, double dYield,bool bPlus) const
{
  double dIssueTime     = m_cbLike.GetIssueTime();
  double dCmpFreq       = GetCBLike().GetCompoundingFrequency();
  double dDiscount      = 1. + dYield / dCmpFreq;
  double dClaim         = GetCBLike().GetIssuePrice() 
                        * pow(dDiscount, dCmpFreq * (dTime - dIssueTime) );

  const CashFlows *pCashFlows =  GetCBLike().GetCashFlows();

  size_t nNbCashFlows = pCashFlows->GetNbCashFlows();
  size_t nIdxCashFlow;

  // Unlike the computation of the claim for recovery value,
  // this is computed by using issue price, which has no problem 
  // for call notice since the coupons after maturity will be anyway ignored.

  if (bPlus) // coupon at dTime (if any) should also be discounted
    for (nIdxCashFlow = 0; 
            nIdxCashFlow < nNbCashFlows  
        && numeric::IsEqualOrBefore(pCashFlows->GetTime(nIdxCashFlow),dTime); 
         nIdxCashFlow++)
      dClaim -= pCashFlows->GetAmount(nIdxCashFlow) 
              * pow( dDiscount, 
                     dCmpFreq * ( dTime - pCashFlows->GetTime(nIdxCashFlow) )
                   );
  else // don't discount the coupon at dTime (if any)
    for (nIdxCashFlow = 0; 
            nIdxCashFlow < nNbCashFlows  
            && numeric::IsBefore(pCashFlows->GetTime(nIdxCashFlow), dTime); 
         nIdxCashFlow++) 
      dClaim -= pCashFlows->GetAmount(nIdxCashFlow) 
              *  pow( dDiscount, 
                      dCmpFreq * ( dTime - pCashFlows->GetTime(nIdxCashFlow) ) 
                    );

  return dClaim;

}//CBLikeParams::GetClaimFromYield(double dTime, double dYield, bool bPlus)

bool CBLikeParams::HasPathDepCoCo() const
{
  return GetConversions()->HasPathDepCoCo(m_dValuationTime, m_dStoppingTime);
}

bool CBLikeParams::HasPathDepCall() const
{
  return GetCalls()->HasPathDepCall(m_dValuationTime, m_dStoppingTime);
}

AutoPtr<CB> CBLikeParams::GetNewShareCB()
{     
  Date TmpDate = m_cbLike.GetFiscalYearStartDate();

  if ( IsEqualOrBefore(GetDoubleFrom(TmpDate), m_dStoppingTime) )
  {
    while ( IsEqualOrBefore(GetDoubleFrom(TmpDate), m_dStoppingTime) )
      TmpDate.AddYears(1);
  }
  else
  {
    while ( IsAfter(GetDoubleFrom(TmpDate), m_dStoppingTime) )
      TmpDate.AddYears(-1);

    TmpDate.AddYears(1);
  }  
  //Then, TmpDate is the first fiscal year after the maturity time

  AutoPtr<CB> cb( new CB() ); 

  cb->SetMaturityTime( GetDoubleFrom(TmpDate) );

  shared_ptr<CashFlows> pCashFlows( new CashFlows() );
  cb->SetCashFlows(pCashFlows);
  
  cb->SetRedemptionValue(0.);
  cb->SetDerivativeCurve( m_cbLike.GetDerivativeCurve() );

  return cb;
}

shared_ptr<CBLikeParams> 
CBLikeParams::GetNewShareParams(const double *pdSpots, const size_t nNbSpots)
{  
  ASSERT_MSG(HasNewShare(), "calling GetNewShareParams for a non new share cb");

  // Construct a particular cb
  AutoPtr<CB> cb = GetNewShareCB();
      
  // Construct the particular cb like params
  shared_ptr<CBLikeParams> pNewShareParams( new CBParams(cb) );

  pNewShareParams->SetSpotSharePrice( GetSpotSharePrice() );
  pNewShareParams->SetValuationTime(m_dStoppingTime);
  
  pNewShareParams->SetYieldCurve( GetYieldCurve() );
  pNewShareParams->SetForeignCurve( GetForeignCurve() );
  
  pNewShareParams->SetDividends(m_pDividends);
  
  // We have to set the payoff manually.    
  pNewShareParams->GetCBLike().SetPayoff
  (shared_ptr<finance::Payoff>(
    new finance::PayoffDiscrete(pdSpots, pdSpots, nNbSpots)) );

  return pNewShareParams;
}

} // namespace pricing

} // namespace ito33
