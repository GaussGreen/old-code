/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/cboptionparams.cpp
// Purpose:     cb option params class
// Author:      Nabil
// Created:     2005/10/13
// RCS-ID:      $Id: cboptionparams.cpp,v 1.7 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2004-2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/payoffconstant.h"

#include "ito33/pricing/dividendevents.h"
#include "ito33/pricing/eventmanager.h"
#include "ito33/pricing/cboption.h"
#include "ito33/pricing/cboptionparams.h"
#include "ito33/pricing/cbevent.h"

#include "ito33/numeric/predicatetime.h"

using namespace ito33;
using namespace ito33::pricing;
using namespace ito33::numeric;

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(pricing::CBOptionParams);
}

void CBOptionParams::Init()
{
  CBParams::Init();
  
  // CB option event

  CBOptionData& cbOptionData = *(m_cboption.GetCBOptionData());

  // maturity of the cb option
  double dCbOptMaturityTime = cbOptionData.GetMaturityTime();

  if ( !IsBefore(dCbOptMaturityTime, m_dValuationTime)
       && IsEqualOrBefore(dCbOptMaturityTime, m_dStoppingTime) )
  {
    m_CBEventManager.AddEvent
      ( make_ptr( new CBEvent(dCbOptMaturityTime, CBET_CBOptionMaturity) ) );
  }

  // floating payment events
  size_t nIdx;

  CashFlows& floatingcashFlows = *(cbOptionData.GetFloatingCashFlows());
  size_t nNbCashFlows = floatingcashFlows.GetNbCashFlows();

  for (nIdx = 0; nIdx < nNbCashFlows; nIdx++)
  {
    double dTime = floatingcashFlows.GetTime(nIdx);

    // floating coupon at maturity will be taken into account
    if (   !IsBefore(dTime, m_dValuationTime) 
        && !IsAfter(dTime, dCbOptMaturityTime) )
    { 
      m_CBEventManager.AddEvent
        ( make_ptr( new CBEvent(dTime, CBET_CBOptionFloatingCoupon, nIdx) ) );
    }
  }

  // m_bInCBOptionWindow must be initialized to false
  m_bInCBOptionWindow = false;
}

void CBOptionParams::SetInitialState(double dTime)
{  
  CBParams::SetInitialState(dTime);
  
  // maturity date >= CBOptionMaturity, so
  m_nIdxASWNextFloatingCoupon
      = GetCBOptionData()->GetFloatingCashFlows()->GetNbCashFlows();

  // Initialize the cb option data
  m_bInCBOptionWindow = false;

  m_dASWSlopeFloatingCoupon = 0.;

  UpdateCBOptionData(dTime, true);
}

void CBOptionParams::Update(double dTime)
{
  CBParams::Update(dTime);

  // Update cb option data
  UpdateCBOptionData();
}

void CBOptionParams::UpdateCBOptionData()
{
  if ( !InCBOptionWindow() )
    return;

  double
    dNPVBalloonCoupon,
    dFloatingAccrued;

  CBOptionData& cbOptionData = *(m_cboption.GetCBOptionData());

  m_dASWNPVDiff *= m_dDerivativeYieldCuveFDF;

  if ( m_dASWSlopeFloatingCoupon != 0)
    dFloatingAccrued = m_dASWSlopeFloatingCoupon * 
                      (m_dCurrentTime - 
                       cbOptionData.GetFloatingCashFlows()
                            ->GetTime(m_nIdxASWNextFloatingCoupon - 1));
  else
    dFloatingAccrued = 0.;

  dNPVBalloonCoupon = cbOptionData.GetBalloonCoupon()
                    * m_cb.GetDerivativeCurve()->GetForwardDiscountFactor( 
                        m_dCurrentTime, cbOptionData.GetMaturityTime() );

  // The strike of the cb option
  m_dCBOptionStrike = cbOptionData.GetCbOptionFaceValue()
                    + m_dASWNPVDiff
                    + dNPVBalloonCoupon
                    + dFloatingAccrued;
}


void CBOptionParams::UpdateCBOptionData(double dTime, bool bPlus)
{
  // todo: check if comparison function is needed!
  if ( GetCBOptionData()->GetMaturityTime() >= dTime )
  {
    // Updates the strike of the cb option
    double
      dFloatingAccrued,
      dNPVBalloonCoupon,
      dNPVFixedLeg,
      dNPVFloatingLeg;

    CBOptionData& cbOptionData = *(m_cboption.GetCBOptionData());

    // Fixed leg
    if ( GetCBOptionData()->GetFixedCashFlows() )
    {
      CashFlows& fixedCashFlows = *(GetCBOptionData()->GetFixedCashFlows());
      dNPVFixedLeg = ComputeNPV( fixedCashFlows, dTime, m_nIdxNextCoupon);
    }
    else
      dNPVFixedLeg = 0.;

    // Floating leg
    CashFlows& 
      floatingCashFlows = *(GetCBOptionData()->GetFloatingCashFlows());

    dNPVFloatingLeg = ComputeNPV( floatingCashFlows, dTime, 
                                  m_nIdxASWNextFloatingCoupon );
    
    dFloatingAccrued = floatingCashFlows.GetAccruedInterest(dTime, bPlus);

    dNPVBalloonCoupon = cbOptionData.GetBalloonCoupon()
                      * m_cb.GetDerivativeCurve()->GetForwardDiscountFactor( 
                          dTime, cbOptionData.GetMaturityTime() );

    m_dASWNPVDiff = dNPVFixedLeg - dNPVFloatingLeg;

    // The strike of the cb option
    m_dCBOptionStrike = cbOptionData.GetCbOptionFaceValue() 
                      + m_dASWNPVDiff 
                      + dNPVBalloonCoupon
                      + dFloatingAccrued;
  }
}

double CBOptionParams::ComputeNPV
       (const CashFlows& cashflows, double dTime, size_t nNextCashFlow) const
{  
  double dNPV = 0.;

  if ( nNextCashFlow  < cashflows.GetNbCashFlows() )
  {
    ASSERT_MSG( cashflows.GetTime(nNextCashFlow) >= dTime,
                "Next cash flow date must be after the current time." );
      
    size_t 
      nIdx,
      nNbCashFlows = cashflows.GetNbCashFlows();
    
    for (nIdx = nNextCashFlow; nIdx < nNbCashFlows; ++nIdx)
      dNPV += cashflows.GetAmount(nIdx) * m_cb.GetDerivativeCurve()
              ->GetForwardDiscountFactor( dTime, cashflows.GetTime(nIdx) );
  }

  return dNPV;
}

void CBOptionParams::MonoDateEventIndexTreatment(const CBEvent* pCBEvent, 
                                                 bool& bConstraintsNeedUpdate)
{
  CBParams::MonoDateEventIndexTreatment(pCBEvent, bConstraintsNeedUpdate);

  // Treatment for cb option
  switch ( pCBEvent->m_cbeventType )
  {
    // I update m_dCBOptionStrike here by playing with fixed/floating coupon.
    // At the end of grid, we need to update CBOptionStrike to update the 
    // constraint for CBOption, while UpdateCBOptionData() is only called 
    // "before step".
    // Does it exist a better solution?
    case CBET_CBOptionFloatingCoupon:  
      {
        m_nIdxASWNextFloatingCoupon = pCBEvent->m_nIndex;

        const CashFlows* 
          pFloatingCashFlows = 
            GetCBOptionData()->GetFloatingCashFlows().get();
        
        if ( m_nIdxASWNextFloatingCoupon > 0  && 
             m_nIdxASWNextFloatingCoupon < pFloatingCashFlows->GetNbCashFlows() 
           )
        {
          double dCoupon = 
            pFloatingCashFlows->GetAmount(m_nIdxASWNextFloatingCoupon);
          
          m_dASWSlopeFloatingCoupon
            = dCoupon
            / 
            ( pFloatingCashFlows->GetTime(m_nIdxASWNextFloatingCoupon)
              - pFloatingCashFlows->GetTime(m_nIdxASWNextFloatingCoupon - 1) 
            );
        
          m_dASWNPVDiff -= dCoupon;
        }
        else
          m_dASWSlopeFloatingCoupon = 0.;
      }
      break;
    // This is necessary to treat correctly the maturity case.
    // The update must be done after all the other events at the same time
    // which modify the cb option data. This is the case if we doesn't
    // change the order of the events in the enum CBEventType.
    case CBET_CBOptionMaturity:
      m_bInCBOptionWindow = true;
      UpdateCBOptionData(pCBEvent->GetTime(), false);
      break;
    case CBET_PayCoupon: // add coupon to ASW_NPV in CBOptionWindow
      // I don't have better solution than do the following check as 
      // m_bInCBOptionWindow can still be false here
      if ( IsBefore( GetCashFlows()->GetTime(m_nIdxNextCoupon),
                            GetCBOptionData()->GetMaturityTime() ) )
      {
        m_dASWNPVDiff += m_dCouponAmount;
        m_dCBOptionStrike += m_dCouponAmount;
      }
      break;
  }
}
