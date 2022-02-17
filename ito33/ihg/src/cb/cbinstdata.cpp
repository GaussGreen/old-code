/////////////////////////////////////////////////////////////////////////////
// Name:        cb/cbinstdata.cpp
// Purpose:     implementations for cb instdata class in the ihg model
// Author:      Nabil
// Created:     2004/03/30
// RCS-ID:      $Id: cbinstdata.cpp,v 1.119 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#if ( defined(_MSC_VER) && _MSC_VER >= 1400 )
// still want to use fill_n but don't know how to avoid the compile warning
// warning C4996: 'std::_Fill_n' was declared deprecated
// so just disable the warning
#pragma warning(disable:4996)
#endif

#include "ito33/sharedptr.h"
#include "ito33/array.h"

#include "ito33/finance/payoffdiscrete.h"
#include "ito33/finance/computationalflags.h"

#include "ito33/numeric/interpolation.h"

#include "ito33/pricing/event.h"
#include "ito33/pricing/cbevent.h"
#include "ito33/pricing/cbconstraints.h"
#include "ito33/pricing/cbmeshmanager.h"
#include "ito33/pricing/resetparams.h"
#include "ito33/pricing/attachedwarrantcbparams.h"

#include "ihg/cbinstdata.h"
#include "ihg/model.h"
#include "ihg/cbpricer.h"
#include "ihg/resetpricer.h"
#include "ihg/attachedwarrantcbpricer.h"
#include "ihg/cbnumoutput.h"

namespace ito33
{

using namespace numeric;
using namespace pricing;

namespace ihg
{

CBInstData::CBInstData( pricing::CBLikeParams& cbparams, 
                        Model& model, 
                        pricing::CBMeshManager& cbmeshes)
                      : InstDataMulti(cbparams, model, cbmeshes),
                        m_cbparams(cbparams), m_cbmeshes(cbmeshes)
{
  m_bHasNewShare = m_cbparams.HasNewShare();
}

void CBInstData::Alloc(size_t nNbS)
{
  InstDataMulti::Alloc(nNbS);

  // Allocate the new share arrays
  if( m_bHasNewShare )
  {
    m_pdNewSharePrices = Array<double>(nNbS);
    m_pdOldNewSharePrices = Array<double>(nNbS);
    m_pdOldOldNewSharePrices = Array<double>(nNbS); 
  }

  m_pdRecoveryValues.resize(nNbS);
}

void CBInstData::ApplyEvent(const pricing::Event *pEvent)
{
  InstDataMulti::ApplyEvent(pEvent);

  if( m_bHasNewShare && pEvent->GetType() == ET_Dividend )
  {
    pEvent->ApplyToPrice( m_pdS, m_pdNewSharePrices.Get(), m_nNbS );
  }
}

void CBInstData::Swap()
{
  InstDataMulti::Swap();

  if (m_bHasNewShare)
  {
    swap(m_pdOldOldNewSharePrices, m_pdNewSharePrices);
    swap(m_pdOldOldNewSharePrices, m_pdOldNewSharePrices);
  }
}

void CBInstData::InterpWithPassageOfSpaceMesh
                 (const double *pdOldSpots, size_t nNbOldSpots, 
                  const double *pdSpots, size_t nNbSpots)
{
  InstDataMulti::InterpWithPassageOfSpaceMesh
                 ( pdOldSpots, nNbOldSpots, pdSpots, nNbSpots );
  
  // Interpolation of new share prices vector on the new grid.
  if( m_bHasNewShare )
    Interpolate(pdOldSpots, m_pdOldNewSharePrices.Get(), nNbOldSpots, 
                pdSpots, m_pdNewSharePrices.Get(), nNbSpots, 
                ExtrapolationMode_Linear, ExtrapolationMode_Linear);
}

void CBInstData::Init()
{
  m_nNbSpotsMax   = m_cbmeshes.GetNbSMax();
 
  Alloc(m_nNbSpotsMax);
  
  // Some initializations.
  m_pdS    = m_cbmeshes.GetS();
  m_pdLogS = m_cbmeshes.GetLogS();  

  bool bTmp;
  m_bDerivativeHasOwnHR = m_cbparams.GetCBLike().IsExchangeable(bTmp);
}

/*
  In this function, we treat the maturity case. At this special time, the price 
  must be computed as follow:
    1. We initialize the price to the Redemption value.
    2. We apply the events (coupon, dividend...) to this price.
    3. We apply constraints (Call, conversion...) to this new price.
*/
void CBInstData::SetInitialValue()
{
  m_nNbS = m_cbmeshes.GetNbS();
     
  if ( m_pPayoff )
    m_pPayoff->Get(m_pdS, m_pdPrices.Get(), m_nNbS);
  else
    m_cbparams.GetInitialValues(m_pdS, m_nNbS, m_pdPrices.Get());

  std::fill_n(m_piFrozenFlags.Get(), m_nNbS, 0); 
  
  if ( m_bComputeVega )
    std::fill_n(m_pdVegas.Get(), m_nNbS, 0.);
  
  if ( m_bComputeFugit )
    std::fill_n(m_pdFugits.Get(), m_nNbS, 0.);

  if ( m_bHasNewShare )
  {
    ComputeNewSharePricesAtMaturity( m_pdS, m_nNbS, m_pdNewSharePrices.Get() );
  }
  
  DoEvents();
}

double CBInstData::GetSpeedCorrection() const
{
  return m_cbmeshes.GetSpeedCorrection();
}

void CBInstData::UpdateBeforeStep()
{
  
  if ( m_cbmeshes.IsEndOfGrid() )
  {
    /*
      Get initial condition for the new subgrid from the followed subgrid.
      Time index is not changed, and is not the one at which the system will 
      be solved at this step.
    */
    UpdateAtEndOfGrid(); 

    /*
      Really go ahead, and change the time index, so constraints need to be 
      updated again.
    */
    m_cbmeshes.TryGoAhead();
  }

  InstDataMulti::UpdateBeforeStep();

  if ( m_cbparams.GetCBLike().IsCrossCurrency() )
    m_cbmeshes.UpdateFXRate();

  if (m_bDerivativeHasOwnHR)
    m_dHROfDerivative = m_model.GetHROfDerivative();

  m_cbparams.GetRecoveryValues(m_pdS, m_nNbS, &m_pdRecoveryValues[0]); 
  
  /*
    Update the constraints, the result will be used to solve the non linear 
    system at current time step.
  */
  UpdateConstraints();
}

double CBInstData::GetValueAfterDefaultAtInitialSpot()
{
  double dResult;
  double dInitialSpot = GetInitialSpot();
  Interpolate(m_pdS, &m_pdRecoveryValues[0], m_nNbS,
              &dInitialSpot, &dResult, 1, 
              ExtrapolationMode_Linear, ExtrapolationMode_Linear);
  return dResult;
}

void CBInstData::UpdateConstraints()
{
  Array<double> pdValues(m_nNbS);
  
  bool bHasConstraints = false;

  // If the calls have a notice period, then a mini cb problem needs to be
  // solved.  The mini cb is model specific, so is solved here.  Without
  // call notice, the params class can determine the constraint
  if (    m_cbparams.GetCalls()->IsActive() 
       && m_cbparams.GetIndexCall() != INVALIDINDEX )
  {    
    bool bHasNoticePeriod = m_cbparams.GetCalls()->HasNoticePeriod();

    if (bHasNoticePeriod)
    {
      if( m_bHasNewShare )
        SolveCallNotice(m_pdS, m_nNbS, m_pdNewSharePrices.Get(), 
                        pdValues.Get());
      else
        SolveCallNotice(m_pdS, m_nNbS, m_pdS, pdValues.Get());
    }

    size_t nIdxStartConversion;
    
    if( m_bHasNewShare )
      m_cbparams.GetCallConstraintValues
                 (  m_pdS, m_nNbS, pdValues.Get(), 
                    nIdxStartConversion, 
                    m_pdNewSharePrices.Get(), 
                    bHasNoticePeriod
                 );
    else
      m_cbparams.GetCallConstraintValues(m_pdS, m_nNbS, pdValues.Get(), 
                                         nIdxStartConversion, 
                                         m_pdS, bHasNoticePeriod);

    m_constraints.UpdateCall(pdValues.Get(), m_nNbS, nIdxStartConversion);

    bHasConstraints = true;
  }
  else
    m_constraints.DisableTheCall();

  // Nothing special about the conversions. Let params determine the constraint
  bool bUpdateConv;

  if( m_bHasNewShare )
    bUpdateConv = m_cbparams.GetConversionConstraintValues
                             ( m_pdS, m_nNbS, 
                               m_pdNewSharePrices.Get(),
                               pdValues.Get() );
  else
    bUpdateConv = m_cbparams.GetConversionConstraintValues
                             ( m_pdS, m_nNbS, m_pdS, pdValues.Get() );
  
  if ( bUpdateConv )
  {
    m_constraints.UpdateConv(pdValues.Get(), m_nNbS);
    bHasConstraints = true;
  }
  else
    m_constraints.DisableTheConv();

  // Nothing special about the puts. Let params determine the constraint
  if ( m_cbparams.GetPutConstraintValues(m_pdS, m_nNbS, pdValues.Get()) )
  {
    m_constraints.UpdatePut(pdValues.Get(), m_nNbS);
    bHasConstraints = true;
  }
  else
    m_constraints.DisableThePut();

  if (bHasConstraints)
    m_pConstraints = &m_constraints;
  else
    m_pConstraints = 0;
}

// --------------------------------------------------------------------------
//  Solve a small cb from t+notice period back to t
//  set a vector of value to the call strike in function of S
// --------------------------------------------------------------------------
void CBInstData::SolveCallNotice(const double *pdS, 
                                 const size_t nNbS,
                                 const double* pdNewSharePrices,
                                 double *pdValues)
{

  // Get the mini cb params from the "big" cb params. Assume
  // it is constructed with the correct, current info
  pricing::CBLikeParams
    *pCallNoticeParams = m_cbparams.GetCallNoticeParams();

  // Now that the mini Cb is not constructed from a financial one
  // we need to set the payoff manually or call a function in CB
  // to let it construct itself the payoff objet from its redemption value.
  // Remark: This should be moved from params to here because now, this
  //         depends on model.
  Array<double> pdInitialValues(nNbS);

  // Get the call strikes at the maturity time of the call notice CB
  // and use it as initial value for the mini cb pricer
  m_cbparams.GetCalls()->GetCallStrikesWithoutCoupon
      (pCallNoticeParams->GetCBLike().GetMaturityTime(),
       pdS, nNbS, pdNewSharePrices, pdInitialValues.Get());

  shared_ptr<finance::Payoff> 
      pPayoff( new finance::PayoffDiscrete(pdS, pdInitialValues.Get(), nNbS) );

  // Do not need to compute anything other than price, so use a default flags
  finance::ComputationalFlags flags;

  // but use the same solver as for the original one
  flags.SetSolverType(m_iSolverType);

  // The pricing output
  AutoPtr<CBNumOutput> pCBNumOutput;

  // Check for reset contract
  pricing::ResetParams* pResetCallNoticeParams 
    = dynamic_cast<pricing::ResetParams*>( pCallNoticeParams );
  
  // Also check for attached warrant contract
  pricing::AttachedWarrantConvertibleBondParams* pAWCallNoticeParams 
    = dynamic_cast<pricing::AttachedWarrantConvertibleBondParams*>
                  ( pCallNoticeParams );

  if ( pResetCallNoticeParams != 0 )
  {
    // We are pricing a reset, so must use reset pricer. If a reset date
    // is not active, the cb pricer is used.
    ResetPricer miniResetPricer(*pResetCallNoticeParams, m_model, flags);

    miniResetPricer.SetInitialValue(pPayoff);

    pCBNumOutput = miniResetPricer.Price(); 
  }
  else if ( pAWCallNoticeParams != 0 )
  {
    // We are pricing an attached warrant, so must use attached warrant
    // pricer.
    AttachedWarrantConvertibleBondPricer
      miniAWPricer(*pAWCallNoticeParams, m_model, flags);

    miniAWPricer.SetInitialValue(pPayoff);

    pCBNumOutput = miniAWPricer.Price(); 
  }
  else
  {
    // Can price as a normal CB
    CBPricer minicbPricer(*pCallNoticeParams, m_model, flags);

    minicbPricer.SetInitialValue(pPayoff);

    pCBNumOutput = minicbPricer.Price(); 
  }

  // Interpolate since meshes could be different
  Interpolate(&pCBNumOutput->GetFinalMesh()[0], 
              &pCBNumOutput->GetFinalPrices()[0], 
              pCBNumOutput->GetFinalPrices().size(),
              pdS, pdValues, nNbS, 
              ExtrapolationMode_Linear, ExtrapolationMode_Linear);
}

void CBInstData::ApplyConstraintsAfterEventsOrConstraintsUpdate()
{
  bool bConstraintsNeedUpdate = m_cbparams.UpdateMonoDateEventIndex(); 

  if (bConstraintsNeedUpdate)
  {
    // @todo only constraints(including the continuous call or conversion 
    // constraint when there is a coupon event) that will be changed by mono
    // date events need to be "updated". Might not be worth it for the moment.
    UpdateConstraints();

    m_bHasEvent = true;
  }

  // Constraints updated because of events, apply them
  if ( m_bHasEvent && m_pConstraints )
    ApplyConstraintsToAll();

}

void CBInstData::UpdateAtEndOfGrid()
{ 
  Swap();

  /*
    Maturity for the moment is considered as an end of grid, although 
    it's not necessary to do the interpolation, it's rather an exception.
    "Real" end of grid can have different space mesh at t+ and t-

    The interpolated values(after having applied the possible constraints)
    will be used as the initial condition of the previous subgrid 
    (we go backward on time).
  */
  InterpWithPassageOfSpaceMesh( m_meshes.GetOldS(), m_meshes.GetOldNbS(), 
                                m_meshes.GetS(), m_meshes.GetNbS() );
  
  m_nNbS = m_meshes.GetNbS();

  if( m_bHasNewShare )
  {
    // Updates the new share price across the fiscal year start date (at t-)
    if( m_cbparams.IsStartOfYear() )
    {
      const double* pdS = m_meshes.GetS();
      
      for(size_t nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
        m_pdNewSharePrices[nIdxS] = pdS[nIdxS];

      m_cbparams.DisableStartOfYear();
    }
  }    

  /*
    Update the constraints at t-, it can be different from constraints at t+
    because of keep accrued or forfeit coupon flags and/or change on index 
    on continuous call or conversion.

    The constraints will not be used for system solving. the constraints used
    to solve the non linear system is done in InstData::UpdateBeforeStep.
  */
  UpdateConstraints();
  
  if( m_pConstraints )
    ApplyConstraintsToAll();
}

void CBInstData::ComputeNewSharePricesAtMaturity(const double *pdSpots, 
                                                 const size_t nNbSpots,
                                                 double *pdValues)
{
  // Get the "new share cb params" from the cb params.
  shared_ptr<pricing::CBLikeParams>
    pNewShareParams = m_cbparams.GetNewShareParams(pdSpots, nNbSpots);

  // Don't need to compute anything other than price, so use a default flags
  finance::ComputationalFlags flags;

  // but use the same solver as for the original one
  flags.SetSolverType(m_iSolverType);

  CBPricer newsharecbPricer(*pNewShareParams.get(), m_model, flags);

  AutoPtr<CBNumOutput> pCBNumOutput = newsharecbPricer.Price(); 

  // Interpolate since meshes could be different
  Interpolate(&pCBNumOutput->GetFinalMesh()[0], 
              &pCBNumOutput->GetFinalPrices()[0], 
              pCBNumOutput->GetFinalPrices().size(),
              pdSpots, pdValues, nNbSpots, 
              ExtrapolationMode_Linear, ExtrapolationMode_Linear);

  //TODO: If Maturity = Fiscal Year, put the new share prices to spot prices
}

} // namespace ihg

} // namespace ito33
