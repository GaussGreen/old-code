/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/cb/cbinstdata.cpp
// Purpose:     implementations for cb instdata class in the HG model
// Created:     2005/04/11
// RCS-ID:      $Id: cbinstdata.cpp,v 1.20 2006/08/19 23:44:26 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#if ( defined(_MSC_VER) && _MSC_VER >= 1400 )
// still want to use fill_n but don't know how to avoid the compile warning
// warning C4996: 'std::_Fill_n' was declared deprecated
// so just disable the warning
#pragma warning(disable:4996)
#endif

#include "ito33/sharedptr.h"
#include "ito33/array.h"

#include "ito33/numeric/interpolation.h"

#include "ito33/pricing/event.h"
#include "ito33/pricing/cbevent.h"
#include "ito33/pricing/cbconstraints.h"
#include "ito33/pricing/cbmeshmanager.h"

#include "ito33/finance/computationalflags.h"

#include "hg/cbinstdata.h"
#include "hg/model.h"
#include "hg/cbpricer.h"
#include "hg/cbnumoutput.h"
#include "hg/payoff.h"

namespace ito33
{

namespace hg
{
   using namespace numeric;


CBInstData::CBInstData( pricing::CBLikeParams& cbparams, 
                        Model& model, 
                        pricing::CBMeshManager& cbmeshes)
                      : InstDataWithConstraints(cbparams, model, cbmeshes),
                        m_cbparams(cbparams), m_cbmeshes(cbmeshes)
                        
{ 
  m_bHasNewShare = m_cbparams.HasNewShare();
}

void CBInstData::Alloc(size_t nNbX)
{
  InstDataWithConstraints::Alloc(nNbX);

  // Unlike option, the recovery values may depend on the spots
  // Recovery values are the same for different regimes, so if we are using the
  // same space mesh for different regimes, it can be allocated with the 
  // maximum size of the space mesh
  m_pdRecoveryValues.resize(m_nNbSpotsMax);

  // Allocate the new share arrays
  if ( m_bHasNewShare )
  {
    m_pdNewSharePrices = Array<double>(nNbX);
    m_pdOldNewSharePrices = Array<double>(nNbX);
    m_pdOldOldNewSharePrices = Array<double>(nNbX); 
  }
}

void CBInstData::SetupFlags(const finance::ComputationalFlags& flags)
{
  InstDataWithConstraints::SetupFlags(flags);

  m_sensitivityMethod = SensitivityMethod_None;

  // m_sensitivityMethod = SensitivityMethod_PDE;

  m_iSolverType = flags.GetSolverType();
}

void CBInstData::Init()
{
  // Allocation of memory by using the maximum size of the space meshes
  m_nNbSpotsMax = m_cbmeshes.GetNbSMax();

  Alloc(m_nNbSpotsMax * m_nNbRegimes);
  
  // Initializes the pointers to the space mesh. Note that they shouldn't 
  // change, otherwise, they need to be updated together with the size
  m_pdS = m_cbmeshes.GetS();
  m_pdLogS = m_cbmeshes.GetLogS();  
}

void 
CBInstData::ComputeNewSharePricesAtMaturity
(const double *pdS, size_t nNbS, double *pdNewSharePrices)
{
  // Get the "new share cb params" from the cb params.
  shared_ptr<pricing::CBLikeParams>
    pNewShareParams = m_cbparams.GetNewShareParams(pdS, nNbS);

  // Don't need to compute anything other than price, so use a default flags
  finance::ComputationalFlags flags;
  
  // but use the same solver as for the original one
  flags.SetSolverType(m_iSolverType);

  CBPricer newsharecbPricer(*pNewShareParams.get(), m_model, flags);

  AutoPtr<CBNumOutput> pCBNumOutput = newsharecbPricer.Price(); 
    
  const std::vector<double>& pdFinalMesh = pCBNumOutput->GetFinalMesh();
  const std::vector<double>& pdFinalPrices = pCBNumOutput->GetFinalPrices();
  const size_t nNbFinalS = pdFinalMesh.size();

  // Interpolate since meshes could be different
  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
    Interpolate(&pdFinalMesh[0], 
                &pdFinalPrices[0] + nIdxR * nNbFinalS, 
                nNbFinalS,
                pdS, pdNewSharePrices + nIdxR * nNbS, nNbS);
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
  m_nNbX = m_nNbS * m_nNbRegimes;

  if ( m_bHasNewShare )
    ComputeNewSharePricesAtMaturity( m_pdS, m_nNbS, m_pdNewSharePrices.Get() );

  if ( m_pPayoff )
    m_pPayoff->Get(m_pdS, m_pdPrices.Get(), m_nNbS);
  else
  {
    double* pdPrices = m_pdPrices.Get();
    for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
      m_cbparams.GetInitialValues(m_pdS, m_nNbS, pdPrices + nIdxR * m_nNbS);
  }

  for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
    m_pdOldPrices[nIdx] = m_pdOldOldPrices[nIdx] = m_pdPrices[nIdx];

  std::fill_n(m_piFrozenFlags.Get(), m_nNbX, 0); 
  
  if ( m_bComputeFugit )
    std::fill_n(m_pdFugits.Get(), m_nNbX, 0.);

  m_aData.m_bIsValid = false;

  DoEvents();
}

void CBInstData::Swap()
{
  InstDataWithConstraints::Swap();

  if (m_bHasNewShare)
  {
    swap(m_pdOldOldNewSharePrices, m_pdNewSharePrices);
    swap(m_pdOldOldNewSharePrices, m_pdOldNewSharePrices);
  }
}

void CBInstData::ApplyEvent(const pricing::Event *pEvent)
{
  InstDataWithConstraints::ApplyEvent(pEvent);

  if ( m_bHasNewShare && pEvent->GetType() == pricing::ET_Dividend )
  {
    for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
      pEvent->ApplyToPrice( m_pdS, m_pdNewSharePrices.Get() + nIdxR * m_nNbS,
                            m_nNbS );
  }
}

void CBInstData::InterpWithPassageOfSpaceMesh
     (const double *pdOldSpots, size_t nNbOldSpots, 
      const double *pdSpots, size_t nNbSpots)
{

  // Interpolation of price vector on the new grid. 
  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
    Interpolate(pdOldSpots, m_pdOldPrices.Get() + nIdxR * nNbOldSpots,
                nNbOldSpots, 
                pdSpots, m_pdPrices.Get() + nIdxR * nNbSpots, nNbSpots, 
                ExtrapolationMode_Linear, ExtrapolationMode_Linear);

  // Interpolation of new share prices vector on the new grid.
  if ( m_bHasNewShare )
    for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
      Interpolate(pdOldSpots, 
                  m_pdOldNewSharePrices.Get() + nIdxR * nNbOldSpots, 
                  nNbOldSpots, 
                  pdSpots,
                  m_pdNewSharePrices.Get() + nIdxR * nNbSpots, 
                  nNbSpots);
}

void CBInstData::ApplyConstraintsToAll()
{
  //Apply Constraints to the price.
  m_pConstraints->Apply(m_pdPrices.Get(), m_piFrozenFlags.Get(), m_nNbX);
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
  InterpWithPassageOfSpaceMesh( m_cbmeshes.GetOldS(), m_cbmeshes.GetOldNbS(), 
                                m_cbmeshes.GetS(), m_cbmeshes.GetNbS() );
  
  m_nNbS = m_cbmeshes.GetNbS();
  m_nNbX = m_nNbS * m_nNbRegimes;
  
  if ( m_bHasNewShare )
  {
    // Updates the new share price across the fiscal year start date (at t-)
    if ( m_cbparams.IsStartOfYear() )
    {
      const double* pdS = m_cbmeshes.GetS();
      
      for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
        for (size_t nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
          m_pdNewSharePrices[nIdxS + nIdxR * m_nNbS] = pdS[nIdxS];

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
  
  if ( m_pConstraints )
    ApplyConstraintsToAll();
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

  InstDataWithConstraints::UpdateBeforeStep();

  if ( m_cbparams.GetCBLike().IsCrossCurrency() )
    m_cbmeshes.UpdateFXRate();

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
  Array<double> pdTmp(m_nNbS);
  Array<double> aConstraints(m_nNbX);
  double* pdConstraints = aConstraints.Get();

  bool bHasConstraints = false;
  
  // If the calls have a notice period, then a mini cb problem needs to be
  // solved.  The mini cb is model specific, so is solved here.  Without
  // call notice, the params class can determine the constraint
  if (    m_cbparams.GetCalls()->IsActive() 
       && m_cbparams.GetIndexCall() != INVALIDINDEX )
  {    
    bool bHasNoticePeriod = m_cbparams.GetCalls()->HasNoticePeriod();

    const double* pdNoticeS;

    if ( m_bHasNewShare )
      pdNoticeS = m_pdNewSharePrices.Get();
    else
      pdNoticeS = 0;

    if ( bHasNoticePeriod )    
      SolveCallNotice(m_pdS, m_nNbS, pdNoticeS, pdConstraints);

    size_t nIdxStartConversion;

    for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
    {
      const double* pdNewS = m_pdS;
      if ( m_bHasNewShare )
        pdNewS = m_pdNewSharePrices.Get() + nIdxR * m_nNbS;

      m_cbparams.GetCallConstraintValues
                 (m_pdS, m_nNbS, pdConstraints + nIdxR * m_nNbS, 
                  nIdxStartConversion, m_pdS, bHasNoticePeriod);
    }

    m_constraints.UpdateCall(pdConstraints, m_nNbX, m_nNbX);

    bHasConstraints = true;
  }
  else
    m_constraints.DisableTheCall();

  // Nothing special about the conversions. Let params determine the constraint
  bool bUpdateConv = false;
  
  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
  {
    const double* pdNewS = m_pdS;
    if ( m_bHasNewShare )
      pdNewS = m_pdNewSharePrices.Get() + nIdxR * m_nNbS;

    bUpdateConv = m_cbparams.GetConversionConstraintValues
                             ( m_pdS, m_nNbS, 
                               pdNewS,
                               pdConstraints + nIdxR * m_nNbS);
  }

  if ( bUpdateConv )
  {
    m_constraints.UpdateConv(pdConstraints, m_nNbX);
    bHasConstraints = true;
  }
  else
    m_constraints.DisableTheConv();

  // Nothing special about the puts. Let params determine the constraint
  if ( m_cbparams.GetPutConstraintValues(m_pdS, m_nNbS, pdTmp.Get()) )
  {
    for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
    {
      double* pdConstraintsTmp = aConstraints.Get() + nIdxR * m_nNbS;
      for (size_t nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
        pdConstraintsTmp[nIdxS] = pdTmp[nIdxS];
    }

    m_constraints.UpdatePut(pdConstraints, m_nNbX);

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
                                 size_t nNbS,
                                 const double* pdNewSharePrices,
                                 double *pdValues)
{

  // Get the mini cb params from the "big" cb params. Assume
  // it is constructed with the correct, current info
  pricing::CBLikeParams
    *pCallNoticeParams = m_cbparams.GetCallNoticeParams();

  // Now that the mini Cb is not constructed from a financial one
  // we need to set the payoff manually or call a function in CB
  // to let it construt itself the payoff objet from its redemption value.
  // Remark: This should be moved from params to here because now, this
  //         depends on model.
  Array<double> pdInitialValues(nNbS * m_nNbRegimes);

  // Get the call strikes at the maturity time of the call notice CB
  // and use it as initial value for the mini cb pricer
  if ( pdNewSharePrices )
    for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
      m_cbparams.GetCalls()->GetCallStrikesWithoutCoupon
          (pCallNoticeParams->GetCBLike().GetMaturityTime(),
          pdS, nNbS, pdNewSharePrices + nIdxR * nNbS, 
          pdInitialValues.Get() + nIdxR * m_nNbS);
  else
    for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
      m_cbparams.GetCalls()->GetCallStrikesWithoutCoupon
          (pCallNoticeParams->GetCBLike().GetMaturityTime(),
          pdS, nNbS, pdS, pdInitialValues.Get() + nIdxR * m_nNbS );

  // Do not need to compute anything other than price, so use a default flags
  finance::ComputationalFlags flags;
  
  // but use the same solver as for the original one
  flags.SetSolverType(m_iSolverType);

  CBPricer minicbPricer(*pCallNoticeParams, m_model, flags);
  
  shared_ptr<Payoff> 
    pPayoff( new Payoff(pdS, pdInitialValues.Get(), nNbS, m_nNbRegimes) );
   
  minicbPricer.SetInitialValue(pPayoff);

  AutoPtr<CBNumOutput> pCBNumOutput = minicbPricer.Price(); 

  // Interpolate since meshes could be different
  const std::vector<double>& pdFinalMesh = pCBNumOutput->GetFinalMesh();
  const std::vector<double>& pdFinalPrices = pCBNumOutput->GetFinalPrices();
  const size_t nNbFinalS = pdFinalMesh.size();

  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
    Interpolate(&pdFinalMesh[0], 
                &pdFinalPrices[0] + nIdxR * nNbFinalS, 
                nNbFinalS,
                pdS, pdValues + nIdxR * nNbS, nNbS);
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


} // namespace hg

} // namespace ito33
