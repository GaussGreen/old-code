/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/params.cpp
// Created:     2004/02/10
// RCS-ID:      $Id: params.cpp,v 1.50 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <algorithm>
#include "ito33/afterstd.h"

#include "ito33/useexception.h"
#include "ito33/dateutils.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/yieldcurve.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/numeric/predicatetime.h"
#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/mesh/optionspacemesh.h"

#include "ito33/pricing/pathdepevent.h"
#include "ito33/pricing/dividendevents.h"
#include "ito33/pricing/dividendevents_forward.h"
#include "ito33/pricing/model.h"
#include "ito33/pricing/params.h"

using ito33::GetDoubleFrom;

using ito33::finance::Dividend;
using ito33::finance::Dividends;

extern const ito33::Error ITO33_UNEXPECTED;

namespace ito33
{

namespace pricing
{

using namespace numeric;
using namespace numeric::mesh;

typedef shared_ptr<DividendEvent> DivEventPtr;

Params::Params(Contracts &contracts)
             : m_contracts(contracts),
               m_pNumParams(new NumParams),
               m_pMeshParams(new MeshParams),
               m_dStoppingTime( m_contracts.GetMaturityTime() ),
               m_dAnalysisTime(-100.)
{
}

Params::Params(Contracts &contracts,
               const ito33::finance::SessionData& sessionData,
               const ito33::shared_ptr<NumParams> &pNumParams,
               const ito33::shared_ptr<MeshParams> &pMeshParams)
  : m_contracts(contracts),
    m_pYieldCurve( sessionData.GetYieldCurve() ),
    m_pYieldCurveForMesh( sessionData.GetYieldCurveForMesh() ),
    m_pForeignCurve( sessionData.GetForeignCurve() ),
    m_pDividends( sessionData.GetDividends() ),
    m_pNumParams(pNumParams),
    m_pMeshParams(pMeshParams),
    m_dValuationTime( GetDoubleFrom(sessionData.GetValuationDate()) ),
    m_dStoppingTime( m_contracts.GetMaturityTime() ),
    m_dAnalysisTime(-100.),
    m_dSpot( sessionData.GetSpotSharePrice() )
{  
}

void Params::SetAnalysisTime(double dAnalysisTime)
{ 
  if (  !IsBefore(dAnalysisTime, m_dValuationTime)
       && IsBefore(dAnalysisTime, m_dStoppingTime) )
    m_dAnalysisTime = dAnalysisTime; 
  else 
  {
    // invalidate the analysis time, otherwise it may crash when param reused
    // so there may have a non negative value left but valuation time and/or 
    // maturity time changed.
    m_dAnalysisTime = -100.;
  }
} 

void Params::Init()
{
  // For call notice problems, the call notice params object is re-used.
  // Hence, we need to reset the stopping time.  This implies that the
  // stopping time can only be reset afer Init.
  //m_dStoppingTime = m_contracts.GetMaturityTime();

  m_eventManager.Init();

  m_pathDepEvents.clear();
}

double Params::GetDiffusionSize(double dSquaredTotalVol) const
{
  return sqrt( (m_dStoppingTime - m_dValuationTime) * dSquaredTotalVol ); 
}

double Params::GetBackwardConvectionSize(double dSquaredTotalVol) const
{
  ASSERT_MSG( m_pYieldCurveForMesh, "Yield curve for mesh not defined.");

  ASSERT_MSG( m_pForeignCurve, "Foreign yield curve not defined.");

  return   ( m_dStoppingTime - m_dValuationTime )
         * (   m_pYieldCurveForMesh->GetContinuousRate(m_dStoppingTime)
             - m_pForeignCurve->GetContinuousRate(m_dStoppingTime)
             - 0.5 * dSquaredTotalVol
           );
}


double Params::GetForwardConvectionSize(double dSquaredTotalVol) const
{
  ASSERT_MSG( m_pYieldCurveForMesh, "Yield curve for mesh not defined.");

  ASSERT_MSG( m_pForeignCurve, "Foreign yield curve not defined.");

  return   ( m_dStoppingTime - m_dValuationTime ) 
         * (   m_pForeignCurve->GetContinuousRate(m_dStoppingTime)
             - m_pYieldCurveForMesh->GetContinuousRate(m_dStoppingTime)
             - 0.5 * dSquaredTotalVol
           );
}

std::vector<double> Params::GenerateSpaceMesh(const Model& model) const
{
  numeric::mesh::OptionSpaceMesh osgGrid;

  // Get the diffusion size from the used numerical model
  double dSquaredTotalVol = model.GetSquaredTotalVolForMesh
                            ( GetStoppingTime(), GetSpotSharePrice() );

  osgGrid.SetDiffusionSize( GetDiffusionSize(dSquaredTotalVol) );

  // Compute the convection size for this problem
  osgGrid.SetConvectionSize( GetBackwardConvectionSize(dSquaredTotalVol) );
 
  // Number of requested points, get it from m_pNumParams 
  size_t nNbPoints = GetNumParams()->GetNbSpaceSteps();
  
  // Generate a log space mesh centered at zero
  std::vector<double> vecGrid;
  
  osgGrid.Build(nNbPoints, 0., vecGrid);

  // Re-center the space mesh on the spot
  double dLogSpot = log( GetSpotSharePrice() );
  
  for (size_t nIdx = 0; nIdx < vecGrid.size(); nIdx++)
    vecGrid[nIdx] += dLogSpot;

  return vecGrid;
}

void Params::ConstructDividendEvents(numeric::ExtrapolationMode emLeft,
                                     numeric::ExtrapolationMode emRight,
                       numeric::InterpolationMethod interpolationMethod)
{ 
  if (!m_pDividends)
    return;
 
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
      shared_ptr<DividendEvent> pEvent;

      switch (pDividend->type)
      {
      case Dividend::PseudoCash :
        pEvent = DivEventPtr
            ( new PseudoCashDividendEvent
                  ( dDividendTime, pDividend->value, pDividend->pseudoYield ) );         
        break;

      case Dividend::Cash :
        pEvent = DivEventPtr( new CashDividendEvent
                                  ( dDividendTime, pDividend->value ) );
        break;
      
      case Dividend::Yield :
        pEvent = DivEventPtr( new YieldDividendEvent
                                  ( dDividendTime, pDividend->value ) );
        break;
      
      default:
        throw EXCEPTION_MSG
             (
              ITO33_UNEXPECTED,
              TRANS("Sorry, unknown dividend type.")
             );
      }

      pEvent->SetEMLeft(emLeft);
      pEvent->SetEMRight(emRight);
      pEvent->SetInterpolationMethod(interpolationMethod);

      m_eventManager.AddEvent(pEvent);
    }
  }
}

void Params::ConstructForwardDividendEvents(numeric::ExtrapolationMode emLeft,
                                            numeric::ExtrapolationMode emRight,
                              numeric::InterpolationMethod interpolationMethod)
{ 
  if (!m_pDividends)
    return;

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
      shared_ptr<DividendEvent> pEvent;

      switch (pDividend->type)
      {
      case Dividend::PseudoCash :
        pEvent = DivEventPtr
            ( new PseudoCashDividendEvent_Forward
                  ( dDividendTime, pDividend->value, pDividend->pseudoYield ) );         
        break;
     
      case Dividend::Yield :
        pEvent = DivEventPtr( new YieldDividendEvent_Forward
                                  ( dDividendTime, pDividend->value ) );
        break;

      case Dividend::Cash :
        pEvent = DivEventPtr( new CashDividendEvent_Forward
                                  ( dDividendTime, pDividend->value ) );
        break;

      default:
        throw EXCEPTION_MSG
              (
                ITO33_UNEXPECTED,
                TRANS("Sorry, unknown dividend type.")
              );
      }

      pEvent->SetEMLeft(emLeft);
      pEvent->SetEMRight(emRight);
      pEvent->SetInterpolationMethod(interpolationMethod);

      m_eventManager.AddEvent(pEvent);
    }
  }
}

void Params::GetSpecialTimes(SpecialTimes& specialTimes) const
{
  specialTimes.clear();

  m_eventManager.GetSpecialTimes(specialTimes);

  SpecialTimes specialTimesContract;
  m_contracts.GetSpecialTimes(specialTimesContract);
  SpecialTimes::const_iterator iter;
  for (iter = specialTimesContract.begin(); 
       iter != specialTimesContract.end();
       ++iter)
    specialTimes.push_back(*iter);

  // Add path dependent events
  std::list< shared_ptr<PathDepEvent> >::const_iterator iterEvent;

  for (iterEvent = m_pathDepEvents.begin(); 
       iterEvent != m_pathDepEvents.end();
       ++iterEvent)
    specialTimes.push_back( SpecialTime((*iterEvent)->GetTime()) );

  // Add AnalysisTime if it's in the range
  if (   !IsBefore(m_dAnalysisTime, m_dValuationTime)
       && IsBefore(m_dAnalysisTime, m_dStoppingTime) )
    specialTimes.push_back(m_dAnalysisTime);

  // Add valuation time
  specialTimes.push_front(SpecialTime(m_dValuationTime) ); 

  // Add stopping time, refine 
  specialTimes.push_back( SpecialTime(m_dStoppingTime, RefineLevel_VeryHigh) );
}

} // namespace pricing

} // namespace ito33

