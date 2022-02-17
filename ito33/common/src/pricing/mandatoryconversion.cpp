///////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/mandatoryconversions.cpp
// Purpose:     Mandatory conversion class
// Author:      Wang
// Created:     2004/08/16
// RCS-ID:      $Id: mandatoryconversion.cpp,v 1.21 2006/02/22 15:25:37 yann Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
///////////////////////////////////////////////////////////////////////////

#include "ito33/dateutils.h"
#include "ito33/finance/bondlike/pepslike.h"
#include "ito33/finance/bondlike/generalizedpepslike.h"
#include "ito33/finance/bondlike/percslike.h"
#include "ito33/finance/bondlike/bondliketerms.h"

#include "ito33/pricing/cblikeparams.h"
#include "ito33/pricing/cbcalls.h"
#include "ito33/pricing/callfixedshare.h"
#include "ito33/pricing/mandatoryconversion.h"

namespace ito33
{

namespace pricing
{
   

MandatoryConversion::MandatoryConversion(const finance::GeneralizedPEPSLike& peps)
  : ConversionProvisions(true, false)
{
  MandatoryPayoffStructure
          structure( peps.GetDownsideConversionRatio(),
                      peps.GetUpsideBaseConversionRatio(),
                      peps.GetLowerStrike(),
                      peps.GetHigherStrike()
                    );

  m_dConversionPrice = structure.GetHigherStrike();
  m_dRatioForRootComputation = structure.GetMinConversionRatio();

  m_nNbConversions = 1;
  if(peps.HasOptionalConversion())
    m_nNbConversions = 2;

  m_pdStartTimes.resize(m_nNbConversions);
  m_pdEndTimes.resize(m_nNbConversions);
  m_pdCashs.resize(m_nNbConversions, 0);
  m_pdTriggerRates.resize(m_nNbConversions, 0.);
  m_pCoCoTypes.resize(m_nNbConversions, finance::CoCoType_Max); 
  m_pbIsActive.resize(m_nNbConversions, true);
  m_pPeriods.resize(m_nNbConversions);

  size_t nIdx = 0;
  double dMaturity = GetDoubleFrom(peps.GetBondLikeTerms()->GetMaturityDate());
  if(peps.HasOptionalConversion())
  {
    m_pdStartTimes[nIdx]
        = GetDoubleFrom(peps.GetBondLikeTerms()->GetIssueDate());
    m_pdEndTimes[nIdx] = dMaturity;
    m_pPeriods[nIdx]
      = MandatoryPayoffStructure(structure.GetMinConversionRatio());
    nIdx++;
  }

  m_pdStartTimes[nIdx] = dMaturity;
  m_pdEndTimes[nIdx] = dMaturity;
  m_pPeriods[nIdx] = structure;
}

MandatoryConversion::MandatoryConversion(const finance::PEPSLike& peps)
  : ConversionProvisions(true, false)
{
  double
    dMaxRatio = peps.GetMaxConversionRatio(),
    dMinRatio = peps.GetMinConversionRatio();
  
  double dIssuePrice = peps.GetBondLikeTerms()->GetIssuePrice()
                * peps.GetBondLikeTerms()->GetNominal();
  m_dConversionPrice = dIssuePrice / dMinRatio;
  m_dRatioForRootComputation = dMinRatio;

  m_nNbConversions = 1;
  if(peps.HasOptionalConversion())
    m_nNbConversions = 2;

  m_pdStartTimes.resize(m_nNbConversions);
  m_pdEndTimes.resize(m_nNbConversions);
  m_pdCashs.resize(m_nNbConversions, 0);
  m_pdTriggerRates.resize(m_nNbConversions, 0.);
  m_pCoCoTypes.resize(m_nNbConversions, finance::CoCoType_Max); 
  m_pbIsActive.resize(m_nNbConversions, true);
  m_pPeriods.resize(m_nNbConversions);

  size_t nIdx = 0;
  double dMaturity = GetDoubleFrom(peps.GetBondLikeTerms()->GetMaturityDate());
  if(peps.HasOptionalConversion())
  {
    m_pdStartTimes[nIdx]
        = GetDoubleFrom(peps.GetBondLikeTerms()->GetIssueDate());
    m_pdEndTimes[nIdx] = dMaturity;
    m_pPeriods[nIdx] = MandatoryPayoffStructure(dMinRatio);
    nIdx++;
  }

  m_pdStartTimes[nIdx] = dMaturity;
  m_pdEndTimes[nIdx] = dMaturity;
  m_pPeriods[nIdx] = MandatoryPayoffStructure(dIssuePrice, dMinRatio, dMaxRatio);
}


MandatoryConversion::MandatoryConversion(const finance::PERCSLike& percs)
  : ConversionProvisions(true, false)
{
  double dCapPrice = percs.GetCapPrice()
                   * percs.GetBondLikeTerms()->GetNominal();
  m_dConversionPrice = dCapPrice / percs.GetMaxConversionRatio();
  m_dRatioForRootComputation = percs.GetMaxConversionRatio();

  m_nNbConversions = 1;

  m_pdStartTimes.resize(1);
  m_pdEndTimes.resize(1);
  m_pdCashs.resize(1, 0);
  m_pdTriggerRates.resize(1, 0.);
  m_pCoCoTypes.resize(1, finance::CoCoType_Max); 
  m_pbIsActive.resize(1, true);
  m_pPeriods.resize(1);

  double dMaturity = GetDoubleFrom(percs.GetBondLikeTerms()->GetMaturityDate());

  m_pdStartTimes[0] = dMaturity;
  m_pdEndTimes[0] = dMaturity;
  m_pPeriods[0]
    = MandatoryPayoffStructure(dCapPrice, percs.GetMaxConversionRatio());
}


bool MandatoryConversion::GetGrossParities
     (const double* pdS, size_t nNbS, 
      const double* pdNewSharePrices,
      double* pdValues) const
{
  size_t nIdxConversion = m_pParams->GetIndexConversion();

  if ( nIdxConversion == INVALIDINDEX )
    return false;

  if ( m_pParams->GetCBLike().IsCrossCurrency() )
    m_pPeriods[nIdxConversion].GetValues
        (pdS, pdValues, nNbS, pdNewSharePrices, m_pParams->GetFXRate());
  else
    m_pPeriods[nIdxConversion].GetValues
        (pdS, pdValues, nNbS, pdNewSharePrices);

  return true;
}

void MandatoryConversion::ComputeRoots
              (
                double dTime, size_t& nNbRoots, 
                numeric::mesh::Root* pRoots, bool bPlus
              )
{
  CBCalls* pCalls = dynamic_cast<CBCalls*>(m_pParams->GetCalls());

  if( pCalls != 0)
    ComputeRootsWithCallCash(*pCalls, dTime, nNbRoots, pRoots, bPlus);
  else
  {
    ComputeRootsWithCallTrigger
      (
        m_pParams->GetCalls()->GetTriggerRate(m_pParams->GetIndexCall()),
        m_pParams->GetCalls()->GetTriggerAsPercentageOf(),
        nNbRoots, pRoots, dTime, bPlus
      );
  }
}

void 
MandatoryConversion::ComputeRootsWithCallTrigger
(double dTriggerRate, finance::TriggerAsPercentageOf asPercentagOf,
 size_t& nNbRoots, numeric::mesh::Root* pRoots, 
 double dTime, bool bPlus)
{
  ASSERT_MSG(dTriggerRate >= 1, 
             "Trigger rate of a fixed share call must be greater than 1.");

  nNbRoots = 1;
   
  pRoots[0].m_dRoot = dTriggerRate 
                    * GetConversionPrice(dTime, asPercentagOf, bPlus);
  
  pRoots[0].m_iRootType = numeric::mesh::RootType_Real;
}

void MandatoryConversion::ComputeRootsWithCallCash
              (
                const CBCalls& calls, double dTime, size_t& nNbRoots, 
                numeric::mesh::Root* pRoots, bool bPlus
              )
{
  // note: the index is setup before calling this function
  size_t
    nIdxCall = m_pParams->GetIndexCall(),
    nIdxConversion = m_pParams->GetIndexConversion();

  nNbRoots = 1;
   
  double dTriggerLevel = calls.GetTriggerRate(nIdxCall)
    * GetConversionPrice(dTime, calls.GetTriggerAsPercentageOf(), bPlus);  

  // in case of a virtual call, call constraint is off, but we still want to
  // add call trigger as root.
  if ( !calls.IsActive() && dTriggerLevel > 0.) 
  {
    pRoots[0].m_dRoot = dTriggerLevel;
    pRoots[0].m_iRootType = numeric::mesh::RootType_Real;

    return;
  }

  double dConvAccrued;

  if (m_bKeepAccrued)
    dConvAccrued = m_pParams->GetCashFlows()->GetAccruedInterest(dTime, bPlus);
  else
    dConvAccrued = 0.;
  
  double dRatio = m_dRatioForRootComputation;

  if ( m_pParams->GetCBLike().IsCrossCurrency() )
    dRatio *= m_pParams->GetFXRate(dTime);
  
  pRoots[0].m_dRoot = (   calls.GetStrikeWithoutCoupon(dTime, bPlus)  
                        + m_pdCashs[nIdxConversion] - dConvAccrued ) / dRatio;
  
  if (pRoots[0].m_dRoot < dTriggerLevel)
    pRoots[0].m_dRoot = dTriggerLevel;

  pRoots[0].m_iRootType = numeric::mesh::RootType_Real;
}

double MandatoryConversion::GetConversionPrice(double /* dTime */,
  finance::TriggerAsPercentageOf /*triggerAsPercentageOf*/,
  bool /* bPlus */) const  
{ 
  // it is always higher strike, even for cross currency
  return m_dConversionPrice;
}


double MandatoryConversion::GetConversionPrice(finance::TriggerAsPercentageOf /*of*/) const
{
  // it is always higher strike, even for cross currency
  return m_dConversionPrice;
}

void MandatoryConversion::SetStockAverage(double dStockAverage)
{
  ASSERT_MSG(dStockAverage >= 0.,"Stock Average must be postive.");
 
  //the last period corresponds to the maturity period
  // where the holder is forced to convert
  m_pPeriods[m_nNbConversions-1].SetStockAverage(dStockAverage);
}
   
double MandatoryConversion::GetConversionRatio(double dS) const
{
 return m_pPeriods[m_nNbConversions-1].GetConversionRatio(dS);
}


void MandatoryConversion::SetConversionRatioAverage(double dConversionRatio)
{
  ASSERT_MSG(dConversionRatio >= 0.,"Conversion ratio must be positive.");

  m_pPeriods[m_nNbConversions-1].SetConversionRatioAverage(dConversionRatio);
}

} // namesapce pricing

} // namespace ito33
