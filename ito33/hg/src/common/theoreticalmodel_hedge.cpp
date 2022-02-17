/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/theoreticalmodel_hedge.cpp
// Purpose:     HG Hedge implementation
// Created:     2005/05/10
// RCS-ID:      $Id: theoreticalmodel_hedge.cpp,v 1.22 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <fstream>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"
#include "ito33/dateutils.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/derivative.h"
#include "ito33/finance/derivatives.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/surfacegeneral.h"
#include "ito33/numeric/predicatedouble.h"

#include "ito33/hg/error.h"
#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/hedgeoutput.h"
#include "ito33/hg/hedgeratiodata.h"
#include "ito33/hg/heroflags.h"
#include "ito33/hg/modeloutput.h"
#include "ito33/hg/version.h"

#include "hg/model.h"
#include "hg/hedgingdata.h"
#include "hg/heroparams.h"
#include "hg/heropricer.h"
#include "hg/heronumoutput.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/common.h"
#include "hg/xml/common.h"

extern const ito33::hg::Error ITO33_HG_HERO_MISSINGDATA;

namespace ito33
{

namespace hg
{

static void 
CheckDataForHero(const shared_ptr<finance::ModelOutput>& pMO, double dEndTime)
{
  // Check for existence of NumOutput: we don't force MO to have a Numoutput
  CHECK_COND(pMO->GetNumOutput(), ITO33_HG_HERO_MISSINGDATA);

  shared_ptr<BackwardNumOutput> 
    pNumOutput( static_pointer_cast<BackwardNumOutput>( pMO->GetNumOutput() ) );

  shared_ptr<numeric::SurfaceGeneral> 
    pPriceSurface( pNumOutput->GetPriceSurface(0) );

  CHECK_COND(pPriceSurface, ITO33_HG_HERO_MISSINGDATA);

  ASSERT( !pPriceSurface->GetDomain()->GetTimes().empty() );

  // Assume backward pricing, but also handle forward pricers
  double dLastTime = pPriceSurface->GetDomain()->GetTimes().front();
  double dTmpTime = pPriceSurface->GetDomain()->GetTimes().back();
  if ( dTmpTime > dLastTime )
    dLastTime = dTmpTime;

  CHECK_COND( numeric::IsEqualOrAfter(dLastTime, dEndTime), 
              ITO33_HG_HERO_MISSINGDATA);
}

shared_ptr<HedgeOutput>
TheoreticalModel::DoComputeHERO
                  (const finance::Derivative& target, 
                   const finance::Derivatives& hedgeInstruments,
                   bool bComputeHERO, shared_ptr<HEROFlags> pHEROFlags)
{
  const finance::Derivatives::Elements& elements = hedgeInstruments.GetAll();
  
  if ( IsDebugOutputEnabled() )
  {
    // dump everything we can before calling the hedge/hero functions, in case
    // they fail
    std::ofstream ofs(GetDebugOutputFile().c_str());

    XML::RootTag tagRoot(XML_TAG_HG_ROOT, ofs);
    tagRoot.precision(10);
    tagRoot.Attr(XML_ATTR_ROOT_VERSION, ITO33_HG_SERIES_DOT_STRING);
     
    target.GetSessionData()->Dump(tagRoot);

    // dump the Model itself
    Dump(tagRoot);

    // dump compute hero boolean
    tagRoot.Element(XML_TAG_HG_COMPUTEHERO)(bComputeHERO);

    // dump target
    tagRoot.Element(XML_TAG_HG_TARGET).Element(target);

    // dump the hedge instruments
    {
      ito33::XML::Tag tagDerivatives(XML_TAG_HG_HEDGEINSTRUMENTS, tagRoot);
    
      finance::Derivatives::Elements::const_iterator iter;
    
      for (iter = elements.begin(); iter != elements.end(); ++iter)
        iter->first->Dump(tagDerivatives);
    }
  }  

  // Use a clone of the current TheoreticalModel since we need to use 
  // specific flags
  shared_ptr<TheoreticalModel> pModelClone( Clone() );
  
  // Use the same quality control
  pModelClone->SetQualityControl(m_pQualityControl);

  // Need values at valuation date for hedging.  Use analysis date.
  shared_ptr<finance::ComputationalFlags> 
    pFlags(new finance::ComputationalFlags);
  pFlags->SetAnalysisDate( target.GetSessionData()->GetValuationDate() );

  // If the HERO is computed, surface data is needed
  if ( bComputeHERO )
    pFlags->SetComputeSurface(true);

  pModelClone->SetExternalFlags(pFlags);

  const size_t nNbHedges = elements.size();  

  // Create the return hedge object here, so model outputs can be set
  shared_ptr<HedgeOutput> pHedgeOutput(new HedgeOutput);
  std::vector< shared_ptr<finance::ModelOutput> > ppModelOutputs(nNbHedges);
  std::vector< shared_ptr<HedgeRatioData> > ppHedgeRatioData(nNbHedges);

  // Need to use the real underlying process for hero/hedging
  shared_ptr<UnderlyingProcess> pRealUP( GetRealUnderlyingProcess() );

  // Even if the hero is not being computed, the HedgingData class can be used
  // to compute the matrices and vectors needed by the hedge calculations.
  // This helps to avoid code duplication.  The HedgingData class needs the
  // numerical outputs of all contracts to do its work.  These will be
  // available through the model output pointers.
  HedgingData hedgingData(target, pRealUP, nNbHedges); 
  
  // Compute the model output for the target and each hedge instrument  
  const finance::Derivative* pDerivative;
  pDerivative = &target;
  shared_ptr<finance::ModelOutput> pMO( pModelClone->Compute(*pDerivative) );

  pHedgeOutput->SetTargetModelOutput( pMO );
  hedgingData.SetTargetModelOutput( pMO );

  // The hero needs data from valuation date to maturity date for the target
  Date targetMaturityDate = target.GetMaturityDate();
  double dTargetMaturityTime = GetDoubleFrom(targetMaturityDate);

  if ( bComputeHERO ) // check for full surface info
    CheckDataForHero(pMO, dTargetMaturityTime);

  for (size_t nIdxH = 0; nIdxH < nNbHedges; nIdxH++)
  {
    pDerivative = elements[nIdxH].first.get();

    shared_ptr<finance::ModelOutput> pMO( pModelClone->Compute(*pDerivative) );

    ppModelOutputs[nIdxH] = pMO;
    
    shared_ptr<HedgeRatioData> hedgeRatioData(new HedgeRatioData);
    hedgeRatioData->SetDerivative( elements[nIdxH].first );
    hedgeRatioData->SetModelOutput( pMO );

    ppHedgeRatioData[nIdxH] = hedgeRatioData;

    // Need hedge contract data from valuation date to the lesser of the 
    // hedge maturity and target maturity.  If the hedge maturity is less 
    // than the target maturity, the code should give zero prices from the
    // hedge to target maturity dates.
    if ( bComputeHERO )
    {    
      // Get lesser of target maturity and hedge maturity
      Date hedgeMaturityDate = pDerivative->GetMaturityDate();
      double dHedgeMaturityTime = GetDoubleFrom(hedgeMaturityDate);

      double dTmp = std::min(dTargetMaturityTime, dHedgeMaturityTime);

      CheckDataForHero(pMO, dTmp);
      
    } // check for surface info
  }

  hedgingData.SetHedgeModelOutputs(ppModelOutputs);
 
  // Compute the hedge ratios
  double dUnderlyingHedge;
  std::vector<double> pdHedges(nNbHedges);
  hedgingData.ComputeFinalHedgeRatios(dUnderlyingHedge, pdHedges);

  for (size_t nIdxH = 0; nIdxH < nNbHedges; nIdxH++)
    ppHedgeRatioData[nIdxH]->SetRatio( pdHedges[nIdxH] );

  // Compute the HERO, if requested
  shared_ptr<finance::ModelOutput> pHeroOutput;
  if ( bComputeHERO )
  {
    // Use the real underlying process
    Model model(*pRealUP, pRealUP);
    model.SetSharpeRatio(m_dSharpeRatio);

    // Only compute surface and analysis date data, and only if requested
    shared_ptr<finance::ComputationalFlags> 
      pFlagsTmp(new finance::ComputationalFlags);

    if ( pHEROFlags )
    {
      pFlagsTmp->SetAnalysisDate( pHEROFlags->GetAnalysisDate() );
      pFlagsTmp->SetComputeSurface( pHEROFlags->GetComputeSurface() );
    }

    // Create the params, pricer, and then price to get the HERO
    shared_ptr<numeric::NumParams> pNumParams( GetNumParams(target) );
    shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);
    HeroParams heroParams(hedgingData, *target.GetSessionData(), pNumParams, 
                          pMeshParams);

    HeroPricer pricer(heroParams, model, *pFlagsTmp);

    AutoPtr<HeroNumOutput> pNumOutput( pricer.Price() );

    pHeroOutput = pNumOutput->GetModelOutput();
  }

  // Finish filling the return hedge output. The target output should 
  // have already been set.
  pHedgeOutput->SetUnderlyingHedgeRatio(dUnderlyingHedge);
  pHedgeOutput->SetHedgeRatioData(ppHedgeRatioData);

  if ( bComputeHERO)
    pHedgeOutput->SetHEROModelOutput(pHeroOutput);

  return pHedgeOutput;
}

shared_ptr<HedgeOutput>
TheoreticalModel::Hedge(const finance::Derivative& target, 
                        const finance::Derivatives& hedgeInstruments)
{
  return DoComputeHERO(target, hedgeInstruments, false, shared_ptr<HEROFlags>());
}

shared_ptr<HedgeOutput>
TheoreticalModel::ComputeHERO
                  (const finance::Derivative& target, 
                   const finance::Derivatives& hedgeInstruments,
                   shared_ptr<HEROFlags> flags)
{
  return DoComputeHERO(target, hedgeInstruments, true, flags);
}

} // namespace hg

} // namespace ito33
