/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/theoreticalmodel_pricepepslike.cpp
// Purpose:     Implementation of PEPS-like pricing using HG model
// Created:     2005/04/11
// RCS-ID:      $Id: theoreticalmodel_pricepepslike.cpp,v 1.9 2006/08/19 23:44:26 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/bondlike/pepslike.h"

#include "ito33/pricing/mandatory.h"
#include "ito33/pricing/mandatoryparams.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/surfacegeneral.h"

#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/bondlikeoutput.h"

#include "hg/model.h"
#include "hg/cbnumoutput.h"
#include "hg/cbpricer.h"
#include "hg/license.h"
#include "hg/cbutil.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_THIS_MODULE(HGPricePEPSLike);

namespace ito33
{

namespace hg
{


shared_ptr<finance::ModelOutput> 
TheoreticalModel::PricePEPSLike(const finance::PEPSLike& pepsLike)
{
  double dFact;
  GetHGCBLicense().Go(&dFact);

  CHECK_COND(!pepsLike.IsExchangeable(), ITO33_HG_EXCHANGEABLE);

  const finance::SessionData& sessionData = *pepsLike.GetSessionData();

  pricing::Mandatory mandatory(pepsLike);

  shared_ptr<numeric::NumParams> pNumParams( GetNumParams(pepsLike) );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::MandatoryParams 
    params(mandatory, sessionData, pNumParams, pMeshParams);

  // create HG model
  Model modelParams(*m_pUnderlyingProcess, m_pUnderlyingProcessForMesh);

  GetHGCBLicense().Check();

  CBPricer pricer(params, modelParams, *m_pComputFlags);

  AutoPtr<CBNumOutput> pNumOutput = pricer.Price(); 

  // activate the following code for step by step debug
  // Don't do it when show convergence is enabled.
  /*
  if ( pNumOutput->GetPriceSurface() )
    pNumOutput->GetPriceSurface()->DumpToFiles();
  */

  shared_ptr<BondLikeOutput> pOutput( pNumOutput->GetBondLikeOutput() );

  pOutput->SetNumOutput( make_ptr( pNumOutput.release() ) );

  return pOutput;
}

static RegisterPriceFunction<finance::PEPSLike> 
  regPEPSLike(&TheoreticalModel::PricePEPSLike);


} // namespace hg

} // namespace ito33
