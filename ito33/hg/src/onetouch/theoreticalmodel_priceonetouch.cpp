/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/onetouch/theoriticalmodel_priceonetouch.cpp
// Purpose:     Implementation of OneTouch pricing using HG model
// Created:     2005/07/04
// RCS-ID:      $Id: theoreticalmodel_priceonetouch.cpp,v 1.8 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/exoticoption/onetouch.h"
#include "ito33/finance/exoticoption/fxonetouch.h"

#include "ito33/pricing/onetouch.h"
#include "ito33/pricing/onetouchparams.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/surfacegeneral.h"

#include "hg/model.h"
#include "hg/backwardnumoutput.h"
#include "hg/onetouchpricer.h"
#include "hg/license.h"

#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/modeloutput.h"

#include "ito33/link.h"
ITO33_FORCE_LINK_THIS_MODULE(HGPriceOneTouch);

namespace ito33
{

namespace hg
{


  shared_ptr<finance::ModelOutput>
TheoreticalModel::PriceOneTouch(const finance::OneTouch& oneTouch) 
{
  double dFact;
  
  GetHGLicense().Go(&dFact);

  const finance::SessionData& sessionData = *oneTouch.GetSessionData();

  pricing::OneTouch oneTouchPricing(oneTouch);

  shared_ptr<numeric::NumParams> pNumParams( GetNumParams(oneTouch) );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::OneTouchParams 
    params(oneTouchPricing, sessionData, pNumParams, pMeshParams);

  Model modelParams(*m_pUnderlyingProcess, m_pUnderlyingProcessForMesh);

  GetHGLicense().Check();

  OneTouchPricer pricer(params, modelParams, *m_pComputFlags);

  AutoPtr<BackwardNumOutput> pNumOutput = pricer.Price();

  /*
  if ( pNumOutput->GetPriceSurface() )
    pNumOutput->GetPriceSurface()->DumpToFiles();
  */

  shared_ptr<ModelOutput> pOutput( pNumOutput->GetModelOutput() );
  
  pOutput->SetNumOutput( make_ptr(pNumOutput.release()) );

  return pOutput;
}

static RegisterPriceFunction<finance::OneTouch> 
  regOneTouch(&TheoreticalModel::PriceOneTouch);

shared_ptr<finance::ModelOutput>
TheoreticalModel::PriceFXOneTouch(const finance::FXOneTouch& oneTouch)
{
  return PriceOneTouch(oneTouch);
}

static RegisterPriceFunction<finance::FXOneTouch> 
  regFXOneTouch(&TheoreticalModel::PriceFXOneTouch);

} // namespace hg

} // namespace ito33
