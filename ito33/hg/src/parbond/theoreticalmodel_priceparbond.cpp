/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/theoriticalmodel_priceparbond.cpp
// Purpose:     Implementation of ParBond pricing using HG model
// Created:     2005/06/09
// RCS-ID:      $Id: theoreticalmodel_priceparbond.cpp,v 1.6 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/parbond.h"

#include "ito33/pricing/parbond.h"
#include "ito33/pricing/parbondparams.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/surfacegeneral.h"

#include "hg/model.h"
#include "hg/numoutputtimeonly.h"
#include "hg/parbondpricer.h"
#include "hg/license.h"

#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/modeloutput.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_THIS_MODULE(HGPriceParBond);

namespace ito33
{

namespace hg
{


shared_ptr<finance::ModelOutput>
TheoreticalModel::PriceParBond(const finance::ParBond& parBond) 
{
  double dFact;
  
  GetHGLicense().Go(&dFact);

  const finance::SessionData& sessionData = *parBond.GetSessionData();

  pricing::ParBond contract(parBond);

  shared_ptr<numeric::NumParams> pNumParams( GetNumParams(parBond) );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::ParBondParams 
    params(contract, sessionData, pNumParams, pMeshParams);

  Model modelParams(*m_pUnderlyingProcess, m_pUnderlyingProcessForMesh);

  GetHGLicense().Check();

  ParBondPricer pricer(params, modelParams, *m_pComputFlags);

  AutoPtr<NumOutputTimeOnly> pNumOutput = pricer.Price();

  /*
  if ( pNumOutput->GetPriceSurface() )
    pNumOutput->GetPriceSurface()->DumpToFiles();
  */

  shared_ptr<ModelOutput> pOutput( pNumOutput->GetModelOutput() );
  
  pOutput->SetNumOutput( make_ptr(pNumOutput.release()) );

  return pOutput;
}

static RegisterPriceFunction<finance::ParBond> 
  regParBond(&TheoreticalModel::PriceParBond);


} // namespace hg

} // namespace ito33
