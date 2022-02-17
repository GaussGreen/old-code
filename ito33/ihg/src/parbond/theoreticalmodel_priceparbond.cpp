/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/theoriticalmodel_priceparbond.cpp
// Purpose:     Implementation of parbond pricing using ihg model
// Author:      ZHANG
// Created:     2005/05/20
// RCS-ID:      $Id: theoreticalmodel_priceparbond.cpp,v 1.7 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/parbond.h"

#include "ito33/pricing/parbond.h"
#include "ito33/pricing/parbondparams.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"

#include "ihg/model.h"
#include "ihg/parbondnumoutput.h"
#include "ihg/parbondnumoutputtimeonly.h"
#include "ihg/parbondpricer.h"
#include "ihg/parbondpricertimeonly.h"
#include "ihg/license.h"

#include "ito33/ihg/modeloutput.h"
#include "ito33/ihg/theoreticalmodel.h"

ITO33_FORCE_LINK_THIS_MODULE(IHGPriceParBond);

namespace ito33
{

namespace ihg
{


shared_ptr<finance::ModelOutput>
TheoreticalModel::PriceParBond(const finance::ParBond& parbond) 
{
  double dFact;
  GetIHGLicense().Go(&dFact);

  CheckAll();

  const finance::SessionData& sessionData = *parbond.GetSessionData();

  pricing::ParBond parbondPricing(parbond);

  shared_ptr<numeric::NumParams> pNumParams( GetNumParams(parbond) );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::ParBondParams 
    params(parbondPricing, sessionData, pNumParams, pMeshParams);

  Model modelParams(GetVolatility(), GetVolatilityForMesh(), GetHazardRate());

  GetIHGLicense().Check();
  
  shared_ptr<ModelOutput> pOutput;

  // Use a special pricer for time only problem
  if ( GetHazardRate()->IsTimeOnly() )
  {
    ParBondPricerTimeOnly pricer(params, modelParams, *m_pComputFlags);

    pOutput = pricer.Price()->GetOutput();
  }
  else
  {      
    ParBondPricer pricer(params, modelParams, *m_pComputFlags);

    pOutput = pricer.Price()->GetOutput();
  }
  
  pOutput->SetPrice(pOutput->GetPrice() * dFact);
  pOutput->SetDelta(pOutput->GetDelta() * dFact);

  return pOutput;
}

static RegisterPriceFunction<finance::ParBond> 
  regParBond(&TheoreticalModel::PriceParBond);


} // namespace ihg

} // namespace ito33

