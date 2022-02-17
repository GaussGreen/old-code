/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/theoriticalmodel_pricecds.cpp
// Purpose:     Implementation of cds pricing using ihg model
// Author:      Wang
// Created:     2004/03/22
// RCS-ID:      $Id: theoreticalmodel_pricecds.cpp,v 1.25 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/cdslike.h"
#include "ito33/finance/cds.h"
#include "ito33/finance/referencecds.h"

#include "ito33/pricing/cds.h"
#include "ito33/pricing/cdsparams.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"

#include "ihg/model.h"
#include "ihg/cdsnumoutput.h"
#include "ihg/cdsnumoutputtimeonly.h"
#include "ihg/cdspricer.h"
#include "ihg/cdspricertimeonly.h"
#include "ihg/license.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/modeloutput.h"

ITO33_FORCE_LINK_THIS_MODULE(IHGPriceCDS);

namespace ito33
{

namespace ihg
{

shared_ptr<finance::ModelOutput>
TheoreticalModel::PriceCDS(const finance::CDS& cds) 
{
  return DoPriceCDSLike(cds);
}


shared_ptr<finance::ModelOutput>
TheoreticalModel::PriceReferenceCDS(const finance::ReferenceCDS& refCDS)
{
  return DoPriceCDSLike(refCDS);
}


shared_ptr<finance::ModelOutput>
TheoreticalModel::DoPriceCDSLike(const finance::CDSLike& cds) 
{
  double dFact;
  GetIHGLicense().Go(&dFact);

  CheckAll();

  const finance::SessionData& sessionData = *cds.GetSessionData();

  pricing::CDS cdsPricing(cds);

  shared_ptr<numeric::NumParams> pNumParams( GetNumParams(cds) );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::CDSParams params(cdsPricing, sessionData, pNumParams, pMeshParams);

  Model modelParams(GetVolatility(), GetVolatilityForMesh(), GetHazardRate());

  GetIHGLicense().Check();

  shared_ptr<ModelOutput> pOutput;

  // Use a special pricer for time only problem
  if ( GetHazardRate()->IsTimeOnly() )
  {
    CDSPricerTimeOnly pricer(params, modelParams, *m_pComputFlags);
    pOutput = pricer.Price()->GetOutput();
  }
  else
  {      
    CDSPricer pricer(params, modelParams, *m_pComputFlags);
    pOutput = pricer.Price()->GetOutput();
  }
  
  pOutput->SetPrice(pOutput->GetPrice() * dFact);
  pOutput->SetDelta(pOutput->GetDelta() * dFact);

  return pOutput;
}


static RegisterPriceFunction<finance::CDS> 
  regCDS(&TheoreticalModel::PriceCDS);

static RegisterPriceFunction<finance::ReferenceCDS> 
  regReferenceCDS(&TheoreticalModel::PriceReferenceCDS);


} // namespace ihg

} // namespace ito33

