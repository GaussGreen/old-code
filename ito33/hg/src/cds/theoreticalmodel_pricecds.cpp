/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/theoriticalmodel_pricecds.cpp
// Purpose:     Implementation of CDS pricing using HG model
// Created:     2005/02/16
// RCS-ID:      $Id: theoreticalmodel_pricecds.cpp,v 1.11 2006/08/19 23:44:26 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/cdslike.h"
#include "ito33/finance/cds.h"
#include "ito33/finance/referencecds.h"

#include "ito33/pricing/cds.h"
#include "ito33/pricing/cdsparams.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/surfacegeneral.h"

#include "hg/model.h"
#include "hg/numoutputtimeonly.h"
#include "hg/cdspricer.h"
#include "hg/license.h"

#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/modeloutput.h"

#include "ito33/link.h"
ITO33_FORCE_LINK_THIS_MODULE(HGPriceCDS);

namespace ito33
{

namespace hg
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
  
  GetHGLicense().Go(&dFact);

  const finance::SessionData& sessionData = *cds.GetSessionData();

  pricing::CDS cdsPricing(cds);

  shared_ptr<numeric::NumParams> pNumParams( GetNumParams(cds) );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::CDSParams params(cdsPricing, sessionData, pNumParams, pMeshParams);

  Model modelParams(*m_pUnderlyingProcess, m_pUnderlyingProcessForMesh);

  GetHGLicense().Check();

  CDSPricer pricer(params, modelParams, *m_pComputFlags);

  AutoPtr<NumOutputTimeOnly> pNumOutput = pricer.Price();

  /*
  if ( pNumOutput->GetPriceSurface() )
    pNumOutput->GetPriceSurface()->DumpToFiles();
  */

  shared_ptr<ModelOutput> pOutput( pNumOutput->GetModelOutput() );
  
  pOutput->SetNumOutput( make_ptr(pNumOutput.release())  );

  return pOutput;
}

static RegisterPriceFunction<finance::CDS> 
  regCDS(&TheoreticalModel::PriceCDS);

static RegisterPriceFunction<finance::ReferenceCDS> 
  regReferenceCDS(&TheoreticalModel::PriceReferenceCDS);


} // namespace hg

} // namespace ito33
