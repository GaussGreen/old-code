/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/theoreticalmodel_pricebond.cpp
// Purpose:     Implementation of bond pricing using HG model
// Created:     2005/04/11
// RCS-ID:      $Id: theoreticalmodel_pricebond.cpp,v 1.9 2006/08/19 23:44:26 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/bond.h"

#include "ito33/pricing/cb.h"
#include "ito33/pricing/cbparams.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"

#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/bondlikeoutput.h"

#include "hg/model.h"
#include "hg/cbnumoutput.h"
#include "hg/cbpricer.h"
#include "hg/license.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_THIS_MODULE(HGPriceBond);

namespace ito33
{

namespace hg
{


shared_ptr<finance::ModelOutput>
TheoreticalModel::PriceBond(const finance::Bond& bond)
{
  double dFact;
  GetHGCBLicense().Go(&dFact);

  const finance::SessionData& sessionData = *bond.GetSessionData();

  pricing::CB pricingCB(bond);

  shared_ptr<numeric::NumParams> pNumParams( GetNumParams(bond) );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::CBParams params(pricingCB, sessionData, pNumParams, pMeshParams);

  // create HG model
  Model modelParams(*m_pUnderlyingProcess, m_pUnderlyingProcessForMesh);

  GetHGCBLicense().Check();

  CBPricer pricer(params, modelParams, *m_pComputFlags);

  AutoPtr<CBNumOutput> pNumOutput = pricer.Price();

  /*
  if ( pNumOutput->GetPriceSurface() )
    pNumOutput->GetPriceSurface()->DumpToFiles();
  */

  shared_ptr<BondLikeOutput> pOutput = pNumOutput->GetBondLikeOutput();
  
  pOutput->SetNumOutput( shared_ptr<CBNumOutput>(pNumOutput.release()) );

  return pOutput;
}

static RegisterPriceFunction<finance::Bond> 
  regBond(&TheoreticalModel::PriceBond);


} // namespace hg

} // namespace ito33
