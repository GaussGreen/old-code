/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/theoreticalmodel_pricecb.cpp
// Purpose:     Implementation of cb pricing using HG model
// Created:     2005/04/11
// RCS-ID:      $Id: theoreticalmodel_pricecb.cpp,v 1.9 2006/08/19 23:44:26 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/link.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/moneymarket.h"
#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/convertiblebond.h"

#include "ito33/pricing/cb.h"
#include "ito33/pricing/cbparams.h"

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

ITO33_FORCE_LINK_THIS_MODULE(HGPriceCB);

namespace ito33
{

namespace hg
{


shared_ptr<finance::ModelOutput>
TheoreticalModel::PriceConvertibleBond(const finance::ConvertibleBond& cb)
{
  double dFact;
  GetHGCBLicense().Go(&dFact);

  CHECK_COND(!cb.IsExchangeable(), ITO33_HG_EXCHANGEABLE);

  const finance::SessionData& sessionData = *cb.GetSessionData();

  pricing::CB pricingCB(cb);

  shared_ptr<numeric::NumParams> pNumParams( GetNumParams(cb) );

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

  shared_ptr<BondLikeOutput> pOutput( pNumOutput->GetBondLikeOutput() );
  
  pOutput->SetNumOutput( make_ptr(pNumOutput.release()) );

  return pOutput;
}

static RegisterPriceFunction<finance::ConvertibleBond> 
  regConvertibleBond(&TheoreticalModel::PriceConvertibleBond);


} // namespace hg

} // namespace ito33
