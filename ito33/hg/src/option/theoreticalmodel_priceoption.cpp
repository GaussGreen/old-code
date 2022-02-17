/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/theoriticalmodel_priceoption.cpp
// Purpose:     Implementation of option pricing using HG model
// Created:     2005/01/13
// RCS-ID:      $Id: theoreticalmodel_priceoption.cpp,v 1.10 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/option.h"

#include "ito33/pricing/option.h"
#include "ito33/pricing/optionparams.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/surfacegeneral.h"

#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/modeloutput.h"

#include "hg/model.h"
#include "hg/optionnumoutput.h"
#include "hg/optionpricer.h"
#include "hg/license.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_THIS_MODULE(HGPriceOption);


namespace ito33
{

namespace hg
{


shared_ptr<finance::ModelOutput>
TheoreticalModel::PriceOption(const finance::Option& option)
{
  double dFact;

  GetHGLicense().Go(&dFact);

  const finance::SessionData& sessionData = *option.GetSessionData();

  pricing::Option opt(option);

  shared_ptr<numeric::NumParams> pNumParams( GetNumParams(option) );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::OptionParams 
    params(opt, sessionData, pNumParams, pMeshParams);

  Model modelParams(*m_pUnderlyingProcess, m_pUnderlyingProcessForMesh);
  
  GetHGLicense().Check();

  OptionPricer pricer(params, modelParams, *m_pComputFlags);

  AutoPtr<OptionNumOutput> pNumOutput = pricer.Price();
  
  /*
  if ( pNumOutput->GetPriceSurface() )
    pNumOutput->GetPriceSurface()->DumpToFiles();
  */

  shared_ptr<ModelOutput> pOutput = pNumOutput->GetModelOutput();
  
  pOutput->SetNumOutput( shared_ptr<OptionNumOutput>(pNumOutput.release())  );
  
  pOutput->SetPrice(pOutput->GetPrice() * dFact);
  pOutput->SetDelta(pOutput->GetDelta() * dFact);

  return pOutput;
}

static RegisterPriceFunction<finance::Option>
  regOption(&TheoreticalModel::PriceOption);


} // namespace hg

} // namespace ito33
