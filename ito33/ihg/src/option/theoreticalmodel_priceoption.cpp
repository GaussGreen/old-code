/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/theoreticalmodel_priceoption.cpp
// Purpose:     Implementation of option pricing using ihg model
// Author:      ZHANG Yunzhi
// Created:     2004/02/26
// RCS-ID:      $Id: theoreticalmodel_priceoption.cpp,v 1.31 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/option.h"

#include "ito33/pricing/option.h"
#include "ito33/pricing/optionparams.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"

#include "ihg/model.h"
#include "ihg/optionnumoutput.h"
#include "ihg/optionpricer.h"
#include "ihg/license.h"

#include "ito33/ihg/modeloutput.h"

#include "ito33/ihg/theoreticalmodel.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_THIS_MODULE(IHGPriceOption);


namespace ito33
{

namespace ihg
{


shared_ptr<finance::ModelOutput>
TheoreticalModel::PriceOption(const finance::Option& option)
{
  double dFact;
  GetIHGLicense().Go(&dFact);

  CheckAll();

  const finance::SessionData& sessionData = *option.GetSessionData();

  pricing::Option opt(option);

  shared_ptr<numeric::NumParams> pNumParams( GetNumParams(option) );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::OptionParams 
    params(opt, sessionData, pNumParams, pMeshParams);

  Model modelParams(GetVolatility(), GetVolatilityForMesh(), GetHazardRate());
  
  GetIHGLicense().Check();

  OptionPricer pricer(params, modelParams, *m_pComputFlags);

  AutoPtr<OptionNumOutput> pOptionNumOutput = pricer.Price();

  shared_ptr<ModelOutput> pOutput = pOptionNumOutput->GetModelOutput();
      
  pOutput->SetNumOutput( shared_ptr<OptionNumOutput>(pOptionNumOutput.release()) );

  pOutput->SetPrice(pOutput->GetPrice() * dFact);
  pOutput->SetDelta(pOutput->GetDelta() * dFact);

  return pOutput;
}

static RegisterPriceFunction<finance::Option>
  regOption(&TheoreticalModel::PriceOption);


} // namespace ihg

} // namespace ito33

