/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/theoreticalmodel_priceasianoption.cpp
// Purpose:     Implementation of Asian option pricing using HG model
// Created:     2006/03/02
// RCS-ID:      $Id: theoreticalmodel_priceasianoption.cpp,v 1.3 2006/08/19 23:44:26 wang Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/exoticoption/asianoption.h"

#include "ito33/pricing/asianoption.h"
#include "ito33/pricing/asianoptionparams.h"

#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/numparams.h"

#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/modeloutput.h"

#include "hg/model.h"
#include "hg/optionnumoutput.h"
#include "hg/asianoptionpricer.h"
#include "hg/license.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_THIS_MODULE(HGPriceAsianOption);


namespace ito33
{

namespace hg
{


shared_ptr<finance::ModelOutput>
TheoreticalModel::PriceAsianOption(const finance::AsianOption& asianOption)
{
  double dFact;

  GetHGLicense().Go(&dFact);

  const finance::SessionData& sessionData = *asianOption.GetSessionData();

  pricing::AsianOption asianOpt(asianOption);

  shared_ptr<numeric::NumParams> pNumParams( GetNumParams(asianOption) );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::AsianOptionParams 
    params(asianOpt, sessionData, pNumParams, pMeshParams);

  Model modelParams(*m_pUnderlyingProcess, m_pUnderlyingProcessForMesh);
  
  GetHGLicense().Check();

  AsianOptionPricer pricer(params, modelParams, *m_pComputFlags);

  AutoPtr<OptionNumOutput> pNumOutput = pricer.Price();
  
  shared_ptr<ModelOutput> pOutput = pNumOutput->GetModelOutput();
  
  pOutput->SetNumOutput( shared_ptr<OptionNumOutput>(pNumOutput.release())  );
  
  pOutput->SetPrice(pOutput->GetPrice() * dFact);
  pOutput->SetDelta(pOutput->GetDelta() * dFact);

  return pOutput;
}

static RegisterPriceFunction<finance::AsianOption>
  regAsianOption(&TheoreticalModel::PriceAsianOption);


} // namespace hg

} // namespace ito33
