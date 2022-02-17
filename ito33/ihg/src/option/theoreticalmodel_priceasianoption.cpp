/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/theoriticalmodel_priceoption.cpp
// Purpose:     Implementation of option pricing using ihg model
// Author:      ZHANG Yunzhi
// Created:     2004/02/26
// RCS-ID:      $Id: theoreticalmodel_priceasianoption.cpp,v 1.7 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"
#include "ito33/dateutils.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/exoticoption/asianoption.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/pricing/asianoption.h"
#include "ito33/pricing/asianoptionparams.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"

#include "ihg/model.h"
#include "ihg/optionnumoutput.h"
#include "ihg/asianoptionpricer.h"
#include "ihg/license.h"


#include "ito33/ihg/theoreticalmodel.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_THIS_MODULE(IHGPriceAsianOption);


namespace ito33
{

namespace ihg
{


shared_ptr<finance::ModelOutput>
TheoreticalModel::PriceAsianOption(const finance::AsianOption& asianOption)
{
  double dFact;
  GetIHGLicense().Go(&dFact);

  CheckAll();

  const finance::SessionData& sessionData = *asianOption.GetSessionData();

  pricing::AsianOption opt(asianOption);

  shared_ptr<numeric::NumParams>
    pNumParams( new numeric::NumParams
                    (
                      *m_pQualityControl,
                      GetDoubleFrom( asianOption.GetMaturityDate() ) 
                    - GetDoubleFrom( sessionData.GetValuationDate() )
                    )
              );

  //use less points to make the code run faster
  //accuracy does not seem to be that much affected
  pNumParams->SetNbSpaceSteps( pNumParams->GetNbSpaceSteps()*1/2);

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::AsianOptionParams   
    params(opt, sessionData, pNumParams, pMeshParams);

  Model modelParams(GetVolatility(), GetVolatilityForMesh(), GetHazardRate());
  
  GetIHGLicense().Check();

  AsianOptionPricer pricer(params, modelParams, *m_pComputFlags);

  shared_ptr<ihg::ModelOutput> pOutput = pricer.Price()->GetModelOutput();

  pOutput->SetPrice(pOutput->GetPrice() * dFact);
  pOutput->SetDelta(pOutput->GetDelta() * dFact);

  return pOutput;
}

static RegisterPriceFunction<finance::AsianOption>
  regAsianOption(&ihg::TheoreticalModel::PriceAsianOption);


} // namespace ihg

} // namespace ito33

