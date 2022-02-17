/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/varianceswap/theoreticalmodel_pricelogcontract.cpp
// Purpose:     Implementation of log contract pricing using HG model
// Created:     2006/07/18
// RCS-ID:      $Id: theoreticalmodel_pricelogcontract.cpp,v 1.3 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2006  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/logcontract.h"

#include "ito33/pricing/logcontract.h"
#include "ito33/pricing/logcontractparams.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"

#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/modeloutput.h"

#include "hg/model.h"
#include "hg/numoutputtimeonly.h"
#include "hg/logcontract_closedform.h"
#include "hg/license.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_THIS_MODULE(HGPriceLogContract);


namespace ito33
{

namespace hg
{


shared_ptr<finance::ModelOutput>
TheoreticalModel::PriceLogContract(const finance::LogContract& derivative)
{
  double dFact;

  GetHGLicense().Go(&dFact);

  const finance::SessionData& sessionData = *derivative.GetSessionData();

  pricing::LogContract logContract(derivative);

  shared_ptr<numeric::NumParams> pNumParams( GetNumParams(derivative) );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::LogContractParams 
    params(logContract, sessionData, pNumParams, pMeshParams);

  Model modelParams(*m_pUnderlyingProcess, m_pUnderlyingProcessForMesh);

  GetHGLicense().Check();

  finance::ComputationalFlags flags(*m_pComputFlags);
  if (    flags.GetAnalysisDate().IsValid()
       && flags.GetAnalysisDate() > derivative.GetStartOfReturnPeriod()
       && flags.GetAnalysisDate() != sessionData.GetValuationDate() )
     flags.SetAnalysisDate(Date());

  LogContractClosedForm closedForm(params, modelParams, flags);

  AutoPtr<NumOutputTimeOnly> pNumOutput( closedForm.Price() );
  shared_ptr<ModelOutput> pOutput( pNumOutput->GetModelOutput() );
  pOutput->SetPrice(pNumOutput->GetFinalPrices()[0]);

  pOutput->SetNumOutput( make_ptr( pNumOutput.release() ) );

  return pOutput;
}

static RegisterPriceFunction<finance::LogContract>
  regLogContract(&TheoreticalModel::PriceLogContract);


} // namespace hg

} // namespace ito33
