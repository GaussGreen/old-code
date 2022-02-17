/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/forwardoption/priceforwardoption.cpp
// Purpose:     Implementation of forward option pricing using HG model
// Created:     2005/05/05
// RCS-ID:      $Id: priceforwardoption.cpp,v 1.10 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"
#include "ito33/dateutils.h"

#include "ito33/finance/moneymarket.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/forwardoption.h"

#include "ito33/pricing/options.h"
#include "ito33/pricing/forwardoptionparams.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/surfacegeneral.h"

#include "ito33/hg/multioutput.h"
#include "ito33/hg/theoreticalmodel.h"

#include "hg/license.h"
#include "hg/model.h"
#include "hg/forwardoptionnumoutput.h"
#include "hg/forwardoptionpricer.h"
#include "hg/priceforwardoption.h"

namespace ito33
{

namespace hg
{


shared_ptr<MultiOutput>
PriceForwardOption(const TheoreticalModel& model,
                   const finance::ForwardOption& forwardOption)
{
  double dFact;

  GetHGLicense().Go(&dFact);

  const finance::SessionData& sessionData = *forwardOption.GetSessionData();

  // Sets the yc for mesh at the yc before the compute process.
  sessionData.SetYieldCurveForMesh( sessionData.GetYieldCurve() );

  pricing::Options options(forwardOption);

  shared_ptr<numeric::NumParams>
    pNumParams( new numeric::NumParams
                    (
                      *model.GetQualityControl(),
                      GetDoubleFrom( forwardOption.GetMaturityDate() ) 
                    - GetDoubleFrom( sessionData.GetValuationDate() )
                    )
              );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::ForwardOptionParams 
    params(options, *forwardOption.GetSessionData(), pNumParams, pMeshParams);

  Model modelParams( *model.GetUnderlyingProcess(),
                     model.GetUnderlyingProcessForMesh() );
  
  GetHGLicense().Check();

  ForwardOptionPricer 
    pricer(params, modelParams, *model.GetComputationalFlags());

  AutoPtr<ForwardOptionNumOutput> pNumOutput = pricer.Price();

  shared_ptr<MultiOutput> pOutput = pNumOutput->GetMultiOutput();

  pOutput->SetNumOutput( shared_ptr<ForwardOptionNumOutput>(pNumOutput.release()) );

  return pOutput;
}


} // namespace hg

} // namespace ito33
