/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/onetouch/priceonetouches.cpp
// Purpose:     Implementation of multiple one touch pricing using HG model
// Created:     2006/02/23
// RCS-ID:      $Id: priceonetouches.cpp,v 1.4 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2006  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"
#include "ito33/dateutils.h"

#include "ito33/finance/moneymarket.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/exoticoption/onetouches.h"

#include "ito33/pricing/onetouch.h"
#include "ito33/pricing/onetouchparams.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/surfacegeneral.h"

#include "ito33/hg/multioutput.h"
#include "ito33/hg/theoreticalmodel.h"

#include "hg/license.h"
#include "hg/model.h"
#include "hg/backwardnumoutput.h"
#include "hg/onetouchpricer.h"
#include "hg/priceonetouches.h"

namespace ito33
{

namespace hg
{


shared_ptr<MultiOutput>
PriceOneTouches(const TheoreticalModel& model,
                const finance::OneTouches& oneTouches)
{
  double dFact;

  GetHGLicense().Go(&dFact);

  const finance::OneTouch& oneTouch = *oneTouches.GetOneTouch();

  const finance::SessionData& sessionData = *oneTouch.GetSessionData();

  // Sets the yc for mesh at the yc before the compute process.
  sessionData.SetYieldCurveForMesh( sessionData.GetYieldCurve() );

  pricing::OneTouch contract(oneTouch);

  shared_ptr<numeric::NumParams>
    pNumParams( new numeric::NumParams
                    (
                      *model.GetQualityControl(),
                      GetDoubleFrom( oneTouch.GetMaturityDate() ) 
                    - GetDoubleFrom( sessionData.GetValuationDate() )
                    )
              );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::OneTouchParams 
    params(contract, sessionData, pNumParams, pMeshParams);

  const std::vector<double>&
    pdSpots( oneTouches.GetSpots() ), 
    pdMarkerPrices( oneTouches.GetMarketPrices() );

  params.SetObservations(pdSpots, pdMarkerPrices);

  Model modelParams( *model.GetUnderlyingProcess(),
                     model.GetUnderlyingProcessForMesh() );
  
  GetHGLicense().Check();

  OneTouchPricer 
    pricer(params, modelParams, *model.GetComputationalFlags());

  AutoPtr<BackwardNumOutput> pNumOutput = pricer.Price();

  shared_ptr<MultiOutput> pOutput(pNumOutput->GetMultiOutput());

  pOutput->SetNumOutput( shared_ptr<BackwardNumOutput>(pNumOutput.release()) );

  return pOutput;
}


} // namespace hg

} // namespace ito33
