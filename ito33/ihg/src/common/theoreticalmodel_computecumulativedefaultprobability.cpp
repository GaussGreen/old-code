/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/theoriticalmodel_computecumulativedefaultprobability.cpp
// Purpose:     Implementation of cumulative default probability computation
// Author:      Wang
// Created:     2004/07/30
// RCS-ID:      $Id: theoreticalmodel_computecumulativedefaultprobability.cpp,v 1.18 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"
#include "ito33/dateutils.h"
#include "ito33/arraycheckers.h"

#include "ito33/finance/error.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/option.h"
#include "ito33/finance/forwardoption.h"
#include "ito33/finance/yieldcurve.h"

#include "ito33/pricing/options.h"
#include "ito33/pricing/forwardoptionparams.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/predicatetime.h"
#include "ito33/numeric/surfacegeneral.h"
#include "ito33/numeric/deltagamma.h"
#include "ito33/numeric/interpolation.h"

#include "ihg/model.h"
#include "ihg/forwardoptionnumoutput.h"
#include "ihg/forwardoptionpricer.h"
#include "ihg/license.h"

#include "ito33/ihg/theoreticalmodel.h"

extern const ito33::finance::Error ITO33_MATURITYBEFOREVALUATION;

namespace ito33
{

using namespace numeric;

namespace ihg
{


std::vector<double>
TheoreticalModel::ComputeCumulativeDefaultProbability
(const finance::SessionData& sessionData, const std::vector<Date>& maturityDates)
{
  double dFact;
  GetIHGLicense().Go(&dFact);
            
  CheckAll();

  CheckIncreasingOrder(maturityDates, "maturity dates");

  CHECK_COND(sessionData.GetValuationDate() < maturityDates[0],
             ITO33_MATURITYBEFOREVALUATION);

  std::list< shared_ptr<finance::Option> > optionList;
  
  for (size_t nIdx = 0; nIdx < maturityDates.size(); nIdx++)
  {
    optionList.push_back(shared_ptr<finance::Option>
                         (
                            new finance::Option(sessionData.GetSpotSharePrice(), 
                                                maturityDates[nIdx],
                                                finance::Option_Call, 
                                                finance::ExerciseType_European)
                         ));
  }
  
  finance::ForwardOption forwardOption(optionList);
  pricing::Options options(forwardOption);

  // Sets the yc for mesh to the yc before the compute process.  
  sessionData.SetYieldCurveForMesh( sessionData.GetYieldCurve() );

  shared_ptr<NumParams>
    pNumParams( new NumParams
                    (
                      *m_pQualityControl,
                      GetDoubleFrom( maturityDates.back() ) 
                    - GetDoubleFrom( sessionData.GetValuationDate() )
                    )
              );
  
  // use temporarily meshparams
  shared_ptr<MeshParams> pMeshParams(new MeshParams);
  
  pricing::ForwardOptionParams 
    params(options, sessionData, pNumParams, pMeshParams);

  Model modelParams(GetVolatility(), GetVolatility(), GetHazardRate());
  
  GetIHGLicense().Check();

  // use a temporary flags, leave the m_pFlags untouched
  finance::ComputationalFlags flags;

  // For the moment, compute the whole surface, need to be changed
  // to store only those values at the maturityDates
  flags.SetComputeSurface(true);

  ForwardOptionPricer pricer(params, modelParams, flags);

  AutoPtr<ForwardOptionNumOutput> pNumOutput = pricer.Price();

  shared_ptr<SurfaceGeneral> pSurface = pNumOutput->GetPriceSurface();
  shared_ptr<Domain> pDomain = pSurface->GetDomain();
  
  const std::vector<double> pdTimes = pDomain->GetTimes();
  shared_ptr<finance::YieldCurve> pYieldCurve = sessionData.GetYieldCurve();
  double dValuationDate = GetDoubleFrom( sessionData.GetValuationDate() );

  size_t nIdxT = 0;

  std::vector<double> pdCDPs;
  
  for (size_t nIdx = 0; nIdx < maturityDates.size(); nIdx++)
  {
    double dTime = GetDoubleFrom(maturityDates[nIdx]);
    
    while ( !AreTimesEqual(dTime, pdTimes[nIdxT]) )
      nIdxT++;

    // FIXME:  FirstValues? or LastValues
    const std::vector<double>& values = pSurface->GetFirstValuesAt(nIdxT);
    const std::vector<double>& spots = pDomain->GetFirstSpaceMeshAt(nIdxT);

    // comupute the delta
    Array<double> pdDeltas( values.size() );
    
    ComputeDelta( &spots[0], &values[0], values.size(), pdDeltas.Get() );
    
    double dS = 0, dGamma0;

    Interpolate(&spots[0], pdDeltas.Get(), values.size(), &dS, &dGamma0, 1,
                ExtrapolationMode_Linear, ExtrapolationMode_Linear);
                     
    double 
      dTmp = pYieldCurve->GetForwardDiscountFactor(dValuationDate, dTime);

    double dCDP = 1. + dGamma0 / dTmp;

    // If hazard rate is set to zero, numerical error can make the value
    // tiny and negative.  Artificially correct.
    if ( dCDP < 0.0 )
      dCDP = 0.0;

    pdCDPs.push_back( dCDP );
  }
  
  return pdCDPs;
}


} // namespace ihg

} // namespace ito33

