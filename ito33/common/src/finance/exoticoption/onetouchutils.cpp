/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/exoticoption/onetouchutils.cpp
// Purpose:     Helper functions for one touch
// Created:     2005/12/08
// RCS-ID:      $Id: onetouchutils.cpp,v 1.7 2006/08/19 22:40:30 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/exoticoption/onetouch.h"
#include "ito33/finance/exoticoption/onetouchutils.h"

#include "ito33/numeric/nonlinearsolver.h"

#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardrateflat.h"
#include "ito33/ihg/underlyingprocess.h"
#include "ito33/ihg/theoreticalmodel.h"

#include "ito33/link.h"
ITO33_FORCE_LINK_MODULE(IHGPriceOneTouch);

namespace ito33
{

namespace finance
{


class RegulaFalsiForOneTouch
{
public:

  RegulaFalsiForOneTouch(const ihg::TheoreticalModel& model, OneTouch& oneTouch)
                  : m_model(model), m_oneTouch(oneTouch)
  {
  }

  double operator()(double dBarrier)
  {
    m_oneTouch.SetBarrier(dBarrier);

    return m_model.Compute(m_oneTouch)->GetPrice() - m_oneTouch.GetMarketPrice();
  }

private:

  NO_COPY_CLASS(RegulaFalsiForOneTouch);

  const ihg::TheoreticalModel& m_model;
 
  OneTouch& m_oneTouch;
};

double 
GetBarrierFromDelta
(const shared_ptr<SessionData>& pSessionData, 
 Date maturityDate, double dBSBarrier, BarrierType barrierType, double dVol)
{
  shared_ptr<ihg::Volatility> vol( new ihg::VolatilityFlat(dVol) );
  shared_ptr<ihg::HazardRate> hr( new ihg::HazardRateFlat(0) );
  shared_ptr<ihg::UnderlyingProcess> 
    pUnderlyingProcess( new ihg::UnderlyingProcess(vol, hr) );

  ihg::TheoreticalModel model(pUnderlyingProcess);

  double dMin, dMax;

  if ( barrierType == Barrier_UpAndOut )
  {
    dMin = pSessionData->GetSpotSharePrice();
    dMax = 1.e10;
  }
  else
  {
    dMin = 1.e-10;
    dMax = pSessionData->GetSpotSharePrice();
  }

  OneTouch oneTouch(maturityDate, 1, barrierType, Rebate_Immediate);
  oneTouch.SetSessionData(pSessionData);
  oneTouch.SetMarketPrice(dBSBarrier);

  numeric::RegulaFalsi solver;
    
  RegulaFalsiForOneTouch oneTouchSolving(model, oneTouch);

  try
  {
    return solver(oneTouchSolving, dMin, dMax);
  }
  catch(const ito33::numeric::Exception&)
  {
    return 0.;
  }
}


} // namespace finance

} // namespace ito33
