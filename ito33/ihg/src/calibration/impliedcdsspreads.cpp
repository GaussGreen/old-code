/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/impliedcdsspreads.cpp
// Purpose:     Class for computing implied cds spreads
// Author:      ITO33
// Created:     2004/11/19
// RCS-ID:      $Id: impliedcdsspreads.cpp,v 1.12 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <vector>
#include <list>
#include "ito33/afterstd.h"

#include "ito33/useexception.h"
#include "ito33/debug.h"
#include "ito33/arraycheckers.h"
#include "ito33/sharedptr.h"

#include "ito33/numeric/predicatetime.h"

#include "ito33/ihg/impliedcdsspreads.h"

#include "ito33/ihg/theoreticalmodel.h"

#include "ito33/finance/cds.h"
#include "ito33/finance/cdslike.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/computationalflags.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"

#include "ito33/pricing/cds.h"
#include "ito33/pricing/forwardcdsparams.h"

#include "ihg/forwardcdsnumoutput.h"
#include "ihg/forwardcdspricer.h"
#include "ihg/model.h"

namespace ito33
{

namespace ihg 
{


std::vector<double> 
ImpliedCDSSpreads::Compute(const TheoreticalModel& theoreticalModel, 
                           const shared_ptr<finance::SessionData>& pSessionData, 
                           const std::vector<Date>& pMaturityDates)
{
  // Verify the data
  CheckIncreasingOrder(pMaturityDates);

  if ( pMaturityDates[0] < m_firstPaymentDate ) 
      throw EXCEPTION_MSG
          (
            ITO33_BAD_DATA,
            TRANS("Implied CDS spread calculator: Maturity dates must be " \
                  "after the first payment date.")
          );

  if ( pSessionData->GetValuationDate() < m_contractingDate ) 
      throw EXCEPTION_MSG
          (
            ITO33_BAD_DATA,
            TRANS("Implied CDS spread calculator: Valuation date must be " \
                  "after the contracting date.")
          );

  if ( pSessionData->GetValuationDate() > pMaturityDates[0] ) 
      throw EXCEPTION_MSG
          (
            ITO33_BAD_DATA,
            TRANS("Implied CDS spread calculator: Valuation date must be " \
                  "before all maturity dates.")
          );

  // Compute using default spread value, then scale results to give 0 price
  double dDefaultSpread = 0.001;
  
  // create the return vector of implied spreads
  size_t nNbDates = pMaturityDates.size();
  std::vector<double> pdImpliedSpreads(nNbDates);

  // If nothing is passed in, return an empty vector.  Could also throw
  // an exception
  if (nNbDates == 0)
    return pdImpliedSpreads;

  // Create a list of cds contracts to match
  // TODO: Convert to cds termstructure
  std::list< shared_ptr<finance::CDSLike> > listOfCDS;
  size_t nIdx;

  for (nIdx = 0; nIdx < nNbDates; nIdx++)
  {
    shared_ptr<finance::CashFlowStreamUniform> pSpreadStream( 
      new finance::CashFlowStreamUniform(m_contractingDate,
                              m_firstPaymentDate, 
                              pMaturityDates[nIdx],
                              dDefaultSpread,
                              m_dcc,
                              m_freq
                              ) );
    shared_ptr<finance::CDSLike> 
      pCDS(new finance::CDS(m_dRecoveryRate, pSpreadStream) );

    pCDS->SetSessionData(pSessionData);

    listOfCDS.push_back(pCDS);

  }

  // Price using forward pricer. This allows us to calculate all the spread
  // values in one pass without any iterations.
  pricing::CDSes cdses(listOfCDS);

  // Sets the yc for mesh to the yc before the compute process.  
  pSessionData->SetYieldCurveForMesh( pSessionData->GetYieldCurve() );

  shared_ptr<numeric::NumParams>
    pNumParams( new numeric::NumParams
                    (
                      *(theoreticalModel.GetQualityControl()),
                      GetDoubleFrom( pMaturityDates[nNbDates-1] ) 
                    - GetDoubleFrom( pSessionData->GetValuationDate() )
                    )
              );

  // use temporarily meshparams
  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);

  pricing::ForwardCDSParams 
    forwardCDSParams(cdses, *pSessionData, pNumParams, pMeshParams);

  Model model(theoreticalModel.GetVolatility(), 
              theoreticalModel.GetVolatility(), 
              theoreticalModel.GetHazardRate());

  // just use default flags, we don't let user pass flags to this function
  finance::ComputationalFlags flags;

  ForwardCDSPricer pricer(forwardCDSParams, model, flags);

  AutoPtr<ForwardCDSNumOutput> numoutput = pricer.Price();
   
  pdImpliedSpreads[nNbDates-1] = numoutput->GetPrice();

  std::vector<double> pdTimes = numoutput->GetTimes();
  std::vector<double> pdRecoveryTerms = numoutput->GetRecoveryTerms();
  std::vector<double> pdAccruedTerms = numoutput->GetAccruedTerms();
  std::vector<double> pdSpreadTerms = numoutput->GetSpreadTerms();
  size_t nNbTimes = pdTimes.size();

  // Match the maturity dates to the numerical output dates and compute
  // the implied spreads
  size_t nIdxTime = 0;
  size_t nIdxDate = 0;
  for (nIdxTime = 0; nIdxTime < nNbTimes; nIdxTime++)
  {
    if ( numeric::AreTimesEqual(GetDoubleFrom(pMaturityDates[nIdxDate]), pdTimes[nIdxTime]) )
    {
      // 0 = recovery - x * (spread  + accrued) / defaultspread
      // x = defaultspread * (recovery ) / (spread + accrued)
      double dImplied = dDefaultSpread * pdRecoveryTerms[nIdxTime]
                      / (pdSpreadTerms[nIdxTime] + pdAccruedTerms[nIdxTime]);

      pdImpliedSpreads[nIdxDate] = dImplied;

      nIdxDate++;
    } // if times equal

    if (nIdxDate == nNbDates)
      break;
  } // loop over maturity dates

  ASSERT_MSG(nIdxDate == nNbDates, "Problem matching dates and times in impliedcdspreads");

  return pdImpliedSpreads;
}


} // namespace ihg

} // namespace ito33
