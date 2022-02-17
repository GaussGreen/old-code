/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/forward/forwardoptionnumoutput.cpp
// Purpose:     implementation of NumOutput for forward options
// Created:     2005/05/05
// RCS-ID:      $Id: forwardoptionnumoutput.cpp,v 1.21 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/*
  @todo We suppose implicitly that the options are well ordered
        according to at first the maturity and then the strike
*/

#include "ito33/autoptr.h"
#include "ito33/dateutils.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/yieldcurve.h"
#include "ito33/finance/optiontype.h"

#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/interpolationmatrix.h"

#include "ito33/pricing/forwardoptionparams.h"
#include "ito33/pricing/options.h"

#include "ito33/hg/multioutput.h"

#include "hg/forwardoptionnumoutput.h"
#include "hg/forwardoptioninstdata.h"
#include "hg/model.h"
#include "hg/sensitivitymethod.h"

// implement the AutoPtrDeleter for ForwardOptionNumOutput
namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(hg::ForwardOptionNumOutput);
}

namespace ito33
{

  using namespace numeric;
  using namespace pricing;

namespace hg
{


void ForwardOptionNumOutput::Init(ForwardOptionInstData& instdata)
{  
  BackwardNumOutput::Init(instdata);

  // allocate space for the price vector for the underlying options
  Options& options = m_params.GetOptions();
  size_t nNbOptions = options.GetStrikes().size();
  m_pdPricesAtObservations.resize( nNbOptions, 0.0 );

  if (   m_computationalFlags.HasSensitivityFlags() 
      && instdata.m_sensitivityMethod != SensitivityMethod_None )
  {
    m_ppdSensitivities.resize( nNbOptions );

    size_t nNbParams = instdata.m_pSensitivityData.size();
    
    for (size_t nIdx = 0; nIdx < nNbOptions; nIdx++)
    {
      m_ppdSensitivities[nIdx].resize(nNbParams, 0.0);
    }
  } // if computing sensitivities

  if ( instdata.m_sensitivityMethod == SensitivityMethod_Adjoint )
    m_bSensitivityOnObjectif = true;

  m_dObjectif = 0;
}

void 
ForwardOptionNumOutput::CalculateFinalScalarResult
(BackwardInstData& /* instdata */)
{
}

// Return the requested data to the user
shared_ptr<MultiOutput> ForwardOptionNumOutput::GetMultiOutput()
{
  // Construct the output class
  shared_ptr<MultiOutput> pOutput( BackwardNumOutput::GetMultiOutput() );

  //  Compute sensitivities
  if ( m_ppdSensitivities.size() > 0 && !m_bSensitivityOnObjectif)
    pOutput->SetMultiSensitivities(m_ppdSensitivities);

  return pOutput;
}

void ForwardOptionNumOutput::UpdateMe(BackwardInstData& instdata, double dTime)
{
  // TODO: Make this more efficient

  // Let the base class process the data first
  BackwardNumOutput::UpdateMe(instdata, dTime);

  // Save data if the current time matches one of the maturity dates in the
  // option list  
  Options& options = m_params.GetOptions();

  const std::vector<Date>& pDates           = options.GetMaturityDates();
  const std::vector<double>& pdStrikes      = options.GetStrikes();    
  const std::vector<double>& pdWeights      = options.GetWeights();       
  const std::vector<double>& pdMarketPrices = options.GetMarketPrices(); 
  const std::vector<finance::OptionType>& 
    pOptionTypes = options.GetOptionTypes();

  const size_t nNbOptions = pDates.size();
  const size_t nNbRegimes = instdata.GetModel().GetNbRegimes();

  std::vector<double> pdStrikesNow;
  size_t nIdxStart = INVALIDINDEX;

  // Loop through the option list, checking if any maturity
  // matches the current time.  Also check for puts, which
  // require extra processing
  bool bHasPut = false;
  for (size_t nIdxOp = 0; nIdxOp < nNbOptions; nIdxOp++)
  {
    if ( AreTimesEqual(dTime, GetDoubleFrom(pDates[nIdxOp])) )
    {
      if (nIdxStart == INVALIDINDEX)
        nIdxStart = nIdxOp;

      pdStrikesNow.push_back(pdStrikes[nIdxOp]);

      if ( pOptionTypes[nIdxOp] == finance::Option_Put )
        bHasPut = true;
    }
  }

  // Return if no maturity at current time
  size_t nNbOpNow = pdStrikesNow.size();
  if ( nNbOpNow == 0 )
    return;

  // Interpolation matrix to compute prices at the given strikes
  AutoPtr<InterpolationMatrix> pMatrix(new InterpolationMatrix
    ( &pdStrikesNow[0], nNbOpNow, m_pdS, m_nNbS, nNbRegimes, false ) );

  // If at least one put is included at current time, update the interpolation
  // matrix for the zero-strike call, and compute the strike discount term.
  // Recall that P = C - C(K=0) + K exp(-r*T)
  double dYieldRate = 0.0;
  if ( bHasPut )
  {
    // Compute the discounting factor for exp(-r*T)
    dYieldRate = m_params.GetYieldCurve()->GetForwardDiscountFactor( 
                 m_params.GetValuationTime(), dTime); 

    // Compute the extrapolation weight for zero strike call C(K=0)
    double dK1 = instdata.m_pdS[0];
    double dK2 = instdata.m_pdS[1];
    double dWeight = -1. / (dK2/dK1 - 1.);

    // Update the interpolation matrix to automatically include C(K=0)
    // (ie. subtract the contribution of C(K=0), needed by price and gradient)
    for ( size_t nIdx = 0; nIdx < nNbOpNow; nIdx++ )
    {
      if ( pOptionTypes[nIdx + nIdxStart] == finance::Option_Put )
      {
        pMatrix->Add(nIdx, 0 , -(1 - dWeight) );
        pMatrix->Add(nIdx, 1 , -dWeight);
      }
    }
  } // if has put

  // Set the option prices by interpolating. Automatically subtracts
  // C(K=0) for puts
  pMatrix->ProductMatrixVector(instdata.m_pdPrices.Get(), 
                               &m_pdPricesAtObservations[nIdxStart]);

  Array<double> pdResiduals;
  if ( pdWeights.size() )
  {
    ASSERT( pdWeights.size() == nNbOptions );
    pdResiduals = Array<double>(nNbOpNow);
  }

  // Compute the objective.  Also update put prices with the discount 
  // term (gradient of this term is zero, so the interpolation matrix is OK)
  for ( size_t nIdx = 0; nIdx < nNbOpNow; nIdx++ )
  {
    size_t nIdxOp = nIdx + nIdxStart;

    // update puts
    if ( pOptionTypes[nIdxOp] == finance::Option_Put )
    {
      double dStrike = pdStrikes[nIdxOp];
      m_pdPricesAtObservations[nIdxOp] += dStrike * dYieldRate;
    }

    // block to compute the objective   
    if ( pdWeights.size() )
    {
      double dError = m_pdPricesAtObservations[nIdxOp] 
                    - pdMarketPrices[nIdxOp];
      pdResiduals[nIdx] = pdWeights[nIdxOp] * dError;

      m_dObjectif += 0.5 * pdResiduals[nIdx] * dError; 
    } //end update objective function
    
  } //end loop over options

  size_t nNbSensitivities = instdata.m_pSensitivityData.size();

  if ( nNbSensitivities 
    && instdata.m_sensitivityMethod == SensitivityMethod_PDE )
  {
    Array<double> pdSensitivities(nNbOpNow);

    for (size_t nIdxD = 0; nIdxD < nNbSensitivities; nIdxD++)
    {
      pMatrix->ProductMatrixVector( &instdata.m_ppdSensitivities[nIdxD][0],
                                    pdSensitivities.Get() );

      for (size_t nIdx = 0; nIdx < nNbOpNow; nIdx++)
        m_ppdSensitivities[nIdx + nIdxStart][nIdxD] = pdSensitivities[nIdx];
    }
  } // if sensitivity by pde
      
  if ( instdata.m_sensitivityMethod == SensitivityMethod_Adjoint )
  {
    SensitivityByAdjointData& aData = m_pAdjointDatas.back();

    aData.m_pRMatrix = pMatrix;
      
    aData.m_pdResiduals = pdResiduals;
  } // sensitivity by adjoint
  
}


} // namespace hg

} // namespace ito33
