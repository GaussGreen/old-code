/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/theoreticalmodel_perfecthedge.cpp
// Purpose:     Implementation of the perfect hedge for the ihg model
// Author:      Zhang
// Created:     2005/02/02
// RCS-ID:      $Id: theoreticalmodel_perfecthedge.cpp,v 1.25 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/array.h"
#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/bondlike/bondliketerms.h"
#include "ito33/finance/bondlike/convertiblelike.h"
#include "ito33/finance/bondlike/cboption.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/ihg/error.h"
#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/perfect_hedge_ratios.h"

#include "ito33/numeric/densematrix.h"
#include "ito33/numeric/densesolver_direct.h"

extern const ito33::ihg::Error 
  ITO33_IHG_PERFECT_HEDGE_FAILED,
  ITO33_IHG_PERFECT_HEDGE_EXCHANGEABLE,
  ITO33_IHG_PERFECT_HEDGE_CBOPTION_EXCHANGEABLE,
  ITO33_IHG_PERFECT_HEDGE_BAD_CROSSCURRENCY_HEDGER,
  ITO33_IHG_PERFECT_HEDGE_DIFFERENT_UNDERLYING;

namespace ito33
{

namespace ihg
{


PerfectHedgeRatios TheoreticalModel::ComputePerfectHedgeRatios
                          (
                            const finance::Derivative& target,
                            const finance::Derivative& hedge
                          ) const
{
  // Verify that the derivative used to hedge and the hedged one have the 
  // same session.

  CHECK_COND(target.GetSessionData() == hedge.GetSessionData(),
               ITO33_IHG_PERFECT_HEDGE_DIFFERENT_UNDERLYING);
  
  // verify neither of the two derivatives is exchangeable
  const ito33::finance::ConvertibleLike*
    pTarget = dynamic_cast<const finance::ConvertibleLike*>(&target);

  if ( pTarget )
  {    
    CHECK_COND(!pTarget->IsExchangeable(),
               ITO33_IHG_PERFECT_HEDGE_EXCHANGEABLE);
  }

  const ito33::finance::ConvertibleLike*
    pHedge = dynamic_cast<const finance::ConvertibleLike*>(&hedge);

  if ( pHedge )
  {
    CHECK_COND(!pHedge->IsExchangeable(),
               ITO33_IHG_PERFECT_HEDGE_EXCHANGEABLE);
  }
  
  // verify neither of the two derivatives is a CBOption with an exchangeable CB
  const ito33::finance::CBOption*
    pTargetCBOption = dynamic_cast<const finance::CBOption*>(&target);

  if ( pTargetCBOption )
  {    
    CHECK_COND(!pTargetCBOption->GetConvertibleBond()->IsExchangeable(),
               ITO33_IHG_PERFECT_HEDGE_CBOPTION_EXCHANGEABLE);
  }

  const ito33::finance::CBOption*
    pHedgeCBOption = dynamic_cast<const finance::CBOption*>(&hedge);

  if ( pHedgeCBOption )
  {
    CHECK_COND(!pHedgeCBOption->GetConvertibleBond()->IsExchangeable(),
               ITO33_IHG_PERFECT_HEDGE_CBOPTION_EXCHANGEABLE);
  }

  bool
    bTargetIsCrossCurrency = target.IsCrossCurrency(),
    bHedgeIsCrossCurrency = hedge.IsCrossCurrency();

  // Instrument (H) used to hedge could be a cross-currency one but 
  // IF AND ONLY IF target instrument (V) is also a cross-currency one and if 
  // the currency of V is the same as the currency of H. 
  if ( bHedgeIsCrossCurrency )
    CHECK_COND(bTargetIsCrossCurrency &&
               target.GetNumeraire() == hedge.GetNumeraire(),
               ITO33_IHG_PERFECT_HEDGE_BAD_CROSSCURRENCY_HEDGER);

  shared_ptr<TheoreticalModel> ptm( Clone() );

  ptm->SetQualityControl(m_pQualityControl);

  // don't use flags in derivative, use default flags
  // this may need to be changed if we let user get the model output
  // in the future, as seen in HG
  ptm->SetExternalFlagsToDefaults();

  shared_ptr<finance::ModelOutput> pTargetOutput = ptm->Compute(target);

  shared_ptr<finance::ModelOutput> pHedgeOutput = ptm->Compute(hedge);  

  PerfectHedgeRatios result;

  if ( !bTargetIsCrossCurrency )
  {
    // b1 = dV/dS,  b2 = valueAfterDefault - price    
    double
      dB1 = pTargetOutput->GetDelta(),
      dB2 = pTargetOutput->GetValueAfterDefault() - pTargetOutput->GetPrice();

    double
      dA11 = pHedgeOutput->GetDelta(),
      dA21 = pHedgeOutput->GetValueAfterDefault() - pHedgeOutput->GetPrice(),
      //
      dA12 = 1.,
      dA22 = - target.GetSessionData()->GetSpotSharePrice();

    double dX1, dX2;

    CHECK_COND
    ( 
      numeric::Solve2DLinearSystem(dA11, dA12, dB1, dA21, dA22, dB2, dX1, dX2), 
      ITO33_IHG_PERFECT_HEDGE_FAILED
    );

    result.SetDefaultHedgeRatio( dX1 );
    result.SetUnderlyingHedgeRatio( dX2 );
  }
  else // bTargetIsCrossCurrency == true
  {
    double
      dFXRate = target.GetSessionData()->GetSpotFXRate
                ( 
                  target.GetSessionData()->GetEquity()->GetNumeraire(),
                  target.GetNumeraire()
                );

    size_t nDim = 3;

    // B matrix
    numeric::DenseMatrix matrixB(nDim, 1);
    
    // b0 = dV/dS,  b1 = dV/dX, b2 = valueAfterDefault - price
    matrixB[0][0] = pTargetOutput->GetDelta();
    matrixB[1][0] = pTargetOutput->GetFXDelta();
    matrixB[2][0] = pTargetOutput->GetValueAfterDefault() 
                  - pTargetOutput->GetPrice();

    // A matrix
    numeric::DenseMatrix matrixA(nDim, nDim);

    if ( !bHedgeIsCrossCurrency )
    {
      matrixA[0][0] = pHedgeOutput->GetDelta() * dFXRate;
      matrixA[1][0] = pHedgeOutput->GetPrice();
      matrixA[2][0] = ( pHedgeOutput->GetValueAfterDefault() 
                    - pHedgeOutput->GetPrice() ) * dFXRate;
    }
    else  // bHedgeIsCrossCurrency == true
    {
      matrixA[0][0] = pHedgeOutput->GetDelta();
      matrixA[1][0] = pHedgeOutput->GetFXDelta();
      matrixA[2][0] = pHedgeOutput->GetValueAfterDefault() 
                    - pHedgeOutput->GetPrice();
    }
    //
    matrixA[0][1] = dFXRate;
    matrixA[1][1] = target.GetSessionData()->GetSpotSharePrice();
    matrixA[2][1] = - target.GetSessionData()->GetSpotSharePrice() * dFXRate;
    //
    matrixA[0][2] = 0.;
    matrixA[1][2] = 1.;
    matrixA[2][2] = 0.;

    CHECK_COND
    ( 
      numeric::SolveLinearSystem(matrixA.Get(), matrixB.Get(), nDim, 1), 
      ITO33_IHG_PERFECT_HEDGE_FAILED
    );

    result.SetDefaultHedgeRatio( matrixB[0][0] );
    result.SetUnderlyingHedgeRatio( matrixB[1][0] );
    result.SetFXHedgeRatio( matrixB[2][0] );   
  }
  
  return result;
}


} // namespace ihg

} // namespace ito33

