/////////////////////////////////////////////////////////////////////////////
// Name:        hg/hedgingdata.h
// Purpose:     contracts class for HERO, and hedge ratio calculator
// Created:     2005/09/26
// RCS-ID:      $Id: hedgingdata.h,v 1.5 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/hedgingdata.h
    @brief The declaration of the contracts class for HERO. 

    The base class for HERO contracts, and hedge ration calculator.  It is 
    not a contract class in the usual sense, but was done this way to 
    satisfy requirements of the engine.  It is primarily a container for all 
    the HERO data objects (target, hedge contracts) and the calculations 
    needed by the HERO PDE and hedge ratios.
 */

#ifndef _HG_HEDGINGDATA_H_
#define _HG_HEDGINGDATA_H_

#include "ito33/array.h"
#include "ito33/vector.h"

#include "ito33/numeric/densematrix.h"

#include "ito33/pricing/contract.h"

namespace ito33
{

namespace finance { class ModelOutput; }

namespace hg
{

  class UnderlyingProcess;

/// The declaration of the HERO contract class.
class HedgingData : public pricing::Contract 
{
public:

  /**
     The ctor.

     Needs to be based on a derivative so it can be derived from Contract.
     It must be derived from Contract so it can be used by the standard
     pricing code.
    
     @param derivative a reference to the target derivative object     
     @param pRealProcess the real (not risk neutral) underlying process
     @param nNbHedges the number of hedge contracts
   */
  HedgingData(const finance::Derivative& derivative, 
              shared_ptr<hg::UnderlyingProcess> pRealProcess,
              size_t nNbHedges);

  /**
    Destructor.
   */
  virtual ~HedgingData() { }
  

  /**
     Set the pricing model outputs of the hedging contracts.

     Does not include the pricing output of the target.

     @param ppHedgeOutputs the modeloutputs of the hedge contracts
   */
  void SetHedgeModelOutputs(
    const std::vector< shared_ptr<finance::ModelOutput> >& ppHedgeOutputs) 
  { 
    m_ppHedgeOutputs = ppHedgeOutputs; 
  }

  /**
     Set the model output of the target contract.

     @param pTargetOutput the modeloutput of the target contract
   */
  void SetTargetModelOutput(shared_ptr<finance::ModelOutput> pTargetOutput)
  { 
    m_pTargetOutput = pTargetOutput; 
  }

  /** 
    Get the number of hedging contracts

    Does not include the target

    @return the number of hedging contracts
   */
  size_t GetNbHedges() const
  {
    return m_nNbHedges;
  }


  /**
    Compute the hedge ratios at the valuation time.

    Assumes that the model outputs were constructed correctly, with
    the analysis date set to the valuation time. It is expected that
    this function will be called by hedging code, which only needs data
    at the valuation date.

    @param dTargetRatio (output) target ratio
    @param pdHedgeRatios (output) hedge contract ratios
   */
  void ComputeFinalHedgeRatios(double& dTargetRatio, 
                               std::vector<double>& pdHedgeRatios) const;


  /**
    Compute the PDE terms at the specified time.

    Computes constant PDE terms for HERO.

    @param dTime the time at which to compute the PDE terms
    @param pdSpots the current PDE grid
    @param nNbS the number of grid points
   */
  std::vector<double> 
  ComputePDETermsAt(double dTime, const double* pdSpots, size_t nNbS) const;

  /**
     Force system solving by NAG.

     Hero calculations use SVD by default (see the cpp file).  Set bForceNAG 
     to force NAG to be used.  Should only be set by test code.
   */
  static bool bForceNAG;

protected:
  
  /// Compute helper hedge matrix and vector
  void ComputeHedgeData(double** ppdVMatrix, double* pdA) const;

  /// Compute the PDE terms 
  std::vector<double> ComputePDETerms(const double* pdSpots,size_t nNbS) const;

  /// Solve system
  void Solve(const double* const* ppdA, const double* pdB, double* pdX, 
             bool bUseNag) const;

  /// The target derivative
  const finance::Derivative& m_target;

  /// The pricing model output of the target contract
  shared_ptr<finance::ModelOutput> m_pTargetOutput;

  /// The pricing model output of the hedging contracts
  std::vector< shared_ptr<finance::ModelOutput> > m_ppHedgeOutputs;

  /// The real underlying process
  shared_ptr<hg::UnderlyingProcess> m_pRealUP;

  /// The number of contracts used to hedge (in addition to underlying)
  size_t m_nNbHedges;

  /// The number of regimes
  size_t m_nNbRegimes;

  /// Inverses of the squares of the total volatilities
  std::vector<double> m_pdInvSqrTotalVols;  

  /// matrix used to store data
  mutable numeric::DenseMatrix m_matrix;

  /// matrix used to solve linear systems
  mutable numeric::DenseMatrix m_matrixToSolve;

  /// The SVD method needs a V matrix and W vector
  mutable std::vector<double> pdW;
  mutable numeric::DenseMatrix m_matrixV;

  /// jump amplitudes and intensities: [regime k][jump from regime k]  
  std::vector< std::vector<double> > m_ppdJumpIntensities;
  std::vector< std::vector<double> > m_ppdJumpAmplitudes;
  
  /// Helper vectors for constructing the final hedges
  mutable std::vector< std::vector<double> > m_ppdOrigSpots;
  mutable std::vector< std::vector<double> > m_ppdOrigPrices;
  mutable std::vector<double> m_pdPrices;
  mutable std::vector<double> m_pdDeltas;
  mutable std::vector<double> m_pdValuesAfterDefault;

  /// Helper data for computing hero
  // the delta values: [hedges][nNbS * nNbRegimes]
  mutable std::vector< std::vector<double> > m_ppdDeltas;

  // price differences: [hedge][regime][jump][S]
  mutable std::vector< std::vector< std::vector< std::vector< double> > > > 
    m_ppppdJumpDiffs;


private:

  NO_COPY_CLASS(HedgingData);

}; // class HedgingData;


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_HEDGINGDATA_H_
