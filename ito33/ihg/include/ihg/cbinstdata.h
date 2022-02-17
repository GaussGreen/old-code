/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/cbinstdata.h
// Purpose:     cbinstdata class for ihg project
// Author:      Nabil
// Created:     2004/03/30
// RCS-ID:      $Id: cbinstdata.h,v 1.53 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/cbinstdata.h
    @brief cbinstdata class for ihg project
 */

#ifndef _IHG_CBINSTDATA_H_
#define _IHG_CBINSTDATA_H_

#include "ito33/vector.h"
#include "ito33/array.h"

#include "ito33/pricing/cblikeparams.h"
#include "ito33/pricing/event.h"
#include "ito33/pricing/cbconstraints.h"

#include "ihg/instdatamulti.h"

namespace ito33
{

namespace finance
{
  class Payoff;
}

namespace pricing
{
  class CBMeshManager;
}

namespace ihg
{


class CBInstData : public InstDataMulti
{

public:

  CBInstData( pricing::CBLikeParams& cbparams, 
              Model& model, 
              pricing::CBMeshManager& cbmeshes);

  virtual ~CBInstData() { }
  
  virtual void Init();
  
  virtual void SetInitialValue();

  virtual void UpdateBeforeStep();

  virtual void UpdateConstraints();
  
  void SetSolverType(int iSolverType) { m_iSolverType = iSolverType; }

  /**
      This function updates the data at the end of the grid 
      (except for the last one).

      @attention{This function updates in particular the member m_nNbS because
                 the size of the space mesh can change only if the grid change:
                 i.e For a given grid, the size of the space mesh is unchanged} 
   */
  virtual void UpdateAtEndOfGrid();

  /**
      Initial value defined by outside.
   */
  void SetOutsideInitialValue(const shared_ptr<finance::Payoff>& pPayoff)
  {
    m_pPayoff = pPayoff;
  }
  
  /**
      Gets the recovery value at initial spot.

      This function is reimplemented as the recovery value function
      is not constant for exchangeable.

      @return the recovery value at initial spot.
   */
  virtual double GetValueAfterDefaultAtInitialSpot();

  /** 
      Gets the correction for the convection term due to the change of variable.

      @return The value of this correction.
   */
  double GetSpeedCorrection() const;

  /**
      Checks if the contract is a cross-currency one.

      @return true if it is a cross-currency contract, false otherwise
   */
  bool IsCrossCurrency() const 
  { 
    return m_cbparams.GetCBLike().IsCrossCurrency(); 
  }

  /**
      Checks if the contract is a fixed quanto one.

      @return true if it is a fixed quanto contract, false otherwise
   */
  bool IsFixedQuanto() const
  { 
    return m_cbparams.GetCBLike().IsFixedQuanto(); 
  }
  
  /**
      Gets the volatility of the FX Rate (for a fixed quanto).

      @return The volatility of the FX rate.
   */
  double GetFXRateVolatility() const
  { 
    return m_cbparams.GetCBLike().GetFXRateVolatility(); 
  }
  
  /**
      Gets the correlation between the underlying share and the FX rate 
      (for a fixed quanto).

      @return The correlation between the underlying share and the FX rate. 
   */
  double GetCorrelationBetweenUnderlyingAndFXRate() const
  { 
    return m_cbparams.GetCBLike().GetCorrelationBetweenUnderlyingAndFXRate(); 
  }

  /**
      Checks if we have a new share feature

      @return true if we have a new share feature, false otherwise.
   */
  bool HasNewShare() const { return m_bHasNewShare; }
  
  /**
      Computes new share prices at maturity. 
     
      To do this, it prices a particular cb from the next fiscal year start 
      date after the maturity to this same maturity.

      @param pdSpots  array of grid points
      @param nNbSpots number of grid points
      @param pdValues array of values
   */
  void ComputeNewSharePricesAtMaturity(const double *pdSpots, 
                                       const size_t nNbSpots,
                                       double *pdValues);  

  /// Recovery rate
  double m_dRecoveryRate;

  /// vector of recovery values
  std::vector<double> m_pdRecoveryValues;
    
  /// @name concerning exchangeable bond
  //@{

  /// whether we have second hazard rate (of the derivative)
  bool m_bDerivativeHasOwnHR;

  /// hazard rate value of the derivative at current time step
  double m_dHROfDerivative;
  //@}

  /// @name concerning new share feature  
  //@{

  /// New share prices at current time step (to be calculated).
  Array<double> m_pdNewSharePrices;

  /// New share prices at last time step.
  Array<double> m_pdOldNewSharePrices;

  /// New share prices at the time step before last (for multistep methods).
  Array<double> m_pdOldOldNewSharePrices;

  //@}


protected:

  /**
      Allocates the memory for the price and vega arrays.

      @param nNbS the size of the arrays to be allocated
   */
  virtual void Alloc(size_t nNbS);

  /**
      Applies an event to prices and greeks.
     
      @param pEvent A pointer to the to be applied event
   */
  virtual void ApplyEvent(const pricing::Event *pEvent);

  /** 
      Swaps the price and Greek (if any) arrays.
   */
  virtual void Swap();

  /**
      This function updates the data depending on space mesh at a change of
      grid time.

      @attention{ The function Swap() must be called just before }
    
      @param pdOldSpots Pointer to the space mesh of the old grid
      @param nNbOldSpots size of the space mesh of the old grid
      @param pdSpots Pointer to the space mesh of the new grid
      @param nNbSpots size of the space mesh of the new grid
   */
  virtual void InterpWithPassageOfSpaceMesh(const double *pdOldSpots, 
                  size_t nNbOldSpots, const double *pdSpots, size_t nNbSpots);

  /**
      Solves a small cb from t+notice period to t and returns
      the values interpolated onto pdSpot

      @param pdSpots  array of grid points
      @param nNbSpots number of grid points
      @param pdNewSharePrices new share prices
      @param pdValues array of values

      @return values of the smal cb interpolated
   */
  void SolveCallNotice(const double *pdSpots, const size_t nNbSpots,
                       const double* pdNewSharePrices, double *pdValues);

  virtual void ApplyConstraintsAfterEventsOrConstraintsUpdate();
 
  /// True if new share feature, false otherwise.
  bool m_bHasNewShare;

  /// The cb constraints (call, put, conversion)
  pricing::CBConstraints m_constraints;

  /// The params for cb
  pricing::CBLikeParams& m_cbparams;

  /// The cb mesh manager
  pricing::CBMeshManager& m_cbmeshes;

  /// The outside payoff if any
  shared_ptr<finance::Payoff> m_pPayoff;

  /// Solver type for coupled equation
  int m_iSolverType;

private:

  NO_COPY_CLASS(CBInstData);

  friend class CBOptionInstData;

}; // class CBInstData


} // namespace ihg

} // namespace ito33 

#endif // #ifndef _IHG_CBINSTDATA_H_
