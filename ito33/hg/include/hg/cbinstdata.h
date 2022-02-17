/////////////////////////////////////////////////////////////////////////////
// Name:        hg/cbinstdata.h
// Purpose:     cbinstdata class for HG model
// Created:     2005/04/11
// RCS-ID:      $Id: cbinstdata.h,v 1.9 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/cbinstdata.h
    @brief cbinstdata class for HG model
 */

#ifndef _HG_CBINSTDATA_H_
#define _HG_CBINSTDATA_H_

#include "ito33/vector.h"
#include "ito33/array.h"

#include "ito33/pricing/cblikeparams.h"
#include "ito33/pricing/event.h"
#include "ito33/pricing/cbconstraints.h"

#include "hg/instdatawithconstraints.h"

namespace ito33
{

namespace pricing
{
  class CBMeshManager;
}

namespace hg
{

  class Payoff;


class CBInstData : public InstDataWithConstraints
{

public:

  CBInstData( pricing::CBLikeParams& cbparams, 
              Model& model, 
              pricing::CBMeshManager& cbmeshes);

  /// Dummy virtual dtor for base class
  virtual ~CBInstData() { }
  
  /// @name Functions required by the SpecialEngine
  //@{

  virtual void Init();
  
  virtual void SetInitialValue();

  virtual void UpdateBeforeStep();

  /**
     Updates the data at the end of the grid(except for the last one).

     @attention{This function updates in particular the member m_nNbS because
                the size of the space mesh can change only if the grid change:
                i.e For a given grid, the size of the space mesh is unchanged}
    
   */
  virtual void UpdateAtEndOfGrid();

  // DoEvents(): Same implementation as in base, no need to reimplement

  //@}

  virtual void SetupFlags(const finance::ComputationalFlags& flags);

  /// Initial value defined by outside.
  void SetOutsideInitialValue(const shared_ptr<Payoff>& pPayoff)
  {
    m_pPayoff = pPayoff;
  }

  /// Gets the maximum mesh size in multi-grid case.
  size_t GetNbSpotsMax() const { return m_nNbSpotsMax; }

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

  /// True if new share feature, false otherwise.
  bool m_bHasNewShare;

  /// New share prices at current time step (to be calculated).
  Array<double> m_pdNewSharePrices;

  /// New share prices at last time step.
  Array<double> m_pdOldNewSharePrices;

  /// New share prices at the time step before last (for multistep methods).
  Array<double> m_pdOldOldNewSharePrices;

  /**
     Compute new share prices at maturity. 
     
     To do this, it prices a particular cb from the next fiscal year start 
     date after the maturity to this same maturity.

     @param pdS array of grid points
     @param nNbS number of grid points
     @param pdValues array of computed values
   */
  void 
  ComputeNewSharePricesAtMaturity
  (const double *pdS, size_t nNbS, double *pdValues);

  //@}

  /// Queries params to updates constraints.
  virtual void UpdateConstraints();


protected:

  /**
     Allocates the memory for recovery values etc.

     @param nNbS the size of the arrays to be allocated
   */
  void Alloc(size_t nNbS); 
  
  /** 
     Swaps the prices and Greek/sesnitivities (if any) arrays.
   */
  virtual void Swap();

  /**
     Applies an event to prices and greeks/sensitivities.
     
     @param pEvent A pointer to the to be applied event
   */
  virtual void ApplyEvent(const pricing::Event *pEvent);

  /**
     Solve a small cb from t+notice period to t and returns
     the values interpolated onto pdSpot

     @param pdSpots  array of grid points
     @param nNbSpots number of grid points
     @param pdValues array of values

     @return values of the smal cb interpolated
   */
  void SolveCallNotice(const double *pdSpots, const size_t nNbSpots,
                       const double* pdNewSharePrices, double *pdValues);

  /**
     Updates the data depending on space mesh at a change of grid time.

     @attention{ The function Swap() must be called just before }
    
     @param pdOldSpots Pointer to the space mesh of the old grid
     @param nNbOldSpots size of the space mesh of the old grid
     @param pdSpots Pointer to the space mesh of the new grid
     @param nNbSpots size of the space mesh of the new grid
   */
  virtual void InterpWithPassageOfSpaceMesh
               (const double *pdOldSpots, size_t nNbOldSpots, 
                const double *pdSpots, size_t nNbSpots);

  virtual void ApplyConstraintsAfterEventsOrConstraintsUpdate();

    /// Apply the constraints to all arrays depends on constraints.
  virtual void ApplyConstraintsToAll();

  /// Maximum mesh size in multi-grid case
  size_t m_nNbSpotsMax;

  /// The cb constraints (call, put, conversion)
  pricing::CBConstraints m_constraints;

  /// The params for cb
  pricing::CBLikeParams& m_cbparams;

  /// The cb mesh manager
  pricing::CBMeshManager& m_cbmeshes;

  /// The outside payoff if any
  shared_ptr<Payoff> m_pPayoff;

  /// The type of the non linear solver (if applicable), 0 penalty, 1 forzen
  int m_iSolverType;

private:

  NO_COPY_CLASS(CBInstData);

  friend class CBOptionInstData;

}; // class CBInstData


} // namespace hg

} // namespace ito33 

#endif // #ifndef _HG_CBINSTDATA_H_

