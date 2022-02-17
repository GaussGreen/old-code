/////////////////////////////////////////////////////////////////////////////
// Name:        hg/cboptioninstdata.h
// Purpose:     cboptioninstdata class for HG
// Created:     2006/01/19
// RCS-ID:      $Id: cboptioninstdata.h,v 1.2 2006/03/20 14:54:09 yann Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/cboptioninstdata.h
    @brief cboptioninstdata class for HG
 */

#ifndef _HG_CBOPTIONINSTDATA_H_
#define _HG_CBOPTIONINSTDATA_H_

#include "ito33/vector.h"
#include "ito33/array.h"

#include "ito33/pricing/cboptionparams.h"
#include "ito33/pricing/event.h"
#include "ito33/pricing/cbconstraints.h"

#include "hg/cbinstdata.h"

namespace ito33
{

namespace pricing
{
  class CBMeshManager;
}

namespace hg
{

class CBOptionInstData : public InstDataWithConstraints
{

public:

  CBOptionInstData( pricing::CBOptionParams& cbparams, 
                    Model& model, 
                    pricing::CBMeshManager& cbmeshes);
  
  /// @name implementation of virtual functions
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

  /// Gets the maximum mesh size in multi-grid case.
  size_t GetNbSpotsMax() const { return m_nNbSpotsMax; }

  virtual void SetupFlags(const finance::ComputationalFlags& flags);

  /** 
      Gets the correction for the convection term due to the change of variable.

      @return The value of this correction.
   */
  double GetSpeedCorrection() const;
      
  /// updates CB Option constraint
  void UpdateCBOptionConstraint();
  
  /// Updates the constraint and the results for the cb option
  void UpdateCBOptionResults();

  /// Gets the constraint of CB Option
  const pricing::Constraints* GetCBOptionConstraint() const
  {
    return &m_CBOptionConstraint;
  }  
  
  /// Gets the reference to the CBInstData object
  CBInstData& GetCBInstData()
  {
    return m_cbinstdata;
  }

  /// True if before maturity of cb option
  bool IsCBOptionWindow() const { return m_cboptionparams.InCBOptionWindow(); }
  
protected:

  /**
      Allocates the memory for the price, sensitivity arrays.

      @param nNbS the size of the arrays to be allocated
   */
  virtual void Alloc(size_t nNbS);

  /**
      Applies an event to prices and greeks.
 
      @param pEvent A pointer to the to be applied event
   */
  virtual void ApplyEvent(const pricing::Event *pEvent);

  /// Swaps the price and sensitivites (if any) arrays
  virtual void Swap();

  /**
      This function updates the data depending on space mesh at a change of 
      grid time.

      @attention{ The function Swap() must be called just before }
      
      @param pdOldS Pointer to the space mesh of the old grid
      @param nNbOldS size of the space mesh of the old grid
      @param pdS Pointer to the space mesh of the new grid
      @param nNbS size of the space mesh of the new grid
   */
  virtual void InterpWithPassageOfSpaceMesh
      (const double *pdOldS, size_t nNbOldS, const double *pdS, size_t nNbS);

  virtual void ApplyConstraintsAfterEventsOrConstraintsUpdate();
  
  virtual void UpdateConstraints() { }
 
  /// The params for cb option
  pricing::CBOptionParams& m_cboptionparams;

  /// The cb mesh manager
  pricing::CBMeshManager& m_cbmeshes;

  /// The cb instdata of the cb
  CBInstData m_cbinstdata;

  /// constraint of cb option
  pricing::MinConstraint m_CBOptionConstraint;

  Array<double> m_pdConstriantsTmp;

private:

  /// Maximum mesh size in multi-grid case
  size_t m_nNbSpotsMax;

  NO_COPY_CLASS(CBOptionInstData);

}; // class CBOptionInstData

} // namespace hg

} // namespace ito33 

#endif // #ifndef _HG_CBOPTIONINSTDATA_H_
