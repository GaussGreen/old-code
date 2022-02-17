/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/cboptioninstdata.h
// Purpose:     cboptioninstdata class for ihg project
// Author:      Nabil
// Created:     2004/10/14
// RCS-ID:      $Id: cboptioninstdata.h,v 1.2 2005/12/27 18:07:53 wang Exp $
// Copyright:   (c) 2004-2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/cboptioninstdata.h
    @brief cboptioninstdata class for ihg project
 */

#ifndef _IHG_CBOPTIONINSTDATA_H_
#define _IHG_CBOPTIONINSTDATA_H_

#include "ito33/vector.h"
#include "ito33/array.h"

#include "ito33/pricing/cboptionparams.h"
#include "ito33/pricing/event.h"
#include "ito33/pricing/cbconstraints.h"

#include "ihg/instdatamulti.h"
#include "ihg/cbinstdata.h"

namespace ito33
{

namespace pricing
{
  class CBMeshManager;
}

namespace ihg
{

class CBInstData;

class CBOptionInstData : public InstDataMulti
{

public:

  CBOptionInstData( pricing::CBOptionParams& cbparams, 
                    Model& model, 
                    pricing::CBMeshManager& cbmeshes);

  virtual ~CBOptionInstData() { }
  
  /// @name implementation of virtual functions
  //@{

  virtual void Init();
  
  virtual void SetInitialValue();

  virtual void UpdateBeforeStep();
  
  /**
    This function updates the data at the end of the grid 
    (except for the last one)

    @attention{This function updates in particular the member m_nNbS because
               the size of the space mesh can change only if the grid change:
               i.e For a given grid, the size of the space mesh is unchanged}
    
    */
  virtual void UpdateAtEndOfGrid();

  //@}

  /**
     Help function to turn off all flags
   */
  virtual void TurnAllFlagsOff() 
  { 
    BackwardInstData::TurnAllFlagsOff();
    
    m_cbinstdata.TurnAllFlagsOff();
  }

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
  
  /// Gets the reference on the CBInstData object
  CBInstData* GetCBInstData()
  {
    return &m_cbinstdata;
  }

  /// True if before maturity of cb option
  bool IsCBOptionWindow() const { return m_cboptionparams.InCBOptionWindow(); }
  
  void ApplyConstraintsToCBOptionGreek(double* pdGreeks,
                                       const double* pdCBGreeks, 
                                       const int* piFlagConstraints, 
                                       size_t nNbValues);
  
  /// @name concerning exchangeable bond
  //@{

  /// whether we have second hazard rate (of the derivative)
  bool m_bDerivativeHasOwnHR;

  /// hazard rate value of the derivative at current time step
  double m_dHROfDerivative;

  //@}

protected:

  /**
    Allocates the memory for the price, vega... arrays.

    @param nNbS the size of the arrays to be allocated
  */
  virtual void Alloc(size_t nNbS);

  /**
     Apply an event to prices and greeks
     
     @param pEvent A pointer to the to be applied event
   */
  virtual void ApplyEvent(const pricing::Event *pEvent);

  /** 
     Swap the price and Greek (if any) arrays
   */
  virtual void Swap();

  /**
    This function updates the data depending on space mesh at a change of grid
    time.

    @attention{ The function Swap() must be called just before }
    
    @param pdOldSpots Pointer to the space mesh of the old grid
    @param nNbOldSpots size of the space mesh of the old grid
    @param pdSpots Pointer to the space mesh of the new grid
    @param nNbSpots size of the space mesh of the new grid

   */
  virtual void InterpWithPassageOfSpaceMesh(const double *pdOldSpots, 
                  size_t nNbOldSpots, const double *pdSpots, size_t nNbSpots);

  virtual void ApplyConstraintsAfterEventsOrConstraintsUpdate();
  
  virtual void UpdateConstraints() {};
 
  /// The params for cb option
  pricing::CBOptionParams& m_cboptionparams;

  /// The cb mesh manager
  pricing::CBMeshManager& m_cbmeshes;

  /// The cb instdata of the cb
  CBInstData m_cbinstdata;

  /// constraint of cb option
  pricing::MinConstraint m_CBOptionConstraint;

  Array<double> m_pdTmp1;

private:

  NO_COPY_CLASS(CBOptionInstData);

}; // class CBOptionInstData

} // namespace ihg

} // namespace ito33 

#endif // #ifndef _IHG_CBOPTIONINSTDATA_H_
