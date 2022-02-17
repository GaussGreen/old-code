/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/backwardinstdata.h
// Purpose:     backward instdata class for ihg projet
// Author:      Wang
// Created:     2004/02/13
// RCS-ID:      $Id: backwardinstdata.h,v 1.15 2005/12/30 11:51:14 nabil Exp $
// Copyright:   (c) 2003-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/backwardinstdata.h
    @brief ihg backward instdata class
 */

#ifndef _IHG_BACKWARDINSTDATA_H_
#define _IHG_BACKWARDINSTDATA_H_

#include "ito33/array.h"

#include "ihg/instdata.h"

namespace ito33
{

namespace ihg
{


/// Base instdata class for ihg projets
class BackwardInstData : public InstData
{

public:

  BackwardInstData(pricing::Params& params, 
                   Model& model, 
                   pricing::MeshManager& meshes)
                 : InstData(params, model, meshes)
  {
    // Set all ihg model dependant compute flags as false
    m_bComputeVega = false;
    m_bComputeFugit = false;
    // the flags should be set manually and explicitly later.
  }

  virtual ~BackwardInstData() { }

  /// @name Functions required by the Engine
 
  // Init() : No implementation at this level, so no need to redeclare

  virtual void UpdateBeforeStep();

  // DoEvents(): Same implementation as in base, no need to reimplement

  // SetInitialValue(): No implementation at this level, so no need to redeclare

  //@}
  
  /**
     Help function to turn off all flags
   */
  virtual void TurnAllFlagsOff() 
  { 
    m_bComputeVega = false;
    m_bComputeFugit = false; 
  }
  
  ///
  virtual const int* GetConstraintFlags() const { return 0; }

  /**
    Gets the recovery value at initial spot.

    Note that for most instrument recovery value is constant with
    respect to spot. So the default implementation just return
    m_dRecoveryValue.

    @return the recovery value at initial spot.
    */
  virtual double GetValueAfterDefaultAtInitialSpot()
  {
    return m_dRecoveryValue;
  }

  /**
    Get initial spot

    @return initial spot in params
    */
  double GetInitialSpot() const;

  /**
    Get space mesh information for numOutput

    @param nNbS (output) number of space points
    @return space points
    */
  const double* GetSpaceMesh(size_t& nNbS)
  {
    nNbS = m_nNbS;
    return m_pdS;
  }

public:  

  /// Default value at current time step
  double m_dRecoveryValue;

  /// @name vega related data
  //@{

  /// Flag indicates if we compute vega
  bool m_bComputeVega;

  /// Current volatilities for vega computation 
  Array<double> m_pdVols;

  /// Vega at current time step (to be calculated). 
  Array<double> m_pdVegas;

  /// Vega at last time step.  
  Array<double> m_pdOldVegas;

  /// Vega at the time step before last (for multistep methods).
  Array<double> m_pdOldOldVegas; 
  
  //@}
 
  /// @name fugit related data
  //@{
   
  /// Flag indicates if we compute fugit
  bool m_bComputeFugit;

  /// Fugit at current time step (to be calculated).
  Array<double> m_pdFugits;

  /// Fugit at last time step.
  Array<double> m_pdOldFugits;

  /// Fugit at the time step before last (for multistep methods).
  Array<double> m_pdOldOldFugits;
  
  //@}


protected:



  /**
     Allocate the memory for the price, vega  and fugit arrays

     @param nNbS the size of the arrays to be allocated
   */
  void Alloc(size_t nNbS);

  /**
     Apply an event to prices and Greeks
     
     @param pEvent A pointer to a to be applied event
   */
  virtual void ApplyEvent(const pricing::Event *pEvent);


  /** 
     Swap the price and Greek (if any) arrays
   */
   virtual void Swap();


private:

  NO_COPY_CLASS(BackwardInstData);

}; // class BackwardInstData


} // namespace ihg

} // namespace ito33 

#endif // #ifndef _IHG_BACKWARDINSTDATA_H_

