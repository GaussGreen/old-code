/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/instdatatimeonly.h
// Purpose:     base instdata class
// Author:      Zhang Yunzhi
// Created:     2003/12/26
// RCS-ID:      $Id: instdatatimeonly.h,v 1.12 2005/03/31 15:47:20 wang Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/instdatatimeonly.h
    @brief base instdata class for time only PDE problem
 */

#ifndef _ITO33_PRICING_INSTDATATIMEONLY_H_
#define _ITO33_PRICING_INSTDATATIMEONLY_H_

#include "ito33/common.h"
#include "ito33/numeric/schemetype.h"

namespace ito33
{

namespace pricing
{

class Event;
class Params;
class MeshManager;

/// base instdata class for all kinds of models and instruments
class InstDataTimeOnly
{
public:

  /**
     ctor
     
     @param params the Params
     @param meshes the mesh manager
   */  
  InstDataTimeOnly(Params &params, MeshManager &meshes)
                 : m_params(params), m_meshes(meshes),
                   m_dTimeStep(0.),
                   m_dRate(0.), m_dDerivativeRate(0.), m_dForeignRate(0.),
                   m_dOldTimeStep(0.), 
                   m_dInverseTimeStep(0.), 
                   m_dOldInverseTimeStep(0.),
                   m_dTimeWeight(0.), 
                   m_dOldTimeWeight(0.), 
                   m_dOldOldTimeWeight(0.)
  { }

  /// dtor
  virtual ~InstDataTimeOnly() { }


  /// @name Functions required by the Engine
 
  //@{

  /** 
    Allocates memory when the maximum size of the space mesh is known.
    It does necessary initialization work as well.

    m_meshes.SetupMe() must be called before calling this function.
     
    */
  virtual void Init() = 0;

  /**
     Update myself at the beginning of each time step
    
     Required by Engine class. At each time step, 
     this function update the instanous parameters from the
     meshes manager. It works also on the Data part.
   */
  virtual void UpdateBeforeStep();

  /**
     Get events from mesh manager and do events

     At the end of each time step, 
     this function check out from the mesh manager all events
     happen at current period, and apply the events' method to the Data part.
   */
  virtual void DoEvents();

  /**
    Set the initial values of myself from the mesh manager.
    
    At the start period, the mesh manager tell me what are the start conditions
    for the Data part. It also works on the InstParam part.
   */
  virtual void SetInitialValue()
  {
    // we set HasEvent as true, as the initial time should be considered
    // as an event.
    m_bHasEvent = true;
  }

  //@}
 
  /// The current timestep. 
  double m_dTimeStep;
  
  /// Current interest rate (yield)
  double m_dRate;
  
  /// Current interest rate (derivative)
  double m_dDerivativeRate;

  /// Current foreign rate
  double m_dForeignRate;

  /// The timestepping scheme at current timestep
  numeric::SchemeType m_schemeType;


  /// @name Helpers for time stepping
  //@{
  
  /// The previous timestep (needed for multistep methods). 
  double m_dOldTimeStep;

  /// The inverse of m_dTimeStep
  double m_dInverseTimeStep;

  /// The inverse of m_dOldTimeStep
  double m_dOldInverseTimeStep;
  
  //@}


  /// @name Helpers for weighting
  //@{

  /// coefficient before unknown for the descretization of dC/dt
  double m_dTimeWeight;

  /**
    coefficient before the "C" obtained by the last time step for
    the descretization of dC/dt.
    */
  double m_dOldTimeWeight;

  /**
    coefficient before the "C" obtained by the time step before last for
    the descretization of dC/dt. It is only available for Three level
    scheme
   */
  double m_dOldOldTimeWeight;
  
  //@}

 /// A boolean indicates if there are events applied at current time
  bool m_bHasEvent;


protected:
  
  /// Helper function for scheme weights computation
  void ComputeSchemeWeights();

  /**
     Apply an event to pricers and/or Greaks
     
     Needs to be overloaded and called by derived class to add model specific
     items
   */
  virtual void ApplyEvent(const Event *pEvent) = 0; 

  /**
     helper function for applying events

     @return true if some events got applied, false otherwise
   */
  void ApplyEvents();

  /// The params
  Params &m_params;

  /// The mesh manager
  MeshManager &m_meshes;

private:

  NO_COPY_CLASS(InstDataTimeOnly);

}; // class InstDataTimeOnly


} // namespace pricing

} // namespace ito33 

#endif // #ifndef _ITO33_PRICING_INSTDATATIMEONLY_H_

