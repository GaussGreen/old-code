/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/forwardmeshmanager.h
// Purpose:     base mesh manager for forward PDE problems
// Author:      Wang
// Created:     2004/02/11
// RCS-ID:      $Id: forwardmeshmanager.h,v 1.12 2006/01/03 17:18:48 zhang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/pricing/forwardmeshmanager.h
   @brief base mesh manager for forward PDE problems
 */

#ifndef _ITO33_PRICING_FORWARDMESHMANAGER_H_
#define _ITO33_PRICING_FORWARDMESHMANAGER_H_

#include "ito33/pricing/meshmanager.h"

namespace ito33
{

namespace pricing
{


/// Forward time mesh manager
class ForwardMeshManager : public MeshManager
{

public:
  
  /**
     Ctor, call the ctor of the base class MeshManager

     @param params The underlying contract paramaters
     @param model The model used for the pricing
   */  
  ForwardMeshManager(Params &params, Model &model) 
                   : MeshManager(params, model) { }

  /// dummy virtual dtor
  virtual ~ForwardMeshManager() { }

  // SetupMe function is the same as in the base class 

  /// set the initial state of the mesh manager
  virtual void SetInitialState() 
  { 
    m_nIdx = 0;

    MeshManager::SetInitialState();
  }

  /**
     Test if we can still go ahead in time.

     By default, check the global time index.
   */
  virtual bool CanGoAhead() const
  {
    return m_nIdx < m_nNbTimes - 1;
  }

  /**
     Modify the rates so that constant can be an exact solution for 
     the discretized PDE.
   */
  virtual void SetupRates();


protected:

  /// Construct the time mesh
  virtual void ConstructTimeMesh(numeric::mesh::SpecialTimes &pdSpecialTimes);
 
  /// Setup the scheme type for each time step
  virtual void SetupSchemeTypes(numeric::mesh::SpecialTimes &pdSpecialTimes);
  
  /// Go ahead in time, by default, change the global time index
  virtual void GoAhead() { m_nIdx++; }

  /// Get the current time step  
  virtual double GetTimeStep() const
  {
    return m_pdTimes[m_nIdx] - m_pdTimes[m_nIdx - 1];
  }


private:

  NO_COPY_CLASS(ForwardMeshManager);

}; // class ForwardMeshManager


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_FORWARDMESHMANAGER_H_

