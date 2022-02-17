/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/backwardmeshmanager.h
// Purpose:     base mesh manager for backward PDE problems
// Author:      Wang
// Created:     2004/02/11
// RCS-ID:      $Id: backwardmeshmanager.h,v 1.13 2006/01/03 17:18:48 zhang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_PRICING_BACKWARDMESHMANAGER_H_
#define _ITO33_PRICING_BACKWARDMESHMANAGER_H_

#include "ito33/pricing/meshmanager.h"

namespace ito33
{

namespace pricing
{


/// Backward time mesh manager
class BackwardMeshManager : public MeshManager
{

public:
  
  /**
     Ctor, call the ctor of the base class MeshManager

     @param params The underlying contract paramaters
     @param model The model used for the pricing
   */  
  BackwardMeshManager(Params &params, Model &model) 
                    : MeshManager(params, model) { }

  /// dummy virtual dtor
  virtual ~BackwardMeshManager() { }

  // SetupMe function is the same as in the base class 

  /// set the initial state of the mesh manager
  virtual void SetInitialState() 
  { 
    m_nIdx = m_nNbTimes - 1;  

    MeshManager::SetInitialState();
  }

  /**
     Test if we can still go ahead in time.
     
     By default, check the global time index.
   */
  virtual bool CanGoAhead() const
  {
    return m_nIdx > 0;
  }

  /**
     Modify the rates so that constant can be an exact solution for 
     the discretized PDE.
   */
  virtual void SetupRates();


protected:

  /// Construct the time mesh
  virtual void ConstructTimeMesh(numeric::mesh::SpecialTimes& specialTimes);
 
  /// Setup the scheme type for each time step
  virtual void SetupSchemeTypes(numeric::mesh::SpecialTimes& specialTimes);
  
  /// Go ahead in time. By default, change the global time index
  virtual void GoAhead() { m_nIdx--; }

  /// Get the current time step  
  virtual double GetTimeStep() const
  {
    return m_pdTimes[m_nIdx + 1] - m_pdTimes[m_nIdx];
  }


private:

  NO_COPY_CLASS(BackwardMeshManager);

}; // class BackwardMeshManager


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_BACKWARDMESHMANAGER_H_

