/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/meshmanager.h
// Purpose:     base mesh manager for PDE problems
// Author:      David Pooley
// Created:     2004/01/05
// RCS-ID:      $Id: meshmanager.h,v 1.48 2006/07/06 14:47:19 dave Exp $
// Copyright:   (c) 2003-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_PRICING_MESHMANAGER_H_
#define _ITO33_PRICING_MESHMANAGER_H_

#include "ito33/common.h"
#include "ito33/array.h"

#include "ito33/numeric/schemetype.h"

namespace ito33
{

namespace numeric
{
  namespace mesh
  {
    class SpecialTimes;
  }
}

namespace pricing
{

  class Params;
  class Model;

/// Base class for space and time grid management
class MeshManager
{
public:
  
  /** 
    Creates MeshManager object by reference to Params and Model objects.

    @param params reference to Params object
    @param model reference to Model object
    */
  MeshManager(Params& params, Model& model)
            : m_params(params),
              m_model(model),
              m_bUseOutsideTimeMesh(false)
  {
  }

  /// Dummy virtual dtor 
  virtual ~MeshManager() { }

  /**  
     Gets some instantaneous data needed by instdata.

     @param dTimeStep the current time step
     @param dRate the current domestic rate
     @param dDerivativeRate the rate of the derivative curve
     @param dForeignRate the current foreign rate
     @param schemeType the scheme type for current time step
   */
  void GetInstValues(double& dTimeStep,
                     double& dRate, double& dDerivativeRate,
                     double& dForeignRate,
                     numeric::SchemeType& schemeType) const;
  
  /**
    Gets the current time.

    @return current time
    */
  double GetTime() const { return m_pdTimes[m_nIdx]; }

  /**
    Test if we can still go ahead in time.
    
    numoutput may need it.

    @return true if we can still go ahead in time, false if not.
   */
  virtual bool CanGoAhead() const = 0;

  /// @name Functions for creating time mesh
  //@{

  /**
    Constructs time mesh and related information by copying. 

    This function is often used for path dependent problem, where more than
    one problem use same time mesh. Note that we also need to copy scheme types

    @param pdTimeMesh given array of time mesh points
    @param pSchemeTypes given array of scheme types
    @param nNbTimesteps size of time mesh
    */
  void SetupTimeMesh(const double* pdTimeMesh, 
                     const numeric::SchemeType* pSchemeTypes, 
                     size_t nNbTimesteps);

  /**
    Constructs time mesh by given special time points 
    and setup related information, such as scheme types at each time point.

    @param specialTimes special time points to be put in time mesh.
    */
  void SetupTimeMesh(numeric::mesh::SpecialTimes& specialTimes);

  //@} // end of name Functions for creating time mesh

  /// @name Functions required by Engine
  //@{
 
  /// Set up the meshes and the data. At this level, time mesh is constructed
  virtual void SetupMe();

  /// Try if another timestep can be made
  bool TryGoAhead();

  /// Reset the manager to the initial state
  virtual void SetInitialState();
  
  //@}

  /**
     Modify the rates so that constant can be an exact solution for 
     the discretized PDE.
   */
  virtual void SetupRates() = 0;

 
  /**
    Gets the time mesh.

    @return pointer to time mesh array
    */
  const double* GetTimes() const
  {
    return m_pdTimes.Get();
  }

  /**
    Gets the number of times.

    @return the number of points in time mesh
    */
  size_t GetNbTimes() const
  {
    return m_nNbTimes;
  }

  /**
    Gets the scheme types.

    @return pointer to scheme types array
    */
  const numeric::SchemeType* GetSchemeTypes() const
  {
    return m_pSchemeTypes.Get();
  }

  /**
    Gets the number of times requested by params. 

    @return the number of times requested by params
    */
  virtual size_t GetNbRequestedTimes() const;


protected:
 
  /// reference to associated Params object
  Params& m_params;

  /// reference to associated Model object
  Model& m_model;

  /// @name internal functions for creating time mesh
  //@{

  /**
    Constructs the time mesh.

    Implementation of this function in sepcific class can use
    ConstructUniformTimeMesh() or ConstructNonuniformTimeMesh()
    already defined in base MeshManager class.

    @param specialTimes special time points to be put in time mesh.
    */
  virtual void ConstructTimeMesh(numeric::mesh::SpecialTimes& specialTimes) = 0;
 
  /**
    Construct a uniform time mesh (or as close to uniform as possible).
    
    @param specialTimes special time points to be put in time mesh.
    */
  void ConstructUniformTimeMesh(numeric::mesh::SpecialTimes& specialTimes);

  /**
    Construct a non-uniform time mesh in requested direction.
    
    @param specialTimes special time points to be put in time mesh.
    @param direction 0 is backward, forward for other value
    */
  void ConstructNonuniformTimeMesh(numeric::mesh::SpecialTimes& specialTimes,
                                   int direction);

  //@} // end of name internal functions for creating time mesh

  /// Setup the scheme type for each time step
  virtual void SetupSchemeTypes(numeric::mesh::SpecialTimes& specialTimes) = 0;

  /// Go ahead in time
  virtual void GoAhead() = 0;
  
  /// Gets the time step  
  virtual double GetTimeStep() const = 0;

  /// time array
  Array<double> m_pdTimes;

  /// number of points in time mesh
  size_t m_nNbTimes;

  /// precalculated zero rates for all time points
  Array<double> m_pdRates;

  /**
    precalculated zero rates related to currency of derivative
    for all time points.
    */
  Array<double> m_pdDerivativeRates;

  /// precalculated zero borrow rates for all time points
  Array<double> m_pdForeignRates;

  /// precalculated scheme types used for all time steps
  Array<numeric::SchemeType> m_pSchemeTypes;

  /// The current index at the time mesh
  size_t m_nIdx;

  /// if the time mesh is set up using outside datas
  bool m_bUseOutsideTimeMesh;


private:

  NO_COPY_CLASS(MeshManager);

}; // class MeshManager


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_MESHMANAGER_H_
