/////////////////////////////////////////////////////////////////////////////
// Name:        hg/instdata.h
// Purpose:     HG base instdata class
// Created:     2005/01/13
// RCS-ID:      $Id: instdata.h,v 1.5 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/instdata.h
   @brief HG base instdata class
 */

#ifndef _HG_INSTDATA_H_
#define _HG_INSTDATA_H_

#include "ito33/array.h"

#include "ito33/pricing/instdata.h"

#include "ito33/numeric/boundary1d.h"

namespace ito33
{

namespace hg
{

  class Model;
  class Payoff;
  
/// HG base instdata class
class InstData : public pricing::InstData
{

public:

  InstData(pricing::Params& params,
           Model& model,
           pricing::MeshManager& meshes);

  virtual ~InstData() { }

  /**
     Gets the boundary condition
   */
  const numeric::Boundary1D &GetBoundaryCondition() const
  {
    return m_BoundaryCondition;
  }

  virtual void UpdateBeforeStep();

  /**
     Apply the boundary condition(primary Dirichlet condition)
     to the right hand side.

     @param pdRHS The right hand side 
   */
  void ApplyBoundaryConditionToSensitivityRHS(double* pdRHS) const;

  /**
      Sets an external payoff.

      Default to doing nothing.  Should be re-implemented by cb-like 
      instdata classes to support call notice. 
  */
  virtual void SetOutsideInitialValue(const shared_ptr<Payoff>&)
  {
  }

  const Model& GetModel() const { return m_model; }

  /// The number of regimes
  size_t m_nNbRegimes;

  /// The dimension of the (non) linear system
  size_t m_nNbX;


protected:
  
  /**
     Apply an event to prices
     
     @param pEvent A pointer to the to be applied event
   */
  virtual void ApplyEvent(const pricing::Event *pEvent);

  /// A reference to a HG model
  Model& m_model;

  /// the boundary condition
  numeric::Boundary1D m_BoundaryCondition;
  // before, we got boundary condition from param, but now it is just
  // set to default in instdata ctor


private:

  NO_COPY_CLASS(InstData);

}; // class InstData


} // namespace hg

} // namespace ito33 

#endif // #ifndef _HG_INSTDATA_H_

