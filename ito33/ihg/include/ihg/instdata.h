/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/instdata.h
// Purpose:     option instdata class
// Author:      David Pooley
// Created:     2003/12/10
// RCS-ID:      $Id: instdata.h,v 1.13 2005/12/30 11:51:14 nabil Exp $
// Copyright:   (c) 2003-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/instdata.h
    @brief ihg instdata class
 */

#ifndef _IHG_INSTDATA_H_
#define _IHG_INSTDATA_H_

#include "ito33/array.h"

#include "ito33/pricing/instdata.h"

#include "ito33/numeric/boundary1d.h"

namespace ito33
{

namespace ihg
{

  class Model;
  
/// Base instdata class for ihg projets
class InstData : public pricing::InstData
{

public:

  InstData(pricing::Params& params,
           Model& model,
           pricing::MeshManager& meshes)
         : pricing::InstData(params, meshes), m_model(model)
  { }

  virtual ~InstData() { }

  /// @name Functions required by the Engine
 
  //@{

  // Init() : No implementation at this level, so no need to redeclare

  virtual void UpdateBeforeStep();

  // DoEvents(): Same implementation as in base, no need to reimplement

  // SetInitialValue(): No implementation at this level, so no need to redeclare

  //@}
  
  /// @name IHG model specific arrays
  //@{

  /// Current volatilities
  Array<double> m_pdVolsSquared;
 
  /// Current hazard rate values
  Array<double> m_pdHazardRates;
  
  //@}

  /**
     Gets the boundary condition
   */
  const numeric::Boundary1D &GetBoundaryCondition() const
  {
    return m_BoundaryCondition;
  }

protected:
 

  /**
     Allocate the memory for the price arrays

     @param nNbS the size of the arrays to be allocated
   */
  virtual void Alloc(size_t nNbS);
  
  /**
     Apply an event to prices
     
     @param pEvent A pointer to the to be applied event
   */
  virtual void ApplyEvent(const pricing::Event *pEvent);

  /// A reference to a ihg model
  Model& m_model;

  /// the boundary condition
  numeric::Boundary1D m_BoundaryCondition;
  // before, we got boundary condition from param, but now it is just
  // set to default in instdata ctor


private:

  typedef ito33::pricing::InstData BaseClass;

  NO_COPY_CLASS(InstData);

}; // class InstData


} // namespace ihg

} // namespace ito33 

#endif // #ifndef _IHG_INSTDATA_H_

