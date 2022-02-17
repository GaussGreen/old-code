/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/instdatawithconstraints.h
// Purpose:     backward instdata class for ihg projet
// Author:      Wang
// Created:     2004/02/13
// RCS-ID:      $Id: instdatawithconstraints.h,v 1.13 2005/12/30 11:51:14 nabil Exp $
// Copyright:   (c) 2003-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/instdatawithconstraints.h
    @brief ihg backward instdata (with constraints) class
 */

#ifndef _IHG_INSTDATAWITHCONSTRAINTS_H_
#define _IHG_INSTDATAWITHCONSTRAINTS_H_

#include "ito33/array.h"

#include "ihg/backwardinstdata.h"

namespace ito33
{
  
namespace pricing
{
  class Params;
  class Constraints;
}

namespace ihg
{


/// Base instdata class for ihg projets
class InstDataWithConstraints : public BackwardInstData
{

public:

  InstDataWithConstraints(pricing::Params& params, 
                          Model& model, 
                          pricing::MeshManager& meshes);

  virtual ~InstDataWithConstraints() { }

  /// @name Functions required by the Engine
  //@{
   
  // Init() : No implementation at this level, so no need to redeclare

  /**
     Update myself at the beginning of each time step
    
     Required by Engine class. At each time step, 
     this function update the instanous parameters from the
     meshes manager. It works also on the Data part.

     It is not same as BackwardInstData::UpdateBeforeStep(). In fact,
     when an event occurs at the last time step (or maturity), we'd
     better initialize flags to 0
   */
  virtual void UpdateBeforeStep();

  /**
     Get events from mesh manager and do events

     At the end of each time step, 
     this function check out from the mesh manager all events
     happen at current period, and apply the events' method to the Data part.
   */
  virtual void DoEvents();

  // SetInitialValue(): No implementation at this level, so no need to redeclare

  //@}

  virtual const int* GetConstraintFlags() const {return m_piFrozenFlags.Get();}

  /**
     Get a pointer to the internal constraint objet
   */
  const pricing::Constraints* GetConstraints() const;

  /// The constraints pointer
  const pricing::Constraints* m_pConstraints;

  /// Array indicates if there is constraints on a point
  Array<int> m_piFrozenFlags;

    
protected:

  /**
     Allocate the memory for the price and vega arrays

     @param nNbS the size of the arrays to be allocated
   */
  virtual void Alloc(size_t nNbS);


  /**
    virtual function called in DoEvents(). When we meet events
    (so prices will be changed) or constraints update, either
    prices or constraints are changed. Therefore, we have to
    apply the constraints.

    The default implmentation of this function considers that
    there constraints stay the same. Specific InstData class
    must implment the function otherly if constraints can be modified.
    @sa CBInstData
    */
  virtual void ApplyConstraintsAfterEventsOrConstraintsUpdate();

private:

  NO_COPY_CLASS(InstDataWithConstraints);

}; // class InstDataWithConstraints


} // namespace ihg

} // namespace ito33 

#endif // #ifndef _IHG_INSTDATAWITHCONSTRAINTS_H_

