/////////////////////////////////////////////////////////////////////////////
// Name:        hg/cbpathdepstructure.h
// Purpose:     cb path dependent structure class with HG model
// Created:     2005/04/11
// RCS-ID:      $Id: cbpathdepstructure.h,v 1.6 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/cbpathdepstructure.h
    @brief path dependent structure class.
 */

#ifndef _HG_CBPATHDEPSTRUCTURE_H_
#define _HG_CBPATHDEPSTRUCTURE_H_

#include "ito33/beforestd.h"
#include <vector>
#include <list>
#include "ito33/afterstd.h"

#include "ito33/common.h"
#include "ito33/sharedptr.h"

#include "ito33/pricing/cb.h"
#include "ito33/pricing/cbmeshmanager.h"
#include "ito33/pricing/pathdepstructure.h"

namespace ito33
{

namespace pricing
{
  class CBLikeParams;
}

namespace numeric
{
 class MeshParams;
 class NumParams;
}

namespace hg
{

  class Model;
  class CBNumOutput;
  class CBInstData;

/// CB Path dependent structure class.
class CBPathDepStructure :public pricing::PathDepStructure
{

public:

  /*
    All these classes are needed to construct the stepper, instdata, etc.

  */
  CBPathDepStructure(
    std::vector<double>& pdGridY,
    std::vector< AutoPtr<pricing::CBLikeParams> >& ppParams,
    Model& model,
    const finance::ComputationalFlags& flags,
    std::list< shared_ptr<pricing::PathDepEvent> >& pathDepEvents,
    size_t nPathToSave=0);

  virtual ~CBPathDepStructure() {};

  /**
     Virtual functions
  */
  /// @name virtual functions from base class
  //@{
  void UpdateOutput(size_t nIdx);

  void SetPathToSave(size_t nIdx);

  void Finalize();
  //@}

  /**
     Get the numerical output.  

     @return AutoPtr to cbnumoutput for the requested path to save
  */
  AutoPtr<CBNumOutput> GetNumOutput();

  /** 
    Prepare for timestepping by initializing the objects for each path

    Used when computing the Greeks to re-initialize the structure

    @param nPathToSave the path to save for numerical output
   */
  void PrepareForTimestepping(size_t nPathToSave);


protected:

  /** 
     Construct identical time meshes in all the paths. 

     Note: this function doesn't use anything specific to CB, so it can
     probably be implemented by using only base objets.

     @param ppMeshes (return) storage for the created mesh objects
     @param ppParams the param objects needed by mesh constructor
     @param pathDepEvents the path dependent events, needed for the times
     @param nNbPaths the number of paths
   */
  void ConstructIdenticalTimeMeshes( 
    shared_ptr< pricing::CBMeshManager >* ppMeshes,    
    std::vector< AutoPtr<pricing::CBLikeParams> >& ppParams, 
    const pricing::Model& model,  
    std::list< shared_ptr<pricing::PathDepEvent> >& pathDepEvents,
    size_t nNbPaths);

  /// The numerical output
  AutoPtr<CBNumOutput> m_pNumOutput;

  /// instdata to update the output. Need the exact type, so can't use m_path.  
  CBInstData* m_pInstDataToSave;

  /// Path to save for output
  size_t m_nPathToSave;

  /// The computational flags (ie. what data to compute)
  const finance::ComputationalFlags& m_flags;

private:

  NO_COPY_CLASS(CBPathDepStructure);
};

} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_CBPATHDEPSTRUCTURE_H_
