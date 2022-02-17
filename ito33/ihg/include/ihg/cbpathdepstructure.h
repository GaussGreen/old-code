/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/cbpathdepstructure.h
// Purpose:     cb path dependent structure class
// Author:      Yann and David
// Created:     18/08/2004
// RCS-ID:      $Id: cbpathdepstructure.h,v 1.20 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/cbpathdepstructure.h
    @brief path dependent structure class.

*/

#ifndef _ITO33_IHG_CBPATHDEPSTRUCTURE_H_
#define _ITO33_IHG_CBPATHDEPSTRUCTURE_H_

#include "ito33/beforestd.h"
#include <vector>
#include <list>
#include "ito33/afterstd.h"
#include "ito33/dlldecl.h"

#include "ito33/common.h"
#include "ito33/sharedptr.h"

#include "ito33/pricing/cb.h"
#include "ito33/pricing/cbmeshmanager.h"
#include "ito33/pricing/pathdepstructure.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL ComputationalFlags;
}

namespace pricing
{
  class CBLikeParams;
}

namespace numeric
{
 class MeshParams;
 class NumParams;
}

namespace ihg
{
  class Model;
  class CBNumOutput;
  class CBInstData;

/// CB Path dependent structure class
class CBPathDepStructure :public pricing::PathDepStructure
{

public:

  /*
    All these classes are needed to construct the stepper, instdata, etc.
  */
  CBPathDepStructure(
    std::vector<double>& pdGridY,
    std::vector< AutoPtr<pricing::CBLikeParams> >& ppParams,
    ihg::Model& model,
    const finance::ComputationalFlags& flags,
    std::list< shared_ptr<pricing::PathDepEvent> >& pathDepEvents,
    size_t nPathToSave=0);

  virtual ~CBPathDepStructure() {};

  /// @name virtual functions from base class
  //@{
  void InitPathToSave();

  void UpdateOutput(size_t nIdx);

  void UpdateOutputEndOfGrid(size_t nIdx);
  
  void Finalize();
  //@}

  /**
      Gets the numerical output.  

      @return AutoPtr to cbnumoutput for the requested path to save
   */
  AutoPtr<CBNumOutput> GetOutput();

  /** 
      Prepares for timestepping by initializing the objects for each path.

      Used when computing the Greeks to re-initialize the structure

      @param nPathToSave the path to save for numerical output
   */
  void PrepareForTimestepping(size_t nPathToSave);

  /**
      Sets an external initial condition/payoff.

      @param pPayoff The external initial condition
   */
  void SetInitialValue(const shared_ptr<finance::Payoff>& pPayoff)
  {
    m_pPayoff = pPayoff;
  }

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
  ihg::CBInstData* m_pInstDataToSave;

  /// The computational flags (ie. what data to compute)
  const finance::ComputationalFlags& m_flags;

  /// External payoff (eg. used for call notice)
  shared_ptr<finance::Payoff> m_pPayoff;

private:

  NO_COPY_CLASS(CBPathDepStructure);
};

} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_IHG_CBPATHDEPSTRUCTURE_H_
