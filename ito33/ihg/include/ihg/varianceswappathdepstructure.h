/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/varianceswappathdepstructure.h
// Purpose:     variance swap path dependent structure class
// Created:     2006/02/21
// RCS-ID:      $Id: varianceswappathdepstructure.h,v 1.5 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/varianceswappathdepstructure.h
    @brief variance swap path dependent structure class.
 */

#ifndef _ITO33_IHG_VARIANCESWAPPATHDEPSTRUCTURE_H_
#define _ITO33_IHG_VARIANCESWAPPATHDEPSTRUCTURE_H_

#include "ito33/beforestd.h"
#include <vector>
#include <list>
#include "ito33/afterstd.h"

#include "ito33/common.h"
#include "ito33/sharedptr.h"
#include "ito33/array.h"

#include "ito33/pricing/pathdepstructure.h"

namespace ito33
{

namespace finance { class ITO33_DLLDECL ComputationalFlags; }

namespace pricing
{
  class VarianceSwapParams;
  class VarianceSwapMeshManager;
}

namespace ihg
{

  class Model;
  class VarianceSwapNumOutput;
  class VarianceSwapInstData;

/// variance swap path dependent structure class
class VarianceSwapPathDepStructure : public pricing::PathDepStructure
{

public:

  /**
      All these classes are needed to construct the stepper, instdata, etc.
   */
  VarianceSwapPathDepStructure(
    const std::vector<double>& pdAvgSqrReturnGrid, 
    const std::vector<double>& pdPreviousSpotGrid,
    const std::vector<double>& pdMasterGrid,
    std::vector< AutoPtr<pricing::VarianceSwapParams> >& ppParams,
    Model& model,
    const finance::ComputationalFlags& flags,
    std::list< shared_ptr<pricing::PathDepEvent> >& pathDepEvents,
    size_t nPathToSave=0);

  virtual ~VarianceSwapPathDepStructure() {}

  /// @name virtual functions from base class
  //@{

  void InitPathToSave();

  void UpdateOutput(size_t nIdx);

  void UpdateOutputEndOfGrid(size_t ) {}

  void Finalize();

  //@}

  /**
      Gets the numerical output.  

      @return AutoPtr to varianceswapnumoutput for the requested path to save
   */
  AutoPtr<VarianceSwapNumOutput> GetOutput();

  /** 
      Prepares for timestepping by initializing the objects for each path.

      Used when computing the Greeks to re-initialize the structure

      @param nPathToSave the path to save for numerical output
   */
  void PrepareForTimestepping(size_t nPathToSave);

protected:

  /** 
      Constructs identical time meshes in all the paths. 

      @param ppMeshes (return) storage for the created mesh objects
      @param ppParams the param objects needed by mesh constructor
      @param pathDepEvents the path dependent events, needed for the times
      @param nNbPaths the number of paths
   */
  void ConstructIdenticalTimeMeshes( 
    shared_ptr< pricing::VarianceSwapMeshManager >* ppMeshes,    
    std::vector< AutoPtr<pricing::VarianceSwapParams> >& ppParams, 
    const pricing::Model& model,  
    std::list< shared_ptr<pricing::PathDepEvent> >& pathDepEvents,
    size_t nNbPaths);

  /// The numerical output
  AutoPtr<VarianceSwapNumOutput> m_pNumOutput;

  /// instdata to update the output. Need the exact type, so can't use m_path.  
  VarianceSwapInstData* m_pInstDataToSave;

  /// The computational flags (ie. what data to compute)
  const finance::ComputationalFlags& m_flags;

private:
  
  /// Initializes the path dependent structure
  void InitPathDepStructure(
    Array< shared_ptr<pricing::VarianceSwapMeshManager> >& ppMeshes,
    const std::vector<double>& pdAvgSqrReturnGrid,
    const std::vector<double>& pdPreviousSpotGrid,
    std::vector< AutoPtr<pricing::VarianceSwapParams> >& ppParams,
    Model& model,
    std::list< shared_ptr<pricing::PathDepEvent> >& pathDepEvents,
    size_t nPathToSave=0);

  NO_COPY_CLASS(VarianceSwapPathDepStructure);

};

} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_IHG_VARIANCESWAPPATHDEPSTRUCTURE_H_
