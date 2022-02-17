/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/optionpathdepstructure.h
// Purpose:     equity path dependent structure class
// Author:      ITO 33 Canada
// Created:     April 7, 2005
// RCS-ID:      $Id: optionpathdepstructure.h,v 1.7 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/optionpathdepstructure.h
    @brief path dependent structure class.

*/

#ifndef _ITO33_IHG_OPTIONPATHDEPSTRUCTURE_H_
#define _ITO33_IHG_OPTIONPATHDEPSTRUCTURE_H_

#include "ito33/beforestd.h"
#include <vector>
#include <list>
#include "ito33/afterstd.h"

#include "ito33/common.h"
#include "ito33/sharedptr.h"

#include "ito33/pricing/option.h"
#include "ito33/pricing/optionmeshmanager.h"
#include "ito33/pricing/pathdepstructure.h"

namespace ito33
{

namespace pricing
{
  class OptionParams;
}

namespace numeric
{
 class MeshParams;
 class NumParams;
}

namespace ihg
{

  class Model;
  class OptionNumOutput;
  class OptionInstData;

/**
   Option Path dependent structure class
*/
  class OptionPathDepStructure :public pricing::PathDepStructure
{

public:

  /*
    All these classes are needed to construct the stepper, instdata, etc.

  */
  OptionPathDepStructure(
    std::vector<double>& pdGridY,
    std::vector< AutoPtr<pricing::OptionParams> >& ppParams,
    ihg::Model& model,
    const finance::ComputationalFlags& flags,
    std::list< shared_ptr<pricing::PathDepEvent> >& pathDepEvents,
    size_t nPathToSave=0);

  
  OptionPathDepStructure(
    Array< shared_ptr<pricing::OptionMeshManager> > &ppMeshes,
    std::vector<double>& pdGridY,
    std::vector< AutoPtr<pricing::OptionParams> >& ppParams,
    ihg::Model& model,
    const finance::ComputationalFlags& flags,
    std::list< shared_ptr<pricing::PathDepEvent> >& pathDepEvents,
    size_t nPathToSave=0);

  virtual ~OptionPathDepStructure() {};

  /**
     Virtual functions
  */
  /// @name virtual functions from base class
  //@{
  void InitPathToSave();

  void UpdateOutput(size_t nIdx);

  /**
    Update the numerical output at the end of grid. 
    Must be done in derived class since type information is needed.
  */
  void UpdateOutputEndOfGrid(size_t ) {}

  void Finalize();
  //@}

  /**
     Get the numerical output.  

     @return AutoPtr to cbnumoutput for the requested path to save
  */
  AutoPtr<OptionNumOutput> GetOutput();

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
    shared_ptr< pricing::OptionMeshManager >* ppMeshes,    
    std::vector< AutoPtr<pricing::OptionParams> >& ppParams, 
    const pricing::Model& model,  
    std::list< shared_ptr<pricing::PathDepEvent> >& pathDepEvents,
    size_t nNbPaths);


  /// The numerical output
  AutoPtr<OptionNumOutput> m_pNumOutput;

  /// instdata to update the output. Need the exact type, so can't use m_path.  
  ihg::OptionInstData* m_pInstDataToSave;

  /// The computational flags (ie. what data to compute)
  const finance::ComputationalFlags& m_flags;

private:
  NO_COPY_CLASS(OptionPathDepStructure);


  /**
    Initialize the path dependent structure
  */
  void InitPathDepStructure(
    Array< shared_ptr<pricing::OptionMeshManager> > &ppMeshes,
    std::vector<double>& pdGridY,
    std::vector< AutoPtr<pricing::OptionParams> >& ppParams,
    ihg::Model& model,
    std::list< shared_ptr<pricing::PathDepEvent> >& pathDepEvents,
    size_t nPathToSave=0);


};

} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_IHG_OPTIONPATHDEPSTRUCTURE_H_

