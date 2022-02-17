/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/resetparams.h
// Purpose:     reset param class
// Author:      Yann and David
// Created:     2004/11/03
// RCS-ID:      $Id: resetparams.h,v 1.9 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/resetparams.h
    @brief reset params class

    Implementation of the reset params class
 */

#ifndef _ITO33_PRICING_RESETPARAMS_H_
#define _ITO33_PRICING_RESETPARAMS_H_

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"
#include "ito33/list.h"
#include "ito33/date.h"

#include "ito33/pricing/reset.h"
#include "ito33/pricing/cblikeparams.h"
#include "ito33/pricing/cbparams.h"

namespace ito33 
{

namespace pricing
{


/// Pricing parameters for resets
class ResetParams : public CBParams
{
public: 

  /**
      Ctor by Reset contract object. 

      Other members must be initialized by the SetXXX() functions.

      @param reset reference to Reset contract object
   */
  ResetParams(Reset& reset) : 
    CBParams(reset), m_reset(reset), m_clonedReset(0)
  {  
  }

  /**
      Ctor by Reset contract object. This Ctor is essentailly the same
      as above, except that the memory for the contract class will
      be managed internally through the autoptr.

      Other members must be initialized by the SetXXX() functions

      @param reset autoptr to Reset contract objet
   */
  ResetParams(AutoPtr<Reset> reset) : 
    CBParams(*reset), m_reset(*reset), m_clonedReset(reset)
  {  
  }


  /**
      Ctor by Reset contract and common financial and numerical datas. 

      @param reset reference to Reset contract objet
      @param sessionData reference to financial SessionData
      @param pNumParams numeric parameters
      @param pMeshParams parameters for mesh builder
   */
  ResetParams(Reset& reset,
              const finance::SessionData& sessionData,
              const shared_ptr<numeric::NumParams>& pNumParams,
              const shared_ptr<numeric::MeshParams>& pMeshParams)
    : CBParams(reset, sessionData, pNumParams, pMeshParams),
      m_reset(reset), m_clonedReset(0)
  {
  }

  /// virtual dtor for base class
  virtual ~ResetParams() { }

  /// @name implmentation of virtual functions
  //@{
  virtual void Init();

  virtual ResetParams* Clone() const;

  //@}

  /**
      Gets a reference to the Reset contract.

      @return reset
   */
  Reset& GetReset() const { return m_reset; }
      
  /**
      Constructs the list of reset events and adds them to the 
       basic event manager.
   */
  void ConstructResetEvents();

  /**
      Indicates if a reset date is specified.

      @return true/false if reset date exists/ does not exist
   */
  bool HasActiveResetDate() const;
  
  /**
      Indicates if the reset is path dependent.

      @return true/false if problem is path dependent
   */
  virtual bool IsPathDependent() const;
  
  /**
      Gets the first active reset time.

      @return the first active reset time
   */
  double GetFirstActiveResetTime() const;

  /// Construct path dependent grid
  virtual std::vector<double> 
  ConstructPathDepGrid(Model& model) const;
  
  /// Construct the path based on the grid
  virtual std::vector< AutoPtr<CBLikeParams> >
  ConstructPaths(const std::vector<double>& grid) const;

  /// Construct the path dependent events
  virtual std::list< shared_ptr<PathDepEvent> >
  ConstructPathDepEvents
  (const std::vector< AutoPtr<CBLikeParams> >& paths) const;

  /// Determine the path to save
  virtual size_t GetPathToSave(const std::vector<double>& grid) const;

  /// initialize the different paths
  virtual void InitPaths(PathDepStructure& pathDepStruct);

protected:
  
  /// reference to the underlying reset contract
  Reset& m_reset;

  /// If the object is cloned, then the clone needs to manage memory
  AutoPtr<Reset> m_clonedReset;


private:

  NO_COPY_CLASS(ResetParams);

}; // class ResetParams


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_RESETPARAMS_H_
