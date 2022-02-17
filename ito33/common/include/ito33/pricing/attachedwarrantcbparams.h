/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/attachedwarrantcbparams.h
// Purpose:     attached warrant convertible bond
// Author:      Ito33
// Created:     2005/01/13
// RCS-ID:      $Id: attachedwarrantcbparams.h,v 1.10 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/attachedwarrantcbparams.h
    @brief sharedependent warrant params class

    Implementation of the attached warrant convertible bond params class
 */

#ifndef _ITO33_PRICING_ATTACHEDWARRANTCBPARAMS_H_
#define _ITO33_PRICING_ATTACHEDWARRANTCBPARAMS_H_

#include "ito33/autoptr.h"

#include "ito33/pricing/attachedwarrantcb.h"
#include "ito33/pricing/cblikeparams.h"
#include "ito33/pricing/cbparams.h"

namespace ito33 
{

namespace pricing
{


/**
    Attached warrant CB  params class

   Pricing parameters for attached warrant cb
 */
class AttachedWarrantConvertibleBondParams : public CBParams
{
public: 

  /**
     Ctor by attached warrant cb contract object. 

     Other members must be initialized by the SetXXX() functions

     @param warrant reference to the attached warrant  cb contract object
   */
  AttachedWarrantConvertibleBondParams(AttachedWarrantConvertibleBond& warrant) 
     : CBParams(warrant), m_warrant(warrant), m_clonedWarrant(0)
  {  
  }

  /**
     Ctor by attached warrant cb contract object. This Ctor is essentially the same
     as above, except that the memory for the contract class will
     be managed internally through the autoptr.

     Other members must be initialized by the SetXXX() functions

     @param warrant reference to attached warrant cb contract object
   */
  AttachedWarrantConvertibleBondParams
    ( AutoPtr<AttachedWarrantConvertibleBond> warrant)  
    : CBParams(*warrant), m_warrant(*warrant), m_clonedWarrant(warrant)
  {  
  }

  /**
     Ctor  contract and common financial and numerical datas. 

     @param warrant reference to attached warrant cb contract object
     @param sessionData reference to financial SessionData
     @param pNumParams numeric parameters
     @param pMeshParams parameters for mesh builder
   */
  AttachedWarrantConvertibleBondParams(
              AttachedWarrantConvertibleBond& warrant,
              const finance::SessionData& sessionData,
              const shared_ptr<numeric::NumParams>& pNumParams,
              const shared_ptr<numeric::MeshParams>& pMeshParams)
    : CBParams(warrant, sessionData, pNumParams, pMeshParams),
      m_warrant(warrant), m_clonedWarrant(0)
  {
  }

  /// virtual dtor for base class
  virtual ~AttachedWarrantConvertibleBondParams() { }

  /// @name implmentation of virtual functions
  //@{

  virtual AttachedWarrantConvertibleBondParams* Clone() const;

  //@}

  /**
     Gets a reference to the contract

     @return reference to contract
   */
  AttachedWarrantConvertibleBond& GetAttachedWarrantCB() const 
  { 
    return m_warrant; 
  }

  /**
     Gets the reset time at which the conversion ratio is fixed.

     @return the reset time
   */
  double GetResetTime() const 
  { 
    return m_warrant.GetResetTime(); 
  }

  /**
     Gets the current conversion ratio.

     @return the current conversion ratio
  */
  double GetCurrentConversionRatio() const
  {
    return m_warrant.GetCurrentConversionRatio();
  }

  /**
     Indicate whether or not a reset time has been set

     @return true/false if a reset time has been set/not set
   */
  bool HasResetTime() const
  {
    return m_warrant.HasResetTime();
  }

 
  /// Does it have a put and call date on the reset date
  bool HasPutCallReset() const;

  virtual bool IsPathDependent() const;
  
  virtual std::vector<double> 
  ConstructPathDepGrid(Model& model) const;
  
  virtual std::vector< AutoPtr<CBLikeParams> >
  ConstructPaths(const std::vector<double>& grid) const;

  virtual std::list< shared_ptr<PathDepEvent> >
  ConstructPathDepEvents
  (const std::vector< AutoPtr<CBLikeParams> >& paths) const;

  virtual size_t GetPathToSave(const std::vector<double>& grid) const;

  virtual void InitPaths(PathDepStructure& pathDepStruct);


protected:
  
  /// reference to the underlying contract
  AttachedWarrantConvertibleBond& m_warrant;

  /// If the object is cloned, then the clone needs to manage memory
  AutoPtr<AttachedWarrantConvertibleBond> m_clonedWarrant;


private:

  NO_COPY_CLASS(AttachedWarrantConvertibleBondParams);

}; // class AttachedWarrantConvertibleBondParams


} // namespace pricing

} // namespace ito33

#endif // #define _ITO33_PRICING_ATTACHEDWARRANTCBPARAMS_H_
