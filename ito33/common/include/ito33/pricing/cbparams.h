/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/cbparams.h
// Purpose:     cb params class
// Author:      Laurence
// Created:     2004/03/15
// RCS-ID:      $Id: cbparams.h,v 1.88 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/cbparams.h
    @brief cb params class

    Implementation of the params class for CB.
 */

#ifndef _ITO33_PRICING_CBPARAMS_H_
#define _ITO33_PRICING_CBPARAMS_H_

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"
#include "ito33/list.h"
#include "ito33/date.h"

#include "ito33/pricing/cb.h"
#include "ito33/pricing/cblikeparams.h"


namespace ito33 
{

namespace pricing
{

/// Pricing parameters for Convertible Bonds
class CBParams : public CBLikeParams
{
public: 

  /**
      Ctor by CB contract objet. 

      Other members must be initialized by the SetXXX() functions.

      @param cb reference to CB contract objet
   */
  CBParams(CB& cb) : CBLikeParams(cb), m_cb(cb), m_clonedCB(0)
  {  
  }

  /**
      Ctor by CB contract objet. This Ctor is essentailly the same
      as above, except that the memory for the contract class will
      be managed internally through the autoptr.

      Other members must be initialized by the SetXXX() functions

      @param cb autoptr to CB contract objet
   */
  CBParams( AutoPtr<CB> cb) : CBLikeParams(*cb), m_cb(*cb), m_clonedCB(cb)
  {  
  }


  /**
      Ctor by CB contract and common financial and numerical datas. 

      @param cb reference to CB contract objet
      @param sessionData reference to financial Session
      @param pNumParams numeric parameters
      @param pMeshParams parameters for mesh builder
   */
  CBParams(CB& cb,
           const finance::SessionData& sessionData,
           const shared_ptr<numeric::NumParams>& pNumParams,
           const shared_ptr<numeric::MeshParams>& pMeshParams)
    : CBLikeParams(cb, sessionData, pNumParams, pMeshParams),
      m_cb(cb), m_clonedCB(0)
  {
  }

  /// virtual dtor for base class
  virtual ~CBParams() { }

  /// @name implmentation of virtual functions
  //@{
  virtual void Init();

  virtual CBParams* Clone() const;

  //@}

  /// Gets a reference to the CB contract
  CB& GetCB() const { return m_cb; }
      
  virtual CBParams* GetCallNoticeParams();


protected:
  
  CB& m_cb;

  /// If the object is cloned, then the clone needs to manage memory
  AutoPtr<CB> m_clonedCB;


private:

  NO_COPY_CLASS(CBParams);

}; // class CBParams


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_CBPARAMS_H_
