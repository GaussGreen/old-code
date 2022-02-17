/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/underlyingprocess.h
// Purpose:     inhomogeneous underlying process class
// Created:     2006/06/01
// RCS-ID:      $Id: underlyingprocess.h,v 1.2 2006/08/20 09:47:41 wang Exp $
// Copyright:   (c) 2004-2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/underlyingprocess.h
    @brief class for inhomogeneous underlying process
 */

#ifndef _ITO33_IHG_UNDERLYINGPROCESS_H_
#define _ITO33_IHG_UNDERLYINGPROCESS_H_

#include "ito33/sharedptr.h"

#include "ito33/ihg/common.h"

#include "ito33/finance/underlyingprocess.h"

namespace ito33
{

namespace XML
{
  class Tag;
}

namespace ihg
{

  class ITO33_IHG_DLLDECL Volatility;
  class ITO33_IHG_DLLDECL HazardRate;

/**
   Class describes an inhomogeneous underlying process.
 */
class ITO33_IHG_DLLDECL UnderlyingProcess : public finance::UnderlyingProcess
{
  
public:

  /// Ctor
  UnderlyingProcess() {}

  /**
     Ctor defines the UnderlyingProcess with necessary parameters.

     @param pVolatility shared pointer to volatility
     @param pHazardRate shared pointer to hazard rate
  */
  UnderlyingProcess(const shared_ptr<Volatility>& pVolatility, 
    const shared_ptr<HazardRate>& pHazardRate):
    m_pVolatility(pVolatility), m_pHazardRate(pHazardRate)
  {
    CheckAll();
  }
  
  /**
      @name Setters.
   */
  //@{

  /**
     The volatility of the underlying process.

     @param pVolatility shared pointer to a Volatility
   */
  void SetVolatility(const shared_ptr<Volatility>& pVolatility)
  {
    m_pVolatility = pVolatility;

    CheckVolatility();
  }

  /**
     The hazard rate of the underlying process.

     @param pHazardRate shared pointer to a HazardRate
   */
  void SetHazardRate(const shared_ptr<HazardRate>& pHazardRate)
  {
    m_pHazardRate = pHazardRate;

    CheckHazardRate();
  } 

  //@} // name Setters

  
  /**
      @name Accessors.
   */
  //@{

  /**
     The volatility of the underlying process.

     @return shared pointer to the volatility
   */
  const shared_ptr<Volatility>& GetVolatility() const
  {
    CheckVolatility();

    return m_pVolatility;
  }

  /**
     The hazard rate of the underlying process.

     @return shared pointer to the hazard rate
   */
  const shared_ptr<HazardRate>& GetHazardRate() const
  {
    CheckHazardRate();

    return m_pHazardRate;
  }

  //@} // name Accessors

  /**
    @internal

    @noexport
   */
  void Dump(ito33::XML::Tag& tagParent) const;

  /**
     @internal
     
     @brief checks if model parameters are valid.

     @noexport
   */
  void CheckAll() const
  {
    CheckVolatility();
    CheckHazardRate();
  }

protected:

  /// checks if Volatility is valid
  void CheckVolatility() const;

  /// checks if HazardRate is valid
  void CheckHazardRate() const;

private:

  /// The volatility of the underlying 
  shared_ptr<Volatility> m_pVolatility;

  /// The hazard rate of the underlying 
  shared_ptr<HazardRate> m_pHazardRate;
  
}; // class UnderlyingProcess


} // namespace ihg

} // namespace ito33

#endif // _ITO33_IHG_UNDERLYINGPROCESS_H_
