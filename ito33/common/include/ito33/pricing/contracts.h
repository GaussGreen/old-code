/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/contracts.h
// Purpose:     a generic contracts class
// Created:     2004/02/11
// RCS-ID:      $Id: contracts.h,v 1.18 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/contracts.h
    @brief The declaration of the contracts class.

    The base class for contracts.  
 */

#ifndef _ITO33_PRICING_CONTRACTS_H_
#define _ITO33_PRICING_CONTRACTS_H_

#include "ito33/sharedptr.h"
#include "ito33/numeric/mesh/specialtimes.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL YieldCurve;
}

namespace numeric
{
  namespace mesh
  {
    class SpecialTimes;
  }
}

namespace pricing
{


/// The base contract class.
class Contracts 
{
public:

  /// Dummy empty ctor
  Contracts() : m_bIsCrossCurrency(false) { }

  /// Dummy virtual dtor for base class
  virtual ~Contracts() { }

  /**
      Gets the maturity (as a fraction of year) of the contracts.
    
      @return the maturity (as a fraction of year) of the contract
   */
  virtual double GetMaturityTime() const = 0;

  /**
      Gets the special times of the contracts.

      @param specialTimes (output) the resulted special times of the contract
   */
  virtual void GetSpecialTimes(numeric::mesh::SpecialTimes& specialTimes) const
  {
    specialTimes.clear();
  }

  /**
      Gets shared ptr to derivative curve object.

      @return derivative curve
   */
  const shared_ptr<finance::YieldCurve>& GetDerivativeCurve() const 
  { 
    return m_pDerivativeCurve; 
  }  

  /**
      Sets the derivative curve.

      @param pDerivativeCurve derivative curve
     
      @remark{This function is required for cross currency instrument, when a
              contract is created from scratch}
   */
  void 
  SetDerivativeCurve(const shared_ptr<finance::YieldCurve>& pDerivativeCurve)
  {
    m_pDerivativeCurve = pDerivativeCurve;
  }
  
  /**
      Indicates if the instrument is a cross currency one.

      @return true if the instrument is a cross currency one, false otherwise
   */
  bool IsCrossCurrency() const { return m_bIsCrossCurrency; }
  

protected:

  /// if the instrument is a cross-currency one
  bool m_bIsCrossCurrency; 

  /// Derivative curve
  shared_ptr<finance::YieldCurve> m_pDerivativeCurve;

}; // class Contracts;


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_CONTRACTS_H_
