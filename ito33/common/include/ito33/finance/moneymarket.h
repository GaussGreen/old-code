/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/moneymarket.h
// Purpose:     Contingency class declaration
// Created:     11.12.2003
// RCS-ID:      $Id: moneymarket.h,v 1.29 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2003 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/moneymarket.h
    @brief Declaration of the MoneyMarket class.

    This class represents the money market, including a currency and a curve.
 */

#ifndef _ITO33_FINANCE_MONEYMARKET_H_
#define _ITO33_FINANCE_MONEYMARKET_H_

#include "ito33/sharedptr.h"
#include "ito33/dlldecl.h"

namespace ito33
{

namespace XML
{
  class Tag;
}


namespace finance
{

class ITO33_DLLDECL YieldCurve;
class ITO33_DLLDECL Numeraire;

/**
    The Money Market.
    
    It includes a currency and a curve.
 */
class ITO33_DLLDECL MoneyMarket
{

public:
  
  /**
      Constructs a MoneyMarket object by a currency and its yield curve.
   */
  MoneyMarket(const shared_ptr<Numeraire>& pCurrency,
              const shared_ptr<YieldCurve>& pYieldCurve);


  /**
      @name Modifiers.
   */
  //@{

  /**
      The yield curve of this MoneyMarket.

      @param pYieldCurve the yield curve for this MoneyMarket
   */
  void SetYieldCurve(const shared_ptr<YieldCurve>& pYieldCurve); 

  //@}

  /**
      @name Accessors.
   */
  //@{

  /**
      The currency of this MoneyMarket.

      @return the currency of this MoneyMarket 
   */
  const shared_ptr<Numeraire>& GetNumeraire() const 
  { 
    return m_pNumeraire; 
  }

  /**
      The YieldCurve of this MoneyMarket.

      @return the yield curve of this MoneyMarket
   */
  const shared_ptr<YieldCurve>& GetYieldCurve() const 
  { 
    return m_pYieldCurve; 
  }

  /**
      @internal
      @brief Gets the yield curve for mesh.

      @noexport
   */
  const shared_ptr<YieldCurve>& GetYieldCurveForMesh() const;
  
  //@}

  /**
      @internal

      @brief Sets the yield curve for mesh.

      @noexport
   */
  void SetYieldCurveForMesh(const shared_ptr<YieldCurve>& yc);

  /**
      @internal

      @noexport
   */
  void Dump(XML::Tag& tagParent) const;


private:
  
  /// the currency of the market
  shared_ptr<Numeraire> m_pNumeraire;

  /// the yield curve for the currency
  shared_ptr<YieldCurve> m_pYieldCurve;

  /// the yield curve used for mesh construction for the currency
  shared_ptr<YieldCurve> m_pYieldCurveForMesh;

}; // class MoneyMarket


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_MONEYMARKET_H_

