/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/equity.h
// Purpose:     Equity class declaration
// Author:      Wang
// Created:     11.12.2003
// RCS-ID:      $Id: equity.h,v 1.32 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/equity.h
    @brief Declaration of the Equity class.

    This class represents the equity, a kind of underlying.
 */

#ifndef _ITO33_FINANCE_EQUITY_H_
#define _ITO33_FINANCE_EQUITY_H_

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

  class ITO33_DLLDECL Issuer;
  class ITO33_DLLDECL YieldCurve;
  class ITO33_DLLDECL Dividends;
  class ITO33_DLLDECL Numeraire;
  

/**
    Class defining an equity.  

    It has the spot as a continuous state variable, as well as the issuer, 
    a currency, a borrow (foreign, continuous dividend) curve and dividends.
 */
class ITO33_DLLDECL Equity
{

public:
  /**
      Constructs an Equity object.

      @param dSpotSharePrice spot share price
      @param pCurrency a valid currency to use
   */
  Equity(double dSpotSharePrice, 
         const shared_ptr<Numeraire>& pCurrency);

  /**
      Constructs an Equity object with its currency.

      @param pCurrency a valid currency to use
   */
  Equity(const shared_ptr<Numeraire>& pCurrency);

  /**
      @name Modifiers.
   */
  //@{

  /**
      The issuer of this equity.

      @param pIssuer the issuer of this equity
   */
  void SetIssuer(const shared_ptr<Issuer>& pIssuer);

  /**
      The spot price of this equity.

      This method may be called multiple times, for example to periodically
      update the spot with market data.

      @param dSpotSharePrice the equity price in the market
   */
  void SetSpotSharePrice(double dSpotSharePrice);

  /**
      The borrow (foreign, continuous dividend) curve of this equity.

      The borrow curve can be entered by the final user.       

      @param pBorrowCurve a shared pointer to the borrow curve of the equity
   */
  void SetBorrowCurve(const shared_ptr<YieldCurve>& pBorrowCurve);

  /**
      The discrete dividends of this equity.

      @param pDividends a shared pointer of the Dividends objet
   */
  void SetDividends(const shared_ptr<Dividends>& pDividends); 
    
  /**
      Previous share price.

      This value is only relevant when the valuation date
      is greater or equal to the start of the sampling period.
      It is only used for pricing some path dep instruments 
      such as variance swap for example.
    
      @param dPreviousSharePrice previous share price
   */
  void SetPreviousSharePrice(double dPreviousSharePrice);

  //@}

  /**
      @name Accessors.
   */
  //@{

  /**
      The borrow (continuous dividend) curve of this equity.

      @return borrow curve of the equity
   */
  const shared_ptr<YieldCurve>& GetBorrowCurve() const 
  { 
    return m_pBorrowCurve;
  }

  /**
      The discrete dividends of this equity.

      @return shared pointer to discrete dividends of this equity
   */
  const shared_ptr<Dividends>& GetDividends()  const
  { 
    return m_pDividends; 
  }

  /**
      The issuer of this equity.

      @return the issuer of this equity
   */
  const shared_ptr<Issuer>& GetIssuer() const 
  { 
    return m_pIssuer; 
  }

  /**
      The currency of this equity, non NULL as required by the ctor.

      @return the currency of this equity
   */
  const shared_ptr<Numeraire>& GetNumeraire() const
  {
    return m_pNumeraire;
  }

  /**
      The spot price of this equity. The function throws if the spot is not set.

      @return spot share price.
   */
  double GetSpotSharePrice() const;
   
  /**
      Previous share price.

      @return the previous share price
   */
  double GetPreviousSharePrice() const
  {
    return m_dPreviousSharePrice;
  }

  //@}


  /**
      @internal
      @brief Dumps all data stored in this object in XML format.

      @param tagParent the parent tag under which our tag(s) should be created
      
      @noexport
   */
  void Dump(XML::Tag& tagParent) const;


private:
  
  /// the current market price
  double m_dSpot;

  /// the currency
  shared_ptr<Numeraire> m_pNumeraire;

  /// the issuer
  shared_ptr<Issuer> m_pIssuer;
  
  /// borrow curve (or foreign curve, or continous dividends)
  shared_ptr<YieldCurve> m_pBorrowCurve;

  /// the associated dividends, NULL if there are no associated dividends
  shared_ptr<Dividends> m_pDividends;

  /// the previous share price 
  double m_dPreviousSharePrice;

}; // class Equity


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_EQUITY_H_
