/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/cblike.h
// Purpose:     contract class for CB-like
// Created:     2004/08/19
// RCS-ID:      $Id: cblike.h,v 1.43 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/cblike.h
    @brief The declaration of the CB-like contract class.

    The base class for CB-like contracts.  

    @todo call and conversion are both stored as objet, but we have to store 
          also a pointer to base class for them, which is annoying.
 */

#ifndef _ITO33_PRICING_CBLIKE_H_
#define _ITO33_PRICING_CBLIKE_H_

#include "ito33/sharedptr.h"

#include "ito33/pricing/callprovisions.h"
#include "ito33/pricing/cbputs.h"
#include "ito33/pricing/conversionprovisions.h"
#include "ito33/pricing/contract.h"


namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL SessionData;
  class ITO33_DLLDECL YieldCurve;
  class ITO33_DLLDECL BondLikeTerms;
  class ITO33_DLLDECL ConvertibleLike;
  class ITO33_DLLDECL Numeraire;
  
  class Payoff;
}

namespace pricing
{
  class CashFlows;


/**
    The declaration of the CB-like contract class.

    IMPORTANT: 
    call provision, put provision and conversion provision object are presented
    as objects instead of pointers in CBLike class and its sub-classes. This is
    because that we often need to copy a contract while Call(put/conversion) 
    provision object can't be shared. @sa CallProvisions.
    Of couse, this can be remedied by add clone() function to these provision 
    classes. But decision has been made that we don't write clone().
 */
class CBLike : public Contract 
{

public:

  /// Default ctor sets recovery rate to 0, and exchangeable property to false.
  CBLike() : m_dRecoveryRate(0.), m_bNewShare(false), m_bExchangeable(false),
             m_bIsExchangeableUponDefault(false),
             m_bIsFixedQuanto(false), m_dFXRateVolatility(0), m_dCorrelation(0),
             m_dSpotFXRate(1.), m_dCmpFreq(1.0), 
             m_dNominal(0.)
  { }

  /// virtual dtor for base class
  virtual ~CBLike() { }
  
  /// @name Some "sets" functions
  //@{
    
  /**
      Defines an exchangeable and sets its ExchangeableUponDefault property.

      @param bExchangeableUponDefault whether Holder can exchange the bond
                                      when default occurs.
   */
  void SetExchangeable(bool bExchangeableUponDefault)
  {
    m_bIsExchangeableUponDefault = bExchangeableUponDefault; 
    m_bExchangeable = true;
  }

  /**
      Sets an external payoff.
   */      
  void SetPayoff(const shared_ptr<finance::Payoff>& pPayoff) 
  {
    m_pPayoff = pPayoff;
  }

  /// Sets the spot FX rate, called for call notice problem
  void SetSpotFXRate(double dSpotFXRate) 
  { 
    m_dSpotFXRate = dSpotFXRate; 
  }
  
  /**
      Sets the cash flows.

      @param pCashFlows the cash flows.
   */
  void SetCashFlows(const shared_ptr<CashFlows>& pCashFlows)
  {
    m_pCashFlows = pCashFlows;
  }

  /**
      Sets the new share.

      @param bNewShare True if there are new share, false otherwise.
   */
  void SetNewShare(bool bNewShare) { m_bNewShare = bNewShare; }
  

  //@}

  /**
      Gets the recovery rate.

      @return the recovery rate (a double)
   */
  double GetRecoveryRate() const { return m_dRecoveryRate; }
  
  /**
      Gets the issue price.

      @return the issue price
   */
  double GetIssuePrice() const { return m_dIssuePrice; }

  /**
      Gets the issue time.

      @return the issue time
   */
  double GetIssueTime() const { return m_dIssueTime; }

  /**
      Gets the nominal.

      @return the nominal
   */
  double GetNominal() const { return m_dNominal; }

  /**
      Gets the FX rate at the pricing date.

      @return the FX rate at the pricing date
   */
  double GetSpotFXRate() const { return m_dSpotFXRate; }
  
  /**
      Gets the volatility of the FX Rate (for a fixed quanto).

      @return The volatility of the FX rate.
   */
  double GetFXRateVolatility() const;
  
  /**
      Gets the correlation between the underlying share and the FX rate 
      (for a fixed quanto).

      @return The correlation between the underlying share and the FX rate. 
   */
  double GetCorrelationBetweenUnderlyingAndFXRate() const;

  /**
      Gets the fiscal year start date.
   
      @return the fiscal year start date.
   */
  Date GetFiscalYearStartDate() const { return m_fiscalYearStartDate; }
  
  /**
      Checks if the instrument is exchangeable in case of default, if yes, gets
      also its ExchangeableUponDefault property.

      @param bExchangeableUponDefault (output) whether Holder can exchange
                                      the bond  when default occurs.
      @return whether it is a exchangeable bond
   */
  bool IsExchangeable(bool& bExchangeableUponDefault) const
  {
    bExchangeableUponDefault = m_bIsExchangeableUponDefault; 
    return m_bExchangeable;
  }

  /**
      Checks if the CB-like is excheangeable.

      @return true it is exchangeable, false otherwise
   */
  bool IsExchangeable() const { return m_bExchangeable; }   

  bool HasNewShare() const { return m_bNewShare; }

  bool IsFixedQuanto() const { return m_bIsFixedQuanto; }

  /**
      Gets the call.

      @return a pointer to the call.
   */
  virtual CallProvisions* GetCalls() = 0;

  /**
      Gets the conversion provision.

      @return a pointer to the conversion provision.
   */
  virtual ConversionProvisions* GetConversions() = 0;

  /**
      Gets the puts.

      @return a pointer to the puts.
   */
  CBPuts* GetPuts() { return &m_puts; }

  /**
      Gets the CashFlows.

      @return a pointer to the CashFlows.
   */
  CashFlows* GetCashFlows() const { return m_pCashFlows.get(); }

  /**
      Gets the payoff.
   */      
  finance::Payoff* GetPayoff() const { return m_pPayoff.get(); }

  /**
      Gets the compounding frequency.
   */
  double GetCompoundingFrequency() 
  {
    return m_dCmpFreq;
  }
  
  /**
      Gets the claim value at the given time dTime.

      @param dTime The time to which we seek to compute the claim. 
      @param bPlus at dTime+ or dTime-
      @return the claim value at the time dTime    
   */
  virtual double GetClaim(double dTime, bool bPlus = true) const = 0;
 
protected:

  /**
      Gets the information from a financial bond-like object. It is normally 
      called by the constructor using financial instrument.

      @param blt financial BondLikeTerms object
      @param sessionData session data of the underlying derivative
      @param pNumeraire the derivative currency
   */
  void GetBondLikeTermsData(const finance::BondLikeTerms& blt, 
                            const finance::SessionData& sessionData,
                            const shared_ptr<finance::Numeraire>& pNumeraire);
  
  /**
      Gets the cross currency information.

      Assumes that currency in contracts has been set.

      @param sessionData session data of the underlying derivative
      @param pNumeraire the derivative currency
   */
  void GetCrossCurrencyData(const finance::SessionData& sessionData,
                            const shared_ptr<finance::Numeraire>& pNumeraire);


  /** 
      Gets the information from a financial convertible like object, conversion
      must have been set, and bondliketerms data already get.

      @param cbLike convertible-like object
   */
  void GetConvertibleLikeData(const finance::ConvertibleLike& cbLike);

  /// recovery rate
  double m_dRecoveryRate;

  /// issue price as percentage rate
  double m_dIssuePrice;

  /// issue time
  double m_dIssueTime;

  /// nominal
  double m_dNominal;

  /// put provision. 
  CBPuts m_puts;
  
  /// coupon cash flow
  shared_ptr<CashFlows> m_pCashFlows;

  /// payoff of the contract
  shared_ptr<finance::Payoff> m_pPayoff; 

  /// new share
  bool m_bNewShare;

  /// FX rate at valuation date for a cross currency.
  double m_dSpotFXRate;  

  /// if the CB-like is exchangeable
  bool m_bExchangeable;  
  
  /// Exchangeable upon default flag for exchangeable CB-like
  bool m_bIsExchangeableUponDefault;

  /// compounding frequency
  double m_dCmpFreq;

  /// The fiscal year start date
  Date m_fiscalYearStartDate;  

  /// Fixed quanto flag.
  bool m_bIsFixedQuanto;

  /// The FX rate volatility
  double m_dFXRateVolatility;

  /** 
      The correlation between the underlying share and the FX rate 
      (for fixed quanto).
   */
  double m_dCorrelation;

}; // class CBLike;


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_CBLIKE_H_
