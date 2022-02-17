/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/bondterms.h
// Purpose:     Common terms of a bond
// Author:      ZHANG Yunzhi
// Created:     2004 may 3
// RCS-ID:      $Id: bondterms.h,v 1.59 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/bondterms.h
    @brief declaration of BondTerms class
 */

#ifndef _ITO33_FINANCE_BONDLIKE_BONDTERMS_H_
#define _ITO33_FINANCE_BONDLIKE_BONDTERMS_H_

#include "ito33/sharedptr.h"
#include "ito33/date.h"

#include "ito33/finance/bondlike/bondliketerms.h"
#include "ito33/finance/frequency.h"

namespace ito33
{

namespace XML
{
  class Tag;
}

namespace finance
{
 
class ITO33_DLLDECL FloatingRates;
class ITO33_DLLDECL SessionData;

/**
   Common characteristics of a bond instrument.
 */
class ITO33_DLLDECL BondTerms : public BondLikeTerms
{

public:
  
  /**
     ctor.

     @param issueDate the issue date
     @param dIssuePrice the isssue price, as a percentage of nominal
     @param dNominal the nominal
     @param maturityDate the maturity date
     @param dRedemptionPrice the redemption price, as a percentage of nominal
     @param dRecoveryRate the recovery rate of the bond
   */
  BondTerms(Date issueDate, 
            double dIssuePrice,
            Date maturityDate,
            double dNominal,                       
            double dRedemptionPrice,
            double dRecoveryRate);

  /// virtual dtor for base class
  virtual ~BondTerms() { }

  // see base class
  virtual void 
  SetCashDistribution(const shared_ptr<CashFlowStream>& pCashDistribution);

  /**
     Gets redemption price at maturity.

     @return redemption price
   */
  double GetRedemptionPrice() const;
  
  /**
     The yield compounding frequency.
  
     @param compoundingFrequency the yield compound frequency
   */
  void SetYieldCompoundingFrequency(Frequency compoundingFrequency);

  /**
     The yield day count convention.
  
     @param dcc the yield day count convention
   */
  void SetYieldDayCountConvention(Date::DayCountConvention dcc);

  /**
     The yield compounding frequency.

     @return the yield compound frequency
   */
  Frequency GetYieldCompoundingFrequency() const;

  /**
     The yield day count convention.

     @return the yield day count convention
   */
  Date::DayCountConvention GetYieldDayCountConvention() const;

  ///@name floating rates property
  //@{
  
  /**
     The floating rates.

     @param pFloatingRates the floating rates.
   */
  void SetFloatingRates(const shared_ptr<FloatingRates>& pFloatingRates);
  
  /**
     The floating rates.
    
     @return the floating rates if any, null otherwise.
   */
  const shared_ptr<FloatingRates>& GetFloatingRates() const 
  {
    return m_pFloatingRates;
  }

  //@} // floating rates property

  ///@name accretion and OID properties
  //@{

  
  /**
     Defines an accreting bond (including OID and premium redemption)
     by setting the gross yield to maturity. Note that the yield can be
     negative.

     @param dYieldToMaturity the gross yield to maturity 

     @method
   */
  void SetAccretingBond(double dYieldToMaturity);

  
  /**
     Whether it is an accreting bond.

     @return true if is consists of an accreting bond, otherwise false
   */
  bool IsAccretingBond() const
  {
    return m_bIsAccretingBond;
  }

  /** 
     The yield to maturity of accreting bonds.

     @return the accreting bond yield
   */
  double GetYieldToMaturityOfAccretingBond() const;


  /**
     Specifies that the bond accretes at the accretion rate
     (as for an OID zero bond) after the last coupon payment.

     @param dAccretionRate the rate at which the bond value
             accretes after cash pay.

     @method
   */
  void SetCashPayToZero(double dAccretionRate);

  /**
     Whether it is a cash pay to zero bond.

     @return true if is consists of a cash pay to zero bond, otherwise false.
   */
  bool IsCashPayToZeroBond() const
  {
    return m_bIsCashPayToZeroBond;
  }

  /**
     Gets the yield rate for cash pay to zero bonds.

     @return the yield rate of a cash pay to zero bond.
   */
  double GetAccretionRateOfCashPayToZero() const;

  //@} // name accretion and OID properties

  /**
     @internal

     @noexport
   */
  virtual XML::Tag Dump(ito33::XML::Tag& tagParent) const;

  
  /**
     Validates the terms: some data should be consistant between them.

     Throws exception if validation fails.
   */
  void Validate() const;

  /**
     Validates with a given session data.

     Throws exception if validation fails.

     @param sessionData The session data of the bond
   */
  void ValidateWith(const SessionData& sessionData) const;

protected:

  /// Redemption price
  double m_dRedemptionPrice;

  /// Floating rates
  shared_ptr<FloatingRates> m_pFloatingRates;

  /// The yield to maturity for OID bonds, can also be negative
  double m_dYield;

  /// Is it an accreting bonds?
  bool m_bIsAccretingBond;

  /// The yield for partial OID bonds, can also be negative
  double m_dAccretionRateOfCashPayToZero;

  /// Is it an partial OID bond
  bool m_bIsCashPayToZeroBond;

  /// Compounding payment frequency
  Frequency m_cmpFreq;

  /// Yield day count convention
  Date::DayCountConvention m_yieldDCC;

}; // class BondTerms


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_BONDTERMS_H_
