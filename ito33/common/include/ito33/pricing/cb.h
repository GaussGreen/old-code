/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/cb.h
// Purpose:     contract class for CB (backward)
// Author:      Nabil
// Created:     2004/03/15
// RCS-ID:      $Id: cb.h,v 1.55 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/cb.h
    @brief The declaration of the CB contract class.

    The base class for CB contracts.  
 */

#ifndef _ITO33_PRICING_CB_H_
#define _ITO33_PRICING_CB_H_

#include "ito33/autoptr.h"
#include "ito33/vector.h"
#include "ito33/sharedptr.h"

#include "ito33/numeric/predicatetime.h"

#include "ito33/pricing/cashflows.h"
#include "ito33/pricing/cbcalls.h"
#include "ito33/pricing/cbputs.h"
#include "ito33/pricing/cbconversions.h"

#include "ito33/pricing/cblike.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL SessionData;
  class ITO33_DLLDECL ConvertibleBond;
  class ITO33_DLLDECL CBBase;
  class ITO33_DLLDECL Bond;
  class ITO33_DLLDECL BondTerms;
  class ITO33_DLLDECL CallSchedule;
  class ITO33_DLLDECL PutSchedule;
  class ITO33_DLLDECL FloatingRates;
  class ITO33_DLLDECL Numeraire;
}

namespace pricing
{

/// The declaration of the (backward) cb contract class.
class CB : public CBLike 
{

public:

  /**
     Creates a CB by financial ConvertibleBond object.

     @param cb financial ConvertibleBond
   */
  CB(const finance::ConvertibleBond& cb);

  /**
     Creates a CB by financial Bond object.

     @param bond financial straight Bond
   */
  CB(const finance::Bond& bond);

  // default ctor
  CB() : CBLike(), 
         m_dAccretionStartTime(-100.),
         m_dYield(-1),
         m_bHasYield(false),
         m_bRedeemedOnce(true) { }

  /// virtual dtor for base class
  virtual ~CB() { }
  
  /// @name Some "sets" functions
  //@{
   
  /**
     Sets the value of the redemption
   */
  void SetRedemptionValue(double dRedemptionValue)
  {
    m_dRedemptionValue = dRedemptionValue; 
  }
  
  /**
     Sets the calls.
   */
  void SetCalls(const CBCalls& calls) 
  { 
    m_calls = calls;  
  }

  /**
     Sets the conversions.
   */
 void SetConversions(const CBConversions& conversions) 
  { 
    m_conversions = conversions;
  }

  //@}

  CallProvisions* GetCalls() { return &m_calls; }

  CBCalls* GetCBCalls() { return &m_calls; }

  virtual ConversionProvisions* GetConversions() 
  { 
    return &m_conversions; 
  }

  /**
     Gets the redemption value

     @return the redemption value
   */
  double GetRedemptionValue() const { return m_dRedemptionValue; }

  /**
     Gets the yield for OID bonds (if any).

     @return the accretion rate
   */
  double GetYield() const { return m_dYield; }
 
  /**
     Gets the accretion start time for OID bonds
 
     @return The accretion start date
   */
  double GetAccretionStartTime() const  { return m_dAccretionStartTime; }
  
  /**
     Checks if the CB has accretion rate at given time.

     @param dTime given time
     @param bPlus dTime+ or dTime-?
     @return Whether the bond has has accretion rate at given time   
   */
  bool IsAccreting(double dTime, bool bPlus = true) const
  {
    if ( !m_bHasYield )
      return false;
    
    if (bPlus)
    {
      if( numeric::IsBefore(dTime, m_dAccretionStartTime) )
        return false;
      else
        return true;
    }
    else
    {
      if( numeric::IsAfter(dTime, m_dAccretionStartTime) )
        return true;
      else
        return false;
    }
  }

  double GetClaim(double dTime, bool bPlus = true) const;


protected:

  /**
     Gets the information from a financial bond object. It is normally called
     by the constructor using financial instrument.

     @param sessionData session data
     @param bondTerms bond terms data
     @param pCalls call data
     @param pPuts put data
     @param pNumeraire the derivative (bond) currency
   */
  void GetBondData(const finance::SessionData& sessionData, 
                   const finance::BondTerms& bondTerms,
                   const shared_ptr<finance::CallSchedule>& pCalls,
                   const shared_ptr<finance::PutSchedule>& pPuts,
                   const shared_ptr<finance::Numeraire>& pNumeraire);

  /**
     Gets the information from a financial CBBase object. It is 
     normally called by constructors accepting derived classes of
     CBBase objects.

     @param cb financial convertible-like object
   */
  void GetCBBaseData(const finance::CBBase& cb);

  /// call provision. Cannot be shared due to cbparams ptr
  CBCalls m_calls;
  
  /// conversion provision. Cannot be shared due to cbparams ptr
  CBConversions m_conversions;
 
  /// redemption value
  double m_dRedemptionValue;

  /// if it is accreting or partial accreting bond
  bool m_bHasYield;

  /// accretion rate for OID, or partial OID, bonds, if any
  double m_dYield;

  /// Start time for accretion rate
  double m_dAccretionStartTime;

  /// floating rates
  shared_ptr<finance::FloatingRates> m_pFloatingRates;

  /// whether the capital is redeemed once
  bool m_bRedeemedOnce;

  //--- for the bond which is redeemed several times
  /// repayment dates for thoses bonds which are redeemed more than once
  std::vector<double> m_pdRepaymentTimes;

  /// repayment rates for thoses bonds which are redeemed more than once
  std::vector<double> m_pdRepaymentRates;

}; // class CB;


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_CB_H_

