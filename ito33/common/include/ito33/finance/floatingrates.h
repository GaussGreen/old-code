///////////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/floatingrates.h
// Purpose:     Class for the description of floating rates defined
// Author:      Nabil
// Created:     2005/07/13
// RCS-ID:      $Id: floatingrates.h,v 1.13 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
///////////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/floatingrates.h
    @brief Class for the description of floating rates.
 */

#ifndef _ITO33_FINANCE_FLOATINGRATES_H_
#define _ITO33_FINANCE_FLOATINGRATES_H_

#include "ito33/vector.h"
#include "ito33/date.h"
#include "ito33/debug.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/frequency.h"
#include "ito33/finance/lastpaymenttype.h"
#include "ito33/finance/yieldcurve.h"

namespace ito33
{

namespace XML
{
  class Tag;
}

namespace finance
{

class ITO33_DLLDECL CashFlowStream;
  
/**
    This class models periodic floating rates (monthly, quaterly, etc).

    The stream is defined by a series of dates and a frequency.   
 */
class ITO33_DLLDECL FloatingRates
{

public:
   
  /**
     Constructs a floating rates uniform stream.

     @param dMargin the margin
     @param startOfAccruedDate the start of accrued date
     @param firstUnknownPaymentDate the date of the first unknown payment
     @param lastUnknownPaymentDate the date of the last unknown payment
     @param frequency the payment frequency between the dates of the first 
            and the before last unknown payments.
     @param lastPaymentType Indicates the type (short or long) of the last 
                            payment
   */
  FloatingRates(double dMargin, 
                Date startOfAccruedDate, 
                Date firstUnknownPaymentDate,
                Date lastUnknownPaymentDate,
                Frequency frequency,
                LastPaymentType lastPaymentType);
   
  /**
     Constructs a floating rates uniform stream.

     @param dMargin the margin
     @param startOfAccruedDate the start of accrued date
     @param firstUnknownPaymentDate the date of the first unknown payment
     @param lastUnknownPaymentDate the date of the last unknown payment
     @param frequency the payment frequency between the dates of the first 
            and the before last unknown payments.
   */
  FloatingRates(double dMargin, 
                Date startOfAccruedDate, 
                Date firstUnknownPaymentDate,
                Date lastUnknownPaymentDate,
                Frequency frequency);
  
  /**
     Constructs a floating rates uniform stream with no unknown payment.

     @param startOfAccruedDate the start of accrued date.
     @param frequency the payment frequency.
   */
  FloatingRates(Date startOfAccruedDate, Frequency frequency);
  
  /** 
     Gets the margin: This is the additional amount which will be paid 
     above the reference rate (LIBOR, EURIBOR).
   
     @return the margin.
   */
  double GetMargin() const;

  /**
     The multiplier coefficient: In deleveraged floaters, the payment rate is
     computed as a fraction of the reference rate plus the margin: 
     Payment rate = b x (reference rate) + margin
     b is called the multiplier and its value is between 0 and 1.

     @return the multiplier coefficient   
   */
  double GetMultiplier() const;
   
  /** 
     The multiplier coefficient.
   
     @param dMultiplier the multiplier coefficient
   */
  void SetMultiplier(double dMultiplier);

  /**
     The cap for the floating payment rates. It is the maximum payment that
     will be paid at any reset date, entered as an annual rate.

     @return the cap for the floating payment rates  
   */
  double GetCap() const;
   
  /** 
     The cap for the floating payment rates. It is the maximum payment that
     will be paid at any reset date, entered as an annual rate.
   
     @param dCap the cap for the floating payment rates 
   */
  void SetCap(double dCap);

  /**
     The floor for the floating payment rates. It is the minimum payment that
     will be paid at any reset date, entered as an annual rate.

     @return the floor for the floating payment rates  
   */
  double GetFloor() const;
   
  /** 
     The floor for the floating payment rates.
   
     @param dFloor the floor for the floating payment rates
   */
  void SetFloor(double dFloor);
  
  /**
     The fixing delay, the number of days between the fixing date of the 
     previous payment date and the date of the relevant LIBOR rate.

     @return the fixing delay   
   */
  int GetFixingDelay() const;
  
  /**
     The fixing delay, the number of days between the fixing date of the 
     previous payment date and the date of the relevant LIBOR rate.

     @param iFixingDelay the fixing delay   
   */
  void SetFixingDelay(int iFixingDelay);  
    
  /** 
     Gets the start of accrued date.
      
     @return the start of accrued date.
   */
  Date GetStartOfAccruedDate() const;
  
  /**
     The day count convention.
   
     @param dcc the day count convention
   */
  void SetDayCountConvention(Date::DayCountConvention dcc);

  /**
     The day count convention.
   
     @return the day count convention
   */
  Date::DayCountConvention GetDayCountConvention() const;

  /** 
     Gets the date of the first unknown payment (if it is available). 
   
     @return the date of the first unknown payment (if it is available).
   */
  Date GetFirstUnknownPaymentDate() const;

  /** 
     Gets the date of the last unknown payment (if it is available). 
   
     @return the date of the last unknown payment (if it is available).
   */
  Date GetLastUnknownPaymentDate() const;

  /**
     Sets the known payment rates.

     @param paymentDates the payment dates (in increasing order)
     @param paymentRates the payment rates.
   */
  void SetKnownPaymentStream(const std::vector<Date>& paymentDates,
                            const std::vector<double>& paymentRates);
  
  /**
     Gets the dates of the known payments.

     @return The dates of the known payments.
   */
  const std::vector<Date>& GetKnownPaymentDates() const;

  /**
     Gets the known payment rates.

     @return The known payment rates.
   */
  const std::vector<double>& GetKnownPaymentRates() const;
  
  /** 
     Gets the payment frequency in order to generate the payments.
   
     @return the payment frequency
   */ 
  Frequency GetPaymentFrequency() const;
  
  /**
     Computes all the floating payments (the known as the unknown).

     @param pYieldcurve the derivative yield curve ie the yield curve used to 
     discount the cash flow of the derivative intrument.
     
     @return all the floating payments (the known as the unknown)
   */
  shared_ptr<CashFlowStream> 
  ComputeFloatingPayments(const shared_ptr<YieldCurve>& pYieldcurve) const;

  /**
     Indicates if tyhe floating rates contains unknown payments. 

     @return true if the floating rates contains unknown payments, 
             false otherwise.
   */
  bool HasUnknownPayments() const { return m_bHasUnknownPayments; }
  
  /**
     @internal
     @brief Generates the dates of the unknown payments.

     @return The dates of the unknown payments.

     @noexport
   */
  std::vector<Date> GenerateUnknownPaymentDates() const;

  /**
     @internal
     @brief Sets the margin, the additional amount which will be paid 
            above the reference rate (LIBOR, EURIBOR).
   
     @param the margin.
     @noexport
   */
  void SetMargin(double dMargin) { m_dMargin = dMargin; }

  /// @noexport
  XML::Tag Dump(XML::Tag& tagParent) const;

protected:

  /// The margin
  double m_dMargin;

  /// The multiplier
  double m_dMultiplier;

  /// The cap for the floating payment rates
  double m_dCap;

  /// The floor for the floating payment rates
  double m_dFloor;

  /// Day count convention (used for the accrued)
  Date::DayCountConvention m_dcc;
  
  /// Start of accrued date
  Date m_startOfAccruedDate;   
  
  /// Date of the first unknown payment
  Date m_firstUnknownPaymentDate;     

  /// Date of the last unknown payment
  Date m_lastUnknownPaymentDate; 

  /// The fixing delay
  int m_iFixingDelay;

  /// boolean indicates if the floating rates contains unknown payments.
  bool m_bHasUnknownPayments;

  /// Frequency payment
  Frequency m_paymentfrequency;

  /// Indicates the type (short or long) of the last unknown payment
  LastPaymentType m_lastPaymentType;

  /// The dates of the Known payment rates
  std::vector<Date> m_pKnownPaymentRatesDates;

  /// The known payment rates
  std::vector<double> m_pKnownPaymentRates;

private:

  void Validate() const;

  /**
     Initialize some data with default values.

     Must be called by all the ctors.
   */
  inline void Init()
  {
    m_dMultiplier = 1.;
    m_dCap = 1.;
    m_dFloor = 0.;
    m_iFixingDelay = 0;

    //Default dcc for the floating rates
    m_dcc = Date::DayCountConvention_Act360;
  };

}; // class FloatingRates

// Inline functions

inline double FloatingRates::GetMargin() const
{
  return m_dMargin;
}

inline double FloatingRates::GetMultiplier() const
{
  return m_dMultiplier;
}

inline double FloatingRates::GetCap() const
{
  return m_dCap;
}

inline double FloatingRates::GetFloor() const
{
  return m_dFloor;
}

inline int FloatingRates::GetFixingDelay() const
{
  return m_iFixingDelay;
}

inline Date FloatingRates::GetStartOfAccruedDate() const
{
  return m_startOfAccruedDate;
}

inline Date::DayCountConvention FloatingRates::GetDayCountConvention() const
{
  return m_dcc;
}

inline Frequency FloatingRates::GetPaymentFrequency() const
{
  return m_paymentfrequency;
}

inline const std::vector<Date>& FloatingRates::GetKnownPaymentDates() const
{
  return m_pKnownPaymentRatesDates;
}

inline const std::vector<double>& FloatingRates::GetKnownPaymentRates() const
{
  return m_pKnownPaymentRates;
}

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_FLOATINGRATES_H_
