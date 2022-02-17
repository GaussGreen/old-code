///////////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/cashflowstream_uniform.h
// Purpose:     Class for the description of cash flows streams defined
//              by a payment frequency
// Author:      Pedro Ferreira
// Created:     Mar 25, 2004
// RCS-ID:      $Id: cashflowstream_uniform.h,v 1.39 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
///////////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/cashflowstream_uniform.h
    @brief  Class for the description of cash flows streams defined
            by a payment frequency
 */

#ifndef _ITO33_FINANCE_CASHFLOWSTREAM_UNIFORM_H_
#define _ITO33_FINANCE_CASHFLOWSTREAM_UNIFORM_H_

#include "ito33/date.h"
#include "ito33/debug.h"

#include "ito33/finance/cashflowstream.h"
#include "ito33/finance/lastpaymenttype.h"

#include "ito33/dlldecl.h"


namespace ito33
{

namespace XML
{
  class Tag;
}

namespace finance
{


/**
    This class models cash flows streams with payments occuring at constant
    intervals (monthly, quaterly, etc).
 */
class ITO33_DLLDECL CashFlowStreamUniform : public CashFlowStream
{
public:

  /**
     Ctor constructs a uniform cash flow stream.

     @param contractingDate The first settlement date of the cash flow stream
     @param firstPaymentDate The first payment date
     @param lastPaymentDate The last payment date
     @param annualAmt Annual payment amount. Expressed as a percentage of the
            nominal as well as absolute value
     @param dcc The day count method used for the calculation of long/short
            payments and for accrued cash payments.
     @param freq The payment frequency between the firstPaymentDate and
                 lastPaymentDate
   */
  CashFlowStreamUniform( Date contractingDate,
                         Date firstPaymentDate,
                         Date lastPaymentDate,
                         double annualAmt,
                         Date::DayCountConvention dcc,
                         Frequency freq )
    : CashFlowStream(contractingDate, dcc, freq),
      m_lastPaymentDate(lastPaymentDate),
      m_firstPaymentDate(firstPaymentDate),
      m_dAnnualPr(annualAmt),
      m_bAdjusted(false),
      m_lastPaymentType(LastPaymentType_Max),
      m_nMonths(0)
  {
    Validate();

    Init();
  }

  /**
     Ctor constructs a uniform cash flow stream which is regular starting from
     contracting date.

     @param contractingDate The first settlement date of the cash flow stream
     @param nMonths The number of months between the contracting date and the
                     last payment date.
     @param freq The payment frequency between the firstPaymentDate and
                 lastPaymentDate
     @param annualAmt Annual payment amount. Expressed as a percentage of the
            nominal as well as absolute value
     @param dcc The day count method used for the calculation of long/short
            payments and for accrued cash payments.
   */
  CashFlowStreamUniform( Date contractingDate,
                         size_t nMonths,
                         Frequency freq,
                         double annualAmt,
                         Date::DayCountConvention dcc);

  /**
     Ctor constructs a uniform cash flow stream where the user may adjust the
     payment dates and amounts to account for week-ends and holidays.
     The user has to ensure the consistency between the rates and the annual 
     amount specified if the sensitivity of the theoretical value of the 
     derivative with respect to a small change in the annual amount is required 
     (for example the sensitivity of the CDS theoretical value to a shift in 
      the CDS spread).

     @param contractingDate The issue date of the cash flow stream
     @param paymentDates The vector of payment dates
     @param paymentAmounts The vector of payment amounts
     @param annualAmt Annual payment amount. Expressed as a percentage of the
            nominal as well as an absolute value.
     @param dcc The day count method used for the calculation of long/short
            payments and for accrued cash payments.
     @param freq The payment frequency.
   */
  CashFlowStreamUniform( Date contractingDate,
                         const std::vector<Date>& paymentDates,
                         const std::vector<double>& paymentAmounts,
                         double annualAmt,
                         Date::DayCountConvention dcc,
                         Frequency freq );

  /**
     Ctor constructs a uniform cash flow stream where the user should
     specify the kind of the last payment (short or long). 
     
     If the last payment is uniform, lastPaymentType will be ignored.

     @param contractingDate The first settlement date of the cash flow stream
     @param firstPaymentDate The first payment date
     @param lastPaymentDate The last payment date
     @param annualAmt Annual payment amount. Expressed as a percentage of the
            nominal as well as absolute value
     @param dcc The day count method used for the calculation of long/short
            payments and for accrued cash payments.
     @param freq The payment frequency between the firstPaymentDate and
                 lastButOnePaymentDate
     @param lastPaymentType Indicates the type (short or long) of the last 
                            payment
   */
  CashFlowStreamUniform( Date contractingDate,
                         Date firstPaymentDate,
                         Date lastPaymentDate,
                         double annualAmt,
                         Date::DayCountConvention dcc,
                         Frequency freq, 
                         LastPaymentType lastPaymentType );

  /**
     Gets the date of the first payment in this stream.

     @return The first payment date
   */
  Date GetFirstPaymentDate() const;

  /**
     Gets the annual cash payment amount.

     @return The annual payment amount.
   */
  double GetAnnualPaymentAmount() const;

  /**
     Checks if the cash flow stream is adjusted.

     @return true if the cash flow stream is adjusted
   */
  bool IsAdjusted() const { return m_bAdjusted; }

  /**
      Creates a new uniform spread stream identical to this one except with a
      different annual spread.

      This function is useful when we need to dynamically adjust the spread
      value of a cashflow stream for an existing object.

      Notice that it doesn't modify this object but creates a new one.

      @noexport
   */
  shared_ptr<CashFlowStreamUniform> ChangeSpread(double annualAmt) const;

  // see base class
  virtual double GetAccrued(Date date) const;
  
  XML::Tag Dump(const char *name, XML::Tag& tagParent) const;


private:

  void Validate() const;

  // generate all the coupons between the first and last dates, called from all
  // ctors
  void Init();

  /**
     Checks if a given period is regular.

     @param startDate the start date of the period.
     @param endDate the end date of the period. 

     @return True if the period is regular, false otherwise.
   */
  bool IsRegularPeriod(Date startDate, Date endDate) const;

  /**
     Generates the coupon amount payed at the date \a currentDate which 
     should be either the first or the last payment date.

     \remark This function should be used only to generate the first 
              and last payment.

     @param previousDate the start date of the coupon period. 
            (i.e: last but one coupon date or contracting date)
     @param currentDate the date of the payment: It should be either 
            the first or the last payment date. 

     @return The generated coupon amount.
   */
  double GenerateCouponAmount(Date previousDate, Date currentDate) const;

  /// Last payment date
  Date m_lastPaymentDate;
  Date m_firstPaymentDate;           // First payment date
  double m_dAnnualPr;                // The annual payment amount

  /// boolean indicates if the cash flow stream is user adjusted
  bool m_bAdjusted;
  
  /// Indicates the type (short or long) of the last payment
  LastPaymentType m_lastPaymentType;

  /**
     Duration in unit of month if the cash flow stream is given by contracting
     date and duration.
   */
  size_t m_nMonths;

}; // class CashFlowStreamUniform

// Inline functions

inline Date CashFlowStreamUniform::GetFirstPaymentDate() const
{
  return m_firstPaymentDate;
}

inline double CashFlowStreamUniform::GetAnnualPaymentAmount() const
{
  return m_dAnnualPr;
}

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_CASHFLOWSTREAM_UNIFORM_H_
