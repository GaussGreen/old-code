/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/referencecds.h
// Purpose:     financial reference cds class 
// Created:     2006/05/17
// RCS-ID:      $Id: referencecds.h,v 1.9 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/referencecds.h
    @brief declaration of the financial reference cds class.
 */

#ifndef _ITO33_FINANCE_REFERENCECDS_H_
#define _ITO33_FINANCE_REFERENCECDS_H_

#include "ito33/common.h"
#include "ito33/date.h"

#include "ito33/finance/frequency.h"
#include "ito33/finance/cdslike.h"


namespace ito33
{

namespace finance
{

  class DerivateVisitor;
  class DerivativeModifyingVisitor;

/**
   This class describes the Credit Default Swap (CDS) contract related to 
   a reference CDS.

   A reference CDS is primarily used for calibration.  It has a price
   of zero by default, and assumes that the issue date is equal
   to the valuation date. It is also called a sliding CDS.
 */
class ITO33_DLLDECL ReferenceCDS : public CDSLike
{

public:

  /**
     Ctor creates a reference CDS (used mainly for calibration).

     The market price defaults to 0.  The function SetSpread() must
     be called before an object of this class can be used.
 
     @param maturity The approximated maturity (number of months) of the CDS
     @param freq The payment frequency
     @param dcc The day count convention
     @param dRecoveryRate The recovery rate
   */
  ReferenceCDS(size_t maturity,
               Frequency freq, Date::DayCountConvention dcc,
               double dRecoveryRate);  

  // Default dtor is ok. ReferenceCDS won't be derived

  /**
     The spread.

     Once the spread is defined, the underlying spread stream can be 
     created, thereby allowing the contract to be priced or calibrated.

     @param dSpread the spread (payment amount)
   */
  void SetSpread(double dSpread);
  
  /**
     Checks if the spread has been set.

     @return true if the spread has been set, false otherwise
   */
  bool HasSpread() const;

  /**
     Gets the maturity of the reference CDS as number of months.

     @return the maturity of the reference CDS as number of months
   */
  size_t GetMonthsToMaturity() const
  {
    return m_nMonthsToMaturity;
  }

  /**
      Gets the day count convention of the cash flow payment.

      @return the day count convention of the cash flow payment.
   */
  Date::DayCountConvention GetDayCountConvention() const
  {
    return m_dcc;
  }

  /**
      Gets the spread payment frequency.

      @return the spread payment frequency.
   */
  Frequency GetPaymentFrequency() const
  {
    return m_freq;
  }

  /**
     Gets the maturity date of the CDS contract which is determined by
     the valuation date of the session data and the maturity of the 
     reference CDS (as number of months). Throws if the session data is
     not set.

     @return the maturity date 
   */
  virtual Date GetMaturityDate() const;

  /**
     Gets the first payment date of the CDS contract which is determined by
     the valuation date of the session data. Throws if the session data is
     not set.

     @return the first payment date
   */
  virtual Date GetFirstPaymentDate() const;

  /**
     The spread.

     Throws if the spread is not defined.

     @return the spread
   */
  virtual double GetSpread() const;

  /**
     Gets the stream of spread payments.

     Throws if the spread or the session data has not been defined.

     @return the stream of spread payments
   */
  virtual shared_ptr<CashFlowStreamUniform> GetSpreadStream() const;

  /**
     Sets the market price for the reference CDS.

     Throws if the market price is not zero.

     @param dPrice market price
   */
  void SetMarketPrice(double dPrice);

  /**
      Sets the session data of the reference CDS. This makes
      maturity date of the reference CDS available.

      @pSessionData session data of the reference CDS.
   */
  virtual void SetSessionData(const shared_ptr<SessionData>& pSessionData);

  // implement base class pure virtuals
  void Visit(DerivativeVisitor& visitor) const;
  void Visit(DerivativeModifyingVisitor& visitor);
  XML::Tag Dump(XML::Tag& tagParent) const;

private:
  /// makes the spread Stream
  void MakeSpreadStream();

private:
  /// The approximated maturity (number of months) of the CDS
  size_t m_nMonthsToMaturity;

  /// The CDS spread
  double m_dSpread;
  
  /// The payment frequency
  Frequency m_freq;

  /// The day count convention
  Date::DayCountConvention m_dcc;
}; // class ReferenceCDS


} //namespace finance

} //namespace ito33

#endif // #ifndef _ITO33_FINANCE_REFERENCECDS_H_

