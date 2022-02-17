/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/parbond.h
// Purpose:     financial ParBond class
// Author:      zhang
// Created:     2003/05/20
// RCS-ID:      $Id: parbond.h,v 1.6 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/parbond.h
    @brief declaration of the financial ParBond class

    ParBond is an abstract financial class in the sense that no contigency is
    directly associated to it.
 */

#ifndef _ITO33_FINANCE_PARBOND_H_
#define _ITO33_FINANCE_PARBOND_H_

#include "ito33/useexception.h"
#include "ito33/date.h"

#include "ito33/finance/frequency.h"
#include "ito33/finance/derivative.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace finance
{

class ITO33_DLLDECL CashFlowStream;

/**
    ParBond represents the financial aspects of a par bond.
 */
class ITO33_DLLDECL ParBond : public Derivative
{
public:
  /**
     Creates a par bond object.

     ParBond object is immutable once it was created, i.e. its description
     can't be changed.

     @param contractingDate valuation date. same as the one in Session.
     @param nMaturity Maturity of the par bond expressed as a number of monthes
     @param dYTM Risk-free yield to maturity
     @param dSpread The spread of the par bond over risk-free yield to maturity
     @param paymentFrequency The payment frequency of the par bond
     @param dcc The day count convention
     @param dRecoveryRate The recovery rate of the par bond
   */
  ParBond(Date contractingDate,
          size_t nMaturity,
          double dYTM,
          double dSpread,
          Frequency paymentFrequency,
          Date::DayCountConvention dcc,
          double dRecoveryRate)
          : m_dYTM(dYTM),
            m_dSpread(dSpread),
            m_nMaturity(nMaturity),
            m_paymentFrequency(paymentFrequency),
            m_dcc(dcc),
            m_dRecoveryRate(dRecoveryRate),
            m_fakeContractingDate(contractingDate)
  {
    SelfValidate();

    m_dPrice = 1;
  }

  /**
     Empty virtual destructor.
   */
  virtual ~ParBond() { }

  /**
      @name Methods for accessing the ParBond.
   */
  //@{

  /**
     the spread over risk-free yield to maturity.

     @return the spread over risk-free yield to maturity.
   */
  double GetSpread() const
  {
    return m_dSpread;
  }

  /**
     risk-free yield to maturity

     @return risk-free yield to maturity
   */
  double GetRiskFreeYTM() const
  {
    return m_dYTM;
  }

  /**
     maturity expressed as a number of monthes.

     @return maturity
   */
  size_t GetMaturity() const
  {
    return m_nMaturity;
  }

  /**
     payment frequency.

     @return payment frequency
   */
  Frequency GetPaymentFrequency() const
  {
    return m_paymentFrequency;
  }

  /**
     day count convention

     @return day count convention
   */
  Date::DayCountConvention GetDayCountConvention() const
  {
    return m_dcc;
  }

  /**
     recovery rate.

     @return recovery rate
   */
  double GetRecoveryRate() const
  {
    return m_dRecoveryRate;
  }

  /**
     Contracting date
 
     @return contracting date
   */
  Date GetContractingDate() const;

  //@}

  /**
     Annual coupon rate.

     @return Annual coupon rate
   */
  double GetCouponRate() const;

  // implement base class pure virtuals
  virtual Date GetMaturityDate() const;

  virtual void Visit(DerivativeVisitor& visitor) const;
  virtual XML::Tag Dump(XML::Tag& tagParent) const;

  /**
     @internal

     @noexport
   */
  shared_ptr<CashFlowStream> GetCashFlowStream() const;


private:

    /// validate data in constructor
    void SelfValidate() const;

    /// session must be defined. maturity date depend on it.
    void CheckSession() const;

    /// risk-free yield-to-maturity
    double m_dYTM;

    /// spead over the ytm
    double m_dSpread;

    /// maturity expressed as a number of months
    size_t m_nMaturity;

    /// payment frequency
    Frequency m_paymentFrequency;
          
    /// day count convention
    Date::DayCountConvention m_dcc;

    /// recovery rate
    double m_dRecoveryRate;

    /// contracting date
    Date m_fakeContractingDate;

}; // class ParBond


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_PARBOND_H_

