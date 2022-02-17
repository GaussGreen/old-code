///////////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/cashflowstream.h
// Purpose:     Class for description of the cash flows streams used in several 
//              instruments
// Author:      Pedro Ferreira
// Created:     Mar 25, 2004
// RCS-ID:      $Id: cashflowstream.h,v 1.50 2006/08/16 11:04:05 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
///////////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/cashflowstream.h
    @brief Class for description of the cash flows used in several instruments
 */

#ifndef _ITO33_FINANCE_CASHFLOWSTREAM_H_
#define _ITO33_FINANCE_CASHFLOWSTREAM_H_

#include "ito33/beforestd.h"
#include <map>
#include <vector>
#include "ito33/afterstd.h"

#include "ito33/debug.h"
#include "ito33/date.h"

#include "ito33/finance/cashflow.h"
#include "ito33/finance/frequency.h"

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

/**
    This class models cash flow streams used by derivatives like credit 
    default swaps or convertible bonds.

    This is the base class that acts like a container for the cash flows. 
    The derived classes can implement simpler (global) descriptions of the 
    stream of cash flows.
    
    @iterator GetAll
    @nocreate
 */
class ITO33_DLLDECL CashFlowStream
{
public:

  /// type of the container
  typedef std::map<Date, double> Elements;

  /// vector of cash flow
  typedef std::vector<CashFlow> CollectionElements;

  /// Embedded type for iteration on the container
  typedef Elements::const_iterator const_iterator;

  // Assignment operator is ok

  /**
     Ctor.

     @param contractingDate The date on which the cash amount paid starts 
            to accrue
     @param dcc The day count convention for accrued calculation
     @param freq The payment frequency 
   */
  CashFlowStream(Date contractingDate, Date::DayCountConvention dcc,
                 Frequency freq);
    
  /// Dummy virtual dtor for base class 
  virtual ~CashFlowStream() { }

  /** 
     Gets the contracting date of this cash flow stream.
   
     This date is usually not the same as the date for the first cash flow in 
     the stream. This is the date at which the accrued cash amounts starts to 
     accrue
   
     @return The contracting date of this cash flow.
   */
  Date GetContractingDate() const;

  /** 
     Gets the date of the last payment of this cash flow stream.
   
     @return The date for the last payment.
   */
  Date GetLastPaymentDate() const;

  /** 
     Gets the frequency of cash payments.
   
     @return The payment frequency
   */ 
  Frequency GetPaymentFrequency() const;

  /**
     Gets the day count convention for the accrued calculation.
   
     @return The day count method
   */
  Date::DayCountConvention GetDayCountConvention() const;

  /** 
     Gets the accrued interest at a given date.

     @param date The date at which the accrued cash amount is calculated
     @return the accrued interest at the given date
   */
  virtual double GetAccrued(Date date) const;

  /**
     Gets the discount of all the coupons starting from the given date.

     @param yc The yield curve used for discounting
     @param date The start date

     @noexport
   */
  double GetDiscount(const YieldCurve& yc, Date date) const;

  /**
     @internal

     @brief Gets the iterator points to the first cash flow.

     @noexport
   */
  const_iterator begin() const;
 
  /**
     @internal

     @brief Returns an iterator that addresses the location succeeding
            the last cash flow.

     @noexport
   */
  const_iterator end() const;

  /**
     @internal

     @brief Returns the number of cash flows.

     @noexport
   */
  size_t size() const;

  /**
     Gets all the cash flows.

     @noexport COM
   */
  const CollectionElements& GetAll() const;

  /**
     @internal

     @noexport
   */
  virtual XML::Tag Dump(const char *name, XML::Tag& tagParent) const = 0;


protected:
  
  /**
     Writes common data to tag root of specific CashFlowStream.

     @param tagRoot CallSchedule etc.
   */
  void DumpMe(XML::Tag& tagRoot) const;

  /**
     Adds one cash flow provision to this container.

     @param date The date of the to be added cash flow provision.
     @param amount The amount of the to be added cash flow provision.
   */
  inline void AddCashFlow(Date date, double amount); 

  /// Get all
  void DoGetAll() const;
  
  /// Start date for the payment
  Date m_contractingDate;
  
  /// The payment frequency
  Frequency m_freq;          

  /// The day count convention for the accrued calculation
  Date::DayCountConvention m_dcc;     

  /// The list of the elements
  Elements m_elements;

  /// vector of the cash flow
  CollectionElements m_collectionElements;

};

// Inline functions

inline Date CashFlowStream::GetContractingDate() const
{
  return m_contractingDate;
}

inline Date CashFlowStream::GetLastPaymentDate() const
{
  return m_elements.rbegin()->first;
}

inline Frequency CashFlowStream::GetPaymentFrequency() const
{
  return m_freq;
}

inline Date::DayCountConvention CashFlowStream::GetDayCountConvention() const
{
  return m_dcc;
}

inline CashFlowStream::const_iterator CashFlowStream::begin() const
{
  return m_elements.begin();
}

inline CashFlowStream::const_iterator CashFlowStream::end() const
{
  return m_elements.end();
}

inline size_t CashFlowStream::size() const
{
  return m_elements.size();
}

inline void CashFlowStream::AddCashFlow(Date dt, double amount) 
{
  ASSERT_MSG
  (
    m_elements.find(dt) == m_elements.end(),
    "Cash flow stream definition: duplicated payment dates"
  );

  /// protected functions, developper should check the following code himeself
  ASSERT_MSG
  (
    amount >= 0,
    "Cash flow definition: Setting negative amount."
  );

  m_elements[dt] = amount;
}

} // namespace finance

} // namespace ito33

#endif // _ITO33_FINANCE_CASHFLOWSTREAM_H_
