/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/putperiod.h
// Purpose:     put period
// Author:      ZHANG Yunzhi
// Created:     2005-04-20
// RCS-ID:      $Id: putperiod.h,v 1.5 2006/05/03 10:12:32 nabil Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/finance/bondlike/putperiod.h
   @brief declaration of the PutPeriod class for CB-like    
 */

#ifndef _ITO33_FINANCE_BONDLIKE_PUTPERIOD_H_
#define _ITO33_FINANCE_BONDLIKE_PUTPERIOD_H_

#include "ito33/date.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace finance
{

/**
   Defines a put period.

   @nocreate
 */
class ITO33_DLLDECL PutPeriod
{
public:
  
  /**
     Default constructor required by std::list.

     @noexport
   */
  PutPeriod() {}

  /**
     Creates a PutPeriod.

     @noexport
   */
  PutPeriod(Date putDate, double dPutValue, bool bHasYield)
  {
    m_date = putDate;
    m_bHasYield = bHasYield;

    if(m_bHasYield)
      m_dYield = dPutValue;
    else
      m_dStrike = dPutValue;
  }

  // copy constructor is ok

  /**
     Gets the date of the put

     @return The put date
   */
  Date GetDate() const { return m_date; }

  /**
     Whether put payment is defined by a guaranteed yield

     @return true if a yield is guaranteed upon put, false if a
             fraction of principal is paied
   */
  bool HasYield() const
  {
    return m_bHasYield;
  }

  /**
     Gets the put strike expressed as a percentage of the principal.

     @return The put strike expressed as a percentage of principal.
   */
  double GetStrike() const;

  /**
     Gets the guaranteed yield upon put.

     @return the guaranteed yield upon put
   */
  double GetGuaranteedYield() const;

private:

  /// put date
  Date m_date;

  /// Strike of the put, in percentage
  double m_dStrike;

  /// Guaranteed yield to put
  double m_dYield;

  /// Has yield
  bool m_bHasYield;

}; // class PutPeriod


} // namespace finance

} // namespace ito33

#endif // _ITO33_FINANCE_BONDLIKE_PUTPERIOD_H_
