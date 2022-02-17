///////////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/cashflowstream_general.h
// Purpose:     Class for the description of cash flows streams defined
//              by a series of dates and amounts
// Author:      Pedro Ferreira
// Created:     Mar 25, 2004
// RCS-ID:      $Id: cashflowstream_general.h,v 1.15 2005/03/31 10:24:34 pedro Exp $
// Copyright:   (c) 2004 Trilemma LLP
///////////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/cashflowstream_general.h
    @brief  Class for the description of cash flows streams defined
            by a series of dates and amounts.
 */

#ifndef _ITO33_FINANCE_CASHFLOWSTREAM_GENERAL_H_
#define _ITO33_FINANCE_CASHFLOWSTREAM_GENERAL_H_

#include "ito33/vector.h"
#include "ito33/date.h"

#include "ito33/finance/cashflowstream.h"

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
    This class models general cash flows streams with payments which are not 
    necessarily uniform.

    The stream is defined by a series of dates and values
    (expressed as percentage rates as well as absolute values). 
 */
class ITO33_DLLDECL CashFlowStreamGeneral : public CashFlowStream
{
public:
  /** 
     Ctor constructs a general cash flow stream. 

     @param contractingDate The date on which the cash amount paid starts 
            to accrue
     @param paymentDates The cash payment dates
     @param paymentAmounts The cash payment amounts. they can be expressed
              as percentage rates of nominal or as absolute values.
     @param dcc The day count convention for accrued calculation
     @param freq The payment frequency
   */
  CashFlowStreamGeneral( Date contractingDate,
                         const std::vector<Date>& paymentDates,
                         const std::vector<double>& paymentAmounts,
                         Date::DayCountConvention dcc,                         
                         Frequency freq );

  XML::Tag Dump(const char *name, XML::Tag& tagParent) const;

}; // class CashFlowStreamGeneral


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_CASHFLOWSTREAM_GENERAL_H_
