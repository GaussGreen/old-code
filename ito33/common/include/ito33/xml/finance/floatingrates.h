/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/floatingrates.h
// Purpose:     Names of elements and attributes used in XML for 
//              FloatingRates class
// Author:      Nabil
// Created:     2005/07/15
// RCS-ID:      $Id: floatingrates.h,v 1.9 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/floatingrates.h
    @brief Contains the names of the elements used in the XML description of a
           floating rates stream and also reader functions.
 */

#ifndef _ITO33_XML_FINANCE_FLOATINGRATES_H_
#define _ITO33_XML_FINANCE_FLOATINGRATES_H_

#include "ito33/common.h"

namespace xml
{
  class node;
}

//=============================================================================
//=              readers                                                      =
//=============================================================================
namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL  FloatingRates;
}

namespace XML
{

/**
  Gets Floatingrates from given node. It throws exception
  when no floating rates is found

  @param node given node
  @return shared pointer to FloatingRates
  */
shared_ptr<finance::FloatingRates> GetFloatingRatesFromNode
      (
        const xml::node& node
      );

} // namespace XML

} // namespace ito33


/**
   @name Tag name macros
*/
//@{

#define XML_TAG_FLOATINGRATES_ROOT "floating_rates"

#define XML_TAG_FLOATINGRATES_DAYCOUNTCONVENTION "day_count_convention"
#define XML_TAG_FLOATINGRATES_MULTIPLIER "multiplier"
#define XML_TAG_FLOATINGRATES_MARGIN "margin"
#define XML_TAG_FLOATINGRATES_CAP "cap"
#define XML_TAG_FLOATINGRATES_FLOOR "floor"
#define XML_TAG_FLOATINGRATES_FIXINGDELAY "fixing_delay"
#define XML_TAG_FLOATINGRATES_KNOWNPAYMENTS "known_payments"
#define XML_TAG_FLOATINGRATES_KNOWNPAYMENT "known_payment"

#define XML_TAG_FLOATINGRATES_STARTOFACCRUEDDATE "start_of_accrued_date"

#define XML_TAG_FLOATINGRATES_FIRSTPAYMENTDATE \
  "first_unknown_payment_date"
#define XML_TAG_FLOATINGRATES_PRELASTDATE \
  "last_but_one_unknown_payment_date"
#define XML_TAG_FLOATINGRATES_LASTDATE \
  "last_unknown_payment_date"
//@}

#endif // _ITO33_XML_FINANCE_FLOATINGRATES_H_

