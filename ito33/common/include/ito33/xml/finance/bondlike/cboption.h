/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/bondlike/cboption.h
// Purpose:     Names of elements and attributes used in XML for a cb option
// Author:      Nabil
// Created:     2005/06/23
// RCS-ID:      $Id: cboption.h,v 1.6 2006/04/04 16:29:47 wang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/bondlike/cboption.h
    @brief Contains the names of the elements used in the XML description of a
           cb option

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_BONDLIKE_CBOPTION_H_
#define _ITO33_XML_FINANCE_BONDLIKE_CBOPTION_H_

#include "ito33/enum_values_names.h"

/**
   @name Tag name macros
*/
//@{
#define XML_TAG_CBOPTION_ROOT "cb_option"

#define XML_TAG_MATURITYDATE "maturity_date"

#define XML_TAG_ASWNOTIONALTYPE "aswnotional_type"

#define XML_VALUE_ASWNOTIONALIS_ISSUEPRICE "issue_price"
#define XML_VALUE_ASWNOTIONALIS_PUTPRICE "put_price"

//__________________floating leg related______________________________
//
#define XML_TAG_FLOATING_DAYCOUNTCONVENTION "day_count_convention"
#define XML_TAG_FLOATING_STARTOFACCRUEDDATE "start_of_accrued_date"
#define XML_TAG_FLOATING_MULTIPLIER "multiplier"
#define XML_TAG_FLOATING_RECALLSPREAD "recall_spread"
#define XML_TAG_FLOATING_CAP "cap"
#define XML_TAG_FLOATING_FLOOR "floor"
#define XML_TAG_FLOATING_FIXINGDELAY "fixing_delay"
#define XML_TAG_FLOATING_KNOWNCOUPONS "known_coupons"
#define XML_TAG_FLOATING_KNOWNCOUPON "known_coupon"

#define XML_TAG_FLOATING_FIRSTCOUPONDATE  "first_floating_coupon_date"
#define XML_TAG_FLOATING_PRELASTDATE      "last_but_one_floating_coupon_date"
#define XML_TAG_FLOATING_LASTDATE         "last_floating_coupon_date"
#define XML_TAG_FLOATING_PAYMENTFREQUENCY "floating_frequency"
//@}

namespace ito33
{
/// Mapping between finance variables and xml macro names
const EnumValuesNames<finance::ASWNotionalIs>
  g_ASWNotionalIs[] =
    {
      { XML_VALUE_ASWNOTIONALIS_ISSUEPRICE,
        finance::ASWNotionalIs_IssuePrice },
      { XML_VALUE_ASWNOTIONALIS_PUTPRICE,
        finance::ASWNotionalIs_PutPrice }
    };

}

#endif // _ITO33_XML_FINANCE_BONDLIKE_CBOPTION_H_

