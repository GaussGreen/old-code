/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/lastpaymenttype.h
// Purpose:     Names of elements and attributes in XML for last payment type
// Created:     2006/06/07
// RCS-ID:      $Id: lastpaymenttype.h,v 1.1 2006/06/07 21:43:22 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/lastpaymenttype.h
    @brief Contains the names of the elements used in the XML description
           of last payment type
 */

#ifndef _ITO33_XML_FINANCE_LASTPAYMENTTYPE_H_
#define _ITO33_XML_FINANCE_LASTPAYMENTTYPE_H_

#include "ito33/common.h"
#include "ito33/enum_values_names.h"

#include "ito33/finance/lastpaymenttype.h"

/**
    @name Tag name macros
 */
//@{

#define XML_TAG_LASTPAYMENTTYPE   "last_payment_type"
#define XML_VALUE_LASTPAYMENTTYPE_SHORT "short"
#define XML_VALUE_LASTPAYMENTTYPE_LONG "long"

//@}

namespace ito33
{

/// Mapping between macro name and finance type
const EnumValuesNames<finance::LastPaymentType>
  g_LastPaymentType[] =
    {
      { XML_VALUE_LASTPAYMENTTYPE_SHORT,
        finance::LastPaymentType_Short },
      { XML_VALUE_LASTPAYMENTTYPE_LONG,
        finance::LastPaymentType_Long }
    };

}

#endif // _ITO33_XML_FINANCE_LASTPAYMENTTYPE_H_
