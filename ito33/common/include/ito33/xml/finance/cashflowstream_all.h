/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/cashflowstream_all.h
// Purpose:     Names of elements and attributes used in XML for cashflowstream
// Author:      Yann d'Halluin, ZHANG
// Created:     2004-06-25
// RCS-ID:      $Id: cashflowstream_all.h,v 1.20 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/cashflowstream_all.h
    @brief Contains the names of the elements used in the XML description of a
           cash flow stream and also reader functions.
 */

#ifndef _ITO33_XML_FINANCE_CASHFLOWSTREAM_ALL_H_
#define _ITO33_XML_FINANCE_CASHFLOWSTREAM_ALL_H_

#include "ito33/common.h"
#include "ito33/enum_values_names.h"

#include "ito33/xml/finance/lastpaymenttype.h"

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
  class ITO33_DLLDECL CashFlowStream;
  class ITO33_DLLDECL CashFlowStreamUniform;
  class ITO33_DLLDECL CashFlowStreamGeneral;
}

namespace XML
{

/**
    Gets CashFlowStream from given node. It throws exception
    when no cash flow stream is found

    @param node given node
    @return shared pointer to CashFlowStream
 */
shared_ptr<finance::CashFlowStream> 
GetCashFlowStreamInNode(const xml::node& node);

/**
    Gets CashFlowStream from given node. It throws exception
    when no cash flow stream is found

    @param node given node
    @return shared pointer to CashFlowStream
 */
shared_ptr<finance::CashFlowStreamUniform>
GetCashFlowStreamUniformInNode(const xml::node& node);

/**
    Gets CashFlowStream from given node. It throws exception
    when no cash flow stream is found

    @param node given node
    @return shared pointer to CashFlowStream
 */
shared_ptr<finance::CashFlowStreamGeneral>
GetCashFlowStreamGeneralInNode(const xml::node& node);

} // namespace XML

} // namespace ito33


/**
   @name Tag name macros
*/
//@{

#define XML_TAG_CASHFLOWSTREAM_DAYCOUNTCONVENTION "day_count_convention"
#define XML_TAG_CASHFLOWSTREAM_CONTRACTINGDATE   "contracting_date"
#define XML_TAG_CASHFLOWSTREAM_PAYMENTFREQUENCY  "payment_frequency"

//__________________cash flow streams uniform related__________________________
//
#define XML_TAG_CASHFLOWSTREAMUNIFORM_ROOT              "cash_flow_stream_uniform"
#define XML_TAG_CASHFLOWSTREAMUNIFORM_FIRSTDATE         "first_date"
#define XML_TAG_CASHFLOWSTREAMUNIFORM_LASTDATE          "last_date"
#define XML_TAG_CASHFLOWSTREAMUNIFORM_ANNUALAMOUNT      "annual_amount"

#define XML_TAG_CASHFLOWSTREAMUNIFORM_DURATION          "duration_in_months"   

//__________________cash flow streams general related__________________________
//
#define XML_TAG_CASHFLOWSTREAMGENERAL_ROOT              "cash_flow_stream_general"
#define XML_TAG_CASHFLOWSTREAMGENERAL_CASHFLOW          "cash_flow"

//@}

#endif // _ITO33_XML_FINANCE_CASHFLOWSTREAM_ALL_H_
