/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/bondlike/common.h
// Purpose:     Tag name that are common
// Author:      Yann d'Halluin
// Created:     2004-10-01
// RCS-ID:      $Id: common.h,v 1.8 2006/04/04 16:29:47 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/bondlike/common.h
    @brief Contains the names of the elements used in the XML description of 
     the different bond characteristic

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_BONDLIKE_COMMON_H_
#define _ITO33_XML_FINANCE_BONDLIKE_COMMON_H_

/**
   @name common tag name macros.
*/
//@{
#define XML_TAG_BONDLIKE_STARTDATE                   "start_date"

#define XML_TAG_BONDLIKE_ENDDATE                     "end_date"

#define XML_TAG_BONDLIKE_KEEPACCRUED                 "keep_accrued"

#define XML_TAG_BONDLIKE_FORFEITCOUPON               "forfeit_coupon"

#define XML_TAG_BONDLIKE_NOTICEPERIOD                "notice_period"

#define XML_TAG_BONDLIKE_COUPON                      "coupon"

#define XML_TAG_BONDLIKE_TRIGGERRATE                 "trigger_rate"

#define  XML_TAG_BONDLIKE_RATIO                      "ratio"

#define XML_TAG_BONDLIKE_KEEPACCRUED                 "keep_accrued"

#define XML_TAG_BONDLIKE_FORFEITCOUPON               "forfeit_coupon"

#define XML_TAG_BONDLIKE_ISSUEDATE                   "issue_date"

#define XML_TAG_BONDLIKE_ISSUEPRICE                  "issue_price"

#define XML_TAG_BONDLIKE_NOMINAL                     "nominal"

#define XML_VALUE_MANDATORY_FIXED_CASH               "fixed_cash"

#define XML_VALUE_MANDATORY_FIXED_SHARE              "fixed_share"

#define XML_VALUE_MANDATORY_VARIABLE_SHARE           "variable_share"
//@}

#endif // _ITO33_XML_FINANCE_BONDLIKE_COMMON_H_

