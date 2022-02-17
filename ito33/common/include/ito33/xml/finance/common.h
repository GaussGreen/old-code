/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/common.h
// Purpose:     Tag name that are common
// Author:      Yann d'Halluin
// Created:     2004-12-06
// RCS-ID:      $Id: common.h,v 1.17 2006/05/31 10:14:45 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/common.h
    @brief Contains the names of the elements used in the XML description of 
     the different contracts

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_COMMON_H_
#define _ITO33_XML_FINANCE_COMMON_H_

/// @name Common tag name macros. 
//@{

#define XML_TAG_FINANCE_DATES                       "dates"
#define XML_TAG_FINANCE_VALUES                      "values"

#define XML_TAG_FINANCE_DATE                        "date"
#define XML_TAG_FINANCE_VALUE                       "value"

#define XML_TAG_FINANCE_RATE                        "rate"

#define XML_TAG_FINANCE_RECOVERYRATE                "recovery_rate"

#define XML_TAG_FINANCE_STRIKE                      "strike"

#define XML_TAG_FINANCE_YIELD                       "yield"

#define XML_TAG_FINANCE_MATURITY                    "maturity"

#define XML_TAG_FINANCE_FLAT                        "flat"

#define XML_TAG_FINANCE_DIVIDENDS                   "dividends"

#define XML_ATTR_ROOT_VERSION                       "version"

#define XML_TAG_DERIVATIVES                         "derivatives"

#define XML_TAG_DERIVATIVEWEIGHTS                   "derivative_weights"

#define XML_TAG_DERIVATIVEWEIGHT                    "weight"

#define XML_TAG_MODEL                               "model"

#define XML_TAG_OUTPUT                              "output"

#define XML_ATTR_ROOT_VERSION                       "version"

#define XML_TAG_MODEL                               "model"

#define XML_TAG_FINANCE_START_DATE                  "start_date"

#define XML_TAG_FINANCE_END_DATE                    "end_date"

#define XML_TAG_FINANCE_TYPE                        "type"

#define XML_TAG_CURRENCY                            "currency"

#define XML_TAG_POST_DEFAULT_VOLATILITY             "post_default_volatility"
//@}

#endif // _ITO33_XML_FINANCE_COMMON_H_
