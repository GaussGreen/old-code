/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/bondlike/percslike.h
// Purpose:     Names of elements and attributes used in XML for 
//              a convertible bond
// Author:      ZHANG Yunzhi
// Created:     2004-09-03
// RCS-ID:      $Id: percslike.h,v 1.5 2006/04/04 16:29:47 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/bondlike/percslike.h
    @brief Contains the names of the elements used in the XML description of a
           percslike instrument

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_BONDLIKE_PERCSLIKE_H_
#define _ITO33_XML_FINANCE_BONDLIKE_PERCSLIKE_H_


/**
   @name Tag name macros
*/
//@{
#define XML_TAG_PERCSLIKE_ROOT                "percs_like"

#define XML_TAG_PERCSLIKE_CONVERSION_RATIO_AT_MATURITY                        \
                                              "conversion_ratio_at_maturity"
                                              
#define XML_TAG_PERCSLIKE_CAP_PRICE           "cap_price"
//@}

#endif // _ITO33_XML_FINANCE_BONDLIKE_percslike_H_

