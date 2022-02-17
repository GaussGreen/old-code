/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/spotsfxrates.h
// Purpose:     Names of elements and attributes used in XML for 
//              the SpotsFXRates
// Author:      ZHANG Yunzhi
// Created:     2004-09-21
// RCS-ID:      $Id: spotfxrates.h,v 1.8 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/spotfxrates.h
    @brief Contains the names of the elements used in the XML description of a
           spotsfxrates.

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_SPOTFXRATES_H_
#define _ITO33_XML_FINANCE_SPOTFXRATES_H_

#include "ito33/sharedptr.h"

/**
    @name Names of the elements and values in the XML description of the
          SpotFXRates class.
 */
//@{

#define XML_TAG_SPOT_FX_RATES_ROOT         "spot_fx_rates"

#define XML_TAG_SPOT_FX_ROOT               "fx_rate"
#define XML_TAG_SPOT_FX_FOREIGN_CURRENCY   "foreign_currency"
#define XML_TAG_SPOT_FX_BASE_CURRENCY      "base_currency"

//@}

namespace xml { class node; }

namespace ito33
{

namespace XML
{

/**
    Gets a SpotFXRates object from XML. 
    
    Exception is thrown if format is wrong.

    @param node the SpotFXRate tag in DOM tree
    @return completed SpotFXRate object with information we read from node
 */
shared_ptr<finance::SpotFXRates> GetSpotFXRatesFromNode(const xml::node& node);

} // namespace XML

} // namespace ito33

#endif // _ITO33_XML_FINANCE_SPOTFXRATES_H_
