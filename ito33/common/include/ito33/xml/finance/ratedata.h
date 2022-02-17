/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/ratedata.h
// Purpose:     Names of elements and attributes used in XML for RateData
// Created:     2006/03/17
// RCS-ID:      $Id: ratedata.h,v 1.3 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/ratedata.h
    @brief Contains the names of the elements used in the XML description of a
           RateData object.

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_RATEDATA_H_
#define _ITO33_XML_FINANCE_RATEDATA_H_

#include "ito33/sharedptr.h"

/**
    @name Names of the elements and values in the XML description of the
          Equity class.
 */
//@{

#define XML_TAG_RATEDATA_ROOT                  "rate_data"

//@}


namespace xml { class node; }

namespace ito33
{

namespace XML
{

/**
   Get a RateData object from XML. 
  
   Exception is thrown if format is wrong.
 
   @param node RateData node in DOM tree
   @return completed RateData object from data read from the node
 */
shared_ptr<finance::RateData> GetRateDataFromNode(const xml::node& node);


} // namespace XML

} // namespace ito33


#endif // _ITO33_XML_FINANCE_SPOTFXRATES_H_

