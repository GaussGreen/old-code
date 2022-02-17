/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/sessiondata.h
// Purpose:     Names of elements and attributes used in XML for a SessionData
// Created:     2006/03/16
// RCS-ID:      $Id: sessiondata.h,v 1.3 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/sessiondata.h
    @brief Contains the names of the elements used in the XML description of
           the SessionData.
 */

#ifndef _ITO33_XML_FINANCE_SESSIONDATA_H_
#define _ITO33_XML_FINANCE_SESSIONDATA_H_

#include "ito33/sharedptr.h"

/**
    @name Names of the elements and values in the XML description of the
          SessionData.
 */
//@{

/// The name of the root tag containing the session data description
#define XML_TAG_SESSIONDATA                  "session_data"

/// Valuation, or pricing, date 
#define XML_TAG_SESSIONDATA_VALUATION_DATE   "valuation_date"

//@}

namespace xml { class node; }

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL SessionData;
}

namespace XML
{

/**
    Get a SessionData object from XML. 
    
    Exception is thrown if format is wrong.

    @param node the session data node in DOM tree
    @return completed SessionData object from data read from the node
 */
shared_ptr<finance::SessionData> GetSessionDataFromNode(const xml::node& node);

} // namespace XML

} // namespace ito33

#endif // _ITO33_XML_FINANCE_SESSIONDATA_H_

