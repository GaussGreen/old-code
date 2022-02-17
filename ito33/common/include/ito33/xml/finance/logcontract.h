/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/logcontract.h
// Purpose:     Names of elements and attributes used in XML for a LogContract.
// Created:     2006/07/18
// RCS-ID:      $Id: logcontract.h,v 1.2 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/logcontract.h
    @brief Contains the names of the elements used in the XML description of a
           log contract.
 */

#ifndef _ITO33_XML_FINANCE_LOGCONTRACT_H_
#define _ITO33_XML_FINANCE_LOGCONTRACT_H_

#include "ito33/sharedptr.h"

/**
   @name Tag name macros
*/
//@{

#define XML_TAG_LOGCONTRACT_ROOT          "log_contract"

#define XML_TAG_LOGCONTRACT_T0            "start_of_return_period"
 
#define XML_TAG_LOGCONTRACT_S0            "start_share_price"

//@}

namespace xml
{
  class node;
}

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL LogContract;
}

namespace XML
{

/// Restore a LogContract from xml node
bool
Restore(const xml::node& node, shared_ptr<finance::LogContract>& pLogContract);

} // namespace XML

} // namespace ito33


#endif // _ITO33_XML_FINANCE_LOGCONTRACT_H_
