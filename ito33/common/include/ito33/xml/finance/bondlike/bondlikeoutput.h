/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/bondlike/bondlikeoutput.h
// Purpose:     Restoring BondLikeOutput from XML
// Author:      ITO33 Canada
// Created:     September 22, 2005
// RCS-ID:      $Id: bondlikeoutput.h,v 1.2 2006/02/28 16:13:48 wang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/bondlike/bondlikeoutput.h
    @brief Support for restoring finance::BondLikeOutput from XML.
 */

#ifndef _ITO33_XML_FINANCE_BONDLIKE_BONDLIKEOUTPUT_H_
#define _ITO33_XML_FINANCE_BONDLIKE_BONDLIKEOUTPUT_H_


/// All bond-specific output values are inside this tag
#define XML_TAG_OUTPUT_BONDLIKE "bond"

/// Bond floor value
#define XML_TAG_OUTPUT_BONDLIKE_FLOOR "floor"


namespace xml { class node; }

namespace ito33
{

namespace finance
{
  class BondLikeOutput;
}

namespace XML
{

/**
    Restore a BondLikeOutput object from XML.

    If the node tag doesn't contain subtag corresponding to bond-specific
    object details, an exception is thrown.

    @param node the "output" tag in DOM tree
    @param output to fill in with information we read from node
    @return true if we found all the data we expected, false if some expected
            tags or values were not found
 */
bool Restore(const xml::node& node, finance::BondLikeOutput& output);

} // namespace XML


} // namespace ito33

#endif //_ITO33_XML_FINANCE_BONDLIKE_BONDLIKEOUTPUT_H_

