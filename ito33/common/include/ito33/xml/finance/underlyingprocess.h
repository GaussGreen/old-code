/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/underlyingprocess.h
// Purpose:     Restore underlying process from xml document
// Created:     May 18, 2006
// RCS-ID:      $Id: underlyingprocess.h,v 1.2 2006/06/22 10:11:49 nabil Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/underlyingprocess.h
    @brief Helper function to read underlying process
 */

#ifndef _ITO33_XML_FINANCE_UNDERLYING_PROCESS_H_
#define _ITO33_XML_FINANCE_UNDERLYING_PROCESS_H_

#define XML_TAG_UNDERLYING_PROCESS "underlying_process"

namespace xml { class node; }

namespace ito33
{

namespace XML
{

/**
    Gets the volatility.

    @param node in DOM tree

    @return the volatility after default
 */
double GetPostDefaultVolatilityFromNode(const xml::node &node);


} // namespace XML

} // namespace ito33

#endif //#ifndef _ITO33_XML_FINANCE_UNDERLYING_PROCESS_H_
