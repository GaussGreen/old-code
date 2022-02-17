/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/bondlike/conversion.h
// Purpose:     function to restore common conversion datas
// Created:     2005/01/06
// RCS-ID:      $Id: conversion.h,v 1.3 2006/01/10 16:53:16 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/bondlike/conversion.h
    @brief function to restore common conversion datas
 */

#ifndef _ITO33_XML_FINANCE_BONDLIKE_CONVERSION_H_
#define _ITO33_XML_FINANCE_BONDLIKE_CONVERSION_H_

#include "ito33/common.h"

namespace xml
{
  class node;
}

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL Conversion;
}

namespace XML
{

/**
   Restore common conversion terms in given node.

   @param node given node, normally root node of conversion
   @param conversion (output)
 */
void RestoreCommonConversionData
     (const xml::node& node, finance::Conversion& conversion);

} // namespace XML

} // namespace ito33

#endif // _ITO33_XML_FINANCE_BONDLIKE_CONVERSION_H_

