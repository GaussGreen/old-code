/////////////////////////////////////////////////////////////////////////////
// Name:        conversion_xml.cpp
// Purpose:     Restore common conversion datas from XML document
// Created:     2005/01/06
// RCS-ID:      $Id: conversion_xml.cpp,v 1.3 2006/06/03 19:39:04 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/bondlike/conversion.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/bondlike/common.h"
#include "ito33/xml/finance/bondlike/conversion.h"

namespace ito33
{

namespace XML
{ 
  

void RestoreCommonConversionData
     (const xml::node& node, finance::Conversion& conversion)
{
  conversion.SetKeepAccrued
             ( GetBoolFromName(node, XML_TAG_BONDLIKE_KEEPACCRUED) );

  conversion.SetForfeitCoupon
             ( GetBoolFromName(node, XML_TAG_BONDLIKE_FORFEITCOUPON) );
}


} // namespace XML

} // namespace ito33
