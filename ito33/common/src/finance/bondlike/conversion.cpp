/////////////////////////////////////////////////////////////////////////////
// Name:        src/finance/bondlike/conversion.cpp
// Purpose:     class for common conversion terms
// Created:     2005/01/06 
// RCS-ID:      $Id: conversion.cpp,v 1.4 2006/08/19 22:43:51 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file src/finance/bondlike/conversion.cpp
 */

#include "ito33/finance/bondlike/conversion.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/bondlike/common.h"

namespace ito33
{

namespace finance
{


void Conversion::DumpMe(XML::Tag& tagRoot) const
{
  tagRoot.Element(XML_TAG_BONDLIKE_KEEPACCRUED)(m_bKeepAccrued);
  tagRoot.Element(XML_TAG_BONDLIKE_FORFEITCOUPON)(m_bForfeitCoupon);
}


} // namespace finance

} // namespace ito33
