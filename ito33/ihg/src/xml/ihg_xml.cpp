/////////////////////////////////////////////////////////////////////////////
// Name:        ihg_xml.cpp
// Purpose:     force linking of ihg xml files
// Author:      Yann d'Halluin
// Created:     2004-07-17
// RCS-ID:      $Id: ihg_xml.cpp,v 1.2 2004/10/04 18:04:08 pedro Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// this file must be included in the main project when using XML parsing code
// in a (static) library because otherwise the definition of the static object
// by the macro below would be discarded by the overzealous MSVC linker

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(volatility_xml);
ITO33_FORCE_LINK_MODULE(hazardrate_xml);
