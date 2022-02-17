/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/xml/link_option_xml.cpp
// Purpose:     force linking of xml files related to option
// Created:     2004/11/22
// RCS-ID:      $Id: link_option_xml.cpp,v 1.2 2005/06/27 14:26:59 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// this file must be included in the main project when using XML parsing code
// in a (static) library because otherwise the definition of the static object
// by the macro below would be discarded by the overzealous MSVC linker

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(option_xml);
ITO33_FORCE_LINK_MODULE(asianoption_xml);

ITO33_FORCE_LINK_MODULE(yieldcurve_xml);
