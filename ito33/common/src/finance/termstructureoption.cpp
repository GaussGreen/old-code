/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/termstructureoption.cpp
// Purpose:     Implementation for TermStructureOption
// Created:     2005/03/04
// RCS-ID:      $Id: termstructureoption.cpp,v 1.5 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/termstructure.h"

#include "ito33/finance/termstructureoption.h"

namespace ito33
{

namespace finance
{


void TermStructureOption::Dump(XML::Tag& tagParent) const
{
  XML::Tag tagTS(XML_TAG_TERMSTRUCTURE_OPTION, tagParent);
  Elements::const_iterator i;
  for (i = m_elements.begin(); i != m_elements.end(); ++i)
    (*i)->Dump(tagTS);
}


} // namespace finance

} // namespace ito33
