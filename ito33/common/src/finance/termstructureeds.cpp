/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/termstructureeds.cpp
// Purpose:     do the necessary for TermStructureEDS class
// Author:      ITO33
// Created:     2005/02/10
// RCS-ID:      $Id: termstructureeds.cpp,v 1.4 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2005- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/termstructure.h"

#include "ito33/finance/termstructureeds.h"

namespace ito33
{

namespace finance
{


void TermStructureEDS::Dump(XML::Tag& tagParent) const
{
  XML::Tag tagTS(XML_TAG_TERMSTRUCTURE_EDS, tagParent);
  Elements::const_iterator i;
  for (i = m_elements.begin(); i != m_elements.end(); ++i)
    (*i)->Dump(tagTS);
}


#ifndef __CPP2ANY__

void TermStructureEDS::Add(const shared_ptr<EDS>& eds)
{
  TermStructure<EDS>::Add(eds);
}

#else

void TermStructureEDS::Add(const shared_ptr<EDS>& eds)
{
  TermStructure_EDS::Add(eds);
}

#endif // #ifndef __CPP2ANY__


} // namespace finance

} // namespace ito33
