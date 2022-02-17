/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/termstructurecds.cpp
// Purpose:     do the necessary for TermStructureCDS class
// Author:      ZHANG Yunzhi
// Created:     June 05, 2004
// RCS-ID:      $Id: termstructurecds.cpp,v 1.9 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/termstructure.h"

#include "ito33/finance/termstructurecds.h"

namespace ito33
{

namespace finance
{


void TermStructureCDS::Dump(XML::Tag& tagParent) const
{
  XML::Tag tagTS(XML_TAG_TERMSTRUCTURE_CDS, tagParent);
  Elements::const_iterator i;
  for (i = m_elements.begin(); i != m_elements.end(); ++i)
    (*i)->Dump(tagTS);
}

#ifndef __CPP2ANY__

void TermStructureCDS::Add(const shared_ptr<CDSLike>& cds)
{
  TermStructure<CDSLike>::Add(cds);
}

#else

void TermStructureCDS::Add(const shared_ptr<CDSLike>& cds)
{
  TermStructure_CDS::Add(cds);
}

#endif // #ifndef __CPP2ANY__


} // namespace finance

} // namespace ito33
