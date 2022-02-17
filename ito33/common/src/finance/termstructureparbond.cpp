/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/termstructureparbond.cpp
// Purpose:     do the necessary for TermStructureParBond class
// Author:      ZHANG Yunzhi
// Created:     June 05, 2004
// RCS-ID:      $Id: termstructureparbond.cpp,v 1.3 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/termstructure.h"

#include "ito33/finance/termstructureparbond.h"

namespace ito33
{

namespace finance
{


void TermStructureParBond::Dump(XML::Tag& tagParent) const
{
  XML::Tag tagTS(XML_TAG_TERMSTRUCTURE_PARBOND, tagParent);
  Elements::const_iterator i;
  for (i = m_elements.begin(); i != m_elements.end(); ++i)
    (*i)->Dump(tagTS);
}

#ifndef __CPP2ANY__

void TermStructureParBond::Add(const shared_ptr<ParBond>& parbond)
{
  TermStructure<ParBond>::Add(parbond);
}

#else

void TermStructureParBond::Add(const shared_ptr<ParBond>& parbond)
{
  TermStructure_PARBOND::Add(parbond);
}

#endif // #ifndef __CPP2ANY__


} // namespace finance

} // namespace ito33
