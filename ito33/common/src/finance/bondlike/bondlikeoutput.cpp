/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/bondlike/bondlikeoutput.cpp
// Purpose:     implementation of finance::BondLikeOutput
// Created:     September 21, 2005
// RCS-ID:      $Id: bondlikeoutput.cpp,v 1.3 2006/08/19 22:43:51 wang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/bondlike/bondlikeoutput.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/bondlike/bondlikeoutput.h"

namespace ito33
{

namespace finance
{

void BondLikeOutput::Dump(ito33::XML::Tag& tagParent) const
{
  ModelOutput::Dump(tagParent);

  ito33::XML::Tag tagBond(XML_TAG_OUTPUT_BONDLIKE, tagParent);
  tagBond.Element(XML_TAG_OUTPUT_BONDLIKE_FLOOR)(m_dBondFloor);
}

} // namespace finance

} // namespace ito33
