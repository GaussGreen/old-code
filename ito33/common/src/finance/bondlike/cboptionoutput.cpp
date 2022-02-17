/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/bondlike/cboptionoutput.cpp
// Purpose:     CBOptionOutput class implementation
// Created:     2006/01/03
// RCS-ID:      $Id: cboptionoutput.cpp,v 1.2 2006/08/19 22:43:51 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/sharedptr.h"

#include "ito33/finance/bondlike/cboptionoutput.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/bondlike/convertiblebond.h"

// ============================================================================
// implementation
// ============================================================================

namespace ito33
{

namespace finance
{

void CBOptionOutput::Dump(ito33::XML::Tag& tagParent) const
{
  ModelOutput::Dump(tagParent);
  
  ito33::XML::Tag tagCB(XML_TAG_CONVERTIBLEBOND_ROOT, tagParent);
  
  m_pCBOutput->Dump(tagCB);
}

} // namespace finance

} // namespace ito33
