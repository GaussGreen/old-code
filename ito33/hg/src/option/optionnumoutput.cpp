/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/option/optionnumoutput.cpp
// Purpose:     implementation of NumOutput for Option
// Created:     2005/01/13
// RCS-ID:      $Id: optionnumoutput.cpp,v 1.5 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"

#include "ito33/pricing/optionparams.h"

#include "hg/backwardnumoutput.h"
#include "hg/optionnumoutput.h"
#include "hg/optioninstdata.h"

// implement the AutoPtrDeleter for OptionNumOutput
namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(hg::OptionNumOutput);
}

namespace ito33
{

namespace hg
{

OptionNumOutput::OptionNumOutput(pricing::OptionParams& params) 
                               : BackwardNumOutput(params), m_params(params)
{ 
}

void OptionNumOutput::Init(OptionInstData& instdata)
{  
  if ( instdata.m_pConstraints )
    m_bHasConstraintFlags = true;

  BackwardNumOutput::Init(instdata);
}

} // namespace hg

} // namespace ito33
