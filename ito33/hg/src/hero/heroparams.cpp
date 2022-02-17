/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/hero/heroparams.cpp
// Purpose:     Implementation of HG HeroParams class
// Created:     2005/09/26
// RCS-ID:      $Id: heroparams.cpp,v 1.3 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/autoptr.h"

#include "ito33/finance/payoffconstant.h"

#include "hg/heroparams.h"

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(hg::HeroParams);
}

namespace ito33
{

namespace hg
{

void HeroParams::Init()
{
  Params::Init();

  ConstructDividendEvents();

  m_pPayoff = make_ptr( new finance::PayoffConstant(0.0) );
}


} // namespace hg

} // namespace ito33
