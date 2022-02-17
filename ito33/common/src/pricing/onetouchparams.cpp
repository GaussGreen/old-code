/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/onetouchparams.cpp
// Created:     2005/01/26
// RCS-ID:      $Id: onetouchparams.cpp,v 1.2 2005/07/05 18:34:57 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/extrapolationmode.h"

#include "ito33/pricing/onetouchparams.h"

namespace ito33
{

namespace pricing
{
  using namespace numeric;


void OneTouchParams::Init()
{
  Params::Init();

  ExtrapolationMode emLeft, emRight;

  if ( m_oneTouch.GetBarrierType() == finance::Barrier_UpAndOut )
  {
    emLeft = ExtrapolationMode_Linear;
    emRight = ExtrapolationMode_Constant;
  }
  else
  {
    emLeft = ExtrapolationMode_Constant;
    emRight = ExtrapolationMode_Linear;
  }

  ConstructDividendEvents(emLeft, emRight, InterpolationMethod_Quadratic);
}


} // namespace pricing

} // namespace ito33
