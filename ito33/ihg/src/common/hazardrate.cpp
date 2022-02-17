/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/hazardrate.cpp
// Purpose:     implement sharedptr for HazardRate class
// Author:      (z)
// Created:     03/11/04
// RCS-ID:      $Id: hazardrate.cpp,v 1.17 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/numeric/mesh/specialtimes.h"

#include "ito33/ihg/hazardrate.h"

namespace ito33
{

namespace ihg
{


void HazardRate::GetSpecialTimes(numeric::mesh::SpecialTimes& specialTimes) const
{
  specialTimes.clear();
}

void HazardRate::GetModelParameters
                 (finance::ModelParametersConsumer& /*visitor*/) const
{
}


} // namespace ihg

} // namespace ito33
