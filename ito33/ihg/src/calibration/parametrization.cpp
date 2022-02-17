/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/parametrization.cpp
// Purpose:     implementation for parametrization class
// Author:      ITO 33
// Created:     2004/12/15
// RCS-ID:      $Id: parametrization.cpp,v 1.8 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/debugparameters.h"

#include "ito33/ihg/parametrization.h"
#include "ito33/ihg/theoreticalmodel.h"

namespace ito33
{

namespace ihg
{

std::string Parametrization::GetDebugOutputFile() const
{
  return m_pDebug->GetDebugOutputFile("ihg_calibration.xml");
}

// Default implementation does nothing so derived class is not forced
// to implement this function.
void Parametrization::Calibrate(const finance::BasketGoodType& /*basket*/)
{
}

// Default implementation returns null pointer so derived class is not forced
// to implement this function.
shared_ptr<finance::TheoreticalModel> Parametrization::GetTheoreticalModel()
{
  return shared_ptr<TheoreticalModel>();
}


} // namespace ihg

} // namespace ito33
