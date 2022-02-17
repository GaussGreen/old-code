/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/steppertridiag.cpp
// Created:     2004/09/03
// RCS-ID:      $Id: steppertridiag.cpp,v 1.9 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file common/src/pricing/steppertridiag.cpp
    @brief implementation of the class StepperTriDiag
 */

#include "ito33/finance/computationalflags.h"

#include "ito33/pricing/steppertridiag.h"

namespace ito33
{

namespace pricing
{

StepperTriDiag::StepperTriDiag(const finance::ComputationalFlags& flags)
{
  m_iSolverType = flags.GetSolverType();
}

} // namespace pricing

} // namespace ito33
