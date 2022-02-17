/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/logcontractparams.cpp
// Created:     2006/07/18
// RCS-ID:      $Id: logcontractparams.cpp,v 1.1 2006/07/19 17:38:55 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/pricing/logcontract.h"
#include "ito33/pricing/logcontractparams.h"

namespace ito33
{

namespace pricing
{

void LogContractParams::Init()
{
  Params::Init();

  // Discrete dividends disabled, since log contract is used to valuate
  // variance swap which is usually protected from dividend jump.
}

} // namespace pricing

} // namespace ito33
