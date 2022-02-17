/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/contract.cpp
// Purpose:     Implementation for the base Contract class
// Created:     2005/01/26
// RCS-ID:      $Id: contract.cpp,v 1.1 2005/01/26 18:32:58 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/dateutils.h"

#include "ito33/finance/derivative.h"

#include "ito33/pricing/contract.h"

namespace ito33
{

namespace pricing
{


Contract::Contract(const finance::Derivative& derivative)
                 : Contracts()
{
  m_dMaturityTime = GetDoubleFrom( derivative.GetMaturityDate() );
}


} // namespace pricing

} // namespace ito33
