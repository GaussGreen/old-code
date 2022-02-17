/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/cdslike.cpp
// Purpose:     implement base financial cds-like class
// Created:     2006/05/17
// RCS-ID:      $Id: cdslike.cpp,v 1.6 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/cdslike.h"

extern const ito33::finance::Error ITO33_INVALID_RECOVERYRATE_1;

namespace ito33
{

namespace finance
{

CDSLike::CDSLike(double dRecoveryRate)
       : Derivative(), m_dRecoveryRate(dRecoveryRate)
{
  CHECK_COND_1(dRecoveryRate >= 0. && dRecoveryRate <= 1., 
               ITO33_INVALID_RECOVERYRATE_1,
               dRecoveryRate);
}

Date CDSLike::GetIssueDate() const
{
  return GetSpreadStream()->GetContractingDate();
}


} // namespace finance

} // namespace ito33
