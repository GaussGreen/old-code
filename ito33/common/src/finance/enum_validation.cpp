/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/enum_validation.cpp
// Purpose:     validate enum in finance namespace
// Created:     2005/05/20
// author:      ZHANG Yunzhi
// RCS-ID:      $Id: enum_validation.cpp,v 1.3 2006/06/15 18:57:04 zhang Exp $
// Copyright:   (c) 2003 - 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/error.h"
#include "ito33/useexception.h"
#include "ito33/finance/frequency.h"
#include "ito33/finance/lastpaymenttype.h"

extern const ito33::finance::Error
  ITO33_INVALID_FREQUENCY,
  ITO33_INVALID_LASTPAYMENTTYPE;

namespace ito33
{

namespace finance
{

void ThrowInvalidFrequency()
{
  throw EXCEPTION(ITO33_INVALID_FREQUENCY);
}

void ThrowInvalidLastPaymentType()
{
  throw EXCEPTION(ITO33_INVALID_LASTPAYMENTTYPE);
}

} // namespace finance

} // namespace ito33
