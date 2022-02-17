/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/option.cpp
// Purpose:     implementation of the contract for option
// Created:     2004/02/10
// RCS-ID:      $Id: option.cpp,v 1.14 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2004 - 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"
#include "ito33/autoptr.h"
#include "ito33/dateutils.h"

#include "ito33/finance/sessiondata.h"

#include "ito33/pricing/option.h"

extern const ito33::Error ITO33_UNEXPECTED;

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(pricing::Option); 

namespace pricing
{


Option::Option(const finance::OptionLike& optionLike)
{
  GetOptionLikeData(optionLike);
}

void Option::GetOptionLikeData(const finance::OptionLike &optionLike)
{
  m_exerciseType  = optionLike.GetExerciseType();
  m_dMaturityTime = GetDoubleFrom( optionLike.GetMaturityDate() );

    // TODO: The cross currency case.
  if( m_bIsCrossCurrency )
  {
    throw EXCEPTION_MSG
        (
          ITO33_UNEXPECTED,
          TRANS("Cross currency option not treated yet")
        );      
  }
  else
    m_pDerivativeCurve = optionLike.GetSessionData()->GetYieldCurve();
  
  // m_dStrik must be initialized by default to initial spot
  // prevent crash in mesh builder
  m_dStrike = optionLike.GetSessionData()->GetSpotSharePrice();
}

Option::Option(const finance::Option& option) 
{
  GetOptionLikeData(option);

  switch ( option.GetOptionType() )
  {
  case finance::Option_Call:
    m_optionType = pricing::Option_Call;
    break;

  case finance::Option_Put:
    m_optionType = pricing::Option_Put;
    break;

  case finance::Option_Digital:
    m_optionType = pricing::Option_Digital;
    break;

  default:
    FAIL("No such option type.");      
  }
  
  m_dStrike = option.GetStrike();

} //Option::Option


} // namespace pricing

} // namespace ito33
