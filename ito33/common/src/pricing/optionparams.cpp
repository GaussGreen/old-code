/////////////////////////////////////////////////////////////////////////////
// Name:        pricing/optionparams.cpp
// Author:      Wang
// Created:     2004/02/24
// RCS-ID:      $Id: optionparams.cpp,v 1.14 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/autoptr.h"

#include "ito33/finance/optiontype.h"
#include "ito33/finance/payoffoption.h"

#include "ito33/pricing/optionparams.h"
#include "ito33/pricing/option.h"

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::pricing;

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(pricing::OptionParams);
}

void OptionParams::Init()
{
  Params::Init();

  ConstructDividendEvents();

  ConstructPayoff();
}

AutoPtr<OptionParams> OptionParams::Clone()
{
  // Copy the underlying contract. The new params class will manage the
  // memory
  AutoPtr<pricing::Option> clonedOption( new Option(m_option) );

  // Construct and setup the cloned cb params
  AutoPtr<pricing::OptionParams> 
    clonedParams(new OptionParams(clonedOption) );

  clonedParams->SetNumParams(      m_pNumParams );
  clonedParams->SetMeshParams(     m_pMeshParams );
  clonedParams->SetYieldCurve(     GetYieldCurve() );
  clonedParams->SetYieldCurveForMesh( GetYieldCurveForMesh() );
  clonedParams->SetForeignCurve(   GetForeignCurve() );
  clonedParams->SetDividends(      GetDividends() );
  clonedParams->SetValuationTime(  GetValuationTime() );
  clonedParams->SetSpotSharePrice( GetSpotSharePrice() );  

  return clonedParams;
}

void OptionParams::ConstructPayoff()
{
  
  double dStrike = m_option.GetStrike();

  switch ( m_option.GetOptionType() )
  {
  case Option_Call: 
   
    m_pPayoff = shared_ptr<Payoff>( new PayoffCall(dStrike) );
    
    break;
  
  case Option_Put:
  
    m_pPayoff = shared_ptr<Payoff>( new PayoffPut(dStrike) );
    
    break;
 
  case Option_Digital:
  
    m_pPayoff = shared_ptr<Payoff>( new PayoffDigital(dStrike) );
  }
}

