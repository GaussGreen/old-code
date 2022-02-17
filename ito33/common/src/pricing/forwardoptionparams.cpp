/////////////////////////////////////////////////////////////////////////////
// Name:        pricing/forwardoptionparams.cpp
// Author:      Wang
// Created:     2004/03/04
// RCS-ID:      $Id: forwardoptionparams.cpp,v 1.11 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/autoptr.h"
#include "ito33/array.h"

#include "ito33/finance/optiontype.h"
#include "ito33/finance/payoffoption.h"

#include "ito33/pricing/problemtype.h"
#include "ito33/pricing/eventmanager.h"
#include "ito33/pricing/forwardoptionparams.h"

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::pricing;

void ForwardOptionParams::Init()
{
  Params::Init();

  // construct the event manager
  m_eventManager.SetProblemType(ProblemType_Forward);

  ConstructForwardDividendEvents();

  ConstructPayoff();
}

void ForwardOptionParams::ConstructPayoff()
{
   
  //everything is priced to be calls
  m_pPayoff = shared_ptr<Payoff>( new PayoffPut(m_dSpot) );

}

