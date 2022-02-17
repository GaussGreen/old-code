/////////////////////////////////////////////////////////////////////////////
// Name:        pricing/forwardcdsparams.cpp
// Author:      David
// Created:     2004/03/31
// RCS-ID:      $Id: forwardcdsparams.cpp,v 1.7 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

//#include "ito33/sharedptr.h"
//#include "ito33/autoptr.h"
//#include "ito33/array.h"

#include "ito33/pricing/eventmanager.h"
#include "ito33/pricing/forwardcdsparams.h"

#include "ito33/finance/payoffoption.h"

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::pricing;

void ForwardCDSParams::Init()
{
  Params::Init();

  // cosntruct the event manager
  m_eventManager.SetProblemType(ProblemType_Forward);

  ConstructForwardDividendEvents();

  ConstructPayoff();
}

void ForwardCDSParams::ConstructPayoff()
{
  m_pPayoff = shared_ptr<Payoff>( new PayoffPut(m_dSpot) ); 
}

