/////////////////////////////////////////////////////////////////////////////
// Name:        pricing/eventmanager.cpp
// Purpose:     Manage events (dividends, barriers, etc)
// Author:      
// Created:     2003/10/03
// RCS-ID:      $Id: eventmanager.cpp,v 1.17 2004/11/12 17:02:19 zhang Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/debug.h"
#include "ito33/autoptr.h"

#include "ito33/numeric/mesh/specialtimes.h"

#include "ito33/pricing/eventmanager.h"

using ito33::pricing::Event;
using ito33::pricing::EventManager;

using namespace ito33::numeric::mesh;

// implement the autoptrdeleter of the event class
namespace ito33
{

ITO33_IMPLEMENT_AUTOPTR(pricing::EventManager);

} // namespace ito33

void EventManager::GetSpecialTimes(SpecialTimes& specialTimes) const
{
  specialTimes.clear();

  Events::const_iterator iter;

  for (iter = m_events.begin(); iter != m_events.end(); ++iter)
    specialTimes.push_back( SpecialTime((*iter)->GetTime()) );
}
