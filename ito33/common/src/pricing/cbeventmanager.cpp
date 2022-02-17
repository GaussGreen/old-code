/////////////////////////////////////////////////////////////////////////////
// Name:        pricing/cbeventmanager.cpp
// Purpose:     Manage cb events
// Author:      Nabil
// Created:     2004/04/19
// RCS-ID:      $Id: cbeventmanager.cpp,v 1.6 2004/10/15 16:38:44 wang Exp $
// Copyright:   (c) 1999-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"

#include "ito33/numeric/mesh/specialtimes.h"

#include "ito33/pricing/cbeventmanager.h"

// implement the autoptrdeleter for the CBEventManager class
namespace ito33
{

ITO33_IMPLEMENT_AUTOPTR(pricing::CBEventManager);

} // namespace ito33

namespace ito33
{

namespace pricing
{
  using namespace numeric::mesh;


void CBEventManager::GetSpecialTimes(SpecialTimes& specialTimes) const
{
  specialTimes.clear();

  Events::const_iterator iter;

  for (iter = m_events.begin(); iter != m_events.end(); ++iter)
  {
    RefineLevel refineLevel;

    switch ((*iter)->m_cbeventType)
    {
    case CBET_EndConversion:
    case CBET_EndCall:
    case CBET_StartCall:
    case CBET_StartConversion:
    case CBET_StartAccretion:
      refineLevel = RefineLevel_VeryHigh;
      break;
    default:
      refineLevel = RefineLevel_Standard;
    }

    specialTimes.push_back( SpecialTime((*iter)->GetTime(), refineLevel) );
  }
}


} // namespace pricing

} // namespace ito33
