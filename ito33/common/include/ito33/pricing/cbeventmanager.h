/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/cbeventmanager.h
// Purpose:     Manage cb events
// Author:      Nabil
// Created:     2004/04/19
// RCS-ID:      $Id: cbeventmanager.h,v 1.6 2006/03/22 13:10:23 yann Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/cbeventmanager.h
    @brief Class to manage cb events
 */

#ifndef _ITO33_PRICING_CBEVENTMANAGER_H_
#define _ITO33_PRICING_CBEVENTMANAGER_H_

#include "ito33/pricing/cbevent.h"
#include "ito33/pricing/baseeventmanager.h"

namespace ito33
{

namespace numeric
{
  namespace mesh
  {
    class SpecialTimes;
  }
}

namespace pricing
{

/// Declaration of the CBEventManager
class CBEventManager: public BaseEventManager<CBEvent>
{

public:

  /** 
     Get a list of event times as special times for time mesh.

     @param specialTimes special times to be filled
   */
  void GetSpecialTimes(numeric::mesh::SpecialTimes& specialTimes) const;

}; // class CBEventManager


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_CBEVENTMANAGER_H_

