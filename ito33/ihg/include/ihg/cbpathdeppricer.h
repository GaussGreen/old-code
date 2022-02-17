/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/cbpathdeppricer.h
// Purpose:     path dependent pricer class for ihg cb problems
// Author:      Yann and David
// Created:     8/03/2005
// RCS-ID:      $Id: cbpathdeppricer.h,v 1.3 2006/08/03 21:25:18 dave Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/cbpathdeppricer.h
    @brief Path dependent pricer class for ihg cb problems.

    Implementation of the pricer class for ihg cb path dependent contracts.
    This class is identical to the common pathdeppricer except for the
    AdvanceToTime function, which takes the cb multi-mesh into account
    (same difference as for engine and specialengine).
*/

#ifndef _ITO33_IHG_CBPATHDEPPRICER_H_ 
#define _ITO33_IHG_CBPATHDEPPRICER_H_


#include "ito33/pricing/pathdeppricer.h"


namespace ito33
{

namespace pricing
{
  class PathDepStructure;
}

namespace ihg
{

/**
   IHG cb path dependent pricer class
*/

class CBPathDepPricer : public pricing::PathDepPricer
{
public:

  /// Constructor.  Nothing to do.
  CBPathDepPricer() : PathDepPricer()
  { 
  }
  

protected:

  /**
      Advance all paths to the specified time.

      Check for the end of subgrids.

      @param dStopTime is the time to stop timestepping at
      @param pathDepStruct contains all information about the paths
  */
  void AdvanceToTime(double dStopTime, pricing::PathDepStructure& pathDepStruct);

private:

  NO_COPY_CLASS(CBPathDepPricer);

}; // class CBPathDepPricer


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_CBPATHDEPPRICER_H_

