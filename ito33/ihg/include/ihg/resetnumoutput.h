/////////////////////////////////////////////////////////////////////////////
// Name:      ihg/resetnumoutput.h
// Purpose:   implementation of reset class 
// Author:    David and Yann
// Created:   2004/11/08
// RCS-ID:    $Id: resetnumoutput.h,v 1.5 2006/05/02 14:57:47 dave Exp $
// Copyright: (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/resetnumoutput.h
    @brief cb numerical output class

    Implementation of the NumOutput class for resettable.
*/

#ifndef _IHG_RESETNUMOUTPUT_H_
#define _IHG_RESETNUMOUTPUT_H_

#include "ihg/cbnumoutput.h"

namespace ito33
{

namespace pricing
{
  class ResetParams;
}

namespace ihg
{

/**
    This class stores all pricing information for 1D reset pricing (i.e. when
    a similarity reduction is possible). It is responsible for constructing
    the model output class returned to the user.  Full 2D reset pricing
    uses the cbnumoutput class.  The primary differences between the two
    classes are:
    - resetnumoutput checks the current time against the first reset time
      and only saves data between valuation and the first reset. For 2D
      pricing, the path dep pricer only saves data after the last path dep
      event has been applied, so automatically gets this behaviour with
      cbnumoutput
    - resetnumoutput undoes the similarity transformation which is 
      embedded in the price data
 */
class ResetNumOutput : public CBNumOutput
{
public:

  /**
     Constructor by params.

     @param params reference to ResetParams
   */
  ResetNumOutput(pricing::ResetParams& params);

  /// virtual dtor
  virtual ~ResetNumOutput() { }
 
  /// @name re-implement base class functions 
  //@{

  virtual void UpdateMe(BackwardInstData& instdata, double dTime);

  virtual void UpdateMeAtEndOfGrid(BackwardInstData& instdata, double dTime);

  virtual void Finalize(BackwardInstData& instdata);

  virtual void CalculateFinalScalarResult(BackwardInstData& instdata);

  virtual void SaveAnalysisData(BackwardInstData& instdata);
  
  virtual void SaveSurface(BackwardInstData& instdata, double dTime);

  virtual void SaveSurfaceAtEndOfGrid(BackwardInstData& instdata, double dTime);
  //@}

protected:
  
  /// reset params
  pricing::ResetParams& m_resetParams;

  /// The inverse of the current conversion ratio
  double m_dInverseCurrentRatio;

  /// Do not store data past the first (forward time) reset event
  double m_dFirstResetTime;

private:

  // Explicitly forbid copy since copying a numoutput doesn't make sense
  NO_COPY_CLASS(ResetNumOutput);

}; // class ResetNumOutput


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_RESETNUMOUTPUT_H_

