/////////////////////////////////////////////////////////////////////////////
// Name:      hg/resetnumoutput.h
// Purpose:   implementation of reset class 
// Created:   2006/04/17
// RCS-ID:    $Id: resetnumoutput.h,v 1.1 2006/05/01 20:26:09 dave Exp $
// Copyright: (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/resetnumoutput.h
    @brief reset numerical output class

    Implementation of the NumOutput class for resets.
*/

#ifndef _HG_RESETNUMOUTPUT_H_
#define _HG_RESETNUMOUTPUT_H_

#include "hg/cbnumoutput.h"

namespace ito33
{

namespace pricing
{
  class ResetParams;
}

namespace hg
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


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_RESETNUMOUTPUT_H_

