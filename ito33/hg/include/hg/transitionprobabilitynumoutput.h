/////////////////////////////////////////////////////////////////////////////
// Name:      hg/transitionprobabilitynumoutput.h
// Purpose:   implementation of HG numoutput class for transition probability
// Created:   2006/03/31 
// RCS-ID:    $Id: transitionprobabilitynumoutput.h,v 1.2 2006/08/19 23:46:52 wang Exp $
// Copyright: (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/transitionprobabilitynumoutput.h
    @brief implementation of HG numoutput class for transition probability
 */

#ifndef _HG_TRANSITIONPROBABILITYNUMOUTPUT_H_
#define _HG_TRANSITIONPROBABILITYNUMOUTPUT_H_

#include "ito33/sharedptr.h"
#include "ito33/vector.h"

#include "ito33/finance/modeloutput.h"

#include "ito33/numeric/surfacedouble.h"

#include "hg/backwardnumoutput.h"
#include "hg/sensitivitybyadjointdata.h"

namespace ito33
{

namespace numeric
{
  class SurfaceGeneral;
  class SurfaceFlag;
}

namespace pricing
{
  class Params;
}

namespace hg
{

class TransitionProbabilityInstData;
class TransitionProbabilityOutput;

/// This class stores all pricing information for transition probability.
class TransitionProbabilityNumOutput : public BackwardNumOutput
{
public:

  /**
      Constructor by params.
 
      @param params reference to pricing params
   */
  TransitionProbabilityNumOutput(pricing::Params& params) 
                               : BackwardNumOutput(params)
  {
  }

  /**
      Solves the transpose system of the backward to get the transition
      probability.
   */
  shared_ptr<TransitionProbabilityOutput>
  Solve(TransitionProbabilityInstData& instdata);

private:

  // Explicitly forbid copy since copying a numoutput doesn't make sense
  NO_COPY_CLASS(TransitionProbabilityNumOutput);

}; // class TransitionProbabilityNumOutput


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_TRANSITIONPROBABILITYNUMOUTPUT_H_
