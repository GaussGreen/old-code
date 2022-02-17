/////////////////////////////////////////////////////////////////////////////
// Name:      hg/varianceswapnumoutput.h
// Purpose:   HG NumOutput class for variance swaps
// Created:   2006/03/05
// RCS-ID:    $Id: varianceswapnumoutput.h,v 1.3 2006/08/02 21:02:38 wang Exp $
// Copyright: (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/varianceswapnumoutput.h
    @brief HG NumOutput class for variance swaps
 */

#ifndef _HG_VARIANCESWAPNUMOUTPUT_H_
#define _HG_VARIANCESWAPNUMOUTPUT_H_

#include "ito33/sharedptr.h"
#include "ito33/vector.h"

#include "ito33/pricing/varianceswapparams.h"

#include "hg/varianceswapinstdata.h"
#include "hg/backwardnumoutput.h"

namespace ito33
{

namespace hg
{

/**
    This class stores all pricing information for variance swaps, and is
    responsible for constructing the model output class returned
    to the user.
 */
class VarianceSwapNumOutput : public BackwardNumOutput
{
public:
  
  VarianceSwapNumOutput(pricing::VarianceSwapParams& params) 
    : BackwardNumOutput(params), m_params(params)
  { 
  }

  virtual ~VarianceSwapNumOutput() { }

  /**
      Initializes the class variables.
 
      instdata::SetInitialState must have been called
   */
  void Init(VarianceSwapInstData& instdata);


protected:

  /// The params of the PDE
  pricing::VarianceSwapParams& m_params;


private:

  // Explicitly forbid copy since copying a numoutput doesn't make sense
  NO_COPY_CLASS(VarianceSwapNumOutput);

}; // class VarianceSwapNumOutput


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_VARIANCESWAPNUMOUTPUT_H_
