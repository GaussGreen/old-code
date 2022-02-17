/////////////////////////////////////////////////////////////////////////////
// Name:      ihg/varianceswapnumoutput.h
// Purpose:   declaration of VarianceSwapNumOutput class 
// Created:   2006/02/21
// RCS-ID:    $Id: varianceswapnumoutput.h,v 1.3 2006/08/20 09:36:16 wang Exp $
// Copyright: (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/varianceswapnumoutput.h
    @brief variance swap numericial output class

    Implementation of the numoutput class for variance swaps.
 */

#ifndef _IHG_VARIANCESWAPNUMOUTPUT_H_
#define _IHG_VARIANCESWAPNUMOUTPUT_H_

#include "ito33/sharedptr.h"
#include "ito33/vector.h"

#include "ito33/ihg/modeloutput.h"

#include "ihg/backwardnumoutput.h"

namespace ito33
{

namespace numeric
{
  class DomainFixedSpaceMesh;
  class SurfaceGeneral;
  class SurfaceFlag;
}

namespace pricing
{
  class VarianceSwapParams;
}

namespace ihg
{

  class VarianceSwapInstData;

/**
    This class stores all pricing information for variance swap, and 
    is responsible for constructing the model output class returned
    to the user.
 */
class VarianceSwapNumOutput : public BackwardNumOutput
{
public:
  
  VarianceSwapNumOutput(pricing::VarianceSwapParams& params) 
    : BackwardNumOutput(), m_params(params)
  { 
  }

  virtual ~VarianceSwapNumOutput() { }

  /**
      Returns the output structure containing the price information.

      @return the model output class containing requested price information
   */
  shared_ptr<ModelOutput> GetModelOutput();

  /**
      Initializes the class variables.

      instdata::SetInitialState must have been called.
   */
  void Init(VarianceSwapInstData& instdata);


protected:
  
  virtual void SaveSurface(BackwardInstData& instdata, double dTime);

  /// The params of the PDE
  pricing::VarianceSwapParams& m_params;


private:

  // Explicitly forbid copy since copying a numoutput doesn't make sense
  NO_COPY_CLASS(VarianceSwapNumOutput);

}; // class VarianceSwapNumOutput


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_VARIANCESWAPNUMOUTPUT_H_
