/////////////////////////////////////////////////////////////////////////////
// Name:      ihg/forwardoptionnumoutput.h
// Purpose:   NumOutput class for European option pricing with forward PDE
// Author:    Wang
// RCS-ID:    $Id: forwardoptionnumoutput.h,v 1.12 2006/08/20 09:36:16 wang Exp $
// Copyright: (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/forwardoptionnumoutput.h
    @brief NumOutput class for European option pricing with forward PDE
 */

#ifndef _IHG_FORWARDOPTIONNUMOUTPUT_H_
#define _IHG_FORWARDOPTIONNUMOUTPUT_H_

#include "ito33/autoptr.h"
#include "ito33/vector.h"

#include "ihg/forwardoptioninstdata.h"
#include "ihg/numoutput.h"

namespace ito33
{

namespace numeric
{
  class DomainFixedSpaceMesh;
  class SurfaceGeneral;
}

namespace finance
{
  class ModelOutput;
}

namespace pricing
{
  class ForwardOptionParams;
}

namespace ihg
{

/**
    This class stores all pricing information for options, and is
    responsible for constructing the model output class returned
    to the user.
 */
class ForwardOptionNumOutput : public NumOutput
{
public:
  
  ForwardOptionNumOutput(pricing::ForwardOptionParams& params) 
    : NumOutput(), m_params(params)
  { 
    m_nNbS = 0;
  }

  // Default dtor is ok

  /**
     Return the output structure containing the price information

     @return the model output class containing requested price information
   */
  shared_ptr<finance::ModelOutput> GetModelOutput();

  /**
     Initialize the class variables

     instdata::SetInitialState must have been called
   */
  void Init(ForwardOptionInstData& instdata);

  /**
     If required, strore the pricing data at the current timestep
   */
  void UpdateMe(ForwardOptionInstData& instdata, double dTime);

  /**
     Finalize the pricing data.

     Called after pricing is complete.  If required, any pointers to external 
     data objects (such as grids, data at the last timestep) must be copied.
   */
  void Finalize(ForwardOptionInstData& instdata);

  shared_ptr<numeric::SurfaceGeneral> GetPriceSurface() 
  {
    return m_pPriceSurface;
  }


protected:

  /// The params of the PDE
  pricing::ForwardOptionParams& m_params;

  /// Temporary price array used at each timestep
  double *m_pdTmpPrices;

  /// The size of the space mesh
  size_t m_nNbS;

  /// The actual space mesh
  std::vector<double> m_pdS;
  
  shared_ptr<numeric::DomainFixedSpaceMesh> m_pDomain;
  
  /// The price surface and array 
  shared_ptr<numeric::SurfaceGeneral> m_pPriceSurface;
  
  CountedArray<double> m_pdValues;


private:

  NO_COPY_CLASS(ForwardOptionNumOutput);

}; // class ForwardOptionNumOutput


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_FORWARDOPTIONNUMOUTPUT_H_

