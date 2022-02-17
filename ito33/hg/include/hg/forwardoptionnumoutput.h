/////////////////////////////////////////////////////////////////////////////
// Name:      hg/forwardoptionnumoutput.h
// Purpose:   HG NumOutput class for forward options
// Created:   2005/05/05
// RCS-ID:    $Id: forwardoptionnumoutput.h,v 1.9 2006/08/19 23:46:52 wang Exp $
// Copyright: (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/forwardoptionnumoutput.h
    @brief HG NumOutput class for forward options
 */

#ifndef _HG_FORWARDOPTIONNUMOUTPUT_H_
#define _HG_FORWARDOPTIONNUMOUTPUT_H_

#include "ito33/vector.h"
#include "ito33/sharedptr.h"

#include "hg/backwardnumoutput.h"

namespace ito33
{

namespace numeric
{
  class DomainFixedSpaceMesh;
  class SurfaceGeneral;
  class SurfaceFlag;
}

namespace hg
{

class ForwardOptionInstData;
class MultiOutput;

/**
    This class stores all pricing information for forward options, and 
    is responsible for constructing the model output class returned
    to the user.
 */
class ForwardOptionNumOutput : public BackwardNumOutput
{
public:
  
  ForwardOptionNumOutput(pricing::ForwardOptionParams& params) 
    : BackwardNumOutput(params), m_params(params)
  { 
  }

  virtual ~ForwardOptionNumOutput() { }

  /**
     Return the output structure containing the price information

     @return the model output class containing requested price information
   */
  virtual shared_ptr<MultiOutput> GetMultiOutput();

  /**
     Initialize the class variables

     instdata::SetInitialState must have been called
  */
  void Init(ForwardOptionInstData& instdata);

  /**
    Update the output
  */
  void UpdateMe(BackwardInstData& instdata, double dTime);


protected:
  
  virtual void CalculateFinalScalarResult(BackwardInstData& instdata);

  /// The params of the PDE
  pricing::ForwardOptionParams& m_params;

  /// The sensitivities for each of the underlying options
  std::vector< std::vector<double> > m_ppdSensitivities;

private:

  // Explicitly forbid copy since copying a numoutput doesn't make sense
  NO_COPY_CLASS(ForwardOptionNumOutput);

}; // class ForwardOptionNumOutput


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_FORWARDOPTIONNUMOUTPUT_H_
