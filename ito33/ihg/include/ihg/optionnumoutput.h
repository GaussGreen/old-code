/////////////////////////////////////////////////////////////////////////////
// Name:      ihg/optionnumoutput.h
// Purpose:   implementation of OptionNumOutput class 
// Author:    ZHANG Yunzhi
// RCS-ID:    $Id: optionnumoutput.h,v 1.25 2006/08/20 09:36:16 wang Exp $
// Copyright: (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/optionnumoutput.h
    @brief option numericial output class

    Implementation of the NumOutput class for options.
*/

#ifndef _IHG_OPTIONNUMOUTPUT_H_
#define _IHG_OPTIONNUMOUTPUT_H_

#include "ito33/sharedptr.h"
#include "ito33/vector.h"

#include "ito33/ihg/modeloutput.h"

#include "ihg/optioninstdata.h"
#include "ihg/backwardnumoutput.h"

namespace ito33
{

namespace numeric
{
  class DomainFixedSpaceMesh;
  class SurfaceGeneral;
  class SurfaceFlag;
}

namespace ihg
{


/**
    This class stores all pricing information for options, and is
    responsible for constructing the model output class returned
    to the user.
 */
class OptionNumOutput : public BackwardNumOutput
{
public:
  
  OptionNumOutput(pricing::OptionParams& params) 
    : BackwardNumOutput(), m_params(params)
  { 
  }

  virtual ~OptionNumOutput() { }

  /**
     Return the output structure containing the price information

     @return the model output class containing requested price information
   */
  shared_ptr<ModelOutput> GetModelOutput();

  /**
     Initialize the class variables

     instdata::SetInitialState must have been called
  */
  void Init(OptionInstData& instdata);


protected:
  
  virtual void SaveSurface(BackwardInstData& instdata, double dTime);

  /// The params of the PDE
  pricing::OptionParams& m_params;


private:

  // Explicitly forbid copy since copying a numoutput doesn't make sense
  NO_COPY_CLASS(OptionNumOutput);

}; // class OptionNumOutput


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_OPTIONNUMOUTPUT_H_

