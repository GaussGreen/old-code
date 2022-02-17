/////////////////////////////////////////////////////////////////////////////
// Name:      hg/optionnumoutput.h
// Purpose:   HG NumOutput class for regular options
// Created:   2005/01/13
// RCS-ID:    $Id: optionnumoutput.h,v 1.5 2006/08/19 23:46:52 wang Exp $
// Copyright: (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/optionnumoutput.h
    @brief HG NumOutput class for regular options
 */

#ifndef _HG_OPTIONNUMOUTPUT_H_
#define _HG_OPTIONNUMOUTPUT_H_

#include "ito33/sharedptr.h"
#include "ito33/vector.h"

#include "hg/optioninstdata.h"
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

class ModelOutput;

/**
    This class stores all pricing information for options, and is
    responsible for constructing the model output class returned
    to the user.
 */
class OptionNumOutput : public BackwardNumOutput
{
public:
  
  OptionNumOutput(pricing::OptionParams& params);

  virtual ~OptionNumOutput() { }

  /**
      Initialize the class variables.

      instdata::SetInitialState must have been called.
   */
  void Init(OptionInstData& instdata);


protected:

  /// The params of the PDE
  pricing::OptionParams& m_params;


private:

  // Explicitly forbid copy since copying a numoutput doesn't make sense
  NO_COPY_CLASS(OptionNumOutput);

}; // class OptionNumOutput


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_OPTIONNUMOUTPUT_H_
