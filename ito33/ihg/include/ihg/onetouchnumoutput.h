/////////////////////////////////////////////////////////////////////////////
// Name:      ihg/onetouchumoutput.h
// Purpose:   implementation NumOutput for OneTouch
// Created:   2006/08/11
// RCS-ID:    $Id: onetouchnumoutput.h,v 1.2 2006/08/20 09:36:16 wang Exp $
// Copyright: (c) 2006  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/onetouchumoutput.h
    @brief OneTouch numericial output class

    @todo This, and EDSNumoutPut can all reduced to a base BackwardNumOutput,
          and returns only a base model output. See HG.
 */

#ifndef _IHG_ONETOUCHNUMOUTPUT_H_
#define _IHG_ONETOUCHNUMOUTPUT_H_

#include "ito33/vector.h"
#include "ito33/autoptr.h"
#include "ito33/array.h"

#include "ihg/onetouchinstdata.h"
#include "ihg/backwardnumoutput.h"


namespace ito33
{

namespace ihg
{

  class ModelOutput;

/**
    This class stores all pricing information for OneTouch contracts, and is
    responsible for constructing the model output class returned
    to the user.
 */
class OneTouchNumOutput : public BackwardNumOutput
{
public:
  
  /**
      Constructor by params and meshes.
 
      @param params reference to OneTouchParams
   */
  OneTouchNumOutput(pricing::OneTouchParams& params) 
                  : BackwardNumOutput(), m_params(params)
  { 
  }

  /// virtual dtor
  virtual ~OneTouchNumOutput() { }
  
  /**
      Initializes the class variables
 
      instdata::SetInitialState must have been called.
   */
  void Init(OneTouchInstData& instdata);

  /**
      Returns the output structure containing the price information

      @return the model output class containing requested price information
   */
  shared_ptr<ModelOutput> GetOutput();

protected:

  virtual void SaveSurface(BackwardInstData& instdata, double dTime);

  /// The params of the PDE
  pricing::OneTouchParams& m_params;


private:

  NO_COPY_CLASS(OneTouchNumOutput);

}; // class OneTouchNumOutput


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_ONETOUCHNUMOUTPUT_H_

