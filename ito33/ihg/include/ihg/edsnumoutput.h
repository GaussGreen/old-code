/////////////////////////////////////////////////////////////////////////////
// Name:      ihg/edsumoutput.h
// Purpose:   implementation NumOutput for EDS
// RCS-ID:    $Id: edsnumoutput.h,v 1.3 2006/08/20 09:36:16 wang Exp $
// Copyright: (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/edsnumoutput.h
    @brief EDS numericial output class

    @todo This, and CDSNumoutPut can all reduced to a base BackwardNumOutput,
          and returns only a base model output.
 */

#ifndef _IHG_EDSNUMOUTPUT_H_
#define _IHG_EDSNUMOUTPUT_H_

#include "ito33/vector.h"
#include "ito33/autoptr.h"
#include "ito33/array.h"

#include "ihg/edsinstdata.h"
#include "ihg/backwardnumoutput.h"


namespace ito33
{

namespace ihg
{

  class ModelOutput;

/**
    This class stores all pricing information for EDS contracts, and is
    responsible for constructing the model output class returned
    to the user.
 */
class EDSNumOutput : public BackwardNumOutput
{
public:
  
  /**
     Constructor by params and meshes
 
     @param params reference to EDSParams
   */
  EDSNumOutput(pricing::EDSParams& params) 
    : BackwardNumOutput(), m_params(params)
  { 
  }

  /// virtual dtor
  virtual ~EDSNumOutput() { }

  /**
     Return the output structure containing the price information

     @return the model output class containing requested price information
  */
  shared_ptr<ModelOutput> GetOutput();

  /**
     Initialize the class variables

     instdata::SetInitialState must have been called
  */
  void Init(EDSInstData& instdata);


protected:

  virtual void SaveSurface(BackwardInstData& instdata, double dTime);

  /// The params of the PDE
  pricing::EDSParams& m_params;


private:

  NO_COPY_CLASS(EDSNumOutput);

}; // class EDSNumOutput


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_EDSNUMOUTPUT_H_

