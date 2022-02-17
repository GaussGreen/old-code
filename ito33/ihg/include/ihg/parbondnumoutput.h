/////////////////////////////////////////////////////////////////////////////
// Name:      ihg/parbondumoutput.h
// Purpose:   implementation of parbond NumOutput class 
// Author:    ZHANG
// RCS-ID:    $Id: parbondnumoutput.h,v 1.3 2006/08/20 09:36:16 wang Exp $
// Copyright: (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/parbondnumoutput.h
    @brief parbond numericial output class

    Implementation of the NumOutput class for ParBond contracts.
 */

#ifndef _IHG_PARBONDNUMOUTPUT_H_
#define _IHG_PARBONDNUMOUTPUT_H_

#include "ito33/vector.h"
#include "ito33/autoptr.h"
#include "ito33/array.h"

#include "ihg/parbondinstdata.h"
#include "ihg/backwardnumoutput.h"


namespace ito33
{

namespace ihg
{

  class ModelOutput;

/**
    This class stores all pricing information for parbond contracts, and is
    responsible for constructing the model output class returned
    to the user.
 */
class ParBondNumOutput : public BackwardNumOutput
{
public:
  
  /**
    Constructor by params and meshes

    @param params reference to ParBondParams
  */
  ParBondNumOutput(pricing::ParBondParams& params) 
    : BackwardNumOutput(), m_params(params)
  { 
  }

  /// virtual dtor
  virtual ~ParBondNumOutput() { }

  /**
     Return the output structure containing the price information

     @return the model output class containing requested price information
  */
  shared_ptr<ModelOutput> GetOutput();

  /**
     Initialize the class variables

     instdata::SetInitialState must have been called
  */
  void Init(ParBondInstData &instdata);


protected:

  virtual void SaveSurface(BackwardInstData& instdata, double dTime);

  /// The params of the PDE
  pricing::ParBondParams& m_params;


private:

  NO_COPY_CLASS(ParBondNumOutput);

}; // class ParBondNumOutput


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_PARBONDNUMOUTPUT_H_

