/////////////////////////////////////////////////////////////////////////////
// Name:      ihg/cdsumoutput.h
// Purpose:   implementation of cds NumOutput class 
// Author:    Wang
// RCS-ID:    $Id: cdsnumoutput.h,v 1.18 2006/08/20 09:36:16 wang Exp $
// Copyright: (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/cdsnumoutput.h
    @brief cds numericial output class

    Implementation of the NumOutput class for CDS contracts.
 */

#ifndef _IHG_CDSNUMOUTPUT_H_
#define _IHG_CDSNUMOUTPUT_H_

#include "ito33/vector.h"
#include "ito33/autoptr.h"
#include "ito33/array.h"

#include "ito33/ihg/modeloutput.h"

#include "ihg/cdsinstdata.h"
#include "ihg/backwardnumoutput.h"


namespace ito33
{

namespace ihg
{

/**
    This class stores all pricing information for cds contracts, and is
    responsible for constructing the model output class returned
    to the user.
 */
class CDSNumOutput : public BackwardNumOutput
{
public:
  
  /**
    Constructor by params and meshes

    @param params reference to CDSParams
  */
  CDSNumOutput(pricing::CDSParams& params) 
    : BackwardNumOutput(), m_params(params)
  { 
  }

  /// virtual dtor
  virtual ~CDSNumOutput() { }

  /**
     Return the output structure containing the price information

     @return the model output class containing requested price information
  */
  shared_ptr<ModelOutput> GetOutput();

  /**
     Initialize the class variables

     instdata::SetInitialState must have been called
  */
  void Init(CDSInstData &instdata);


protected:

  virtual void SaveSurface(BackwardInstData& instdata, double dTime);

  /// The params of the PDE
  pricing::CDSParams& m_params;


private:

  NO_COPY_CLASS(CDSNumOutput);

}; // class CDSNumOutput


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_CDSNUMOUTPUT_H_

