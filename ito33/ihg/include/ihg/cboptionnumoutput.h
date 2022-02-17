/////////////////////////////////////////////////////////////////////////////
// Name:      ihg/cboptionnumoutput.h
// Purpose:   implementation of CBOptionNumOutput class 
// Author:    Nabil
// Created:   2005/10/18
// RCS-ID:    $Id: cboptionnumoutput.h,v 1.6 2006/08/20 09:36:16 wang Exp $
// Copyright: (c) 2004-2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/cboptionnumoutput.h
    @brief cb option numerical output class

    Implementation of the NumOutput class for cb options.
*/

#ifndef _IHG_CBOPTIONNUMOUTPUT_H_
#define _IHG_CBOPTIONNUMOUTPUT_H_

#include "ito33/array.h"
#include "ito33/sharedptr.h"
#include "ito33/vector.h"

#include "ito33/ihg/cboptionoutput.h"
#include "ito33/ihg/bondlikeoutput.h"

#include "ihg/backwardnumoutput.h"
#include "ihg/cbnumoutput.h"
#include "ihg/cboptioninstdata.h"


namespace ito33
{

namespace pricing
{
  class CBLikeParams;
}

namespace numeric
{
  class DomainFixedSpaceMesh;
  class SurfaceGeneral;
}

namespace ihg
{

/**
    This class stores all pricing information for CB options, and is
    responsible for constructing the model output class returned
    to the user.
 */
class CBOptionNumOutput : public BackwardNumOutput
{
public:

  /**
    Constructor by params

    @param params reference to CBLikeParams
    */
  CBOptionNumOutput(pricing::CBLikeParams& params);

  /// virtual dtor
  virtual ~CBOptionNumOutput() {}
  
  /**
     Return the output structure containing the price information

     @return the model output class containing requested price information
  */
  virtual shared_ptr<CBOptionOutput> GetOutput(); 

  /**
     Gets the CBNumOutput object

     @return the CBNumOutput object
  */
  shared_ptr<CBNumOutput> GetCBNumOutput() const
  {
    return m_pCBNumOutput;
  }

  /**
     Initialize the class variables

     instdata::SetInitialState must have been called
  */
  void Init(CBOptionInstData &instdata);

  /**
     see base class
   */
  virtual void UpdateMe(BackwardInstData& instdata, double dTime);
 
  /**
     see base class
   */
  virtual void UpdateMeAtEndOfGrid(BackwardInstData& instdata, double dTime);

  /**
     see base class.
  */
  virtual void Finalize(BackwardInstData& instdata);

  
protected:  

  virtual void SaveSurface(BackwardInstData& instdata, double dTime);

  virtual 
    void SaveSurfaceAtEndOfGrid(BackwardInstData& instdata, double dTime);

  /// The params of the PDE
  pricing::CBLikeParams& m_params;

  /// Numerical output for the underlying cb
  shared_ptr<CBNumOutput> m_pCBNumOutput;


private:

  // Explicitly forbid copy since copying a numoutput doesn't make sense
  NO_COPY_CLASS(CBOptionNumOutput);

}; // class CBOptionNumOutput


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_CBOPTIONNUMOUTPUT_H_
