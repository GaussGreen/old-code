/////////////////////////////////////////////////////////////////////////////
// Name:      hg/cboptionnumoutput.h
// Purpose:   implementation of CBOptionNumOutput class 
// Created:   2006/01/19
// RCS-ID:    $Id: cboptionnumoutput.h,v 1.3 2006/08/19 23:46:52 wang Exp $
// Copyright: (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/cboptionnumoutput.h
    @brief cb option numerical output class

    Implementation of the NumOutput class for cb options.
 */

#ifndef _HG_CBOPTIONNUMOUTPUT_H_
#define _HG_CBOPTIONNUMOUTPUT_H_

#include "ito33/array.h"
#include "ito33/autoptr.h"

#include "ito33/hg/common.h"

#include "hg/backwardnumoutput.h"

namespace ito33
{

namespace pricing
{
  class CBLikeParams;
}

namespace hg
{
  class CBOptionInstData;
  class CBNumOutput;
  class ITO33_HG_DLLDECL CBOptionOutput;

/**
    This class stores all pricing information for CB options, and is
    responsible for constructing the model output class returned
    to the user.
 */
class CBOptionNumOutput : public BackwardNumOutput
{
public:

  /**
      Constructor by params.

      @param params reference to CBLikeParams
   */
  CBOptionNumOutput(pricing::CBLikeParams& params);

  /// virtual dtor
  virtual ~CBOptionNumOutput() {}
  
  /**
      Returns the output structure containing the price information.

      @return the model output class containing requested price information
   */
  shared_ptr<CBOptionOutput> GetCBOptionOutput(); 

  /**
      Gets the CBNumOutput object.

      @return the CBNumOutput object
   */
  CBNumOutput* GetCBNumOutput() const
  {
    return m_pCBNumOutput.get();
  }

  /**
      Initializes the class variables.

      instdata::SetInitialState() must have been called.
   */
  void Init(CBOptionInstData &instdata);

  virtual void UpdateMe(BackwardInstData& instdata, double dTime);
 
  virtual void UpdateMeAtEndOfGrid(BackwardInstData& instdata, double dTime);

  virtual void Finalize(BackwardInstData& instdata);
  
protected:  

  virtual void SaveSurface(BackwardInstData& instdata, double dTime);

  virtual void 
  SaveSurfaceAtEndOfGrid(BackwardInstData& instdata, double dTime);

  /// The params of the PDE
  pricing::CBLikeParams& m_params;

  /// Numerical output for the underlying cb
  AutoPtr<CBNumOutput> m_pCBNumOutput;

private:

  // Explicitly forbid copy since copying a numoutput doesn't make sense
  NO_COPY_CLASS(CBOptionNumOutput);

}; // class CBOptionNumOutput


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_CBOPTIONNUMOUTPUT_H_
