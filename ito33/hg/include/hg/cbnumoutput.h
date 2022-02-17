/////////////////////////////////////////////////////////////////////////////
// Name:      hg/cbnumoutput.h
// Purpose:   implementation of CBNumOutput class 
// Created:   2005/04/11
// RCS-ID:    $Id: cbnumoutput.h,v 1.6 2006/08/19 23:46:52 wang Exp $
// Copyright: (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/cbnumoutput.h
    @brief cb numerical output class

    Implementation of the NumOutput class for convertible bonds.
 */

#ifndef _HG_CBNUMOUTPUT_H_
#define _HG_CBNUMOUTPUT_H_

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
  class CBInstData;
  class ITO33_HG_DLLDECL BondLikeOutput;

/**
    This class stores all pricing information for CBs, and is
    responsible for constructing the model output class returned
    to the user.
 */
class CBNumOutput : public BackwardNumOutput
{
public:

  /**
      Constructor by params and meshes.

      @param params reference to CBParams
      @param meshes reference to CBMeshManager
   */
  CBNumOutput(pricing::CBLikeParams& params);

  /// virtual dtor
  virtual ~CBNumOutput() { }
  
  /**
      Returns the output structure containing the price information.

      @return the bond like output class containing requested price information
   */
  shared_ptr<BondLikeOutput> GetBondLikeOutput();

  /**
      Initializes the class variables.

      instdata::SetInitialState must have been called.
   */
  void Init(CBInstData& instdata);


protected:  

  // see base class
  virtual void CalculateFinalScalarResult(BackwardInstData& instdata);

  virtual void SaveSurface(BackwardInstData& instdata, double dTime);

  virtual void SaveSurfaceAtEndOfGrid(BackwardInstData& instdata, double dTime);

  /// The params of the PDE
  pricing::CBLikeParams& m_params;

  /// The bond floor
  double m_dBondFloor;


private:

  // Explicitly forbid copy since copying a numoutput doesn't make sense
  NO_COPY_CLASS(CBNumOutput);

}; // class CBNumOutput


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_CBNUMOUTPUT_H_
