/////////////////////////////////////////////////////////////////////////////
// Name:      ihg/cbnumoutput.h
// Purpose:   implementation of CBNumOutput class 
// Author:    Nabil
// Created:   2004/04/13
// RCS-ID:    $Id: cbnumoutput.h,v 1.22 2006/08/20 09:36:16 wang Exp $
// Copyright: (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/cbnumoutput.h
    @brief cb numerical output class

    Implementation of the NumOutput class for convertible bonds.
*/

#ifndef _IHG_CBNUMOUTPUT_H_
#define _IHG_CBNUMOUTPUT_H_

#include "ito33/beforestd.h"
#include <vector>
#include "ito33/afterstd.h"

#include "ito33/array.h"

#include "ito33/ihg/bondlikeoutput.h"

#include "ihg/cbinstdata.h"
#include "ihg/backwardnumoutput.h"


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
    This class stores all pricing information for CBs, and is
    responsible for constructing the model output class returned
    to the user.
 */
class CBNumOutput : public BackwardNumOutput
{
public:

  /**
    Constructor by params and meshes

    @param params reference to CBParams
    */
  CBNumOutput(pricing::CBLikeParams& params) 
    : BackwardNumOutput(), m_params(params)
  {
    m_bHasConstraintFlags = true;
  }

  /// virtual dtor
  virtual ~CBNumOutput() { }
  
  /**
     Return the output structure containing the price information

     @return the model output class containing requested price information
  */
  virtual shared_ptr<BondLikeOutput> GetOutput();

  /**
     Initialize the class variables

     instdata::SetInitialState must have been called
  */
  void Init(CBInstData &instdata);


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


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_CBNUMOUTPUT_H_

