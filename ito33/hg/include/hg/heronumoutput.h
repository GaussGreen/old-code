/////////////////////////////////////////////////////////////////////////////
// Name:      hg/heronumoutput.h
// Purpose:   HG NumOutput class for HERO
// Created:   2005/09/26
// RCS-ID:    $Id: heronumoutput.h,v 1.6 2006/08/19 23:46:52 wang Exp $
// Copyright: (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/heronumoutput.h
    @brief HG NumOutput class for HERO
 */

#ifndef _HG_HERONUMOUTPUT_H_
#define _HG_HERONUMOUTPUT_H_

#include <cmath>

#include "ito33/sharedptr.h"
#include "ito33/vector.h"

#include "hg/heroinstdata.h"
#include "hg/heroparams.h"
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

/**
    Function object for taking the square root of the computed hero values.
 */
struct HeroSQRT
{

  double operator() (double dValue) const
  {
    if ( dValue > 0.0 )
      return sqrt(dValue);
    else
      return 0.0; 
  }

};


/**
    This class stores all pricing information for HERO
 */
class HeroNumOutput : public BackwardNumOutput
{
public:
  
  HeroNumOutput(HeroParams& params) 
    : BackwardNumOutput(params), m_params(params)
  { 
  }

  virtual ~HeroNumOutput() { }

  /**
      Initializes the class variables.

      instdata::SetInitialState must have been called
   */
  void Init(HeroInstData& instdata);

  virtual shared_ptr<ModelOutput> GetModelOutput();


protected:

  /// The params of the PDE
  HeroParams& m_params;


private:

  // Explicitly forbid copy since copying a numoutput doesn't make sense
  NO_COPY_CLASS(HeroNumOutput);

}; // class HeroNumOutput


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_HERONUMOUTPUT_H_
