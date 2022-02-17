/////////////////////////////////////////////////////////////////////////////
// Name:      ihg/parbondnumoutputtimeonly.h
// Purpose:   implementation of time only parbond NumOutput class 
// Author:    ZHANG
// RCS-ID:    $Id: parbondnumoutputtimeonly.h,v 1.3 2006/08/20 09:36:16 wang Exp $
// Copyright: (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/parbondnumoutputtimeonly.h
    @brief time only parbond numericial output class
 */

#ifndef _IHG_ParBondNUMOUTPUTTIMEONLY_H_
#define _IHG_ParBondNUMOUTPUTTIMEONLY_H_

#include "ito33/vector.h"
#include "ito33/sharedptr.h"
#include "ito33/array.h"
#include "ito33/dateutils.h"

#include "ihg/backwardnumoutput.h"
#include "ihg/parbondinstdatatimeonly.h"


namespace ito33
{

namespace numeric
{
  class DomainFixedSpaceMesh;
  class SurfaceGeneral;
}

namespace ihg
{

  class ModelOutput;

/**
    This class stores all pricing information for options, and is
    responsible for constructing the model output class returned
    to the user.
 */
class ParBondNumOutputTimeOnly : public BackwardNumOutput
{
public:
  
  ParBondNumOutputTimeOnly(pricing::ParBondParams& params);

  // Default dtor is ok

  /**
     Return the output structure containing the price information

     @return the model output class containing requested price information
   */
  shared_ptr<ModelOutput> GetOutput();

  /**
     Initialize the class variables

     instdata::SetInitialState must have been called
  */
  void Init(ParBondInstDataTimeOnly& instdata);

  /**
     If required, strore the pricing data at the current timestep
  */
  void UpdateMe(ParBondInstDataTimeOnly& instdata, double dTime);

  /// finalize
  void Finalize(ParBondInstDataTimeOnly& instdata);


protected:

  /// TOFIX : this function doesn't make sense here.
  virtual void SaveSurface(BackwardInstData& instdata, double dTime);

  /// The params of the PDE
  pricing::ParBondParams& m_params;


private:
 
  /// the price at the analysis date
  double m_dPriceAtAnalysisDate;

  /// the theta at the analysis date
  double m_dThetaAtanalysisDate;

  /// domain/ use correct type.
  shared_ptr<numeric::DomainFixedSpaceMesh> m_pFixedDomain;

  NO_COPY_CLASS(ParBondNumOutputTimeOnly);

}; // class ParBondNumOutputTimeOnly


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_ParBondNUMOUTPUTTIMEONLY_H_

