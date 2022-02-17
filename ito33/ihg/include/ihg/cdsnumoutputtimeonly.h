/////////////////////////////////////////////////////////////////////////////
// Name:      ihg/cdsnumoutputtimeonly.h
// Purpose:   implementation of time only cds NumOutput class 
// RCS-ID:    $Id: cdsnumoutputtimeonly.h,v 1.24 2006/08/21 14:59:57 wang Exp $
// Copyright: (c) 2003-2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/cdsnumoutputtimeonly.h
    @brief time only cds numericial output class
 */

#ifndef _IHG_CDSNUMOUTPUTTIMEONLY_H_
#define _IHG_CDSNUMOUTPUTTIMEONLY_H_

#include "ito33/sharedptr.h"

#include "ihg/backwardnumoutput.h"

namespace ito33
{

namespace numeric
{
  class DomainFixedSpaceMesh;
  class SurfaceGeneral;
}

namespace ihg
{
  class CDSInstDataTimeOnly;
  class ModelOutput;

/**
    This class stores all pricing information for CDSs, and is
    responsible for constructing the model output class returned
    to the user.
 */
class CDSNumOutputTimeOnly : public BackwardNumOutput
{
public:
  
  CDSNumOutputTimeOnly(pricing::CDSParams& params);

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
  void Init(CDSInstDataTimeOnly& instdata);

  /**
     If required, strore the pricing data at the current timestep
  */
  void UpdateMe(CDSInstDataTimeOnly& instdata, double dTime);

  /// finalize
  void Finalize(CDSInstDataTimeOnly& instdata);


protected:

  /// TOFIX : this function doesn't make sense here.
  virtual void SaveSurface(BackwardInstData& instdata, double dTime);

  /// The params of the PDE
  pricing::CDSParams& m_params;


private:
 
  /// the price at the analysis date
  double m_dPriceAtAnalysisDate;

  /// the theta at the analysis date
  double m_dThetaAtanalysisDate;

  /// domain/ use correct type.
  shared_ptr<numeric::DomainFixedSpaceMesh> m_pFixedDomain;

  NO_COPY_CLASS(CDSNumOutputTimeOnly);

}; // class CDSNumOutputTimeOnly


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_CDSNUMOUTPUTTIMEONLY_H_

