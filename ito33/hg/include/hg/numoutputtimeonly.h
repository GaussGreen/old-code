/////////////////////////////////////////////////////////////////////////////
// Name:      hg/numoutputtimeonly.h
// Purpose:   implementation of cds NumOutput class with HG model
// Created:   2005/06/09
// RCS-ID:    $Id: numoutputtimeonly.h,v 1.6 2006/08/19 23:46:52 wang Exp $
// Copyright: (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/numoutputtimeonly.h
    @brief time only numericial output class using HG model
 */

#ifndef _HG_NUMOUTPUTTIMEONLY_H_
#define _HG_NUMOUTPUTTIMEONLY_H_

#include "ito33/vector.h"
#include "ito33/sharedptr.h"
#include "ito33/array.h"

#include "ito33/finance/modeloutput.h"

#include "hg/backwardnumoutput.h"
#include "hg/instdatatimeonly.h"


namespace ito33
{

namespace numeric
{
  class DomainFixedSpaceMesh;
  class SurfaceGeneral;
}

namespace hg
{

  class ModelOutput;

/**
    This class stores all pricing information for options, and is
    responsible for constructing the model output class returned
    to the user.
 */
class NumOutputTimeOnly : public BackwardNumOutput
{
public:
  
  NumOutputTimeOnly(pricing::Params& params);

  // Default dtor is ok

  /**
     Initialize the class variables

     instdata::SetInitialState must have been called
  */
  void Init(InstDataTimeOnly& instdata);

  /**
     If required, strore the pricing data at the current timestep
  */
  void UpdateMe(InstDataTimeOnly& instdata, double dTime);

  /// finalize
  void Finalize(InstDataTimeOnly& instdata);

  /**
      Sets the final prices.

      A temporary work around for the current time only analytical problem,
      will be removed once proper support for the latter is implemented.

      @param pdPrices the final prices
   */
  void SetPrices(const std::vector<double>& pdPrices);

  virtual shared_ptr<ModelOutput> GetModelOutput();

protected:

  /// TOFIX : this function doesn't make sense here.
  virtual void SaveSurface(BackwardInstData& instdata, double dTime);


private:
 
  /// the price at the analysis date
  double m_dPriceAtAnalysisDate;

  /// the theta at the analysis date
  double m_dThetaAtanalysisDate;

  /// domain/ use correct type.
  shared_ptr<numeric::DomainFixedSpaceMesh> m_pFixedDomain;

  NO_COPY_CLASS(NumOutputTimeOnly);

}; // class NumOutputTimeOnly


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_CDSNUMOUTPUT_H_
