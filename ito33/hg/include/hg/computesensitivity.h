/////////////////////////////////////////////////////////////////////////////
// Name:        hg/computesensitivity.h
// Purpose:     Compute sensitivity by finite difference
// Created:     2005/05/17
// RCS-ID:      $Id: computesensitivity.h,v 1.6 2006/08/22 18:04:18 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/computesensitivity.h
   @brief Compute sensitivity by finite difference
 */

#ifndef _HG_COMPUTESENSITIVITY_H_
#define _HG_COMPUTESENSITIVITY_H_

#include "ito33/common.h"
#include "ito33/constants.h"

#include <boost/optional/optional.hpp>

namespace ito33
{

  namespace finance
  {
    class ITO33_DLLDECL Derivative;
    class ITO33_DLLDECL QualityControl;
  }

namespace hg
{

  class ITO33_HG_DLLDECL TheoreticalModel;
  class ITO33_HG_DLLDECL UnderlyingProcess;

/**
   Computes the sensitivities by finite difference.

   Since the computational flags within the theoretical model may
   be empty (eg. if the model was cloned), this function calls
   ComputeSpecifiedSensitivites with all elements of the sensitivity
   vector parameter set to true.

   @param model The model to be used to compute the sensitivity
   @param derivative The derivative of which we need the sensitivity
   @param dShift The shift on the parameters
   @param dUnshiftedPrice The unshifted price, if available

   @return The computed sensitivity in a predefined order
 */
std::vector<double> 
ComputeSensitivity
(const TheoreticalModel& model, const finance::Derivative& derivative,
 double dShift = SHIFT,  
 const boost::optional<double>& dUnshiftedPrice = boost::optional<double>() );


/**
   Computes the specified sensitivities by finite difference.

   This function uses finite differencing to compute sensitivities. As 
   such, it can be called with a cloned theoretical model with empty
   computational flags (no point computing extra Greeks, etc).  Thus, 
   the sensitivity flags vector in the theoretical model may be empty 
   and cannot be trusted.  This is why the sensitivity flags are passed 
   in separately.

   @param pUnderlyingProcess The process for which sensitivities are needed
   @param pQualityControl The quality control to use for pricing
   @param derivative The derivative of which we need the sensitivity
   @param dShift The shift on the parameters
   @param dUnshiftedPrice The unshifted price
   @param pbSensitivityFlags Booleans indicating which sensitivities to compute
 
   @return The computed sensitivity in a predefined order
 */
std::vector<double> 
ComputeSpecifiedSensitivities(
                    const shared_ptr<UnderlyingProcess>& pUnderlyingProcess,
                    const shared_ptr<finance::QualityControl>& pQualityControl,
                    const finance::Derivative& derivative,
                    double dShift, 
                    const boost::optional<double>& dPrice,
                    const std::vector<bool>& pdSensitivityFlags);


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_COMPUTESENSITIVITY_H_
