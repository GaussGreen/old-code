/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/hrspotcomponentpower.h
// Purpose:     spot power hazard rate class
// Created:     2004/06/11
// RCS-ID:      $Id: hrspotcomponentpower.h,v 1.23 2006/03/02 18:04:19 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/ihg/hrspotcomponentpower.h
   @brief spot power function as spot component of hazard rate
 */

#ifndef _ITO33_IHG_HRSPOTCOMPONENTPOWER_H_
#define _ITO33_IHG_HRSPOTCOMPONENTPOWER_H_

#include "ito33/ihg/spotcomponent.h"

namespace ito33
{

namespace ihg
{


/** 
    Class for a spot component as a power function of the spot.

    \f$ \lambda (S,t) = (S_0 / S)^\beta \f$ 

    This is a typical parametrization of a space only hazard rate.
 */
class ITO33_IHG_DLLDECL HRSpotComponentPower : public SpotComponent
{
public:

  /**
     Constructor.

     @param dBeta beta
     @param dS0 s0
   */
  HRSpotComponentPower(double dBeta, double dS0);

  /**
     Constructor takes only beta. S_0 will be determined during calibration.

     @param dBeta beta
     
     @noexport
   */
  HRSpotComponentPower(double dBeta);

  // Default dtor is ok

  /**
     Gets the value of beta.

     @return beta
   */
  double GetBeta() const { return m_dBeta; }
 
  void GetValues(const double *pdS, double *pdValues, size_t nNbS) const;

  virtual ito33::XML::Tag Dump(const char *name, XML::Tag& tagParent) const;

  virtual 
  void GetModelParameters(finance::ModelParametersConsumer& visitor,
                             const char* componentName) const;

private:

  // validate beta 
  void DoValidate();

  double m_dBeta;

}; // class HRSpotComponentPower


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_HRSPOTCOMPONENTPOWER_H_
