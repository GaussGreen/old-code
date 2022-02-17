/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/spotcomponent.h
// Purpose:     Space component of a volatility or a hazard rate
// Author:      Wang
// Created:     2004/06/04
// RCS-ID:      $Id: spotcomponent.h,v 1.24 2005/12/02 15:52:57 zhang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_IHG_SPOTCOMPONENT_H_
#define _ITO33_IHG_SPOTCOMPONENT_H_

#include "ito33/ihg/common.h"
#include "ito33/ihg/parametrization.h"

namespace ito33
{

namespace XML
{
  class Tag;
}


namespace ihg
{

/**
    Spot component of a volatility or hazard rate

    @nocreate
 */
class ITO33_IHG_DLLDECL SpotComponent
{
public:

  /**
     Ctor takes the reference spot.

     @param dS0 The reference spot

     @noexport
   */
  SpotComponent(double dS0) : m_dS0(dS0) { }

  /// Dummy virtual dtor for base class
  virtual ~SpotComponent() { }

  /**
     @internal
     @brief Sets the reference spot if it's not yet, used by calibration.

     @param dS0 The reference spot
     @noexport
   */
  void TrySetS0(double dS0) 
  { 
    if ( m_dS0 < 0 )
      m_dS0 = dS0; 
  }
   
  /**
     Gets the reference spot.

     @return The reference spot
   */
  double GetS0() const { return m_dS0; }

  /**
     @internal
     @brief Get the values of the spot component at a given time for an array of spots

     @param pdS an ascending array of spots
     @param pdValues the values for given spots
     @param nNbS the number of spots

     @noexport
   */
  virtual void 
  GetValues(const double *pdS, double *pdValues, size_t nNbS) const = 0;

  /**
      @internal
      @brief Dump all data of this spot component in XML format

      Note that this method does @b  dump the contents of the associated
      spot component.

      @param tagParent the parent tag under which our tag(s) should be created
      @param name the name of the tag to create

      @noexport
   */
  virtual ito33::XML::Tag
  Dump(const char *name, ito33::XML::Tag& tagParent) const = 0;

  /**
     @internal
     @brief Serializes the parameters used

     @noexport
   */
  virtual 
    void GetModelParameters(finance::ModelParametersConsumer& visitor,
                               const char* componentName) const = 0;

protected:

  /// The reference spot
  double m_dS0;

}; // class SpotComponent;


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_SPOTCOMPONENT_H_
