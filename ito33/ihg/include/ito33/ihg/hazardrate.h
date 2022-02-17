/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/hazardrate.h
// Purpose:     base hazard rate class
// Author:      (z)
// Created:     03/11/04
// RCS-ID:      $Id: hazardrate.h,v 1.45 2005/12/02 15:52:57 zhang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_IHG_HAZARDRATE_H_
#define _ITO33_IHG_HAZARDRATE_H_

#include "ito33/sharedptr.h"
#include "ito33/vector.h"

#include "ito33/ihg/common.h"
#include "ito33/ihg/parametrization.h"

namespace ito33
{

namespace finance
{
  class ModelParametersConsumer;
}

namespace numeric
{
  namespace mesh
  {
    class SpecialTimes;
  }
}

namespace XML
{
  class Tag;
}

namespace ihg
{

class HazardRateVisitor;


/**
    Base class for all kinds of hazard rates (flat, parametrized, etc).

    @nocreate
 */
class ITO33_IHG_DLLDECL HazardRate
{
public:

  /// Dummy virtual dtor for base class
  virtual ~HazardRate() { }

  /**
     @internal
     @brief Gets the values of the hazard rate at a given time for an array of spots.

     @param dTime the time at which the values will be obtained
     @param pdS an ascending array of spots
     @param pdValues the (computed) hazard rate values
     @param nNbS the number of spots

     @noexport
  */
  virtual void GetHazardRates(double dTime,
                              const double *pdS,
                              double *pdValues,
                              size_t nNbS) const = 0;
 
  /*
     By default, it returns false. Otherwise, specific hazard rate
     classes must reimplement this function.
  */
  
  /**
     @internal
     @brief Return whether or not the hazard rate is time only.

     Certain optimizations can be made for time only (space independent)
     hazard rates. 

     @return true if the hazard rate is time only, false otherwise

     @noexport
   */
  virtual bool IsTimeOnly() const
  {
    return false;
  }

  /**
     @internal
     @brief Gets the special times in the hazard rate surface.

     Used by the mesh manager to force these special points in the
     mesh.  The mesh is also refined around these points.
     Helps with the convergence when the hazard rate is non-smooth
     or discontinuous.

     By default, no special time is considered

     @noexport
   */
  virtual void GetSpecialTimes(numeric::mesh::SpecialTimes& specialTimes) const;

  /**
      @internal
      @brief Dumps all data of this instrument in XML format.

      Note that this method does @b  dump the contents of the associated
      hazard rate.

      @param tagParent the parent tag under which our tag(s) should be created
      
      @noexport
   */
  virtual void Dump(ito33::XML::Tag& tagParent) const = 0;

  /**
      @internal
      @brief This method is part of the implementation of the visitor pattern.

      Visitor pattern allows to do different things for different kinds of
      hazard without adding virtual functions to do all of them in the
      base hazard rate class itself nor doing the dreaded switchs on typeid of
      the concrete type. Instead, define a new class deriving from
      HazardRateVisitor and do whatever is required in its methods.

      @noexport
   */ 
  virtual void Visit(HazardRateVisitor& visitor) const = 0;

    
  /**
     @internal
     @brief Serializes the parameters used

     @noexport
  */
  virtual 
    void GetModelParameters(finance::ModelParametersConsumer& visitor) const;

}; // class HazardRate


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_HAZARDRATE_H_

