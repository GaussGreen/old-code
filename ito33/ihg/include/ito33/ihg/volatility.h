/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/volatiliy.h
// Purpose:     interface of volatilities(flat, grid, parametrized etc)
// Author:      z
// Created:     03.09.23
// RCS-ID:      $Id: volatility.h,v 1.46 2006/08/20 09:47:41 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/volatility.h
 */

#ifndef _ITO33_IHG_VOLATILITY_H_
#define _ITO33_IHG_VOLATILITY_H_

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

namespace XML
{
  class Tag;
}

namespace numeric
{
  namespace mesh
  {
    class SpecialTimes;
  }
}

namespace ihg
{

  class VolatilityVisitor;


/**
    Base class for all kinds of volatilities (flat, parametrized, etc)

    @nocreate
 */
class ITO33_IHG_DLLDECL Volatility
{
public:
  
  // default empty ctor is ok
  
  /// dummy virtual dtor
  virtual ~Volatility() { }

  /**
     @internal
     @brief Get the values of the volatility at a given time for an array of spots

     @param dTime the time at which the values will be obtained
     @param pdS an ascending array of spots
     @param pdVols the (computed) values of valitity
     @param nNbS the number of spots
     
     @noexport
   */
  virtual void GetVols(double dTime,
                       const double *pdS,
                       double *pdVols,
                       size_t nNbS) const = 0;

  /**
     @internal
     Get the squares of the volatility at a given time for an array of spot

     By default, it's implemented by calling the function GetVols. Flat 
     volatility can have its own implementation

     @param dTime the time at which the values will be obtained
     @param pdS an ascending array of spots
     @param pdVolsSquared the (computed) values of valitity
     @param nNbS the number of spots
     
     @noexport
   */
  virtual void GetVolsSquared(double dTime,
                              const double *pdS,
                              double *pdVolsSquared,
                              size_t nNbS) const;

  /**
     @internal
     @brief Returns a perturbed volatility

     To calculate the derivatives with respect to the volatility by finite 
     difference, we must be able to shift it (volaitility surface) by some 
     (small) amount. This method must be implemented by the derived classes
     to do this.

     Note that it doesn't modify this object at all but rather returns a new,
     perturbed, copy of it.

     @param dShift the perturbation parameter, its meaning varies for
                   different derived classes

     @return a new perturbed volatility

     @noexport
   */
  virtual shared_ptr<Volatility> Perturb(double dShift) const;

  /**
     @internal
     @brief Gets the special times in the volatility surface

     Used by the mesh manager to force these special points in the
     mesh.  The mesh is also refined around these points.
     Helps with the convergence when the volatility is non-smooth
     or discontinuous.

     By default, no special time is considered

     @noexport
   */
  virtual void GetSpecialTimes(numeric::mesh::SpecialTimes& specialTimes) const;

  /**
      @internal
      @brief Dump all data of this volatility in XML format.


      Note that this method does @b  dump the contents of the associated
      volatility.

      @param tagParent the parent tag under which our tag(s) should be created
      
      @noexport

   */
  virtual void Dump(ito33::XML::Tag& tagParent) const = 0;


  /**
      @internal
      @brief This method is part of the implementation of the visitor pattern.

      Visitor pattern allows to do different things for different kinds of
      volatility without adding virtual functions to do all of them in the
      base volatility class itself nor doing the dreaded switchs on typeid of
      the concrete type. Instead, define a new class deriving from
      VolatilityVisitor and do whatever is required in its methods.

      @noexport
   */
  virtual void Visit(VolatilityVisitor& visitor) const = 0;

  /**
     @internal
     @brief Return whether or not the volatility is time only.

     @return true if time only, false otherwise

     @noexport
   */
  virtual bool IsTimeOnly() const
  {
    return false;
  }

  /**
     @internal
     @brief Serializes the parameters used

     @noexport
   */ 
  virtual void 
    GetModelParameters(finance::ModelParametersConsumer& visitor) const;

}; // class Volatility


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_VOLATILITY_H_
