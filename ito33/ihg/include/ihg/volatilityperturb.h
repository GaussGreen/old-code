/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/volatilityperturb.h
// Purpose:     Helper class for volatility perturbation
// Created:     2005/02/10
// RCS-ID:      $Id: volatilityperturb.h,v 1.7 2005/12/02 15:53:54 zhang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/volatilityperturb.h
 */

#ifndef _IHG_VOLATILITYPERTURB_H_
#define _IHG_VOLATILITYPERTURB_H_

#include "ito33/sharedptr.h"
#include "ito33/debug.h"

#include "ito33/ihg/volatility.h"

namespace ito33
{

namespace numeric
{
  namespace mesh
  {
    class SpecialTimes;
  }
}

namespace ihg
{


/**
    Base class for all kinds of volatilities (flat, parametrized, etc)

    @nocreate
 */
class VolatilityPerturb : public Volatility
{
public:
  
  /**
     Ctor takes a volatility reference and a shift to the volatility. 

     It's obvious that the life of the volatility should be long enough.

     @param volatility the volatility to be perturbed
     @param dShift the shift to the volatility
   */
  VolatilityPerturb(const Volatility& volatility, double dShift)
                  : m_volatility(volatility), m_dShift(dShift)
  {
  }
  
  // Default dtor is ok 

  void 
  GetVols(double dTime, const double *pdS, double *pdVols, size_t nNbS) const
  {
    m_volatility.GetVols(dTime, pdS, pdVols, nNbS);

    for (size_t nIdxS = 0; nIdxS < nNbS; nIdxS++)
      pdVols[nIdxS] += m_dShift;
  }

  void GetSpecialTimes(numeric::mesh::SpecialTimes& specialTimes) const
  {
    m_volatility.GetSpecialTimes(specialTimes);
  }

  void Dump(ito33::XML::Tag& /* tagParent */) const 
  { 
    FAIL("no implementation, something goes really wrong here");
  }

  void Visit(VolatilityVisitor& /* visitor */) const
  { 
    FAIL("no implementation, something goes really wrong here");
  }

  bool IsTimeOnly() const
  {
    return m_volatility.IsTimeOnly();
  }

  virtual void 
    GetModelParameters(finance::ModelParametersConsumer& /*visitor*/) const
  {
    FAIL("no implementation, something goes really wrong here");
  }


private:

  const Volatility& m_volatility;

  double m_dShift;

  NO_COPY_CLASS(VolatilityPerturb);

}; // class VolatilityPerturb


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_VOLATILITYPERTURB_H_

