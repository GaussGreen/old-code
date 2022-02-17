/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/model.h
// Purpose:     base model class
// Author:      Wang
// Created:     2004/02/11
// RCS-ID:      $Id: model.h,v 1.15 2006/05/19 18:49:19 yann Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/model.h
    @brief Base model class 
 */

#ifndef _ITO33_PRICING_MODEL_H_
#define _ITO33_PRICING_MODEL_H_

#include "ito33/numeric/mesh/specialtimes.h"

namespace ito33
{

namespace pricing
{


/// Base class for different models(HG and IHG)
class Model
{

public:
 
  /// Default ctor
  Model(): m_dPostDefaultVolatility(0.0) {}

  /// Dummy virtual destructor
  virtual ~Model() { }

  /**
     Sets the initial state of the model
     
     @param dTime the intial time
   */
  virtual void SetInitialState(double dTime) { m_dTime = dTime; }

  /** 
     Updates the state of the model
     
     @param dTime the current time
   */
  void Update(double dTime) { m_dTime = dTime; }

  /**
     Gets an estimation of the squared total vol

     @param dMaturity the maturity at which the squared total vol is estimated
     @param dSpot the spot at which the squared total vol is estimated

     @return the estimated squared total vol
   */
  double GetSquaredTotalVol(double dMaturity, double dSpot) const
  {
    return GetSquaredTotalVol(dMaturity, dSpot, false);
  }

  /**
     Gets an estimation of the squared total vol using the volatility for mesh.

     @param dMaturity the maturity at which the squared total vol is estimated
     @param dSpot the spot at which the squared total vol is estimated

     @return the estimated squared total vol for mesh
   */
  double GetSquaredTotalVolForMesh(double dMaturity, double dSpot) const
  {
    return GetSquaredTotalVol(dMaturity, dSpot, true);
  }

  /**
     Gets an estimation of the convection size

     @param dMaturity the maturity at which the convection size is estimated
     @param dSpot the spot at which the convection size is estimated

     @return the estimated convection size
   */
  virtual double GetConvection(double dMaturity, 
                                    double dSpot) const = 0;

  /**
     Gets the special times in the model.

     Used by the mesh manager to force these special points in the
     mesh.  The mesh is also refined around these points.
     Helps with the convergence when the model parameter is non-smooth
     or discontinuous.

     By default, no special time is considered.

     @param specialTimes the container that the special times will be filled to.
   */
  virtual void GetSpecialTimes(numeric::mesh::SpecialTimes& specialTimes) const
  {
    specialTimes.clear();
  }
   
  /**
     The post default volatility.

     @return the post default volatility
  */
  double GetPostDefaultVolatility() const
  {
    return m_dPostDefaultVolatility;
  }

  /**
     The post default volatility.

     @param dPostDefaultVolatility post default volatility
  */
  void SetPostDefaultVolatility(double dPostDefaultVolatility)
  {
    m_dPostDefaultVolatility = dPostDefaultVolatility;
  }

protected:
   
  /** 
     Gets the squared total vol with the appropriate vol.
     
     For the moment, hazard rate is ignored.

     @param dMaturity the maturity time
     @param dSpot the spot
     @param bVolForMesh true to use the vol for mesh, 
            false to use the actual vol.
      
     @return the squared total vol
   */
  virtual double 
  GetSquaredTotalVol(double dMaturity, double dSpot, bool bVolForMesh) const 
    = 0;

  /// the current time which indicates the state of the model
  double m_dTime;

  /// The post default volatility
  double m_dPostDefaultVolatility;

}; // class Model


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_MODEL_H_
