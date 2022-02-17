/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/model.h
// Purpose:     ihg model class. It supports second hr of the derivative
// Created:     2004/02/11
// RCS-ID:      $Id: model.h,v 1.42 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/model.h
    @brief model class for ihg project
 */

#ifndef _IHG_MODEL_H_
#define _IHG_MODEL_H_

#include "ito33/sharedptr.h"

#include "ito33/pricing/model.h"

#include "ito33/ihg/volatility.h"
#include "ito33/ihg/hazardrate.h"

namespace ito33
{

namespace ihg
{


/// Base params class for IHG model. It supports second hr of the derivative.
class Model : public pricing::Model
{

public:

  /**
      Ctor sets volatility, hazard rate. By default, we don't have 
      the hr of the derivative.
     
      @param pVolatility the volatility to be set for the model
      @param pVolatilityForMesh the volatility to be used for the mesh
      @param pHazardRate the hazard rate to be set for the model
   */
  Model(const shared_ptr<Volatility>& pVolatility,
        const shared_ptr<Volatility>& pVolatilityForMesh,
        const shared_ptr<HazardRate>& pHazardRate)
      : m_pVolatility(pVolatility), m_pVolatilityForMesh(pVolatilityForMesh),
        m_pHazardRate(pHazardRate)
  {
  }

  // Default dtor is ok

  /**
      Sets the volatility of the model.

      @param pVolatility a shared pointer to volatility
   */
  void SetVolatility(const shared_ptr<Volatility>& pVolatility)
  {
    m_pVolatility = pVolatility;
  }

  /**
      Sets the hazard rate of the model.

      @param pHazardRate a shared pointer to hazard rate
   */
  void SetHazardRate(const shared_ptr<HazardRate>& pHazardRate)
  {
    m_pHazardRate = pHazardRate;
  }
  
  /**
      Sets the hazard rate of the derivative.

      @param pHazardRate a shared pointer to hazard rate

      REQUIRE: the pointer must not be null and the HR should be time only
   */
  void SetHazardRateOfDerivative(const shared_ptr<HazardRate>& pHazardRate)
  {
    ASSERT_MSG(pHazardRate && pHazardRate->IsTimeOnly(),
      "Invalid hazard rate of the derivative: null pointer or not time only.");

    m_pHazardRateOfDerivative = pHazardRate;
  }

  /**
      Gets the volatilities for the pdS array.

      @param pdS an ascending  array of spots
      @param pdVols the array of vol that is returned
      @param nNbS number of spot points
   */
  void GetVols(const double *pdS, double *pdVols, size_t nNbS) const
  {
    m_pVolatility->GetVols(m_dTime, pdS, pdVols, nNbS);
  }

  /**
      Gets the squared volatilities for the pdS array.

      @param pdS an ascending  array of spots
      @param pdVolsSqr returns array of square of volatilities 
      @param nNbS number of spot points
   */
  void GetVolsSquared(const double *pdS, double *pdVolsSqr, size_t nNbS) const
  {
    m_pVolatility->GetVolsSquared(m_dTime, pdS, pdVolsSqr, nNbS);
  }

  /**
      Gets the values of hazard rates on the spot grid.

      @param pdS an ascending  array of spots
      @param pdhrs array of hazard rates 
      @param nNbS number of spot points
   */
  void GetHazardRates(const double *pdS, double *pdhrs, size_t nNbS) const
  {
    m_pHazardRate->GetHazardRates(m_dTime, pdS, pdhrs, nNbS);
  }

  /**
      Checks if the hazard rate is time only.

      @return true if the hazard rate is time only, false otherwise
   */
  bool IsHazardRateTimeOnly() const
  {
    return m_pHazardRate->IsTimeOnly();
  }

  /**
      Checks if volatility is time only.

      @return true if time only, false otherwise
   */
  bool IsVolatilityTimeOnly() const
  {
    return m_pVolatility->IsTimeOnly();
  }

  /**
      Gets the value of hazard rate of derivative at current time step.

      @return value of hazard rate of derivative at this time step
   */
  double GetHROfDerivative() const
  {
    ASSERT_MSG(m_pHazardRateOfDerivative,
      "The second harzard rate is not defined.");

    double
      dS = 0,
      dR;

    m_pHazardRateOfDerivative->GetHazardRates(m_dTime, &dS, &dR, 1);

    return dR;
  }

  void SetInitialState(double dTime)
  {
    Model::BaseClass::SetInitialState(dTime);
  }
  
  /// See base class
  double GetConvection(double dMaturity, double dSpot) const
  {
    double dHR;
    
    m_pHazardRate->GetHazardRates(dMaturity, &dSpot, &dHR, 1);
    
    return GetSquaredTotalVolForMesh(dMaturity, dSpot) - dHR * 2;
  }

  /**
      Gets the special times in the ihg model.

      Used by the mesh manager to force these special points in the
      mesh.  The mesh is also refined around these points.
      Helps with the convergence when the volatility or hazard rate
      is non-smooth or discontinuous.

      @param specialTimes the container that the special times will be filled to
   */
  void GetSpecialTimes(numeric::mesh::SpecialTimes& specialTimes) const
  {
    // Get the special times from volatility
    m_pVolatility->GetSpecialTimes(specialTimes);

    numeric::mesh::SpecialTimes specialTimesTmp;

    // Get the special times from hazard rate
    m_pHazardRate->GetSpecialTimes(specialTimesTmp);

    // join the two lists
    specialTimes.insert(specialTimes.end(), 
                        specialTimesTmp.begin(), specialTimesTmp.end() );

    if (m_pHazardRateOfDerivative)
    {
      // Get the special times from second hazard rate
      m_pHazardRateOfDerivative->GetSpecialTimes(specialTimesTmp);

      // join the two lists
      specialTimes.insert(specialTimes.end(), 
                          specialTimesTmp.begin(), specialTimesTmp.end() );
    }
  }


protected:
  
  /// See base class
  double 
  GetSquaredTotalVol(double dMaturity, double dSpot, bool bForMesh) const
  {
    Volatility* pVolatility;
     
    if( bForMesh )
    {
      ASSERT_MSG( m_pVolatilityForMesh, "Volatility for mesh not defined.");

      pVolatility = m_pVolatilityForMesh.get();
    }
    else
      pVolatility = m_pVolatility.get();

    double dVol;
    
    pVolatility->GetVols(dMaturity, &dSpot, &dVol, 1);

    if (dVol < 0.1)
      dVol = 0.1;
    
    return dVol * dVol;
  }

  /// @name IHG specific parameters.
  //@{

  /// The volatility for mesh 
  shared_ptr<Volatility> m_pVolatilityForMesh;

  /// The volatility of the underlying 
  shared_ptr<Volatility> m_pVolatility;

  /// The hazard rate (default intensity) of the underlying
  shared_ptr<HazardRate> m_pHazardRate;

  /// The hazard rate (default intensity) of the derivative
  shared_ptr<HazardRate> m_pHazardRateOfDerivative;

  //@}

private:
  
  typedef pricing::Model BaseClass;

}; // class Model


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_MODEL_H_
