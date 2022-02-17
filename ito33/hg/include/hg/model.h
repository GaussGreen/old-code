/////////////////////////////////////////////////////////////////////////////
// Name:        hg/model.h
// Purpose:     HG pricing model class
// Created:     2005/01/13
// RCS-ID:      $Id: model.h,v 1.9 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/model.h
   @brief HG pricing model class
 */

#ifndef _HG_MODEL_H_
#define _HG_MODEL_H_

#include "ito33/pricing/model.h"

#include "ito33/hg/jumps.h"
#include "ito33/hg/underlyingprocess.h"

namespace ito33
{

namespace hg
{


/// HG pricing model class. 
class Model : public pricing::Model
{

public:

  /**
     Ctor copies from the UnderlyingProcess.

     @param underlyingProcess The underlying process of the model
     @param pUnderlyingProcessForMesh The process for mesh construction
   */
  Model(const UnderlyingProcess& underlyingProcess,
        const shared_ptr<UnderlyingProcess>& pUnderlyingProcessForMesh);

  // Default dtor is ok

  /// Set the underlying process for mesh construction
  void SetUnderlyingProcessForMesh
    (const shared_ptr<UnderlyingProcess>& pUnderlyingProcessForMesh)
  {
    ASSERT_MSG( pUnderlyingProcessForMesh, 
                "Invalid process for mesh in Model class" );
    
    m_pUnderlyingProcessForMesh = pUnderlyingProcessForMesh;
  }
  
  /**
     Gets the number of parameters determing the underlying process.

     @return The number of parameters determing the underlying process

     @noexport
   */
  size_t GetNbParameters() const
  {
    return m_underlyingProcess.GetNbParameters();
  }  
  
  /**
     The volatilities of each non default regime.

     @return The volatilities of each non default regime
   */
  const std::vector<double>& GetVolatilities() const
  {
    return m_underlyingProcess.GetVolatilities();
  }

  /**
     The intensities of jumps from non default regimes to the default regime.

     @return The default jump intensities
   */
  const std::vector<double>& GetJumpsToDefault() const
  { 
    return m_underlyingProcess.GetJumpsToDefault();
  }

    /**
     The jumps from one non default regime to another non default regime.
     
     @param nRegimeFrom The Regime the jumps jump from
     @param nRegimeTo The Regime the jumps jump to

     @return all the jumps from one non default regime to another non default
   */
  const Jumps& GetJumps(size_t nRegimeFrom, size_t nRegimeTo) const
  {
    return m_underlyingProcess.GetJumps(nRegimeFrom, nRegimeTo);
  }

  /**
     The number of regimes.

     @return The number of regimes
   */
  size_t GetNbRegimes() const
  {
    return m_underlyingProcess.GetNbRegimes();
  }

  /**
     Computes the total volatility at each non default regime.

     @return A vector contains the total vol at each non default regime
   */
  std::vector<double> ComputeTotalVolatilities() const
  {
    return m_underlyingProcess.ComputeTotalVolatilities();
  }

  /**
     Computes the regime transition proba.

     @param dT The time period
     @return The matrix for the regime transition probas
   */
  numeric::DenseMatrix ComputeRegimeTransitionProba(double dT) const;

  /**
     Computes the spot regime transition proba.

     @param dT The time period
     @return The matrix for the spot regime transition probas
   */
  numeric::DenseMatrix ComputeSpotRegimeTransitionProba(double dT) const;

  /// Get the convection size
  double GetConvection(double dMaturity, double dSpot) const;

  /// Set the Sharpe ratio
  void SetSharpeRatio(double dSharpeRatio)
  {
    m_dSharpeRatio = dSharpeRatio;
  }

  /// Get the Sharpe ratio
  double GetSharpeRatio() const
  {
    return m_dSharpeRatio;
  }


protected:

  /// Get the square of the total vol
  double 
    GetSquaredTotalVol(double dMaturity, double dSpot, bool bForMesh) const;

  /// The Sharpe ratio
  double m_dSharpeRatio;

  /// Reference to the underlying process
  const UnderlyingProcess& m_underlyingProcess;
  
  /// The underlying process for mesh construction
  shared_ptr<UnderlyingProcess> m_pUnderlyingProcessForMesh;


private:
  
  NO_COPY_CLASS(Model);

}; // class Model


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_MODEL_H_
