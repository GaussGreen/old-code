/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/hg/underlyingprocess.h
// Purpose:     homogeneous underlying process class
// Created:     2005/04/15
// RCS-ID:      $Id: underlyingprocess.h,v 1.15 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/hg/underlyingprocess.h
    @brief class for homogeneous underlying process
 */

#ifndef _ITO33_HG_UNDERLYINGPROCESS_H_
#define _ITO33_HG_UNDERLYINGPROCESS_H_

#include "ito33/vector.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/underlyingprocess.h"

#include "ito33/hg/jumps.h"

namespace ito33
{

namespace XML
{
  class Tag;
}

namespace numeric
{
  class DenseMatrix;
}

namespace hg
{

  /// Maximum of regimes we will handle
  const size_t NMAXREGIMES = 3;

/**
   Class describes a homogeneous underlying process.
 */
class ITO33_HG_DLLDECL UnderlyingProcess : public finance::UnderlyingProcess
{
  
public:

  /**
     ctor defines the UnderlyingProcess with necessary parameters.
     
     @param nNbRegimes Number of non default regimes
     @param pdVols The volatilities at non default regimes
     @param pdIntensities The intensities of jumps from non default regimes
                          to default regime
   */
  UnderlyingProcess(size_t nNbRegimes, 
                    const std::vector<double>& pdVols,
                    const std::vector<double>& pdIntensities);

  /**
      @name Initialization functions.
   */
  //@{
  
  /**
     The jumps from one non default regime to another non default regime.

     @param nRegimeFrom The Regime the jumps jump from
     @param nRegimeTo The Regime the jumps jump to
     @param jumps The jumps from regime nRegimeFrom to regime nRegimeTo
   */
  void SetJumps(size_t nRegimeFrom, size_t nRegimeTo, const Jumps& jumps);

  /**
     The jumps from one non default regime to another non default regime.

     @param nRegimeFrom The Regime the jumps jump from
     @param nRegimeTo The Regime the jumps jump to
     @param pdIntensities The intensities of the jumps
     @param pdAmplitudes The amplitudes of the jumps

     @noexport
   */
  void SetJumps(size_t nRegimeFrom, size_t nRegimeTo, 
                const std::vector<double>& pdIntensities,
                const std::vector<double>& pdAmplitudes);


  /**
     The intensities of jumps from non default regimes to the default regime.

     @param pdIntensities The intensities of jumps from non default regimes
                          to default regime
   */
  void SetJumpsToDefault(const std::vector<double>& pdIntensities);

  /**
     The volatilities of each non default regime.

     @param pdVols The volatilities of each non default regime
   */
  void SetVolatilities(const std::vector<double>& pdVols);

  /// @noexport
  void SetJumpToDefault(size_t nIdxR, double dIntensity);
  
  /// @noexport
  void SetVolatility(size_t nIdxR, double dVol);
 
  //@} // name Initialization functions

  /**
     The volatilities of each non default regime.

     @return The volatilities of each non default regime
   */
  const std::vector<double>& GetVolatilities() const
  {
    return m_pdVols;
  }

  /**
     The intensities of jumps from non default regimes to the default regime.

     @return The default jump intensities
   */
  const std::vector<double>& GetJumpsToDefault() const
  { 
    return m_pdDefaultIntensities;
  }

  /**
     The jumps from one non default regime to another non default regime.
     
     @param nRegimeFrom The Regime the jumps jump from
     @param nRegimeTo The Regime the jumps jump to

     @return all the jumps from one non default regime to another non default
   */
  const Jumps& GetJumps(size_t nRegimeFrom, size_t nRegimeTo) const
  {
    return m_ppJumps[nRegimeFrom][nRegimeTo];
  }

  /**
     The number of regimes.

     @return The number of regimes
   */
  size_t GetNbRegimes() const
  {
    return m_nNbRegimes;
  }

  /**
     Computes the total volatility at each non default regime.

     @return A vector contains the total vol at each non default regime
   */
  std::vector<double> ComputeTotalVolatilities() const;

  /**
     Computes the real(historic) underlying process corresponding to the
     risk neutral underlying process.

     @param dSharpeRatio The sharpe ratio of the real underlying process
     @return The real underlying process

     @noexport
   */
  shared_ptr<UnderlyingProcess> 
  ComputeUnderlyingProcess(double dSharpeRatio) const;

  /**
     @internal
     @brief Computes the regime transition probabilities.

     @return The probability of switching from one regime to another regime

     @return A dense matrix contains the regime transition proba, element
             (i, j) represents the proba to jump from regime i to regime
             j( i, j non default regimes. when j = NbRegimes, it is the
             the probability to default from regime i.

     @noexport
   */
  numeric::DenseMatrix ComputeRegimeTransitionProba(double dT) const;

  /**
     @internal
     @brief Computes the spot regime transition probabilities.

     @return The probability of switching from one regime to another regime

     @return A dense matrix contains the transition proba, element
             (i, j) represents the proba to jump from regime i to regime
             j( i, j non default regimes. when j = NbRegimes, it is the
             the probability to default from regime i.

     @noexport
   */
  numeric::DenseMatrix ComputeSpotRegimeTransitionProba(double dT) const;

  /**
     Gets the number of parameters determing the underlying process

     @return The number of parameters determing the underlying process

     @noexport
   */
  size_t GetNbParameters() const;

  /**
     @internal

     @noexport
   */
  void Dump(ito33::XML::Tag& tagParent) const;


protected:

  /// The number of non default regimes
  size_t m_nNbRegimes;

  /// The volatility values for each non default regime
  std::vector<double> m_pdVols;

  /// The intensities of the jumps to the default regime
  std::vector<double> m_pdDefaultIntensities;

  /// The jumps between non default regimes
  Jumps m_ppJumps[NMAXREGIMES][NMAXREGIMES];


private:

  /// Check the regime index
  void CheckRegimeIndex(size_t nIdxR) const;
  
  /// Check the regime size
  void CheckRegimeSize(size_t nNbRegimes) const;

  /// Check the jump intensity
  void CheckIntensity(double dIntensity) const;

  /// Check the jump amplitude
  void CheckAmplitude(double dAmplitude) const;
  
}; // class UnderlyingProcess


} // namespace hg

} // namespace ito33

#endif // _ITO33_HG_UNDERLYINGPROCESS_H_
