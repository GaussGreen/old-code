/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/modeloutput.h
// Purpose:     Base output class
// Created:     Dec 8, 2003
// RCS-ID:      $Id: modeloutput.h,v 1.34 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/modeloutput.h
    @brief base output class
 */

#ifndef _ITO33_FINANCE_MODELOUTPUT_H_
#define _ITO33_FINANCE_MODELOUTPUT_H_

#include "ito33/sharedptr.h"

#include "ito33/finance/surfacedouble.h"
#include "ito33/finance/surfaceflag.h"

namespace ito33
{

namespace XML { class Tag; }

namespace finance
{

typedef std::vector<double> Values;
typedef shared_ptr<finance::SurfaceDouble> SharedSurface;
#ifdef __CPP2ANY__
typedef std::vector<int> Flags;
#else
typedef finance::SurfaceFlag::Flags Flags;
#endif

class NumOutput;

/**
   Base output class for all computations. It contains common results 
   for all instruments.

   @nocreate
 */
class ITO33_DLLDECL ModelOutput
{
public:
  
  ModelOutput() : m_dPrice(0), m_dValueAfterDefault(0),
                  m_dDelta(0.), m_dGamma(0.), m_dFXDelta(0.), m_dTheta(0.),
                  m_dVega(0.), m_dRho(0.), m_dUnderlyingRho(0.), m_dFugit(0.),
                  m_bHasDataAtAnalysisDate(false),
                  m_bHasSurface(false),
                  m_bHasFXDelta(false),
                  m_bHasTheta(false),
                  m_bHasVega(false),
                  m_bHasRho(false),
                  m_bHasUnderlyingRho(false),
                  m_bHasFugit(false),
                  m_bHasThetaAtAnalysisDate(false),
                  m_bHasVegaAtAnalysisDate(false),
                  m_bHasRhoAtAnalysisDate(false),
                  m_bHasUnderlyingRhoAtAnalysisDate(false),
                  m_bHasFugitAtAnalysisDate(false),
                  m_bHasThetaSurface(false),
                  m_bHasVegaSurface(false),
                  m_bHasRhoSurface(false),
                  m_bHasUnderlyingRhoSurface(false),
                  m_bHasFugitSurface(false)
  { 
  }

  /// Dummy virtual dtor for base class
  virtual ~ModelOutput() { }
  
  /// @name Modifiers for scalar values
  //@{
  
  /**
     @internal
     @brief Sets price value.

     @param dPrice the price value
     
     @noexport
   */
  void SetPrice(double dPrice)
  {
    m_dPrice = dPrice;
  }
  
  /**
     @internal
     @brief Sets the value after default of the instrument.

     @param dValueAfterDefault the value after default of the instrument
     
     @noexport
   */
  void SetValueAfterDefault(double dValueAfterDefault)
  {
    m_dValueAfterDefault = dValueAfterDefault;
  }

  /**
     @internal

     @noexport
   */
  void SetDelta(double dDelta) { m_dDelta = dDelta; }

  /**
     @internal

     @noexport
   */
  void SetGamma(double dGamma) { m_dGamma = dGamma; }

  /**
     @internal

     @noexport
   */
  void SetFXDelta(double dFXDelta) 
  { 
    m_bHasFXDelta = true;    
    m_dFXDelta = dFXDelta; 
  }
  
  /**
     @internal

     @noexport
   */
  void SetTheta(double dTheta) { m_bHasTheta = true; m_dTheta = dTheta; }

  /**
     @internal

     @noexport
   */
  void SetVega(double dVega) { m_bHasVega = true; m_dVega = dVega; }

  /**
     @internal

     @noexport
   */
  void SetRho(double dRho) { m_bHasRho = true; m_dRho = dRho; }

  /**
     @internal

     @noexport
   */
  void SetUnderlyingRho(double dUnderlyingRho) 
  { 
    m_bHasUnderlyingRho = true; 
    m_dUnderlyingRho = dUnderlyingRho; 
  }

  /**
     @internal

     @noexport
   */
  void SetFugit(double dFugit) { m_bHasFugit = true; m_dFugit = dFugit; }

  //@}
 
  /// @name Accessors for scalar values
  //@{

  /**
     Gets the price value.

     @return the price value
   */
  double GetPrice() const
  {
    return m_dPrice;
  }

  /**
     Gets the value after default of the instrument.

     @return the value after default of the instrument
   */
  double GetValueAfterDefault() const
  {
    return m_dValueAfterDefault;
  }
  
  double GetDelta() const { return m_dDelta; }

  double GetGamma() const { return m_dGamma; }

  double GetFXDelta() const;
  
  double GetTheta() const;

  double GetVega() const;

  double GetRho() const;

  double GetUnderlyingRho() const;

  double GetFugit() const;

  //@}

  /// @name Modifiers and Accessors for numerical output
  //@{
  
  /**
     Sets the numerical output.

     @internal

     @noexport
   */
  void SetNumOutput(const shared_ptr<NumOutput>& pNumOutput)
  {
    m_pNumOutput = pNumOutput;
  }

  /**
     @internal

     @noexport
   */
  const shared_ptr<NumOutput>& GetNumOutput() const
  {
    return m_pNumOutput;
  }

  //@}

  /**
     @name Modifiers for values at the analysis date

   */
  //@{

  /**
     @internal

     @noexport
   */
  void SetThetasAtAnalysisDate(const Values& pdThetas);
     
  /**
     @internal

     @noexport
   */
  void SetVegasAtAnalysisDate(const Values& pdVegas);
  
  /**
     @internal

     @noexport
   */
  void SetRhosAtAnalysisDate(const Values& pdRhos);

  /**
     @internal

     @noexport
   */
  void SetUnderlyingRhosAtAnalysisDate(const Values &pdUnderlyingRhos);

  /**
     @internal

     @noexport
   */
  void SetFugitsAtAnalysisDate(const Values& pdFugits);

  //@}

  /**
     @name Modifiers for values at the analysis date

     If SetSpotsAtAnalysisDate() is called, so should be the setters for 
     price, delta, gamma.
   */
  //@{

  /**
     @internal

     @noexport
   */
  void SetSpotsAtAnalysisDate(const Values& pdSpots);

  /**
     @internal

     @noexport
   */
  void SetPricesAtAnalysisDate(const Values& pdPrices);

  /**
     @internal

     @noexport
   */
  void SetDeltasAtAnalysisDate(const Values& pdDeltas);

  /**
     @internal

     @noexport
   */
  void SetGammasAtAnalysisDate(const Values& pdGammas);

  //@}


  /**
     @name Accessors for values at analysis date

     Price, Delta, Gamma, Theta share the same check.
   */
  //@{

  const Values& GetSpotsAtAnalysisDate() const;

  const Values& GetPricesAtAnalysisDate() const;

  const Values& GetDeltasAtAnalysisDate() const;

  const Values& GetGammasAtAnalysisDate() const;

  const Values& GetThetasAtAnalysisDate() const;

  const Values& GetVegasAtAnalysisDate() const;

  const Values& GetRhosAtAnalysisDate() const;

  const Values& GetUnderlyingRhosAtAnalysisDate() const;

  const Values& GetFugitsAtAnalysisDate() const;

  //@}


  /**
     @name Modifiers for surfaces

     If SetDomain() is called, so should be the setters for price,
     delta, gamma.
   */
  //@{

  /**
     @internal

     @noexport
   */
  void SetDomain(const shared_ptr<finance::Domain>& pDomain);

  /**
     @internal

     @noexport
   */
  void SetPriceSurface(const SharedSurface& pPriceSurface);

  /**
     @internal

     @noexport
   */
  void SetDeltaSurface(const SharedSurface& pDeltaSurface);

  /**
     @internal

     @noexport
   */
  void SetGammaSurface(const SharedSurface& pGammaSurface);

  /**
     @internal

     @noexport
   */
  void SetThetaSurface(const SharedSurface& pThetaSurface);
  
  /**
     @internal

     @noexport
   */
  void SetVegaSurface(const SharedSurface& pVegaSurface);

  /**
     @internal

     @noexport
   */
  void SetRhoSurface(const SharedSurface& pRhoSurface);

  /**
     @internal

     @noexport
   */
  void SetUnderlyingRhoSurface(const SharedSurface & pUnderlyingRhoSurface);

  /**
     @internal

     @noexport
   */
  void SetFugitSurface(const SharedSurface& pFugitSurface);

  //@}



  /// @name Accessors for surfaces
  //@{

  shared_ptr<finance::Domain> GetDomain() const;

  SharedSurface GetPriceSurface() const;

  SharedSurface GetDeltaSurface() const;

  SharedSurface GetGammaSurface() const;

  SharedSurface GetThetaSurface() const;
  
  SharedSurface GetVegaSurface() const;

  SharedSurface GetRhoSurface() const;

  SharedSurface GetUnderlyingRhoSurface() const;

  SharedSurface GetFugitSurface() const;

  //@}



  /// @name Checkers
  //@{

  bool HasSpotAtAnalysisDate() const { return HasDataAtAnalysisDate(); }

  bool HasPriceAtAnalysisDate() const { return HasDataAtAnalysisDate(); }

  bool HasDeltaAtAnalysisDate() const { return HasDataAtAnalysisDate(); }

  bool HasGammaAtAnalysisDate() const { return HasDataAtAnalysisDate(); }

  bool HasDomain() const { return HasSurface(); }

  bool HasPriceSurface() const { return HasSurface(); }

  bool HasDeltaSurface() const { return HasSurface(); }

  bool HasGammaSurface() const { return HasSurface(); }
  
  bool HasFXDelta() const { return m_bHasFXDelta; } 
  
  bool HasTheta() const { return m_bHasTheta; }
  
  bool HasVega() const { return m_bHasVega; }
  
  bool HasRho() const { return m_bHasRho; }

  bool HasUnderlyingRho() const { return m_bHasUnderlyingRho; }

  bool HasFugit() const { return m_bHasFugit; }

  bool HasThetaAtAnalysisDate() const { return m_bHasThetaAtAnalysisDate; }

  bool HasVegaAtAnalysisDate() const { return m_bHasVegaAtAnalysisDate; }

  bool HasRhoAtAnalysisDate() const { return m_bHasRhoAtAnalysisDate; }

  bool HasUnderlyingRhoAtAnalysisDate() const 
  { 
    return m_bHasUnderlyingRhoAtAnalysisDate; 
  }

  bool HasFugitAtAnalysisDate() const { return m_bHasFugitAtAnalysisDate; }

  bool HasThetaSurface() const { return m_bHasThetaSurface; }

  bool HasVegaSurface() const { return m_bHasVegaSurface; }

  bool HasRhoSurface() const { return m_bHasRhoSurface; }

  bool HasUnderlyingRhoSurface() const { return m_bHasUnderlyingRhoSurface; }

  bool HasFugitSurface() const { return m_bHasFugitSurface; }

  //@}

  /**
     @internal
     @brief Sets all the results concerning the FX delta.

     @param pModelOutputNew The modeloutput containing the results of the 
            pricing with the shifted FX 
     @param dInverseShift The inverse of the shift of the FX

     @noexport
   */
  void SetFXDeltaResults( const shared_ptr<ModelOutput>& pModelOutputNew, 
                          double dInverseShift );

  /**
     @internal
     @brief Sets all the results concerning the rho.

     @param pModelOutputNew The modeloutput containing the results of the 
            pricing with the shifted YC 
     @param dInverseShift The inverse of the shift of the YC

     @noexport
   */
  void SetRhoResults( const shared_ptr<ModelOutput>& pModelOutputNew, 
                      double dInverseShift );

  /**
     @internal
     @brief Sets all the results concerning the underlying rho.

     @param pModelOutputNew The modeloutput containing the results of the 
            pricing with the shifted YC 
     @param dInverseShift The inverse of the shift of the YC

     @noexport
   */
  void SetUnderlyingRhoResults
       (const shared_ptr<ModelOutput>& pModelOutputNew, double dInverseYCShift);

  /**
     @internal
     @brief Sets all the results concerning the vega.

     @param pModelOutputNew The modeloutput containing the results of the 
            pricing with the shifted volatility 
     @param dInverseShift The inverse of the shift of the volatility

     @noexport
   */
  void SetVegaResults( const shared_ptr<ModelOutput>& pModelOutputNew, 
                       double dInverseShift );

  /**
     @internal
     @brief Dump all data stored in this object in XML format.
     This method is usually called by the function doing the pricing,
     calibration &c but can also be called "manually" if needed.
     It will normally be overridden in the derived classes but don't forget to
     call the base class version from there.

     @param tagParent the parent tag under which our tag(s) should be created

     @noexport
   */
  virtual void Dump(XML::Tag& tagParent) const;


protected:

  /// Throw exception is value requested not available.
  static void ThrowValueNotAvailable();

  /// Throw exception if analysis data not available.
  static void ThrowDataAtAnalysisDateNotAvailable();

  /// Throw exception is surface not available.
  static void ThrowSurfaceNotAvailable();

  /// @name Helper for Checking availability
  //@{

  void CheckAtAnalysisDate() const;

  void CheckSurface() const;

  void CheckFXDelta() const;

  void CheckTheta() const;

  void CheckVega() const;
  
  void CheckRho() const;

  void CheckUnderlyingRho() const;

  void CheckUnderlyingRhoAtAnalysisDate() const;

  void CheckFugit() const;

  void CheckThetaAtAnalysisDate() const;

  void CheckVegaAtAnalysisDate() const;

  void CheckRhoAtAnalysisDate() const;

  void CheckUnderlyingRhoSurface() const;

  void CheckFugitAtAnalysisDate() const;

  void CheckThetaSurface() const;

  void CheckVegaSurface() const;

  void CheckRhoSurface() const;

  void CheckFugitSurface() const;

  //@}

  /// Check if there is data at the analysis date.
  bool HasDataAtAnalysisDate() const { return m_bHasDataAtAnalysisDate; }

  /// Check if there is a surface.
  bool HasSurface() const { return m_bHasSurface; }

  // The scalar values 
  /// Price
  double  m_dPrice;

  /// Value after default
  double m_dValueAfterDefault;
  
  /// Greek: delta
  double m_dDelta;

  /// Greek: gamma
  double m_dGamma;
  
  /// Greek: FX delta
  double m_dFXDelta;

  /// Greek: theta
  double m_dTheta;

  /// Greek: vega
  double m_dVega;

  /// Greek: Rho
  double m_dRho;

  /// Underlying Rho
  double m_dUnderlyingRho;

  /// Fugit
  double m_dFugit;

  /// The values at the analysis date
  Values m_pdSpotsAtAnalysisDate;

  /// Prices at analysis date
  Values m_pdPricesAtAnalysisDate;

  /// Deltas at analysis date
  Values m_pdDeltasAtAnalysisDate;

  /// Gammas at analysis date
  Values m_pdGammasAtAnalysisDate;

  /// Thetas at analysis date
  Values m_pdThetasAtAnalysisDate;

  /// Vegas at analysis date
  Values m_pdVegasAtAnalysisDate;

  /// Rhos at analysis date
  Values m_pdRhosAtAnalysisDate;

  /// Underlying rhos at analysis date
  Values m_pdUnderlyingRhosAtAnalysisDate;

  /// Fugits at analysis date
  Values m_pdFugitsAtAnalysisDate;

  /// The Surfaces
  shared_ptr<finance::Domain> m_pDomain;

  /// Price surface
  SharedSurface m_pPriceSurface;

  /// Delta surface
  SharedSurface m_pDeltaSurface;

  /// Gamma surface
  SharedSurface m_pGammaSurface;

  /// Theta surface
  SharedSurface m_pThetaSurface;
  
  /// Vega surface
  SharedSurface m_pVegaSurface;

  /// Rho surface
  SharedSurface m_pRhoSurface;

  /// Underlying rho surface
  SharedSurface m_pUnderlyingRhoSurface;

  /// Fugit surface
  SharedSurface m_pFugitSurface;

  // The Flags

  /// Common flags for price, delta, gamma, and theta at analysis date
  bool m_bHasDataAtAnalysisDate;

  /// Common flags for price, delta, gamma, and theta surfaces
  bool m_bHasSurface;  

  /// Common flag for fx delta
  bool m_bHasFXDelta;

  /// Common flag for theta
  bool m_bHasTheta;

  /// Common flag for vega
  bool m_bHasVega;

  /// Common flag for rho
  bool m_bHasRho;

  /// Common flag for underlying rho
  bool m_bHasUnderlyingRho;

  /// Common flag for fugit
  bool m_bHasFugit;

  /// Common flag for theta at analysis date
  bool m_bHasThetaAtAnalysisDate;

  /// Common flag for vega at analysis date
  bool m_bHasVegaAtAnalysisDate;

  /// Common flag for rho at analysis date
  bool m_bHasRhoAtAnalysisDate;

  /// Common flag for underlying rho at analysis date
  bool m_bHasUnderlyingRhoAtAnalysisDate;

  /// Common flag for fugit at anlysis date
  bool m_bHasFugitAtAnalysisDate;

  /// Common flag for theta surface
  bool m_bHasThetaSurface;

  /// Common flag for vega surface
  bool m_bHasVegaSurface;

  /// Common flag for rho surface
  bool m_bHasRhoSurface;

  /// Common flag for underlying rho surface
  bool m_bHasUnderlyingRhoSurface;

  /// Common flag for fugit surface
  bool m_bHasFugitSurface;

  /// Shared pointer to numerical output
  shared_ptr<NumOutput> m_pNumOutput;
};  // class ModelOutput


// Helpers for Checking availability

inline void ModelOutput::CheckAtAnalysisDate() const
{
  if (!m_bHasDataAtAnalysisDate)
    ThrowDataAtAnalysisDateNotAvailable();
}

inline void ModelOutput::CheckSurface() const
{
  if (!m_bHasSurface)
    ThrowSurfaceNotAvailable();
}

// Modifiers for values at analysis date
inline void ModelOutput::SetSpotsAtAnalysisDate(const Values& pdSpots)
{
  m_bHasDataAtAnalysisDate = true;

  m_pdSpotsAtAnalysisDate = pdSpots;
}

inline void ModelOutput::SetPricesAtAnalysisDate(const Values& pdPrices)
{
  m_pdPricesAtAnalysisDate = pdPrices;
}

inline void ModelOutput::SetDeltasAtAnalysisDate(const Values& pdDeltas)
{
  m_pdDeltasAtAnalysisDate = pdDeltas;
}

inline void ModelOutput::SetGammasAtAnalysisDate(const Values& pdGammas)
{
  m_pdGammasAtAnalysisDate = pdGammas;
}

// Accessors for values at analysis date
inline const Values& ModelOutput::GetSpotsAtAnalysisDate() const
{
  CheckAtAnalysisDate();

  return m_pdSpotsAtAnalysisDate;
}

inline const Values& ModelOutput::GetPricesAtAnalysisDate() const
{
  CheckAtAnalysisDate();

  return m_pdPricesAtAnalysisDate;
}

inline const Values& ModelOutput::GetDeltasAtAnalysisDate() const
{
  CheckAtAnalysisDate();

  return m_pdDeltasAtAnalysisDate;
}

inline const Values& ModelOutput::GetGammasAtAnalysisDate() const
{
  CheckAtAnalysisDate();

  return m_pdGammasAtAnalysisDate;
}

// Modifiers for surfaces
inline void ModelOutput::SetDomain(const shared_ptr<finance::Domain>& pDomain)
{
  m_bHasSurface = true;

  m_pDomain = pDomain;
}

inline void ModelOutput::SetPriceSurface(const SharedSurface& pPriceSurface)
{
  m_pPriceSurface = pPriceSurface;
}

inline void ModelOutput::SetDeltaSurface(const SharedSurface& pDeltaSurface)
{
  m_pDeltaSurface = pDeltaSurface;
}

inline void ModelOutput::SetGammaSurface(const SharedSurface& pGammaSurface)
{
  m_pGammaSurface = pGammaSurface;
}

// Accessors for surfaces
inline shared_ptr<Domain> ModelOutput::GetDomain() const
{
  CheckSurface();

  return m_pDomain;
}

inline SharedSurface ModelOutput::GetPriceSurface() const
{
  CheckSurface();

  return m_pPriceSurface;
}

inline SharedSurface ModelOutput::GetDeltaSurface() const
{
  CheckSurface();

  return m_pDeltaSurface;
}

inline SharedSurface ModelOutput::GetGammaSurface() const
{
  CheckSurface();

  return m_pGammaSurface;
}

// Helpers for Checking availability
inline void ModelOutput::CheckFXDelta() const
{
  if (!m_bHasFXDelta)
    ThrowValueNotAvailable();
}

inline void ModelOutput::CheckTheta() const
{
  if (!m_bHasTheta)
    ThrowValueNotAvailable();
}

inline void ModelOutput::CheckVega() const
{
  if (!m_bHasVega)
    ThrowValueNotAvailable();
}

inline void ModelOutput::CheckRho() const
{
  if (!m_bHasRho)
    ThrowValueNotAvailable();
}

inline void ModelOutput::CheckUnderlyingRho() const
{
  if (!m_bHasUnderlyingRho)
    ThrowValueNotAvailable();
}

inline void ModelOutput::CheckFugit() const
{
  if (!m_bHasFugit)
    ThrowValueNotAvailable();
}

inline void ModelOutput::CheckThetaAtAnalysisDate() const
{
  if (!m_bHasThetaAtAnalysisDate)
    ThrowDataAtAnalysisDateNotAvailable();
}

inline void ModelOutput::CheckVegaAtAnalysisDate() const
{
  if (!m_bHasVegaAtAnalysisDate)
    ThrowDataAtAnalysisDateNotAvailable();
}

inline void ModelOutput::CheckRhoAtAnalysisDate() const
{
  if (!m_bHasRhoAtAnalysisDate)
    ThrowDataAtAnalysisDateNotAvailable();
}

inline void ModelOutput::CheckUnderlyingRhoAtAnalysisDate() const
{
  if (!m_bHasUnderlyingRhoAtAnalysisDate)
    ThrowDataAtAnalysisDateNotAvailable();
}

inline void ModelOutput::CheckFugitAtAnalysisDate() const
{
  if (!m_bHasFugitAtAnalysisDate)
    ThrowDataAtAnalysisDateNotAvailable();
}

inline void ModelOutput::CheckThetaSurface() const
{
  if (!m_bHasThetaSurface)
    ThrowSurfaceNotAvailable();
}

inline void ModelOutput::CheckVegaSurface() const
{
  if (!m_bHasVegaSurface)
    ThrowSurfaceNotAvailable();
}

inline void ModelOutput::CheckRhoSurface() const
{
  if (!m_bHasRhoSurface)
    ThrowSurfaceNotAvailable();
}

inline void ModelOutput::CheckUnderlyingRhoSurface() const
{
  if (!m_bHasUnderlyingRhoSurface)
    ThrowSurfaceNotAvailable();
}

inline void ModelOutput::CheckFugitSurface() const
{
  if (!m_bHasFugitSurface)
    ThrowSurfaceNotAvailable();
}

// Accessors for scalar

inline double ModelOutput::GetFXDelta() const
{
  CheckFXDelta();

  return m_dFXDelta;
}

inline double ModelOutput::GetTheta() const
{
  CheckTheta();

  return m_dTheta;
}

inline double ModelOutput::GetVega() const
{
  CheckVega();

  return m_dVega;
}

inline double ModelOutput::GetRho() const
{
  CheckRho();

  return m_dRho;
}

inline double ModelOutput::GetUnderlyingRho() const
{
  CheckUnderlyingRho();

  return m_dUnderlyingRho;
}

inline double ModelOutput::GetFugit() const
{
  CheckFugit();

  return m_dFugit;
}

// Modifiors for values at analysis date

inline void ModelOutput::SetThetasAtAnalysisDate(const Values& pdThetas)
{
  m_bHasThetaAtAnalysisDate = true;

  m_pdThetasAtAnalysisDate = pdThetas;
}

inline void ModelOutput::SetVegasAtAnalysisDate(const Values& pdVegas)
{
  m_bHasVegaAtAnalysisDate = true;

  m_pdVegasAtAnalysisDate = pdVegas;
}

inline void ModelOutput::SetRhosAtAnalysisDate(const Values& pdRhos)
{
  m_bHasRhoAtAnalysisDate = true;

  m_pdRhosAtAnalysisDate = pdRhos;
}

inline void 
ModelOutput::SetUnderlyingRhosAtAnalysisDate(const Values &pdUnderlyingRhos)
{
  m_bHasUnderlyingRhoAtAnalysisDate = true;

  m_pdUnderlyingRhosAtAnalysisDate = pdUnderlyingRhos;
}

inline void ModelOutput::SetFugitsAtAnalysisDate(const Values& pdFugits)
{
  m_bHasFugitAtAnalysisDate = true;

  m_pdFugitsAtAnalysisDate = pdFugits;
}

// Accessors for values at analysis date
inline const Values& ModelOutput::GetThetasAtAnalysisDate() const
{
  CheckThetaAtAnalysisDate();

  return m_pdThetasAtAnalysisDate;
}

inline const Values& ModelOutput::GetVegasAtAnalysisDate() const
{
  CheckVegaAtAnalysisDate();

  return m_pdVegasAtAnalysisDate;
}

inline const Values& ModelOutput::GetRhosAtAnalysisDate() const
{
  CheckRhoAtAnalysisDate();

  return m_pdRhosAtAnalysisDate;
}

inline const Values & ModelOutput::GetUnderlyingRhosAtAnalysisDate() const
{
  CheckUnderlyingRhoAtAnalysisDate();

  return m_pdUnderlyingRhosAtAnalysisDate;
}

inline const Values& ModelOutput::GetFugitsAtAnalysisDate() const
{
  CheckFugitAtAnalysisDate();

  return m_pdFugitsAtAnalysisDate;
}

// Modifiors for surfaces
inline void ModelOutput::SetThetaSurface(const SharedSurface& pThetaSurface)
{
  m_bHasThetaSurface = true;

  m_pThetaSurface = pThetaSurface;
}

inline void ModelOutput::SetVegaSurface(const SharedSurface& pVegaSurface)
{
  m_bHasVegaSurface = true;

  m_pVegaSurface = pVegaSurface;
}

inline void ModelOutput::SetRhoSurface(const SharedSurface& pRhoSurface)
{
  m_bHasRhoSurface = true;

  m_pRhoSurface = pRhoSurface;
}

inline void 
ModelOutput::SetUnderlyingRhoSurface(const SharedSurface & 
                                     pUnderlyingRhoSurface)
{
  m_bHasUnderlyingRhoSurface = true;

  m_pUnderlyingRhoSurface = pUnderlyingRhoSurface;
}

inline void ModelOutput::SetFugitSurface(const SharedSurface& pFugitSurface)
{
  m_bHasFugitSurface = true;

  m_pFugitSurface = pFugitSurface;
}


// Accessors for surfaces
inline SharedSurface ModelOutput::GetThetaSurface() const
{
  CheckThetaSurface();

  return m_pThetaSurface;
}

inline SharedSurface ModelOutput::GetVegaSurface() const
{
  CheckVegaSurface();

  return m_pVegaSurface;
}

inline SharedSurface ModelOutput::GetRhoSurface() const
{
  CheckRhoSurface();

  return m_pRhoSurface;
}

inline SharedSurface ModelOutput::GetUnderlyingRhoSurface() const
{
  CheckUnderlyingRhoSurface();

  return m_pUnderlyingRhoSurface;
}

inline SharedSurface ModelOutput::GetFugitSurface() const
{
  CheckFugitSurface();

  return m_pFugitSurface;
}


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_MODELOUTPUT_H_
