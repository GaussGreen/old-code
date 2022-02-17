/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/computationalflags.h
// Purpose:     Computational flags class
// Created:     04/02/26
// RCS-ID:      $Id: computationalflags.h,v 1.32 2006/07/20 03:17:44 dave Exp $
// Copyright:   (c) 2004 -2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/finance/computationalflags.h
   @brief Class for computational flags.
 */

#ifndef _ITO33_FINANCE_COMPUTATIONALFLAGS_H_
#define _ITO33_FINANCE_COMPUTATIONALFLAGS_H_

#include "ito33/date.h"
#include "ito33/vector.h"
#include "ito33/dlldecl.h"

namespace ito33
{

namespace finance
{


/**
   Class for computational flags.

   Note: delta, gamma and theta are always computed so they do not need flags.
 */
class ITO33_DLLDECL ComputationalFlags
{
public:

  /**
     Default ctor.
    
     All flags are set to false by default. The analysis date will be invalid.
   */
  ComputationalFlags() : m_bComputeSurface(false),
                         m_bComputeVega(false),
                         m_bComputeFugit(false),
                         m_bComputeRho(false),
                         m_bActivateAllSensitivities(false),
                         m_bHasSensitivityFlags(false),
                         m_iSensitivityMethod(0),
                         m_iDiscretizationMethod(0),
                         m_iSolverType(0),
                         m_bUseAnalyticSolvers(true),
                         m_bUseSimilarityReductions(true)
  {  
  }

  // default dtor is ok
 
  ///@name Modifier functions
  //@{

  /**
     The flag for vega computation. 

     @param bComputeVega Flag indicates if we need to calculate vega
   */
  void SetComputeVega(bool bComputeVega)
  {
    m_bComputeVega = bComputeVega;
  }
  
  /**
     The flag for fugit computation. 

     @param bComputeFugit Flag indicates if we need to calculate fugit
   */
  void SetComputeFugit(bool bComputeFugit)
  {
    m_bComputeFugit = bComputeFugit;
  }

  /**
     The flag for rho computation.

     @param bComputeRho Flag indicates if we need to calculate rho
   */
  void SetComputeRho(bool bComputeRho)
  {
    m_bComputeRho = bComputeRho;
  }

  /**
     The flag for the surface computation.

     @param bComputeSurface Flag for the surface computation
   */
  void SetComputeSurface(bool bComputeSurface = true)
  {
    m_bComputeSurface = bComputeSurface;
  }

  /**
     The analysis date at which the user wishes to view the outputs (price and 
     greeks) against stock levels.
     
     As a consequence, the ModelOutput object will hold the array of price 
     values at this date as well as other model dependent values.

     The user can also set an invalid date to de-activate this functionality.

     @param analysisDate The analysis date
   */
  void SetAnalysisDate(Date analysisDate)
  {
    m_analysisDate = analysisDate;
  }

  /**
     Computes all the sensitivities when the flag passed in is set
     to true.

     @param bActivateAllSensitivities flag indicates if we need to calculate 
            all sensitivities
   */
  void ActivateAllSensitivities(bool bActivateAllSensitivities)
  {
    m_bActivateAllSensitivities = bActivateAllSensitivities;
  }

  /**
     @brief The flag for sensitivity computation.

     @param pbComputeSensitivities flags indicate if we need to calculate 
            sensitivity for each variable

   */
  void SetSensitivityFlags(const std::vector<bool>& pbComputeSensitivities)
  {
    m_pbComputeSensitivities = pbComputeSensitivities;

    // Only set Has flag to true if at least one flag is true
    m_bHasSensitivityFlags = false;
    for (size_t nIdx = 0; nIdx < pbComputeSensitivities.size(); nIdx++)
    {
      if ( pbComputeSensitivities[nIdx] == true )      
      {
        m_bHasSensitivityFlags = true;
        break;
      }
    } // check if any flag is true

  }

  /**
     @internal
     @brief The method used to compute the sensitivity.

     The pricing code uses an enum, but use an int here for simplicity.
     Also, since this will be a hidden flag, do no want to include
     other files.

     @param iMethod The sensitivity method: 0 = none, 1 = PDE, 2 = adjoint

     @noexport
   */
  void SetSensitivityMethod(int iMethod)
  {
    m_iSensitivityMethod = iMethod;
  }

  /**
     @internal
     @brief The (space) discretization method used to discretize the PDE.

     @param iDiscretizationMethod The discretization method: 0 = FD, 1 = FE

     @noexport
   */
  void SetDiscretizationMethod(int iDiscretizationMethod)
  {
    m_iDiscretizationMethod = iDiscretizationMethod;
  }

  /**
     @internal
     @brief The type of the solver used to solve the non linear system.

     @param iSolverType Type of the solver used to solve the non linear system:
                        0 = Penalty, 1 = Frozen

     @noexport
   */
  void SetSolverType(int iSolverType)
  {
    m_iSolverType = iSolverType;
  }

  /**
     @internal
     @brief Whether or not analytic solvers are used (when available).

     @param bUseAnalyticSolvers true if analytic solvers are to be used, false
                                otherwise

     @noexport
   */
  void SetUseAnalyticSolvers(bool bUseAnalyticSolvers)
  {
    m_bUseAnalyticSolvers = bUseAnalyticSolvers;
  }

  /**
     @internal
     @brief Whether or not similarity reductions are used (when available).

     @param bUseSimilarityReductions true if similarity reductions are to be 
                                     used, false otherwise

     @noexport
   */
  void SetUseSimilarityReductions(bool bUseSimilarityReductions)
  {
    m_bUseSimilarityReductions = bUseSimilarityReductions;
  }

  //@} // name initialization functions

  
  ///@name Accessor functions
  //@{

  /**
     The flag to compute all all sensitivities.
     Indicates if the computation of all the sensitivities has been requested.

     @return ComputeSensitivities flag 
   */
  bool AreAllSensitivitiesActivated() const
  {
    return m_bActivateAllSensitivities;
  }

  /**     
     The flag for sensitivities computation.

     @return ComputeSensitivities flag
   */
  const std::vector<bool>& GetSensitivityFlags() const
  {
    return m_pbComputeSensitivities;
  }

  /**      
     Whether or not sensitivity flags were set.

     @return true if sensitivity flags were set, false otherwise
   */
  bool HasSensitivityFlags() const
  {
    return m_bHasSensitivityFlags;
  }

  /**
     @internal 
     @brief The sensitivity method.

     @return the sensitivity method, as an int

     @noexport
   */
  int GetSensitivityMethod() const
  {
    return m_iSensitivityMethod;
  }

  /**
     @internal 
     @brief The discretization method, FD = 0, FE = 1.

     @return the discretization method, as an int

     @noexport
   */
  int GetDiscretizationMethod() const
  {
    return m_iDiscretizationMethod;
  }

  /**
     @internal 
     @brief The type of the solver used to solve the non linear system,
            0 = Penalty, 1 = Frozen.

     @return the type of the solver, as an int

     @noexport
   */
  int GetSolverType() const
  {
    return m_iSolverType;
  }

  /**
     @internal
     @brief Whether or not analytic solvers are used when available.

     @return true if analytic solvers are to be used, false otherwise

     @noexport
   */
  bool GetUseAnalyticSolvers() const
  {
    return m_bUseAnalyticSolvers;
  }

  /**
     @internal
     @brief Whether or not similarity reductions are used (when available).

     @param bUseSimilarityReductions true if similarity reductions are to be 
                                     used, false otherwise

     @noexport
   */
  bool GetUseSimilarityReductions() const
  {
    return m_bUseSimilarityReductions;
  }

  /**
     The flag for vega computation.

     @return vega computation flag 
   */
  bool GetComputeVega() const { return m_bComputeVega; }

  /**
     The flag for fugit computation.

     @return fugit computation flag
   */
  bool GetComputeFugit() const { return m_bComputeFugit; } 

  /**
     The flag for rho computation. 

     @return rho computation flag
   */
  bool GetComputeRho() const
  {
    return m_bComputeRho;
  }

  /** 
     The flag for surface computation.

     @return The flag for surface computation
   */
  bool GetComputeSurface() const { return m_bComputeSurface; }

  /** 
     The analysis date at which the user wishes to view the outputs (price and
     greeks) against stock levels.

     @return The analysis date
   */
  Date GetAnalysisDate() const { return m_analysisDate; }

  //@}

  /**
     @internal 
     @brief Turn off all flags to compute the greeks.
     @todo Either remove this function(only used by ihg attachwarrant pricer),
           or do it correctly also for sensitivities

     @noexport
   */
  void TurnOffAllGreeks()
  {
    m_bComputeVega = false;
    m_bComputeFugit = false;
    m_bComputeRho = false;
  }

private:

  /// true if output surfaces (prices and greeks) should be computed
  bool m_bComputeSurface;
  
  /// the analysis date - user wishes to view the outputs against stock levels
  Date m_analysisDate;

  /// boolean indicates if we compute the vega
  bool m_bComputeVega;

  /// boolean indicates if we compute the fugit
  bool m_bComputeFugit;
  
  /// boolean indicates if we compute the rho
  bool m_bComputeRho;

    /// boolean indicates if we compute all sensitivities
  bool m_bActivateAllSensitivities;
  
  /// booleans indicate if we compute partial sensitivities
  std::vector<bool> m_pbComputeSensitivities;

  /// boolean indicates if sensitivity flags were set
  bool m_bHasSensitivityFlags;

  /// The method used to compute sensitivities
  int m_iSensitivityMethod;

  /// The method used to discretize the PDE
  int m_iDiscretizationMethod;

  /// The type of the solver for non linear system
  int m_iSolverType;

  /// Whether or not analytic solvers are to be used
  bool m_bUseAnalyticSolvers;

  /// Whether or not similarity reductions are to be used
  bool m_bUseSimilarityReductions;

}; // class ComputationalFlags


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_COMPUTATIONALFLAGS_H_
