/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/hg/heroflags.h
// Purpose:     Computational flags class for HERO
// Created:     2006/05/19
// RCS-ID:      $Id: heroflags.h,v 1.1 2006/05/19 18:16:30 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/hg/heroflags.h
    @brief Computational flags class for HERO.
 */

#ifndef _ITO33_HG_HEROFLAGS_H_
#define _ITO33_HG_HEROFLAGS_H_

#include "ito33/date.h"

#include "ito33/hg/common.h"

namespace ito33
{

namespace hg
{

/// Computational flags class for HERO.
class ITO33_HG_DLLDECL HEROFlags
{
public:

  /**
      Default ctor.
    
      All flags are set to false by default. The analysis date will be invalid.
   */
  HEROFlags() : m_bComputeSurface(false)
  {  
  }

  // default dtor is ok
 
  ///@name Modifier functions
  //@{

  /**
      The flag for HERO surface computation.

      @param bComputeSurface Flag for the HERO surface computation
   */
  void SetComputeSurface(bool bComputeSurface = true)
  {
    m_bComputeSurface = bComputeSurface;
  }

  /**
      The analysis date at which the user wishes to view the HERO outputs 
      against stock levels.

      The user can also set an invalid date to de-activate this functionality.

      @param analysisDate The analysis date
   */
  void SetAnalysisDate(Date analysisDate)
  {
    m_analysisDate = analysisDate;
  }

  //@} // name initialization functions

  
  ///@name Accessor functions
  //@{
  
  /** 
      The analysis date at which the user wishes to view the HERO outputs 
      against stock levels.

      @return The analysis date
   */
  Date GetAnalysisDate() const { return m_analysisDate; }

  /** 
      The flag for HERO surface computation.

      @return The flag for HERO surface computation
   */
  bool GetComputeSurface() const { return m_bComputeSurface; }

  //@}

private:

  /// the analysis date - user wishes to view the outputs against stock levels
  Date m_analysisDate;  
  
  /// true if output surfaces (prices and greeks) should be computed
  bool m_bComputeSurface;

}; // class HEROFlags


} // namespace hg

} // namespace ito33

#endif // #ifndef _ITO33_HG_HEROFLAGS_H_
