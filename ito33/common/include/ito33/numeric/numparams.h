/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/numparams.h
// Purpose:     numerical parameters class
// Author:      Wang
// Created:     2004/02/09
// RCS-ID:      $Id: numparams.h,v 1.12 2004/10/28 14:34:41 wang Exp $
// Copyright:   (c) 2003 - 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
  @file ito33/numeric/numparams.h
  @brief definition of the class for pure Numerical Parameters 
 */
  
#ifndef _ITO33_NUMERIC_NUMPARAMS_H_
#define _ITO33_NUMERIC_NUMPARAMS_H_

#include <cstddef> // for size_t

#include "ito33/numeric/schemetype.h"

namespace ito33
{

namespace finance
{
  class QualityControl;
}

namespace numeric
{

  /// Numerical Parameters class
class NumParams 
{

public:
  /// ctor
  NumParams() : m_nNbSpaceSteps(100),
                m_nNbTimeSteps(100),
                m_dMaxDeltaLogS(0.15),
                m_SchemeType(SchemeType_ThreeLevel)
  {
  }

  /**
    ctor by the QuatilityControl and the size of the computation domain
    in the time axe.

    @param control requested QualityControl 
    @param dMaturity the size of the time domain.
    */
  NumParams(const finance::QualityControl& control, 
            double dMaturity);

  /** 
    Sets the suggested number of time steps.
   
    The model may use more time steps in order to take in account all 
    event's dates to generate a grid 
    that is adapted to the problem. 

    @param nNbTimeSteps the suggested number of time steps.
   */
  void SetNbTimeSteps(size_t nNbTimeSteps)
  {
    m_nNbTimeSteps = nNbTimeSteps; 
  }

  /** 
    Sets the suggested number of space steps.
   
    The model may use more space steps in order to generate a grid 
    that is adapted to the problem. This remark is only valid if the grid 
    is not uniform.

    @param nNbSpaceSteps the suggested number of space steps.
   */
  void SetNbSpaceSteps(size_t nNbSpaceSteps)
  { 
    m_nNbSpaceSteps = nNbSpaceSteps; 
  }

  /**
    Set the tolerated highest mesh size in the main computation domain 
    of log(s).
    
    The mesh size DX is determined by the following way: 
      DX = min( Z / NbSpaceSteps, MaxDeltaLogS)
    where Z is the size of the main computation domain.

    @param dValue the value of MaxDeltaLogS
    */
  void SetMaxDeltaLogS(double dValue)
  {
    m_dMaxDeltaLogS = dValue;
  }

  /** 
    Sets the suggested schemetype
   
    The model may use implicit timestepping at special time points (time zero,
    dividend dates, etc.), but will otherwise use the suggested
    schemetype (eg. Crank-Nicolson)

    @param eSchemeType the suggested schemetype.
   */
  void SetSchemeType(SchemeType eSchemeType)
  { 
    m_SchemeType = eSchemeType; 
  }


  /**
    Get the suggested number of space steps.

    @return the suggested number of space steps.
   */
  size_t GetNbSpaceSteps() const { return m_nNbSpaceSteps; }  

  /**
    Get the suggested number of time steps.

    @return the suggested number of time steps.
   */
  size_t GetNbTimeSteps() const { return m_nNbTimeSteps; }

  /**
    Get the suggested schemetype

    @return the suggested schemetype.
   */
  SchemeType GetSchemeType() const { return m_SchemeType; }

  /**
    Gets the tolerated highest mesh size in the main computation domain 
    of log(s).
    
    @return the value of MaxDeltaLogS
    */
  double GetMaxDeltaLogS()
  {
    return m_dMaxDeltaLogS;
  }
protected:
  
  /// The suggested number of space steps.
  size_t m_nNbSpaceSteps;
    
  /// The suggested number of time steps.
  size_t m_nNbTimeSteps;  

  /// The tolerated highest mesh size in the main computation domain of log(s)
  double m_dMaxDeltaLogS;

  /// The suggested scheme type (implicit, Crank-Nicolson, etc).
  SchemeType m_SchemeType;  

}; // class NumParams 


} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_NUMPARAMS_H_

