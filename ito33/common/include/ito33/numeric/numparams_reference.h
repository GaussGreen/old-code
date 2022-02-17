/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/numparams_reference.h
// Purpose:     bridge between finance/QualityControl and numeric/NumParams
// Author:      ZHANG Yunzhi
// Created:     2004/05/15
// RCS-ID:      $Id: numparams_reference.h,v 1.6 2004/11/10 16:08:57 afrolov Exp $
// Copyright:   (c) 2003 - 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
  @file ito33/numeric/numparams_reference.h
  @brief the bridge between finance/userparams and numeric/NumParamsReference

  This file offers functions to get values of some numerical variables 
  to translate a finance::QualityControl object to a numeric::NumParams 
  (or MeshParams) object.
 */
  
#ifndef _ITO33_NUMERIC_NUMPARAMS_REFERENCE_H_
#define _ITO33_NUMERIC_NUMPARAMS_REFERENCE_H_

#include <stddef.h>
#include "ito33/numeric/schemetype.h"

namespace ito33
{

namespace numeric
{

namespace NumParamsReference
{

  /**
    Get NbTimeStepsFor5Years which is the minimum number of points we need
    in the time mesh for computation of a instrument whose maturity is exactly
    in 5 years.
    
    The value of NbTimeSteps in numeric::NumParams class could be calculated
    in following way by NbTimeStepsFor5Years "N" and the compution time "T" 
    (maturity of the instrument)
    NbTimeSteps = N / 5        if T <= 1
                = N * T / 5    if T <= 5
                = N            if T <= 15
                = N * T / 15   otherwise

    The general idea is that for the instruments last shortly, we
    prefer fixing the number of the time points instead of fixing the mesh 
    size. That is, more closed is the maturity, more precise should be the 
    result. While for the long maturity, we need less precision even though
    we don't like to see ridiculously large mesh size.

    There is not particular reason for the choice of these values: 1 year, 
    5 years and 15 years. The feeling is that most of bonds that we are 
    instrested in have around 5 years' maturity or less. Most of options
    have short maturities, but we don't need as many time points for 1 year
    and 5 years. 

    The default value is set to 150.

    @return the value of NbTimeStepsFor5Years
    */
  size_t GetNbTimeStepsFor5Years();

  /**
    Get MinNbSpaceSteps which is the minimum number of space points
    we put in the main computation domain of log(s).

    The value of NbSpaceSteps in numeric::NumParams class is determined by 
    maximum of MinNbSpaceSteps and NbTimeSteps.

    Note that NbSpaceSteps should not be too greater than NbTimeSteps. 

    Its default value is set to 100.

    @return the value of MinNbSpaceSteps
    */
  size_t GetMinNbSpaceSteps();

  /**
    Get MaxDeltaLogS which is the highest tolerated mesh size in the main
    computation domain of log(S).

    The real mesh size DX is determined by the following way: 
      DX = min( Z / NumParams::NbSpaceSteps, NumParams::MaxDeltaLogS)
    where Z is the size of the main computation domain.
    
    Its default value is set to 0.15.

    @return the value of MaxDeltaLogS
    */
  double GetMaxDeltaLogS();

  /**
    Get the suggested schemetype

    @return the suggested schemetype.
   */
  SchemeType GetSchemeType();
  
} // namespace NumParamsReference 


} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_NUMPARAMS_REFERENCE_H_

