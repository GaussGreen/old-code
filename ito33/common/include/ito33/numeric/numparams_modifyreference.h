/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/numparams_modifyreference.h
// Purpose:      
// Author:      ZHANG Yunzhi
// Created:     2004/05/15
// RCS-ID:      $Id: numparams_modifyreference.h,v 1.4 2004/10/05 09:13:38 pedro Exp $
// Copyright:   (c) 2003 - 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
  @file ito33/numeric/numparams_modifyreference.h
  @brief the bridge between finance/userparams and numeric/NumParamsReference

  STOP!
  THIS FILE SHOULD ONLY INCLUDED BY PROJECTS FOR TESTS.

  This file offers some functions to modify the numerical variables in 
  NumParamsReference namespace. The functions must be used with many attention.
  Normally, Developement&Testing need them for INTERNAL use, such as 
  convergence quality test etc.
 */
  
#ifndef _ITO33_NUMERIC_NUMPARAMS_MODIFYREFERENCE_H_
#define _ITO33_NUMERIC_NUMPARAMS_MODIFYREFERENCE_H_

#ifdef ITO33_TEST_MODENV

namespace ito33
{

namespace numeric
{

enum SchemeType;

namespace NumParamsReference
{
  
  /**
    Set the value of MinNbSpaceSteps.

    See numparams_reference.h for its description.

    @return nNumber the value of MinNbSpaceSteps.
    */
  void SetMinNbSpaceSteps(size_t nNumber);

  /**
    Set the value of NbTimeStepsFor5Years.

    See numparams_reference.h for its description.

    @return nNumber the value of NbTimeStepsFor5Years.
    */
  void SetNbTimeStepsFor5Years(size_t nNumber);

  /**
    Set the value of MaxDeltaLogS.

    See numparams_reference.h for its description.
    
    @param dValue the value of MaxDeltaLogS
    */
  void SetMaxDeltaLogS(double dValue);

  /** 
    Sets the suggested schemetype
   
    The model may use implicit timestepping at special time points (time zero,
    dividend dates, etc.), but will otherwise use the suggested
    schemetype (eg. Crank-Nicolson)

    @param eSchemeType the suggested schemetype.
   */
  void SetSchemeType(SchemeType schemeType);

} // namespace NumParamsReference 


} // namespace numeric

} // namespace ito33

#endif

#endif // #ifndef _ITO33_NUMERIC_NUMPARAMS_MODIFYREFERENCE_H_

