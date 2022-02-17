/******************************************************************************
 * File.........: ito33/numeric/solvemethod.h
 * Purpose......: Definition of common enumerations
 * Author.......: (z)
 * Created......: 10/08/2003
 * RCS-ID.......: $Id: solvemethod.h,v 1.3 2006/05/26 13:34:31 nabil Exp $
 * Copyright....: (c) 2003 Trilemma LLP
 ******************************************************************************/

/**
  \file ito33/numeric/solvemethod.h
  \brief enumeration:  the method to use for solving PDE with constraints
  */
#ifndef _ITO33_NUMERIC_SOLVEMETHOD_H_
#define _ITO33_NUMERIC_SOLVEMETHOD_H_

namespace ito33
{

namespace numeric
{

/// the method to use for solving PDE with constraints
enum ITO33_SolveMethod
{  
  Solve_Linear    = 1,
  Solve_NonLinear  = 2,
  SolveMethod_Max
};

} // namespace nuemric

} // namespace ito33


#endif // #ifndef _ITO33_COMMONENUM_H
