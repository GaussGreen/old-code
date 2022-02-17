/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/numericerror.h
// Purpose:     error codes for numeric namespace
// Author:      Ito33
// Date:        2004/12/14
// RCS-ID:      $Id: numericerror.h,v 1.3 2005/05/25 14:27:56 yann Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
  @file ito33/numeric/numericerror.h
  @brief This file contains the Error class declaration in numeric namespace.

 */
#ifndef _ITO33_NUMERIC_ERROR_H_
#define _ITO33_NUMERIC_ERROR_H_

namespace ito33
{

namespace numeric
{

/**
   This class represents an error code for the numerical classes.
 */
enum NumericError
     {
       ITO33_NO_ERROR,
       ITO33_NOT_CONVERGED,
       ITO33_TOO_MANY_FUNCTION_CALL,
       ITO33_TOO_MANY_ITERATION,
       ITO33_TEMPERATURE_TOO_SMALL
     };


} // namespace finance

} // namespace ito33

#endif // _ITO33_NUMERIC_ERROR_H_

