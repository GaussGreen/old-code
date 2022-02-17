///////////////////////////////////////////////////////////////////////////////
// Name:        error.cpp
// Purpose:     implementation of Error
// Author:      Vadim Zeitlin
// Created:     03.11.03
// RCS-ID:      $Id: error_base.cpp,v 1.99 2006/06/15 18:57:04 zhang Exp $
// Copyright:   (c) 2003 - 2006 Trilemma LLP
///////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/error.h"

// ----------------------------------------------------------------------------
// Implement Error class
// ----------------------------------------------------------------------------
namespace ito33
{
  ITO33_IMPLEMENT_ERROR_CLASS;
}

// ----------------------------------------------------------------------------
// definitions of the Error objects
// ----------------------------------------------------------------------------
using ito33::Error;

extern const Error ITO33_NO_LICENSE
("No license for the requested feature.");

extern const Error ITO33_UNEXPECTED
("Unexpected error.");

extern const Error ITO33_NOT_IMPLEMENTED
("The function is not implemented yet.");

extern const Error ITO33_SYS_ERROR
("System error occured.");

extern const Error ITO33_UNDEFINED_ENV_VAR
("Environment variable not defined.");

extern const Error ITO33_DLL_ERROR
("Dynamic link error.");

extern const Error ITO33_EXCEL_ERROR
("Excel error.");

extern const Error ITO33_OUT_OF_MEMORY
("Out of memory: new() or malloc() failed.");  

extern const Error ITO33_BUFFER_TOO_SMALL
("Insufficient space in the buffer (increase it and call again).");

extern const Error ITO33_OUTOFBOUND
("The index is out of bound.");

extern const Error ITO33_DIV0
("Divide by zero.");

extern const Error ITO33_NEG_TOL
("Negative tolerance.");

extern const Error ITO33_MAX_ITER
("Maximum number of iterations exceeded in algorithm.");

extern const Error ITO33_BAD_PARAM
("Bad input parameter.");

extern const Error ITO33_DIFF_PARAM
("Parameters should be different.");

extern const Error ITO33_NULL_PARAM
("Null parameter.");

extern const Error ITO33_BAD_DATA
("Bad input data.\n"
"This often happens when we do the posterior input checking.");

extern const Error ITO33_UNDEF_DATA
("Accessing undefined data.");

extern const Error ITO33_NULL_PTR
("Required output parameter is NULL.");

extern const Error ITO33_ARRAY_ZEROSIZE
("Array size is zero.");

extern const Error ITO33_ARRAY_NOTINCREASING
("The array is not strictly increasing.");

extern const Error ITO33_ARRAY_NONPOSITIVE
("Array elements must be strictly positive.");

extern const Error ITO33_ARRAY_NEGATIVE
("Array elements must be non negative.");

extern const Error ITO33_BAD_DATE
("Invalid date specified.");

extern const Error ITO33_INVALID_DAYCOUNTCONVENTION
("Invalid day count convention value.");

//_________________________Optional value_____________________________
extern const Error ITO33_NO_VALUE
("Atempt to access invalid optional value.");   
