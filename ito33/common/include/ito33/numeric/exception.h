/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/exception.h
// Purpose:     numerical Excpetions
// Author:      Wang
// Created:     2004/07/15
// RCS-ID:      $Id: exception.h,v 1.2 2004/10/05 09:13:38 pedro Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/exception.h
    @brief Classes for working with numerical exceptions.

    We would like to hide the numerical exception to the user, so we need to
    be able to catch these exceptions, and reinterpret them.
 */

#ifndef _ITO33_NUMERIC_EXCEPTION_H_
#define _ITO33_NUMERIC_EXCEPTION_H_

#include "ito33/string.h"
#include "ito33/exception.h"

namespace ito33
{

namespace numeric
{


/**
   Class represents numerical exceptions.
 */
class Exception : public ito33::Exception
{
public:
  /**
      Ctor for the Exception object.

      Use the standard EXCEPTION macro to create Exception objects, this
      frees you from having to type __FILE__, __LINE__ and __FUNCTION__
    */
  Exception(int errorCode,
            const std::string& message,
            const char *filename,
            size_t line,
            const char *function)
    : ito33::Exception(errorCode, message, filename, line, function)
    {
    }

}; // class Exception


} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_EXCEPTION_H_

