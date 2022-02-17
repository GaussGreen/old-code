/////////////////////////////////////////////////////////////////////////////
// Name:        test/error/myerror.h
// Purpose:     example of a custom error class for the error test
// Author:      Vadim Zeitlin
// Created:     12.03.04
// RCS-ID:      $Id: myerror.h,v 1.3 2006/06/13 14:46:33 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_TESTS_ERROR_MYERROR_H_
#define _ITO33_TESTS_ERROR_MYERROR_H_

#include "ito33/error_base.h"

namespace MyNameSpace
{
  // 66 is the offset for all error codes for this project, they have values
  // > it (and // warnings < the opposite of it)
  ITO33_DECLARE_ERROR_CLASS(600);
}

#endif // _ITO33_TESTS_ERROR_MYERROR_H_
