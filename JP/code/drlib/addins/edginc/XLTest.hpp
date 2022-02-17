//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLTest.hpp
//
//   Description : Creating/running Tests for Addins and XL Kit
//
//   Author      : Mark A Robson
//
//   Date        : 27 Feb 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_XLTEST_HPP
#define EDG_XLTEST_HPP

#include "edginc/xlapi.hpp"
#include <stdio.h>
#include <stdarg.h>

DRLIB_BEGIN_NAMESPACE

class ADDINS_DLL XLTest{
public:
    /** Run any addin functin where the name is specified by the first OPER - 
        addin function name needs to be without the EDR_ prefix */
    static XL_OPER* generic(
        const XL_OPER*  addinNameOper, /* (I) */
        const XL_OPER*  fileNameOper,  /* (I) may be NULL */
        const XL_OPER*  a0,            /* (I) first real parameter */
        va_list         args);         /* (I) remaining args */

    /* Run any addin functin where the name is specified by the first OPER - 
       addin function name needs to be without the EDG_ prefix */
    static XL_OPER* genericXL(
        bool            calledFromVB,  // (I) true if called from VB
        const XL_OPER*  addinNameOper, /* (I) */
        const XL_OPER*  fileNameOper,  /* (I) */
        va_list         args);         /* (I) remaining args */

};

DRLIB_END_NAMESPACE

#endif
