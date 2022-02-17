//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : SystemException.cpp
//
//   Description : System (structured) exception handler.
//
//   Author      : Anwar E Sidat
//
//   Date        : 07 October 2006
//
//----------------------------------------------------------------------------

#include <iostream>
#include "edginc/coreConfig.hpp"
#include "edginc/ErrorHandler.hpp"
#include "edginc/ModelException.hpp"
#include "edginc/SystemException.hpp"

CORE_BEGIN_NAMESPACE

#if defined(WIN32)

/* static */
_se_translator_function Win32StructuredException::sm_savedTranslator = NULL;


/** Translate function for win32 structured exceptions (converts to C++ exceptions). */
void Win32StructuredException::TranslateException(
    unsigned int        uExceptionCode,
    EXCEPTION_POINTERS* pExp)
{
    switch (uExceptionCode)
    {
    case EXCEPTION_FLT_UNDERFLOW:
    case EXCEPTION_FLT_DIVIDE_BY_ZERO:
    case EXCEPTION_FLT_DENORMAL_OPERAND:
    case EXCEPTION_FLT_INEXACT_RESULT:
    case EXCEPTION_FLT_INVALID_OPERATION:
    case EXCEPTION_FLT_OVERFLOW:
    case EXCEPTION_FLT_STACK_CHECK:
        _clearfp();
        _fpreset();
        throw Win32FloatingPointException(pExp);

    default:
        throw Win32StructuredException(pExp);
    };
}

/**  Constructor for Win32StructuredException. */
Win32StructuredException::Win32StructuredException(EXCEPTION_POINTERS* pExp)
    : SystemError("")
{
    m_dwExceptionCode = pExp->ExceptionRecord->ExceptionCode;
    m_pvAddress = pExp->ExceptionRecord->ExceptionAddress;
    string strMessage;

    switch ( m_dwExceptionCode )
    {
    case EXCEPTION_FLT_DENORMAL_OPERAND:
        strMessage = "Floating Point Denormal Operand";
        break;

    case EXCEPTION_FLT_DIVIDE_BY_ZERO:
        strMessage = "Floating Point Division By Zero";
        break;

    case EXCEPTION_FLT_INEXACT_RESULT:
        strMessage = "Floating Point Inexact Result";
        break;

    case EXCEPTION_FLT_INVALID_OPERATION:
        strMessage = "Floating Point Invalid Operation";
        break;

    case EXCEPTION_FLT_OVERFLOW:
        strMessage = "Floating Point Overflow";
        break;

    case EXCEPTION_FLT_STACK_CHECK:
        strMessage = "Floating Point Stack Over/Underflow";
        break;

    case EXCEPTION_FLT_UNDERFLOW:
        strMessage = "Floating Point Underflow";
        break;

    case EXCEPTION_ACCESS_VIOLATION:
        strMessage = "Access Violation";
        break;

    case EXCEPTION_DATATYPE_MISALIGNMENT:
        strMessage = "Data Type Misalignmenet";
        break;
    
    case EXCEPTION_BREAKPOINT:
        strMessage = "Breakpoint Encountered";
        break;

    case EXCEPTION_SINGLE_STEP:
        strMessage = "Single Step";
        break;

    case EXCEPTION_ARRAY_BOUNDS_EXCEEDED:
        strMessage = "Array Bounds Exceeded";
        break;

    case EXCEPTION_ILLEGAL_INSTRUCTION:
        strMessage = "Illegal Instruction";
        break;

    case EXCEPTION_IN_PAGE_ERROR:
        strMessage = "Page Error";
        break;

    case EXCEPTION_INT_DIVIDE_BY_ZERO:
        strMessage = "Integer Division By Zero";
        break;

    case EXCEPTION_INT_OVERFLOW:
        strMessage = "Integer Overflow";
        break;

    case EXCEPTION_PRIV_INSTRUCTION:
        strMessage = "Privileged Instruction";
        break;

    case EXCEPTION_STACK_OVERFLOW:
        strMessage = "Stack Overflow";
        break;

    default:
        strMessage = "Unknown System Exception";
        break;
    }

    addMsg("Win32StructuredException - " + strMessage + "!");
}

/** Structured exception hook enable. */
void Win32StructuredException::enableHook()
{
    sm_savedTranslator = _set_se_translator(Win32StructuredException::TranslateException);
    _clearfp();
    int nControlWord = _controlfp(0, 0);
    nControlWord &= ~(EM_OVERFLOW | EM_ZERODIVIDE | EM_DENORMAL);   // EM_INEXACT | EM_UNDERFLOW
    _controlfp(nControlWord, _MCW_EM);
}

/** Structured exception hook disable. */
void Win32StructuredException::disableHook()
{
    _set_se_translator(sm_savedTranslator);   // restore translator
    _clearfp();
    int nControlWord = _controlfp(0, 0);
    nControlWord |= (EM_OVERFLOW | EM_ZERODIVIDE | EM_DENORMAL);   // EM_INEXACT | EM_UNDERFLOW
    _controlfp(nControlWord, _MCW_EM);
}

#endif //defined(WIN32)

CORE_END_NAMESPACE
