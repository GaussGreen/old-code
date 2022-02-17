//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SystemException.hpp
//
//   Description : System (structured) exception handler.
//
//                 WARNING: DO NOT INCLUDE THIS FILE UNLESS YOU HAVE TO AS IT
//                          INCLUDES windows.h on WIN32 platform.
//
//   Author      : Anwar E Sidat
//
//   Date        : 07 October 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_SystemException_HPP
#define QLIB_SystemException_HPP

#include <string>
#include "ModelException.hpp"
#if defined(WIN32)
    #include <windows.h>
#endif

CORE_BEGIN_NAMESPACE

/**
* Operating system exceptions such as access violations and floating point exceptions etc.
*/
class CORE_DLL SystemError : public ModelException
{
public:
	/** Initialisation constructor. */
	SystemError(const std::string& msg = "System Error")
		: ModelException(msg) {;}

	/** Copy constructor. */
	SystemError(const SystemError& rhs)
		: ModelException(rhs) {;}
};


#if defined(WIN32)
/**
* Win32 specific structured exception such as access violations or floating point errors.
*/
class CORE_DLL Win32StructuredException	: public SystemError
{
public:
    /** Initialisation constructor */
    Win32StructuredException( EXCEPTION_POINTERS* pExp );

    /** Copy constructor */
    Win32StructuredException( const Win32StructuredException& rhs )
	    : SystemError(rhs), m_dwExceptionCode(rhs.m_dwExceptionCode), m_pvAddress(rhs.m_pvAddress) {;}

    /** Destructor */
    virtual ~Win32StructuredException() throw() {;}

    /** Enables the hook to convert Win32 structured exceptions into C++ exceptions. */
    static void enableHook();

    /** Disables the Win32 structured exception handler hook. */
    static void disableHook();

    /**
    *  HookGuard class for automatically enabling/disabling the Win32 structured exception
    *  handler by declaring an instance of it on the stack at library entry points (eg. at main
    *  entry of Excel functions etc.)
    */
    class CORE_DLL HookGuard
    {
    public:
	    /** Default Constructor - enables the Win32 structured exception handler hook. */
        HookGuard()  { Win32StructuredException::enableHook(); }

	    /** Destructor - disables the Win32 structured exception handler hook. */
        ~HookGuard() { Win32StructuredException::disableHook(); }
    };

protected:
	DWORD m_dwExceptionCode;
	PVOID m_pvAddress;

private:
	static void TranslateException(unsigned int u, EXCEPTION_POINTERS* pExp);
	static _se_translator_function sm_savedTranslator;
};


/**
* Floating Point exceptions (specific types of Win32 structured exceptions.
*/
class CORE_DLL Win32FloatingPointException : public Win32StructuredException
{
public:
    /** Initialisation constructor. */
    Win32FloatingPointException(EXCEPTION_POINTERS* pExp)
        : Win32StructuredException(pExp) {;}

    /** Copy constructor. */
    Win32FloatingPointException(const Win32FloatingPointException& rhs)
        : Win32StructuredException(rhs) {;}

    /** Destructor. */
    virtual ~Win32FloatingPointException() throw() {;}
};

#endif  // defined(WIN32)


CORE_END_NAMESPACE

#endif  // QLIB_SystemException_HPP
