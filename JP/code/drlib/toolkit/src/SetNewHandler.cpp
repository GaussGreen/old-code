//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SetNewHandler.cpp
//
//   Description : Cross-platform implmentation of set_new_handler()
//                 Exists because Microsoft C++ does not comply with ANSI 
//                 standard
//
//   Author      : Andrew J Swain
//
//   Date        : 16 November 2000
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SetNewHandler.hpp"

#if defined(WIN32) && !defined(__CYGWIN__)
#include <new.h>
#endif

using namespace std;

// default handlers
#if defined(WIN32) && !defined(__CYGWIN__)
int newHandler(size_t i)
{
    throw bad_alloc();
}
#else
void newHandler()
{
    throw bad_alloc();
}
#endif

DRLIB_BEGIN_NAMESPACE

// default handler - new throws a bad_alloc exception on failure
void CSetNewHandler::handler() {
#if defined(WIN32) && !defined(__CYGWIN__)
    _set_new_handler(newHandler);  
#else
    set_new_handler(newHandler);   
#endif
}

DRLIB_END_NAMESPACE

