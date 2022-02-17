//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ClientRunnable.hpp
//
//   Description : Interface for client runnable methods via EdrAction
//
//   Author      : Andrew J Swain
//
//   Date        : 14 June 2001
//
//
//----------------------------------------------------------------------------

#ifndef CLIENTRUNNABLE_HPP
#define CLIENTRUNNABLE_HPP

#include "edginc/Object.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interface for client runnable methods via EdrAction */
class TOOLKIT_DLL ClientRunnable{
public:
    static CClassConstSP const TYPE;

    virtual ~ClientRunnable();

    /** do something */
    virtual IObjectSP run() = 0;

protected:
    ClientRunnable();
};


DRLIB_END_NAMESPACE
#endif
