//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Library.hpp
//
//   Description : Starts up & shuts down the DR library
//
//   Author      : Andrew J Swain
//
//   Date        : 29 January 2001
//
//
//----------------------------------------------------------------------------

#ifndef LIBRARY_HPP
#define LIBRARY_HPP

DRLIB_BEGIN_NAMESPACE

/** Starts up & shuts down the DR library */

class ADDINS_DLL Library{
public:
    /** This is invoked at 'start up' immediately after SpreadSheetMode is
        configured */
    static void configureErrorLog();

    /** call before using any library functions (well almost any) */
    static void startup();

    /** call before program exit */
    static void shutdown();
        
    //// The prefix to use for Excel functions. Required by 
    //// XLAddin and xlregister. Set to "EDR_"
    static const char* XL_PREFIX;

    //// the name for the DR Service Interface eg EDG or COGS etc
    static const string SERVICE_NAME;

    //// the name for the thread-safe version of DR Service Interface
    static const string THREADSAFE_SERVICE_NAME;

    //// the version (as encoded by DRLIB_VERSION) for the DR Service
    //// Interface eg EDG or COGS etc
    static const int SERVICE_VERSION;

private:
    Library();
};

DRLIB_END_NAMESPACE

#endif
