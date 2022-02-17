//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SetNewHandler.hpp
//
//   Description : Cross-platform implementation of set_new_handler()
//                 Exists because Microsoft C++ does not comply with ANSI 
//                 standard
//
//   Author      : Andrew J Swain
//
//   Date        : 16 November 2000
//
//
//----------------------------------------------------------------------------

#ifndef SETNEWHANDLER_HPP
#define SETNEWHANDLER_HPP

DRLIB_BEGIN_NAMESPACE

class TOOLKIT_DLL CSetNewHandler{
    public:
    // default handler - new throws a bad_alloc exception on failure
    static void handler();
    
private:
    CSetNewHandler();
};

DRLIB_END_NAMESPACE

#endif
