//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : StubFactory.hpp
//
//   Description : Factory class for building Stubs
//
//   Author      : Andrew J Swain
//
//   Date        : 1 March 2002
//
//
//----------------------------------------------------------------------------

#ifndef STUBFACTORY_HPP
#define STUBFACTORY_HPP

#include "edginc/Stub.hpp"

DRLIB_BEGIN_NAMESPACE

/** Factory class for building Stubs */
class MARKET_DLL StubFactory{
public:
    static Stub* make(const string& stub);
    static Stub* clone(const Stub* stub);

    static const string SIMPLE;
    static const string BOND;
    static const string NONE;
    
private:
    StubFactory();  
};

DRLIB_END_NAMESPACE

#endif
