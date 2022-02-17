//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : StubFactory.cpp
//
//   Description : Factory class for building Stubs
//
//   Author      : Andrew J Swain
//
//   Date        : 1 March 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/StubFactory.hpp"
#include "edginc/StubSimple.hpp"
#include "edginc/StubNone.hpp"
#include "edginc/StubBond.hpp"

DRLIB_BEGIN_NAMESPACE

const string StubFactory::SIMPLE = "Simple";
const string StubFactory::BOND   = "Bond";
const string StubFactory::NONE   = "None";


Stub* StubFactory::make(const string& stub) {
    static string method = "StubFactory::make";
    if (CString::equalsIgnoreCase(stub, SIMPLE) ||
        CString::equalsIgnoreCase(stub, SIMPLE, 1)) {
        return new StubSimple();
    }
    else if (CString::equalsIgnoreCase(stub, NONE) ||
             CString::equalsIgnoreCase(stub, NONE, 1)) {
        return new StubNone();
    }
    else if (CString::equalsIgnoreCase(stub, BOND) || 
             CString::equalsIgnoreCase(stub, BOND, 1)) {
        return new StubBond();
    }
    else {
        throw ModelException(method,
                             "Unknown stub type: " + stub);
    }
}

Stub* StubFactory::clone(const Stub* stub) {
    return make(stub->toString());
}

StubFactory::StubFactory() {
    // empty
}

DRLIB_END_NAMESPACE

