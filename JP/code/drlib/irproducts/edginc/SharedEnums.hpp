#ifndef _SHAREDENUMS_HPP
#define _SHAREDENUMS_HPP

//----------------------------------------------------------------------------
//
//   Description : This file defines enumeration types for QLib.
//
//----------------------------------------------------------------------------

#include "edginc/BoxedEnum.hpp"

DRLIB_BEGIN_NAMESPACE

/******/

struct RateType { 
    enum Enum { SIMPLE, CONTINUOUS };
};
typedef BoxedEnum<RateType::Enum> RateTypeBoxedEnum;

/******/

struct StubPos { 
    enum Enum { STUB_START, STUB_END, STUB_NONE };
};
typedef BoxedEnum<StubPos::Enum> StubPosBoxedEnum;

/******/

struct StubType { 
    enum Enum { STUB_NONE, STUB_BOND, STUB_SIMPLE, STUB_PAR, STUB_NOT_ALLOWED };
};
typedef BoxedEnum<StubType::Enum> StubTypeBoxedEnum;

/******/

DRLIB_END_NAMESPACE

#endif
