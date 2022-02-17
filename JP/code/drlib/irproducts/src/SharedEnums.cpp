#include "edginc/config.hpp"
#include "edginc/SharedEnums.hpp"

DRLIB_BEGIN_NAMESPACE

START_PUBLIC_ENUM_DEFINITION(RateType::Enum, "Types of rates");
ENUM_VALUE_AND_NAME(RateType::SIMPLE, "SIMPLE", "Simple rate: dcf * rate");
ENUM_VALUE_AND_NAME(RateType::CONTINUOUS, "CONTINUOUS", "Continuously compounded rate: (1+rate)^^dcf");
END_ENUM_DEFINITION(RateType::Enum);

START_PUBLIC_ENUM_DEFINITION(StubPos::Enum, "Stub position");
ENUM_VALUE_AND_NAME(StubPos::STUB_START, "STUB_START", "");
ENUM_VALUE_AND_NAME(StubPos::STUB_END, "STUB_END", "");
ENUM_VALUE_AND_NAME(StubPos::STUB_NONE, "STUB_NONE", "");
END_ENUM_DEFINITION(StubPos::Enum);

START_PUBLIC_ENUM_DEFINITION(StubType::Enum, "Stub type");
ENUM_VALUE_AND_NAME(StubType::STUB_NONE, "STUB_NONE", "");
ENUM_VALUE_AND_NAME(StubType::STUB_BOND, "STUB_BOND", "");
ENUM_VALUE_AND_NAME(StubType::STUB_SIMPLE, "STUB_SIMPLE", "");
ENUM_VALUE_AND_NAME(StubType::STUB_PAR, "STUB_PAR", "");
ENUM_VALUE_AND_NAME(StubType::STUB_NOT_ALLOWED, "STUB_NOT_ALLOWED", "");
END_ENUM_DEFINITION(StubType::Enum);

bool SharedEnumsLoad() {
    return RateTypeBoxedEnum::TYPE != 0 
        && StubPosBoxedEnum::TYPE != 0 
        && StubTypeBoxedEnum::TYPE != 0;
}

DRLIB_END_NAMESPACE
