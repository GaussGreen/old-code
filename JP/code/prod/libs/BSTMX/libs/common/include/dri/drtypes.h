#ifndef _dri_drtypes_h
#define _dri_drtypes_h

#include <drtypes.h>

#include <drtraits.h>

namespace drdri {


LG_ENUM_BEGIN (KStubType)
    LG_DESCRIPTION (    "Coupon payment stub types")
    LG_ENUM_ID (NONE,   "NONE", "")
    LG_ENUM_ID (BOND,   "BOND", "")
    LG_ENUM_ID (SIMPLE, "SIMPLE", "")
    LG_ENUM_ID (PAR,    "PAR", "")
LG_END()
 

LG_ENUM_BEGIN (KStubLocation)
    LG_DESCRIPTION ("Stub location")
    LG_ENUM_ID (SHORT_FRONT, "SHORT_FRONT", "")
    LG_ENUM_ID (SHORT_BACK,  "SHORT_BACK", "")
    LG_ENUM_ID (LONG_FRONT,  "LONG_FRONT", "")
    LG_ENUM_ID (LONG_BACK,   "LONG_BACK", "")
LG_END()

LG_ENUM_BEGIN (KDcc)
    LG_DESCRIPTION ("Day count convention")
    LG_ENUM_ID (DCC_30_360,  "30/360", "")
    LG_ENUM_ID (DCC_ACT_360, "ACT/360", "")
    LG_ENUM_ID (DCC_ACT_365, "ACT/365", "")
    LG_ENUM_ID (DCC_ACT_ACT, "ACT/ACT", "")
LG_END()


LG_ENUM_BEGIN (KProbDistType)
    LG_DESCRIPTION ("Probability distribution type")
    LG_ENUM_ID (NORMAL,     "NORMAL", "")
    LG_ENUM_ID (LOGNORMAL,  "LOGNORMAL", "")
LG_END()


LG_ENUM_BEGIN (KFrequency)
    LG_DESCRIPTION ("Frequency")
    LG_ENUM_ID (ANNUAL,         "A", "ANNUAL")
    LG_ENUM_ID (SEMI_ANNUAL,    "S", "SEMI_ANNUAL")
    LG_ENUM_ID (QUARTERLY,      "Q", "QUARTERLY")
    LG_ENUM_ID (IMM,            "I", "IMM")
    LG_ENUM_ID (MONTHLY,        "M", "MONTHLY")
    LG_ENUM_ID (WEEKLY,         "W", "WEEKLY")
    LG_ENUM_ID (DAILY,          "D", "DAILY")
LG_END()

LG_ENUM_BEGIN (KMMBasis)
    LG_DESCRIPTION ("Money market basis")
    LG_ENUM_ID (MMB_360,    "360", "")
    LG_ENUM_ID (MMB_365,    "365", "")
    LG_ENUM_ID (MMB_ACT,    "ACT", "")
LG_END()

}

#endif
