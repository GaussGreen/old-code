#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "drtypes.h"

char const* KStubTypeStrings[] =
{
    "NONE",
    "BOND",
    "SIMPLE",
    "PAR"
};

char const* KStubLocationStrings[] =
{
    "FRONT",
    "BOND",
    "LONG",
    "SHORT"
};

char const* KDccStrings[] =
{
    "30/360",
    "ACT/360",
    "ACT/365",
    "ACT/ACT"
};

char const* KProbDistTypeStrings[] =
{
    "NORMAL",
    "LOGNORMAL"
};

char const* KOptTypeStrings[] =
{
    "CALL",
    "PUT"
};


char const* KFrequencyStrings[] =
{
    "ANNUAL",
    "SEMI_ANNUAL",
    "QUARTERLY",
    "IMM",
    "MONTHLY",
    "WEEKLY",
    "DAIL"
};


KOptType
Char2OptType(char optionType)
{
    switch((char)toupper(optionType))
    {
    case 'C':
        return OPT_CALL;
    case 'P':
        return OPT_PUT;
    default:
        return -1;
    }

}

