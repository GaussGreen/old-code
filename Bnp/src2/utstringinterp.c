/* ===============================================================
   FILE_NAME:	utstringinterp.c

   PURPOSE:     A few translation functions from strings to types
   =============================================================== */

#include "utallhdr.h"

#define UTMAXSTR 30

/* ---------------------------------------------------------------------------
                                        DATES
   --------------------------------------------------------------------------- */
Err interp_month(const char* constStr, SrtMonth* val)
{
    char str[UTMAXSTR + 1];
    strncpy(str, constStr, UTMAXSTR);
    str[UTMAXSTR] = '\0';

    strupper(str);
    strip_white_space(str);

    /* get month */
    if (!strcmp(str, "JAN"))
    {
        *val = SRT_JAN;
        return 0;
    }
    if (!strcmp(str, "JANUARY"))
    {
        *val = SRT_JAN;
        return 0;
    }
    if (!strcmp(str, "FEB"))
    {
        *val = SRT_FEB;
        return 0;
    }
    if (!strcmp(str, "FEBRUARY"))
    {
        *val = SRT_FEB;
        return 0;
    }
    if (!strcmp(str, "MAR"))
    {
        *val = SRT_MAR;
        return 0;
    }
    if (!strcmp(str, "MARCH"))
    {
        *val = SRT_MAR;
        return 0;
    }
    if (!strcmp(str, "APR"))
    {
        *val = SRT_APR;
        return 0;
    }
    if (!strcmp(str, "APRIL"))
    {
        *val = SRT_APR;
        return 0;
    }
    if (!strcmp(str, "MAY"))
    {
        *val = SRT_MAY;
        return 0;
    }
    if (!strcmp(str, "JUN"))
    {
        *val = SRT_JUN;
        return 0;
    }
    if (!strcmp(str, "JUNE"))
    {
        *val = SRT_JUN;
        return 0;
    }

    if (!strcmp(str, "JUL"))
    {
        *val = SRT_JUL;
        return 0;
    }
    if (!strcmp(str, "JULY"))
    {
        *val = SRT_JUL;
        return 0;
    }
    if (!strcmp(str, "AUG"))
    {
        *val = SRT_AUG;
        return 0;
    }
    if (!strcmp(str, "AUGUST"))
    {
        *val = SRT_AUG;
        return 0;
    }
    if (!strcmp(str, "SEP"))
    {
        *val = SRT_SEP;
        return 0;
    }
    if (!strcmp(str, "SEPTEMBER"))
    {
        *val = SRT_SEP;
        return 0;
    }
    if (!strcmp(str, "OCT"))
    {
        *val = SRT_OCT;
        return 0;
    }
    if (!strcmp(str, "OCTOBER"))
    {
        *val = SRT_OCT;
        return 0;
    }
    if (!strcmp(str, "NOV"))
    {
        *val = SRT_NOV;
        return 0;
    }
    if (!strcmp(str, "NOVEMBER"))
    {
        *val = SRT_NOV;
        return 0;
    }
    if (!strcmp(str, "DEC"))
    {
        *val = SRT_DEC;
        return 0;
    }
    if (!strcmp(str, "DECEMBER"))
    {
        *val = SRT_DEC;
        return 0;
    }

    return serror("unknown month %s", str);
}

Err translate_month(char** str, SrtMonth val)
{
    /* get month */
    if (val == SRT_JAN)
    {
        *str = "JAN";
        return 0;
    }
    if (val == SRT_FEB)
    {
        *str = "FEB";
        return 0;
    }
    if (val == SRT_MAR)
    {
        *str = "MAR";
        return 0;
    }
    if (val == SRT_APR)
    {
        *str = "APR";
        return 0;
    }
    if (val == SRT_MAY)
    {
        *str = "MAY";
        return 0;
    }
    if (val == SRT_JUN)
    {
        *str = "JUN";
        return 0;
    }
    if (val == SRT_JUL)
    {
        *str = "JUL";
        return 0;
    }
    if (val == SRT_AUG)
    {
        *str = "AUG";
        return 0;
    }
    if (val == SRT_SEP)
    {
        *str = "SEP";
        return 0;
    }
    if (val == SRT_OCT)
    {
        *str = "OCT";
        return 0;
    }
    if (val == SRT_NOV)
    {
        *str = "NOV";
        return 0;
    }
    if (val == SRT_DEC)
    {
        *str = "DEC";
        return 0;
    }

    return serror("unknown month code %d", val);
}

Err interp_unit(const char* constStr, SrtUnit* val)
{
    char str[UTMAXSTR + 1];
    strncpy(str, constStr, UTMAXSTR);
    str[UTMAXSTR] = '\0';

    strupper(str);
    strip_white_space(str);

    if (!strcmp(str, "D"))
    {
        *val = SRT_DAY;
        return 0;
    }
    if (!strcmp(str, "DAY"))
    {
        *val = SRT_DAY;
        return 0;
    }
    if (!strcmp(str, "DAYS"))
    {
        *val = SRT_DAY;
        return 0;
    }

    if (!strcmp(str, "BD"))
    {
        *val = SRT_BDAY;
        return 0;
    }
    if (!strcmp(str, "B"))
    {
        *val = SRT_BDAY;
        return 0;
    }
    if (!strcmp(str, "BUSINESS DAY"))
    {
        *val = SRT_BDAY;
        return 0;
    }
    if (!strcmp(str, "BUSINESS DAYS"))
    {
        *val = SRT_BDAY;
        return 0;
    }

    if (!strcmp(str, "W"))
    {
        *val = SRT_WEEK;
        return 0;
    }
    if (!strcmp(str, "WEEK"))
    {
        *val = SRT_WEEK;
        return 0;
    }
    if (!strcmp(str, "WEEKS"))
    {
        *val = SRT_WEEK;
        return 0;
    }

    if (!strcmp(str, "M"))
    {
        *val = SRT_MONTH;
        return 0;
    }
    if (!strcmp(str, "MONTH"))
    {
        *val = SRT_MONTH;
        return 0;
    }
    if (!strcmp(str, "MONTHS"))
    {
        *val = SRT_MONTH;
        return 0;
    }

    if (!strcmp(str, "Y"))
    {
        *val = SRT_YEAR;
        return 0;
    }
    if (!strcmp(str, "YEAR"))
    {
        *val = SRT_YEAR;
        return 0;
    }
    if (!strcmp(str, "YEARS"))
    {
        *val = SRT_YEAR;
        return 0;
    }

    return serror("unknown unit: %s", str);
}

Err interp_bus_day_conv(String str, SrtBusDayConv* val)
{
    if (!strcmp(str, "MODSUCCEEDING"))
    {
        *val = MODIFIED_SUCCEEDING;
        return 0;
    }
    if (!strcmp(str, "SUCCEEDING"))
    {
        *val = SUCCEEDING;
        return 0;
    }

    if (!strcmp(str, "MODFOLLOWING"))
    {
        *val = MODIFIED_SUCCEEDING;
        return 0;
    }
    if (!strcmp(str, "FOLLOWING"))
    {
        *val = SUCCEEDING;
        return 0;
    }

    if (!strcmp(str, "NONE"))
    {
        *val = NO_BUSDAY_CONVENTION;
        return 0;
    }
    if (!strcmp(str, "NOBUSDAY"))
    {
        *val = NO_BUSDAY_CONVENTION;
        return 0;
    }

    return serror("unknown busday convention:%s", str);
}

Err translate_bus_day_conv(String* str, SrtBusDayConv m)
{
    switch (m)
    {
    case NO_BUSDAY_CONVENTION:
    {
        *str = "NONE";
        break;
    }
    case SUCCEEDING:
    {
        *str = "SUCCEEDING";
        break;
    }
    case MODIFIED_SUCCEEDING:
    {
        *str = "MODSUCCEEDING";
        break;
    }
    default:
        return "translate bus day conv: unknown bus_day_conv.";
    }
    return 0;
}

/* ---------------------------------------------------------------------------
             SWAP CONVENTIONS: BASIS, COMPOUNDING
   --------------------------------------------------------------------------- */

Err interp_basis(const char* constStr, SrtBasisCode* val)
{
    char str[UTMAXSTR + 1];
    strncpy(str, constStr, UTMAXSTR);
    str[UTMAXSTR] = '\0';

    strupper(str);
    strip_white_space(str);

    if (!strcmp(str, "ACT/ACT"))
    {
        *val = BASIS_ACT_ACT;
        return 0;
    }
    if (!strcmp(str, "ACT/365"))
    {
        *val = BASIS_ACT_365;
        return 0;
    }
    if (!strcmp(str, "ACT/360"))
    {
        *val = BASIS_ACT_360;
        return 0;
    }
    if (!strcmp(str, "30/360"))
    {
        *val = BASIS_30_360;
        return 0;
    }
    if (!strcmp(str, "30/360E"))
    {
        *val = BASIS_30_360E;
        return 0;
    }
    if (!strcmp(str, "30E/360"))
    {
        *val = BASIS_30_360E;
        return 0;
    }
    if (!strcmp(str, "ACT/USD"))
    {
        *val = BASIS_ACT_USD;
        return 0;
    }

    if (!strcmp(str, "AA"))
    {
        *val = BASIS_ACT_ACT;
        return 0;
    }
    if (!strcmp(str, "A5"))
    {
        *val = BASIS_ACT_365;
        return 0;
    }
    if (!strcmp(str, "MM"))
    {
        *val = BASIS_ACT_360;
        return 0;
    }
    if (!strcmp(str, "BB"))
    {
        *val = BASIS_30_360;
        return 0;
    }
    if (!strcmp(str, "BE"))
    {
        *val = BASIS_30_360E;
        return 0;
    }
    if (!strcmp(str, "US"))
    {
        *val = BASIS_ACT_USD;
        return 0;
    }

    return serror("unknown basis: %s", str);
}

Err translate_basis(char** str, SrtBasisCode val)
{
    if (val == BASIS_ACT_ACT)
    {
        *str = "ACT/ACT";
        return 0;
    }
    if (val == BASIS_ACT_365)
    {
        *str = "ACT/365";
        return 0;
    }
    if (val == BASIS_ACT_360)
    {
        *str = "ACT/360";
        return 0;
    }
    if (val == BASIS_30_360)
    {
        *str = "30/360";
        return 0;
    }
    if (val == BASIS_30_360E)
    {
        *str = "30/360E";
        return 0;
    }
    if (val == BASIS_ACT_USD)
    {
        *str = "ACT/USD";
        return 0;
    }

    return "translate_basis: unknown basis";
}

Err interp_compounding(const char* constStr, SrtCompounding* val)
{
    char str[UTMAXSTR + 1];

    if (!constStr)
        return "Empty string in interp_compounding";

    strncpy(str, constStr, UTMAXSTR);
    str[UTMAXSTR] = '\0';

    strupper(str);
    strip_white_space(str);

    if (!strcmp(str, "ANNUAL"))
    {
        *val = SRT_ANNUAL;
        return 0;
    }
    if (!strcmp(str, "ANN"))
    {
        *val = SRT_ANNUAL;
        return 0;
    }
    if (!strcmp(str, "SEMI"))
    {
        *val = SRT_SEMIANNUAL;
        return 0;
    }
    if (!strcmp(str, "SEMIANNUAL"))
    {
        *val = SRT_SEMIANNUAL;
        return 0;
    }
    if (!strcmp(str, "SA"))
    {
        *val = SRT_SEMIANNUAL;
        return 0;
    }
    if (!strcmp(str, "QUARTERLY"))
    {
        *val = SRT_QUARTERLY;
        return 0;
    }
    if (!strcmp(str, "MONTHLY"))
    {
        *val = SRT_MONTHLY;
        return 0;
    }
    if (!strcmp(str, "SIMPLE"))
    {
        *val = SRT_SIMPLE;
        return 0;
    }

    if (!strcmp(str, "A"))
    {
        *val = SRT_ANNUAL;
        return 0;
    }
    if (!strcmp(str, "S"))
    {
        *val = SRT_SEMIANNUAL;
        return 0;
    }
    if (!strcmp(str, "Q"))
    {
        *val = SRT_QUARTERLY;
        return 0;
    }
    if (!strcmp(str, "M"))
    {
        *val = SRT_MONTHLY;
        return 0;
    }
    if (!strcmp(str, "N"))
    {
        *val = SRT_SIMPLE;
        return 0;
    }

    return serror("unknown compounding. %s", str);
}
Err translate_compounding(char** str, SrtCompounding val)
{
    if (val == SRT_SIMPLE)
    {
        *str = "SIMPLE";
        return 0;
    }
    if (val == SRT_ANNUAL)
    {
        *str = "ANNUAL";
        return 0;
    }
    if (val == SRT_SEMIANNUAL)
    {
        *str = "SEMIANNUAL";
        return 0;
    }
    if (val == SRT_QUARTERLY)
    {
        *str = "QUARTERLY";
        return 0;
    }
    if (val == SRT_MONTHLY)
    {
        *str = "MONTHLY";
        return 0;
    }

    return "translate_compounding unknown compounding.";
}

Err interp_interp_method(String str, SrtInterpMethod* val)
{
    if (!strcmp(str, "LIN_RT"))
    {
        *val = LIN_RT;
        return 0;
    }
    if (!strcmp(str, "LIN_R"))
    {
        *val = LIN_R;
        return 0;
    }

    return serror("unknown interp method: %s", str);
}

Err translate_interp_method(String* str, SrtInterpMethod val)
{
    if (val == LIN_RT)
    {
        *str = "LIN_RT";
        return 0;
    }
    if (val == LIN_R)
    {
        *str = "LIN_R";
        return 0;
    }

    return serror("unknown interp method: %d", val);
}

Err interp_rec_pay(const char* constStr, SrtReceiverType* val)
{
    char str[UTMAXSTR + 1];
    strncpy(str, constStr, UTMAXSTR);
    str[UTMAXSTR] = '\0';

    strupper(str);
    strip_white_space(str);

    if (!strcmp(str, "RECEIVER"))
    {
        *val = SRT_RECEIVER;
        return 0;
    }
    if (!strcmp(str, "REC"))
    {
        *val = SRT_RECEIVER;
        return 0;
    }
    if (!strcmp(str, "RECV"))
    {
        *val = SRT_RECEIVER;
        return 0;
    }
    if (!strcmp(str, "RECEIVERS"))
    {
        *val = SRT_RECEIVER;
        return 0;
    }
    if (!strcmp(str, "FLOOR"))
    {
        *val = SRT_RECEIVER;
        return 0;
    }
    if (!strcmp(str, "PAYER"))
    {
        *val = SRT_PAYER;
        return 0;
    }
    if (!strcmp(str, "PAY"))
    {
        *val = SRT_PAYER;
        return 0;
    }
    if (!strcmp(str, "PAYERS"))
    {
        *val = SRT_PAYER;
        return 0;
    }
    if (!strcmp(str, "CAP"))
    {
        *val = SRT_PAYER;
        return 0;
    }
    if (!strcmp(str, "STRADDLE"))
    {
        *val = SRT_STRADDLE;
        return 0;
    }
    if (!strcmp(str, "RATE"))
    {
        *val = SRT_FORWARD;
        return 0;
    }
    if (!strcmp(str, "FORWARD"))
    {
        *val = SRT_FORWARD;
        return 0;
    }
    if (!strcmp(str, "FWD"))
    {
        *val = SRT_FORWARD;
        return 0;
    }

    return "Enter REC (FLOOR) or PAY (CAP) or RATE or STRADDLE for swaption (capfloor).";
}

Err interp_shift(const char* constStr, SrtShiftType* val)
{
    char str[UTMAXSTR + 1];
    strncpy(str, constStr, UTMAXSTR);
    str[UTMAXSTR] = '\0';

    strupper(str);
    strip_white_space(str);

    if (!strcmp(str, "VOL"))
    {
        *val = SRT_SHIFTVOL;
        return 0;
    }
    if (!strcmp(str, "YC"))
    {
        *val = SRT_SHIFTYC;
        return 0;
    }

    return "Enter VOL or YC.";
}
/* ---------------------------------------------------------------------------
             OPTIONS : CALL/PUT, UP/DOWN...
   --------------------------------------------------------------------------- */

Err interp_call_put(const char* constStr, SrtCallPutType* val)
{
    char str[UTMAXSTR + 1];
    strncpy(str, constStr, UTMAXSTR);
    str[UTMAXSTR] = '\0';

    strupper(str);
    strip_white_space(str);

    if (!strcmp(str, "CAP"))
    {
        *val = SRT_CALL;
        return 0;
    }
    if (!strcmp(str, "CALL"))
    {
        *val = SRT_CALL;
        return 0;
    }
    if (!strcmp(str, "C"))
    {
        *val = SRT_CALL;
        return 0;
    }
    if (!strcmp(str, "FLOOR"))
    {
        *val = SRT_PUT;
        return 0;
    }
    if (!strcmp(str, "PUT"))
    {
        *val = SRT_PUT;
        return 0;
    }
    if (!strcmp(str, "P"))
    {
        *val = SRT_PUT;
        return 0;
    }
    if (!strcmp(str, "STRADDLE"))
    {
        *val = SRT_STRADDLE;
        return 0;
    }

    return "Enter C (CALL, CAP) or P (PUT, FLOOR) or STRADDLE.";
}

Err conv_rec_in_call(SrtReceiverType rec_pay, SrtCallPutType* call_put)
{
    if (rec_pay == SRT_RECEIVER)
    {
        *call_put = SRT_PUT;
        return 0;
    }
    if (rec_pay == SRT_PAYER)
    {
        *call_put = SRT_CALL;
        return 0;
    }

    return "Cannot convert SrtReceiverType in SrtCallPutType.";
}

Err interp_diffusion_type(const char* constStr, SrtDiffusionType* val)
{
    char str[UTMAXSTR + 1];
    strncpy(str, constStr, UTMAXSTR);
    str[UTMAXSTR] = '\0';

    strupper(str);
    strip_white_space(str);

    if (!strcmp(str, "LOGNORMAL"))
    {
        *val = SRT_LOGNORMAL;
        return 0;
    }
    if (!strcmp(str, "L"))
    {
        *val = SRT_LOGNORMAL;
        return 0;
    }
    if (!strcmp(str, "LOG"))
    {
        *val = SRT_LOGNORMAL;
        return 0;
    }
    if (!strcmp(str, "0"))
    {
        *val = SRT_LOGNORMAL;
        return 0;
    }
    if (!strcmp(str, "NORMAL"))
    {
        *val = SRT_NORMAL;
        return 0;
    }
    if (!strcmp(str, "N"))
    {
        *val = SRT_NORMAL;
        return 0;
    }
    if (!strcmp(str, "NORM"))
    {
        *val = SRT_NORMAL;
        return 0;
    }
    if (!strcmp(str, "1"))
    {
        *val = SRT_NORMAL;
        return 0;
    }
    if (!strcmp(str, "BETA"))
    {
        *val = SRT_BETAVOL;
        return 0;
    }
    if (!strcmp(str, "BETAVOL"))
    {
        *val = SRT_BETAVOL;
        return 0;
    }
    if (!strcmp(str, "BV"))
    {
        *val = SRT_BETAVOL;
        return 0;
    }
    if (!strcmp(str, "2"))
    {
        *val = SRT_BETAVOL;
        return 0;
    }
    if (!strcmp(str, "HESTON"))
    {
        *val = SRT_HESTONVOL;
        return 0;
    }
    if (!strcmp(str, "H"))
    {
        *val = SRT_HESTONVOL;
        return 0;
    }
    if (!strcmp(str, "HV"))
    {
        *val = SRT_HESTONVOL;
        return 0;
    }
    if (!strcmp(str, "3"))
    {
        *val = SRT_HESTONVOL;
        return 0;
    }

    if (!strcmp(str, "SABR"))
    {
        *val = SRT_SABRVOL;
        return 0;
    }
    if (!strcmp(str, "BVM"))
    {
        *val = SRT_BVMVOL;
        return 0;
    }
    if (!strcmp(str, "BVM2"))
    {
        *val = SRT_BVM2VOL;
        return 0;
    }
    if (!strcmp(str, "BVMH"))
    {
        *val = SRT_BVMHVOL;
        return 0;
    }
    if (!strcmp(str, "BVMH2"))
    {
        *val = SRT_BVMH2VOL;
        return 0;
    }
    if (!strcmp(str, "BVMC"))
    {
        *val = SRT_BVMCVOL;
        return 0;
    }
    if (!strcmp(str, "BMM"))
    {
        *val = SRT_BMMVOL;
        return 0;
    }
    if (!strcmp(str, "BMM2"))
    {
        *val = SRT_BMM2VOL;
        return 0;
    }
    if (!strcmp(str, "BMM3"))
    {
        *val = SRT_BMM3VOL;
        return 0;
    }
    if (!strcmp(str, "BMMBVM"))
    {
        *val = SRT_BMMBVMVOL;
        return 0;
    }
    if (!strcmp(str, "BMMBVM2"))
    {
        *val = SRT_BMMBVM2VOL;
        return 0;
    }
    if (!strcmp(str, "BMMBVMH"))
    {
        *val = SRT_BMMBVMHVOL;
        return 0;
    }
    if (!strcmp(str, "BMMBVMH2"))
    {
        *val = SRT_BMMBVMH2VOL;
        return 0;
    }
    if (!strcmp(str, "BMMBVMC"))
    {
        *val = SRT_BMMBVMCVOL;
        return 0;
    }
    if (!strcmp(str, "BMMCAL"))
    {
        *val = SRT_BMMCALVOL;
        return 0;
    }
    if (!strcmp(str, "BMMCAL2"))
    {
        *val = SRT_BMM2CALVOL;
        return 0;
    }
    if (!strcmp(str, "BMMCAL3"))
    {
        *val = SRT_BMM3CALVOL;
        return 0;
    }

    return "enter L  (LOGNORMAL) or N (NORMAL).";
}

Err interp_barrier_type(const char* constStr, SrtBarrierType* val)
{
    char str[UTMAXSTR + 1];
    strncpy(str, constStr, UTMAXSTR);
    str[UTMAXSTR] = '\0';

    strupper(str);
    strip_white_space(str);

    if (!strcmp(str, "DOWN"))
    {
        *val = SRT_DOWN;
        return 0;
    }
    if (!strcmp(str, "D"))
    {
        *val = SRT_DOWN;
        return 0;
    }
    if (!strcmp(str, "BELOW"))
    {
        *val = SRT_DOWN;
        return 0;
    }
    if (!strcmp(str, "B"))
    {
        *val = SRT_DOWN;
        return 0;
    }
    if (!strcmp(str, "UP"))
    {
        *val = SRT_UP;
        return 0;
    }
    if (!strcmp(str, "U"))
    {
        *val = SRT_UP;
        return 0;
    }
    if (!strcmp(str, "ABOVE"))
    {
        *val = SRT_UP;
        return 0;
    }
    if (!strcmp(str, "A"))
    {
        *val = SRT_UP;
        return 0;
    }

    return "enter D (DOWN),  or U (UP).";
}

Err interp_min_max(const char* constStr, SrtMinmaxType* val)
{
    char str[UTMAXSTR + 1];
    strncpy(str, constStr, UTMAXSTR);
    str[UTMAXSTR] = '\0';

    strupper(str);
    strip_white_space(str);

    if (!strcmp(str, "MINIMUM"))
    {
        *val = SRT_MIN;
        return 0;
    }
    if (!strcmp(str, "MIN"))
    {
        *val = SRT_MIN;
        return 0;
    }
    if (!strcmp(str, "MAXIMUM"))
    {
        *val = SRT_MAX;
        return 0;
    }
    if (!strcmp(str, "MAX"))
    {
        *val = SRT_MAX;
        return 0;
    }

    return "enter MIN or MAX ";
}

Err interp_grad_best(const char* constStr, SrtLadderType* val)
{
    char str[UTMAXSTR + 1];
    strncpy(str, constStr, UTMAXSTR);
    str[UTMAXSTR] = '\0';

    strupper(str);
    strip_white_space(str);

    if (!strcmp(str, "GRADUATED"))
    {
        *val = SRT_GRADUATED;
        return 0;
    }
    if (!strcmp(str, "GRAD"))
    {
        *val = SRT_GRADUATED;
        return 0;
    }
    if (!strcmp(str, "G"))
    {
        *val = SRT_GRADUATED;
        return 0;
    }
    if (!strcmp(str, "BEST"))
    {
        *val = SRT_BEST;
        return 0;
    }
    if (!strcmp(str, "B"))
    {
        *val = SRT_BEST;
        return 0;
    }

    return "enter G (GRAD) or B (BEST).";
}

Err interp_best_worst(const char* constStr, SrtBestWorstType* val)
{
    char str[UTMAXSTR + 1];
    strncpy(str, constStr, UTMAXSTR);
    str[UTMAXSTR] = '\0';

    strupper(str);
    strip_white_space(str);

    if (!strcmp(str, "BEST"))
    {
        *val = SRT_BESTOFF;
        return 0;
    }
    if (!strcmp(str, "WORST"))
    {
        *val = SRT_WORSTOFF;
        return 0;
    }
    if (!strcmp(str, "B"))
    {
        *val = SRT_BESTOFF;
        return 0;
    }
    if (!strcmp(str, "W"))
    {
        *val = SRT_WORSTOFF;
        return 0;
    }

    return "enter B (BEST) or WORST (WORST)";
}

/* ---------------------------------------------------------------------------
                             GREEKS
   --------------------------------------------------------------------------- */

Err interp_greeks(const char* constStr, SrtGreekType* val)
{
    char str[UTMAXSTR + 1];
    strncpy(str, constStr, UTMAXSTR);
    str[UTMAXSTR] = '\0';

    strupper(str);
    strip_white_space(str);

    if (!strcmp(str, "PREMIUM"))
    {
        *val = PREMIUM;
        return 0;
    }
    if (!strcmp(str, "PREM"))
    {
        *val = PREMIUM;
        return 0;
    }
    if (!strcmp(str, "PV"))
    {
        *val = PREMIUM;
        return 0;
    }
    if (!strcmp(str, "DELTA"))
    {
        *val = DELTA;
        return 0;
    }
    if (!strcmp(str, "DELTAX"))
    {
        *val = DELTAX;
        return 0;
    }
    if (!strcmp(str, "DELTAY"))
    {
        *val = DELTAY;
        return 0;
    }
    if (!strcmp(str, "DELTAZ"))
    {
        *val = DELTAZ;
        return 0;
    }
    if (!strcmp(str, "DELTA_FWD"))
    {
        *val = DELTA_FWD;
        return 0;
    }
    if (!strcmp(str, "DELTA_FWD1"))
    {
        *val = DELTA_FWD1;
        return 0;
    }
    if (!strcmp(str, "DELTA_FWD2"))
    {
        *val = DELTA_FWD2;
        return 0;
    }
    if (!strcmp(str, "GAMMA"))
    {
        *val = GAMMA;
        return 0;
    }
    if (!strcmp(str, "GAMMAX"))
    {
        *val = GAMMAX;
        return 0;
    }
    if (!strcmp(str, "GAMMAY"))
    {
        *val = GAMMAY;
        return 0;
    }
    if (!strcmp(str, "GAMMAZ"))
    {
        *val = GAMMAZ;
        return 0;
    }
    if (!strcmp(str, "GAMMAXY"))
    {
        *val = GAMMAXY;
        return 0;
    }
    if (!strcmp(str, "GAMMAXZ"))
    {
        *val = GAMMAXZ;
        return 0;
    }
    if (!strcmp(str, "GAMMAYZ"))
    {
        *val = GAMMAYZ;
        return 0;
    }
    if (!strcmp(str, "GAMMA_FWD"))
    {
        *val = GAMMA_FWD;
        return 0;
    }
    if (!strcmp(str, "GAMMA_FWD1"))
    {
        *val = GAMMA_FWD1;
        return 0;
    }
    if (!strcmp(str, "GAMMA_FWD2"))
    {
        *val = GAMMA_FWD2;
        return 0;
    }
    if (!strcmp(str, "THETA"))
    {
        *val = THETA;
        return 0;
    }
    if (!strcmp(str, "THETA_1D"))
    {
        *val = THETA_1D;
        return 0;
    }
    if (!strcmp(str, "THETAGAMMA"))
    {
        *val = THETAGAMMA;
        return 0;
    }
    if (!strcmp(str, "THETAVOLGA"))
    {
        *val = THETAVOLGA;
        return 0;
    }
    if (!strcmp(str, "THETAVANNA"))
    {
        *val = THETAVANNA;
        return 0;
    }
    if (!strcmp(str, "THETACARRY"))
    {
        *val = THETACARRY;
        return 0;
    }
    if (!strcmp(str, "VEGA"))
    {
        *val = VEGA;
        return 0;
    }
    if (!strcmp(str, "VEGAX"))
    {
        *val = VEGAX;
        return 0;
    }
    if (!strcmp(str, "VEGAY"))
    {
        *val = VEGAY;
        return 0;
    }
    if (!strcmp(str, "VEGAZ"))
    {
        *val = VEGAZ;
        return 0;
    }
    if (!strcmp(str, "VEGAXY"))
    {
        *val = VEGAXY;
        return 0;
    }
    if (!strcmp(str, "VEGA1"))
    {
        *val = VEGA1;
        return 0;
    }
    if (!strcmp(str, "VEGA2"))
    {
        *val = VEGA2;
        return 0;
    }
    if (!strcmp(str, "VANNA"))
    {
        *val = VANNA;
        return 0;
    }
    if (!strcmp(str, "VOLGA"))
    {
        *val = VOLGA;
        return 0;
    }
    if (!strcmp(str, "RHO"))
    {
        *val = RHO;
        return 0;
    }
    if (!strcmp(str, "ALPHA"))
    {
        *val = ALPHA;
        return 0;
    }
    if (!strcmp(str, "BETA"))
    {
        *val = BETA;
        return 0;
    }
    if (!strcmp(str, "DU1"))
    {
        *val = DU1;
        return 0;
    }
    if (!strcmp(str, "DLAMBDA1"))
    {
        *val = DLAMBDA1;
        return 0;
    }
    if (!strcmp(str, "DU2"))
    {
        *val = DU2;
        return 0;
    }
    if (!strcmp(str, "DLAMBDA2"))
    {
        *val = DLAMBDA2;
        return 0;
    }
    if (!strcmp(str, "RHORHO"))
    {
        *val = RHORHO;
        return 0;
    }
    if (!strcmp(str, "RHODELTA"))
    {
        *val = RHODELTA;
        return 0;
    }
    if (!strcmp(str, "DENSITY"))
    {
        *val = DENSITY;
        return 0;
    }
    if (!strcmp(str, "CUMDENSITY"))
    {
        *val = CUMDENSITY;
        return 0;
    }
    return serror("unknown greek: %s", str);
}

/* ---------------------------------------------------------------------------
                              RESET
   --------------------------------------------------------------------------- */
Err interp_reset_optimised(const char* constStr, SrtResOptType* val)
{
    char str[UTMAXSTR + 1];
    strncpy(str, constStr, UTMAXSTR);
    str[UTMAXSTR] = '\0';

    strupper(str);
    strip_white_space(str);

    if (!strcmp(str, "RESET"))
    {
        *val = SRT_AUTORESET;
        return 0;
    }
    if (!strcmp(str, "OPTIMISED"))
    {
        *val = SRT_OPTIMISED;
        return 0;
    }

    return "Enter reset or optimised";
}

Err interp_price_type(const char* constStr, SrtPriceType* val)
{
    char str[UTMAXSTR + 1];
    strncpy(str, constStr, UTMAXSTR);
    str[UTMAXSTR] = '\0';

    strupper(str);
    strip_white_space(str);

    if (!strcmp(str, "PREMIUM"))
    {
        *val = SRT_PREMIUM;
        return 0;
    }
    if (!strcmp(str, "VOLATILITY"))
    {
        *val = SRT_VOLATILITY;
        return 0;
    }

    return "Enter PREMIUM or VOLATILTY";
}

/* ---------------------------------------------------------------------------
                          FOR NUMERICAL METHODS
   --------------------------------------------------------------------------- */
Err interp_derivatives(const char* constStr, SrtDerType* val)
{
    char str[UTMAXSTR + 1];
    strncpy(str, constStr, UTMAXSTR);
    str[UTMAXSTR] = '\0';

    strupper(str);
    strip_white_space(str);

    if (!strcmp(str, "FIRSTDER"))
    {
        *val = SRT_FIRSTDER;
        return 0;
    }
    if (!strcmp(str, "SECONDDER"))
    {
        *val = SRT_SECONDDER;
        return 0;
    }
    if (!strcmp(str, "TERDER"))
    {
        *val = SRT_TERDER;
        return 0;
    }

    return "Enter FIRSTDER or SECONDDER in c code";
}

/* ---------------------------------------------------------------------------
                               UNDERLYINGS
   --------------------------------------------------------------------------- */

Err interp_dom_for_type(const char* constStr, SrtDomForType* val)
{
    char str[UTMAXSTR + 1];
    strncpy(str, constStr, UTMAXSTR);
    str[UTMAXSTR] = '\0';

    strupper(str);
    strip_white_space(str);

    if (!strcmp(str, "DOMESTIC"))
    {
        *val = SRT_DOMESTIC;
        return 0;
    }
    if (!strcmp(str, "D"))
    {
        *val = SRT_DOMESTIC;
        return 0;
    }
    if (!strcmp(str, "DOM"))
    {
        *val = SRT_DOMESTIC;
        return 0;
    }
    if (!strcmp(str, "0"))
    {
        *val = SRT_DOMESTIC;
        return 0;
    }
    if (!strcmp(str, "FOREIGN"))
    {
        *val = SRT_FOREIGN;
        return 0;
    }
    if (!strcmp(str, "F"))
    {
        *val = SRT_FOREIGN;
        return 0;
    }
    if (!strcmp(str, "FOR"))
    {
        *val = SRT_FOREIGN;
        return 0;
    }
    if (!strcmp(str, "1"))
    {
        *val = SRT_FOREIGN;
        return 0;
    }

    return "enter (D)OMESTIC or (F)OREIGN.";
}

/* ---------------------------------------------------------------------------
                               SABRComponent
   --------------------------------------------------------------------------- */

Err interp_SABRVolComponent(const char* constStr, SABRVolComponent* val)
{
    char str[UTMAXSTR + 1];
    strncpy(str, constStr, UTMAXSTR);
    str[UTMAXSTR] = '\0';

    strupper(str);
    strip_white_space(str);

    if (!strcmp(str, "ALPHA"))
    {
        *val = SABR_ALPHA;
        return 0;
    }
    if (!strcmp(str, "WING"))
    {
        *val = SABR_ALPHA;
        return 0;
    }
    if (!strcmp(str, "BETA"))
    {
        *val = SABR_BETA;
        return 0;
    }
    if (!strcmp(str, "SKEW"))
    {
        *val = SABR_BETA;
        return 0;
    }
    if (!strcmp(str, "RHO"))
    {
        *val = SABR_RHO;
        return 0;
    }

    if (!strcmp(str, "LOGVOL"))
    {
        *val = SABR_LOGVOL;
        return 0;
    }
    if (!strcmp(str, "VOLLOG"))
    {
        *val = SABR_LOGVOL;
        return 0;
    }
    if (!strcmp(str, "LOG"))
    {
        *val = SABR_LOGVOL;
        return 0;
    }
    if (!strcmp(str, "LOGNORMAL"))
    {
        *val = SABR_LOGVOL;
        return 0;
    }

    if (!strcmp(str, "NORMVOL"))
    {
        *val = SABR_NORMVOL;
        return 0;
    }
    if (!strcmp(str, "VOLNORM"))
    {
        *val = SABR_NORMVOL;
        return 0;
    }
    if (!strcmp(str, "NORM"))
    {
        *val = SABR_NORMVOL;
        return 0;
    }
    if (!strcmp(str, "NORMAL"))
    {
        *val = SABR_NORMVOL;
        return 0;
    }

    if (!strcmp(str, "BETAVOL"))
    {
        *val = SABR_BETAVOL;
        return 0;
    }
    if (!strcmp(str, "VOLBETA"))
    {
        *val = SABR_BETAVOL;
        return 0;
    }

    if (!strcmp(str, "ATMLOG"))
    {
        *val = SABR_ATMLOG;
        return 0;
    }
    if (!strcmp(str, "ATMLOGNORMAL"))
    {
        *val = SABR_ATMLOG;
        return 0;
    }
    if (!strcmp(str, "ATMLOGVOL"))
    {
        *val = SABR_ATMLOG;
        return 0;
    }
    /* Warning ATMVOL would be understood as ATMLOG */
    if (!strcmp(str, "ATMVOL"))
    {
        *val = SABR_ATMLOG;
        return 0;
    }

    if (!strcmp(str, "ATMNORM"))
    {
        *val = SABR_ATMNORM;
        return 0;
    }
    if (!strcmp(str, "ATMNORMAL"))
    {
        *val = SABR_ATMNORM;
        return 0;
    }
    if (!strcmp(str, "ATMNORMVOL"))
    {
        *val = SABR_ATMNORM;
        return 0;
    }

    return "enter ALPHA, BETA, RHO, ATMLOG, ATMNORM, BETAVOL, LOGVOL, NORMVOL";
}

/* ============================================================================= */

/* ============================================================================= */

/*===========================FOR DELTA REPORT=================================== */

Err interp_HedgeType(const char* constStr, SrtHedgeType* val)
{
    char str[UTMAXSTR + 1];
    strncpy(str, constStr, UTMAXSTR);
    str[UTMAXSTR] = '\0';

    strupper(str);
    strip_white_space(str);

    if (!strcmp(str, "SWAP"))
    {
        *val = SRT_SWAP;
        return 0;
    }
    if (!strcmp(str, "FRA"))
    {
        *val = SRT_FRA;
        return 0;
    }

    return "enter (S)SWAP or (F)FRA";
}

Err interp_UndFRAType(const char* constStr, SrtUndFRAType* val)
{
    char str[UTMAXSTR + 1];
    strncpy(str, constStr, UTMAXSTR);
    str[UTMAXSTR] = '\0';

    strupper(str);
    strip_white_space(str);

    if (!strcmp(str, "SWAP"))
    {
        *val = SRT_UNDFRA_SWAP;
        return 0;
    }
    if (!strcmp(str, "FUT"))
    {
        *val = SRT_UNDFRA_FUT;
        return 0;
    }

    return "enter (S)SWAP or (F)FUT";
}

Err interp_inflation_indexed_bond_asset_swap(
    const char* constStr, SrtInflationIdxBondAssetSwapType* val)
{
    char str[UTMAXSTR + 1];
    strncpy(str, constStr, UTMAXSTR);
    str[UTMAXSTR] = '\0';

    strupper(str);
    strip_white_space(str);

    if (!strcmp(str, "FLAT_NOTIONAL"))
    {
        *val = SRT_FLAT_NOTIONAL;
        return 0;
    }
    if (!strcmp(str, "FLAT_INDEXED_NOTIONAL"))
    {
        *val = SRT_FLAT_INDEXED_NOTIONAL;
        return 0;
    }
    if (!strcmp(str, "INDEXED_NOTIONAL"))
    {
        *val = SRT_INDEXED_NOTIONAL;
        return 0;
    }

    return "enter FLAT_NOTIONAL or  FLAT_INDEXED_NOTIONAL or INDEXED_NOTIONAL";
}

Err interp_BoolType(const char* constStr, SRT_Boolean* val)
{
    char str[UTMAXSTR + 1];
    strncpy(str, constStr, UTMAXSTR);
    str[UTMAXSTR] = '\0';

    strupper(str);
    strip_white_space(str);

    if (!strcmp(str, "YES"))
    {
        *val = SRT_YES;
        return 0;
    }
    if (!strcmp(str, "Y"))
    {
        *val = SRT_YES;
        return 0;
    }
    if (!strcmp(str, "NO"))
    {
        *val = SRT_NO;
        return 0;
    }
    if (!strcmp(str, "N"))
    {
        *val = SRT_NO;
        return 0;
    }

    return "enter (Y)ES or (N)O";
}

Err interp_CalibrationType(const char* constStr, SrtCalibrationType* val)
{
    char str[UTMAXSTR + 1];
    strncpy(str, constStr, UTMAXSTR);
    str[UTMAXSTR] = '\0';

    strupper(str);
    strip_white_space(str);

    if (!strcmp(str, "MATCH_CONV"))
    {
        *val = MATCH_CONV;
        return 0;
    }
    if (!strcmp(str, "CHI2_MIN"))
    {
        *val = CHI2_MIN;
        return 0;
    }
    if (!strcmp(str, "CONV"))
    {
        *val = MATCH_CONV;
        return 0;
    }
    if (!strcmp(str, "CHI2"))
    {
        *val = CHI2_MIN;
        return 0;
    }

    return "enter (CONV) or (CHI2)";
}