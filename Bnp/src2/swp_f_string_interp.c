/* ================================================================================

   FILENAME:      swp_f_string_interp.c

   PURPOSE:       string interpretation functions for swaps parameters

   ================================================================================ */

#include "swp_h_all.h"

Err interp_swap_message(String str, Message* mess)
{
    strupper(str);
    strip_white_space(str);

    if (!strcmp(str, "PV"))
    {
        *mess = COMPUTE_PV;
        return 0;
    }
    if (!strcmp(str, "FWD_PV"))
    {
        *mess = COMPUTE_FWD_PV;
        return 0;
    }
    if (!strcmp(str, "FIXED_RATE"))
    {
        *mess = COMPUTE_FWD_RATE;
        return 0;
    }
    if (!strcmp(str, "FWD_RATE"))
    {
        *mess = COMPUTE_FWD_RATE;
        return 0;
    }
    if (!strcmp(str, "LEVEL"))
    {
        *mess = COMPUTE_LEVEL;
        return 0;
    }
    if (!strcmp(str, "IRR"))
    {
        *mess = COMPUTE_IRR;
        return 0;
    }
    if (!strcmp(str, "FWD_IRR"))
    {
        *mess = COMPUTE_FWD_IRR;
        return 0;
    }
    if (!strcmp(str, "SPREAD"))
    {
        *mess = COMPUTE_FWD_SPREAD;
        return 0;
    }
    if (!strcmp(str, "FWD_SPREAD"))
    {
        *mess = COMPUTE_FWD_SPREAD;
        return 0;
    }
    return serror("unknown swap_message: %s", str);
}

Err interp_struct(String str, Message* val)
{
    strupper(str);
    strip_white_space(str);

    if (!strcmp(str, "BOND_OPTION"))
    {
        *val = BOND_OPTION;
        return 0;
    }
    if (!strcmp(str, "BOND_OPT"))
    {
        *val = BOND_OPTION;
        return 0;
    }
    if (!strcmp(str, "BONDOPTION"))
    {
        *val = BOND_OPTION;
        return 0;
    }
    if (!strcmp(str, "BOND"))
    {
        *val = BOND_OPTION;
        return 0;
    }

    if (!strcmp(str, "SWAPTION"))
    {
        *val = SWAPTION;
        return 0;
    }

    ///// added by Albert Wang 08/25/03 - begin
    if (!strcmp(str, "SIMPLEMIDAT"))
    {
        *val = SIMPLEMIDAT;
        return 0;
    }
    ///// added by Albert Wang 08/25/03 - end

    if (!strcmp(str, "CAPFLOOR"))
    {
        *val = CAPFLOOR;
        return 0;
    }
    if (!strcmp(str, "CAP_FLOOR"))
    {
        *val = CAPFLOOR;
        return 0;
    }
    if (!strcmp(str, "C_F"))
    {
        *val = CAPFLOOR;
        return 0;
    }

    if (!strcmp(str, "RESETCAPFLOOR"))
    {
        *val = RESETCAPFLOOR;
        return 0;
    }
    if (!strcmp(str, "RESETCMSOPTION"))
    {
        *val = RESETCMSOPTION;
        return 0;
    }

    return serror("unknown structure: %s", str);
}

/************ interleaving and toys *********************/

Err interp_no_overwrite_flg(String cs, SRT_Boolean* b)
{
    if (!strcmp(cs, "OVERWRITE"))
    {
        *b = SRT_NO;
        return 0;
    }
    if (!strcmp(cs, "NOOVERWRITE"))
    {
        *b = SRT_YES;
        return 0;
    }
    return serror("toy or ilv no overwrite flag: don't know %s", cs);
}
Err translate_no_overwrite_flg(String* cs, SRT_Boolean b)
{
    if (b == SRT_NO)
    {
        *cs = "OVERWRITE";
        return 0;
    }
    else
    {
        *cs = "NOOVERWRITE";
        return 0;
    }
}
Err interp_insert_flg(String cs, SRT_Boolean* b)
{
    if (!strcmp(cs, "INSERT"))
    {
        *b = SRT_YES;
        return 0;
    }
    if (!strcmp(cs, "NOINSERT"))
    {
        *b = SRT_NO;
        return 0;
    }
    return serror("toy or ilv insert flag: don't know %s", cs);
}
Err translate_insert_flg(String* cs, SRT_Boolean b)
{
    if (b == SRT_NO)
    {
        *cs = "NOINSERT";
        return 0;
    }
    else
    {
        *cs = "INSERT";
        return 0;
    }
}
