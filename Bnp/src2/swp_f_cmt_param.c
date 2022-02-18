/******************************************************************************/
/*                                                                            */
/*      SYSTEM:         SRT     SORT, Fixed Income 2020 Addins                */
/*      SUB_SYSTEM:     SWT     Swap Tools                                    */
/*                                                                            */
/*      MODULE NAME:    SWP_F_CMT_PARAM                                       */
/*                                                                            */
/*      PURPOSE:                                                              */
/*                                                                            */
/*      AUTHORS:        Olivier VAN EYSEREN                   		      */
/*                                                                            */
/*      DATE:           23rd October 1995                                     */
/*                                                                            */
/*      VERSION:        01                                                    */
/*                                                                            */
/*      DESCRIPTION:    Contains default settings for each currency           */
/*                      currently supported                                   */
/*                                                                            */
/*                                                                            */
/******************************************************************************/

#include "math.h"
#include "swp_h_all.h"

CMT_Param_Struct* init_CMT_Param(
    String yc_name,
    String vc_name,
    String mkt_name,
    Err (*GetVol)(long, long, double, double*), /* volatility function for reference swptns */
    CMTCode CMT_code)
{
    CMT_Param_Struct* cps;
    SrtCurvePtr       yc_crv;
    SrtCcyParam*      ccy_param;

    yc_crv = lookup_curve(yc_name);
    if (!yc_crv)
        return NULL;

    ccy_param = (SrtCcyParam*)get_ccyparam_from_yldcrv(yc_crv);
    if (!ccy_param)
        return NULL;

    cps = new_CMT_Param();

    cps->cmt_code = CMT_code;
    cps->cmt_mat  = (long)fabs(CMT_code);

    /* initialise the ref rate with the currency */
    strcpy(cps->swap_ref_rate, get_curve_ccy(yc_crv));

    /* initialise market and vol function */
    strcpy(cps->yc_name, yc_name);
    strcpy(cps->vc_name, vc_name);
    strcpy(cps->mkt_name, mkt_name);
    cps->cmt_getvol = GetVol;
    cps->flatvol    = 0.;
    cps->voltype    = SRT_LOGNORMAL;

    cps->cmt_bus_day_conv = ccy_param->swap_bus_day_conv;
    cps->cms_bus_day_conv = ccy_param->swap_bus_day_conv;

    cps->cmt_basis_code  = BASIS_30_360;
    cps->bond_basis_code = BASIS_30_360;
    cps->cms_basis_code  = ccy_param->swap_basis_code;
    cps->swap_basis_code = ccy_param->swap_basis_code;

    cps->cmt_freq = SRT_QUARTERLY;
    cps->cms_freq = SRT_QUARTERLY;

    cps->swap_compd = ccy_param->compd;
    cps->bond_compd = SRT_SEMIANNUAL;

    cps->interp_method        = ccy_param->interp_method;
    cps->use_prop_vol_flg     = SRT_NO;
    cps->interp_on_spread_flg = SRT_YES;

    cps->spot_lag = ccy_param->spot_lag;

    return cps;
}

/* ------------------------------------------------------------------------- */

CMT_Param_Struct* new_CMT_Param()
{
    CMT_Param_Struct* cps;
    cps = (CMT_Param_Struct*)srt_calloc(1, sizeof(CMT_Param_Struct));
    if (!cps)
        return NULL;
    cps->allocated = 1;
    return cps;
}

/* ------------------------------------------------------------------------- */

int free_CMT_Param(CMT_Param_Struct* cps)
{
    if (cps->allocated)
        srt_free(cps);
    return 0;
}

/* ======================================================================== */

Err interp_vol_calc_method(String str, Message* val)
{
    if (!strcmp(str, "YES"))
    {
        *val = SRT_YES;
        return 0;
    }
    if (!strcmp(str, "NO"))
    {
        *val = SRT_NO;
        return 0;
    }

    return serror("Unknown answer to proportionnal vol %s <> YES || NO", val);
}

Err interp_interp_on_spread(String str, Message* val)
{
    if (!strcmp(str, "YES"))
    {
        *val = SRT_YES;
        return 0;
    }
    if (!strcmp(str, "NO"))
    {
        *val = SRT_NO;
        return 0;
    }

    return serror("Unknown answer to interp on spread %s <> YES || NO", val);
}

/* ------------------------------------------------------------------------- */

Err CMT_string_set_param(CMT_Param_Struct* cmt_param, String param_name, String param_val)
{
    Err     err;
    Message val;

    if (!strcmp(param_name, "CMTBUSDAYCONV"))
    {
        if (err = (interp_bus_day_conv(param_val, &val)))
            return err;
        cmt_param->cmt_bus_day_conv = val;
        return NULL;
    }
    if ((!strcmp(param_name, "CMTBASIS")) || (!strcmp(param_name, "TECBASIS")))
    {
        if (err = (interp_basis(param_val, &val)))
            return err;
        cmt_param->cmt_basis_code = val;
        return NULL;
    }
    if (!strcmp(param_name, "CMSBASIS"))
    {
        if (err = (interp_basis(param_val, &val)))
            return err;
        cmt_param->cms_basis_code = val;
        return NULL;
    }
    if (!strcmp(param_name, "SWAPBASIS"))
    {
        if (err = (interp_basis(param_val, &val)))
            return err;
        cmt_param->swap_basis_code = val;
        return NULL;
    }
    if (!strcmp(param_name, "BONDBASIS"))
    {
        if (err = (interp_basis(param_val, &val)))
            return err;
        cmt_param->bond_basis_code = val;
        return NULL;
    }
    if ((!strcmp(param_name, "CMTFREQUENCY")) || (!strcmp(param_name, "TECFREQUENCY")))
    {
        if (err = (interp_compounding(param_val, &val)))
            return err;
        cmt_param->cmt_freq = val;
        return NULL;
    }
    if (!strcmp(param_name, "CMSFREQUENCY"))
    {
        if (err = (interp_compounding(param_val, &val)))
            return err;
        cmt_param->cms_freq = val;
        return NULL;
    }
    if (!strcmp(param_name, "BONDCOMPOUNDING"))
    {
        if (err = (interp_compounding(param_val, &val)))
            return err;
        cmt_param->bond_compd = val;
        return NULL;
    }
    if (!strcmp(param_name, "SWAPCOMPOUNDING"))
    {
        if (err = (interp_compounding(param_val, &val)))
            return err;
        cmt_param->swap_compd = val;
        return NULL;
    }
    if (!strcmp(param_name, "INTERPMETHOD"))
    {
        if (err = (interp_interp_method(param_val, &val)))
            return err;
        cmt_param->interp_method = val;
        return NULL;
    }
    if (!strcmp(param_name, "USEPROPVOL"))
    {
        if (err = (interp_vol_calc_method(param_val, &val)))
            return err;
        cmt_param->use_prop_vol_flg = val;
        return NULL;
    }
    if (!strcmp(param_name, "INTERPONSPREAD"))
    {
        if (err = (interp_interp_on_spread(param_val, &val)))
            return err;
        cmt_param->interp_on_spread_flg = val;
        return NULL;
    }
    if (!strcmp(param_name, "SWAPREFRATE"))
    {
        strcpy(cmt_param->swap_ref_rate, param_val);
        return NULL;
    }

    return serror("Do not know parameter name %s", param_name);
}

/* ------------------------------------------------------------------------- */

Err interp_cmt_string(String cmt_index_string, CMTCode* cmt_code)
{
    Err err = NULL;

    if (!strcmp(cmt_index_string, "CMT1") || !strcmp(cmt_index_string, "CMT_1"))
    {
        *cmt_code = CMT1;
    }
    else if (!strcmp(cmt_index_string, "CMT2") || !strcmp(cmt_index_string, "CMT_2"))
    {
        *cmt_code = CMT2;
    }
    else if (!strcmp(cmt_index_string, "CMT3") || !strcmp(cmt_index_string, "CMT_3"))
    {
        *cmt_code = CMT3;
    }
    else if (!strcmp(cmt_index_string, "CMT5") || !strcmp(cmt_index_string, "CMT_5"))
    {
        *cmt_code = CMT5;
    }
    else if (!strcmp(cmt_index_string, "CMT7") || !strcmp(cmt_index_string, "CMT_7"))
    {
        *cmt_code = CMT7;
    }
    else if (!strcmp(cmt_index_string, "CMT10") || !strcmp(cmt_index_string, "CMT_10"))
    {
        *cmt_code = CMT10;
    }
    else if (!strcmp(cmt_index_string, "CMT20") || !strcmp(cmt_index_string, "CMT_20"))
    {
        *cmt_code = CMT20;
    }
    else if (!strcmp(cmt_index_string, "CMT30") || !strcmp(cmt_index_string, "CMT_30"))
    {
        *cmt_code = CMT30;
    }
    else if (!strcmp(cmt_index_string, "TEC10") || !strcmp(cmt_index_string, "TEC_10"))
    {
        *cmt_code = TEC10;
    }
    else
        return serror("Unknown CMT type...");

    return err;
}
