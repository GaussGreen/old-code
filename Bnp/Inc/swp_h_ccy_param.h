/* ===========================================================================================
        Filename    swp_h_ccy_param.h

        EA:	declaration and calling functions for ccy_params

        KNL: 	added SrtYCConvTyp      , 7th Sept 1994

        KNL:	renamed to standard filename convention      , 28th Nov 1994
   ============================================================================================
 */

#ifndef SWP_H_CCY_PARAM_H
#define SWP_H_CCY_PARAM_H

/* ---------------------------------- */
/* new enum type for curve convention */
/*  - used in swp_f_market_list.c     */
/* ---------------------------------- */

typedef enum SrtYCConvTyp_ {
  SRT_YC_CONV_SRT,
  SRT_YC_CONV_SPG
} SrtYCConvTyp,
    *SrtYCConvTypPtr;

/* ---------------------------------- */

typedef struct {
  CcyCode currency;
  SrtYCConvTyp curve_conv;
  BusDayConv cash_bus_day_conv;
  BusDayConv swap_bus_day_conv;
  BasisCode cash_basis_code;
  BasisCode m1fut_basis_code;
  BasisCode m3fut_basis_code;
  BasisCode swap_basis_code;
  SrtCompounding compd;
  InterpMethod interp_method;
  int spot_lag;
  Err (*fut_func)(int nm, int spot, int sy, int sm, Date *fut_last_trading,
                  Date *fut_start, Date *fut_end);
  SRT_Boolean ilv_insert_flg;
  SRT_Boolean ilv_no_overwrite_flg;
  SRT_Boolean toy_insert_flg;
  SRT_Boolean toy_no_overwrite_flg;
  SRT_Boolean customized;
} SrtCcyParam, SrtCcyParamStruct;

SrtCcyParam *new_CcyParam();

int free_CcyParam(SrtCcyParam *cps);

Err swp_f_get_CcyParam_from_CcyStr(char *CcyStr, SrtCcyParam **ccy_param);

Err ccy_string_get_param(SrtCcyParam *ccy_param, String param_name,
                         String *param_val);

Err ccy_string_set_param(SrtCcyParam *ccy_param, String param_name,
                         String param_val);

Err ccy_add_definition(
    String ccy_code, SrtYCConvTyp curve_conv, BusDayConv cash_bus_day_conv,
    BusDayConv swap_bus_day_conv, BasisCode cash_basis_code,
    BasisCode m1fut_basis_code, BasisCode m3fut_basis_code,
    BasisCode swap_basis_code, SrtCompounding compd, InterpMethod interp_method,
    int spot_lag, SRT_Boolean ilv_insert_flg, SRT_Boolean ilv_no_overwrite_flg,
    SRT_Boolean toy_insert_flg, SRT_Boolean toy_no_overwrite_flg);

void ccy_clear_dynamic_defs(void);

#endif
