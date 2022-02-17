/* ==================================================================================

   FILENAME :  srt_h_und_fct.h

   PURPOSE  :  macros and functions to access underlying information from the
               underlying object stored in the double linked list (or its name):
               SrtUnd == SrtUndDesc

   ==================================================================================
 */

#ifndef SRT_H_UND_FCT_H
#define SRT_H_UND_FCT_H

/* -------------------------------------------------------------------------------
   Some Macros to access immediately underlying information from the SrtUnd
   -------------------------------------------------------------------------------
 */

/* ----------------------------- Info from any underlying
 * ------------------------ */

#define get_underlying_name(undptr) (undptr->underl_name)

#define get_underlying_lbl(undptr) (undptr->underl_lbl)

#define get_underlying_type(undptr) (undptr->underl_type)

#define get_underlying_ticker(undptr) (undptr->underl_ticker)

#define get_underlying_ccy(undptr) (undptr->underl_ccy)

#define ISUNDTYPE(undptr, type) (undptr->underl_type == (type))

/* ---------------------- Interest Rate Underlyings
 * ------------------------------ */
/* GET */
#define get_ycname_from_irund(undptr)                                          \
  (((SrtIrDesc *)(undptr->spec_desc))->yc_name)

#define get_mdltype_from_irund(undptr)                                         \
  (((SrtIrDesc *)(undptr->spec_desc))->mdl_type)

#define get_mdldim_from_irund(undptr)                                          \
  (((SrtIrDesc *)(undptr->spec_desc))->mdl_dim)

#define get_ts_from_irund(undptr) (((SrtIrDesc *)(undptr->spec_desc))->ts)

/* SET */
#define set_irund_mdltype(undptr, mdltype)                                     \
  (((SrtIrDesc *)(undptr->spec_desc))->mdl_type = mdltype)

#define set_irund_ts(undptr, termstruct)                                       \
  (((SrtIrDesc *)(undptr->spec_desc))->ts = termstruct)

/* ------------------------- Equity Underlyings
 * ---------------------------------- */

#define get_ts_from_eqund(undptr) (((SrtEqDesc *)(undptr->spec_desc))->ts)

#define get_discname_from_eqund(undptr)                                        \
  (((SrtEqDesc *)(undptr->spec_desc))->disc_name)

#define get_dvdname_from_eqund(undptr)                                         \
  (((SrtEqDesc *)(undptr->spec_desc))->dvd_name)

#define get_reponame_from_eqund(undptr)                                        \
  (((SrtEqDesc *)(undptr->spec_desc))->repo_name)

#define get_mdltype_from_eqund(undptr)                                         \
  (((SrtEqDesc *)(undptr->spec_desc))->mdl_type)

#define get_mdldim_from_eqund(undptr)                                          \
  (((SrtEqDesc *)(undptr->spec_desc))->mdl_dim)

#define get_spot_from_eqund(undptr) (((SrtEqDesc *)(undptr->spec_desc))->spot)

/* --------------------- Foreign Exchange Underlyings
 * ---------------------------- */

char *get_discname_from_fxund(SrtUndPtr undptr);

#define get_ts_from_fxund(undptr) (((SrtFxDesc *)(undptr->spec_desc))->ts)

#define get_domname_from_fxund(undptr)                                         \
  (((SrtFxDesc *)(undptr->spec_desc))->dom_name)

#define get_forname_from_fxund(undptr)                                         \
  (((SrtFxDesc *)(undptr->spec_desc))->for_name)

#define get_mdltype_from_fxund(undptr)                                         \
  (((SrtFxDesc *)(undptr->spec_desc))->mdl_type)

#define get_mdldim_from_fxund(undptr)                                          \
  (((SrtFxDesc *)(undptr->spec_desc))->mdl_dim)

#define get_spot_from_fxund(undptr) (((SrtFxDesc *)(undptr->spec_desc))->spot)

/* -------------------------- Bond Underlyings
 * ---------------------------------- */

#define get_ts_from_bndund(undptr) (((SrtBndDesc *)(undptr->spec_desc))->ts)

#define get_discname_from_bndund(undptr)                                       \
  (((SrtBndDesc *)(undptr->spec_desc))->disc_name)

#define get_reponame_from_bndund(undptr)                                       \
  (((SrtBndDesc *)(undptr->spec_desc))->repo_name)

#define get_mdltype_from_bndund(undptr)                                        \
  (((SrtBndDesc *)(undptr->spec_desc))->mdl_type)

#define get_mdldim_from_bndund(undptr)                                         \
  (((SrtBndDesc *)(undptr->spec_desc))->mdl_dim)

#define get_yield_from_bndund(undptr)                                          \
  (((SrtBndDesc *)(undptr->spec_desc))->spot_yield)

#define get_cleanprice_from_bndund(undptr)                                     \
  (((SrtBndDesc *)(undptr->spec_desc))->clean_price)

/* ----------------------------------------------------------------------
   Some functions to access  underlying information from the SrtUndPtr
   ---------------------------------------------------------------------- */

Err get_underlying_mdltype(SrtUndPtr und, SrtMdlType *mdl_type);
#define get_und_mdltype get_underlying_mdltype

Err get_underlying_mdldim(SrtUndPtr und, SrtMdlDim *mdl_dim);
#define get_und_mdldim get_underlying_mdldim

Err get_underlying_ts(SrtUndPtr und, TermStruct **ts);
#define get_und_ts get_underlying_ts

Err free_underlying_ts(SrtUndPtr und);
#define free_und_ts free_underlying_ts

/* Get the Discount Curve Name attached to the underlying */
Err get_underlying_discname(SrtUndPtr und, String *yc_name);
#define get_underlying_ycname get_underlying_disc_name

String get_discname_from_underlying(SrtUndPtr und);
#define get_ycname_from_underlying get_discname_from_underlying

/* Get the Growth Curve Name attached to the underlying */
Err get_underlying_growname(SrtUndPtr und, String *curve_name);

String get_dividend_name_from_underlying(SrtUndPtr und);
String get_repo_name_from_underlying(SrtUndPtr und);

/* Checks whether an interest rate underlying fits into the Cheyette framework
 */
SRT_Boolean is_model_Cheyette_type(SrtMdlType mdl_type);
SRT_Boolean is_irund_Cheyette_type(SrtUndPtr);

/* Checks whether an interest rate underlying fits into the New framework */
SRT_Boolean is_model_New_type(SrtMdlType mdl_type);
SRT_Boolean is_irund_New_type(SrtUndPtr);

/* -------------------------------------------------------------------------
   Functions to extract Yield Curve Information through the Underlying
   ------------------------------------------------------------------------- */

Date get_today_from_underlying(SrtUndPtr und);

Date get_spotdate_from_underlying(SrtUndPtr und);

int get_spotlag_from_underlying(SrtUndPtr und);

/* --------------------------------------------------------------------------------
   Specifically for an FX underlying: gets the currencies through the growth
   curves
   ---------------------------------------------------------------------------------
 */

Err get_fx_underlying_currencies(SrtUndPtr fxund, String *dom_ccy,
                                 String *for_ccy);

#endif
