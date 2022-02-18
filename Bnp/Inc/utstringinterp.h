/* ======================================================
   FILENAME:  utstringinterp.h

   PURPOSE:   A few functions to interpret a string
              into a type
   ====================================================== */

#ifndef UTSTRINGINTERP_H
#define UTSTRINGINTERP_H

/* ---------------------------------------------------------------------------
                                        DATES
   --------------------------------------------------------------------------- */
Err interp_month(const char* constStr, SrtMonth* val);
Err translate_month(char** str, SrtMonth val);

Err interp_bus_day_conv(String str, SrtBusDayConv* val);
Err translate_bus_day_conv(String* str, SrtBusDayConv m);

Err interp_unit(const char* constStr, SrtUnit* val);

/* ---------------------------------------------------------------------------
             SWAP CONVENTIONS: BASIS, COMPOUNDING
   --------------------------------------------------------------------------- */

Err interp_basis(const char* constStr, SrtBasisCode* val);
Err translate_basis(char** str, SrtBasisCode val);

Err interp_compounding(const char* constStr, SrtCompounding* val);
Err translate_compounding(char** str, SrtCompounding val);

Err interp_interp_method(String str, SrtInterpMethod* val);
Err translate_interp_method(String* str, SrtInterpMethod m);

Err interp_rec_pay(const char* constStr, SrtReceiverType* val);
Err interp_shift(const char* constStr, SrtShiftType* val);

/* ---------------------------------------------------------------------------
             OPTIONS : CALL/PUT, UP/DOWN...
   --------------------------------------------------------------------------- */

Err interp_call_put(const char* constStr, SrtCallPutType* val);

Err conv_rec_in_call(SrtReceiverType rec_pay, SrtCallPutType* call_put);

Err interp_diffusion_type(const char* constStr, SrtDiffusionType* val);

Err interp_min_max(const char* constStr, SrtMinmaxType* val);

Err interp_barrier_type(const char* constStr, SrtBarrierType* val);

Err interp_grad_best(const char* constStr, SrtLadderType* val);

Err interp_greeks(const char* constStr, SrtGreekType* val);

Err interp_best_worst(const char* constStr, SrtBestWorstType* val);

/* ---------------------------------------------------------------------------
                         RESET TYPES
   --------------------------------------------------------------------------- */
Err interp_reset_optimised(const char* constStr, SrtResOptType* val);
Err interp_price_type(const char* constStr, SrtPriceType* val);

/* ---------------------------------------------------------------------------
                          FOR NUMERICAL METHODS
   --------------------------------------------------------------------------- */
Err interp_derivatives(const char* constStr, SrtDerType* val);

/* ---------------------------------------------------------------------------
                          FOR SABRGetVol
   --------------------------------------------------------------------------- */
Err interp_SABRVolComponent(const char* constStr, SABRVolComponent* val);

/*-----------------------------------------------------------------------------
                                                FOR DELTA REPORT
------------------------------------------------------------------------------*/
Err interp_HedgeType(const char* constStr, SrtHedgeType* val);
Err interp_UndFRAType(const char* constStr, SrtUndFRAType* val);

/*-----------------------------------------------------------------------------
                                                FOR BOOLEAN
------------------------------------------------------------------------------*/

Err interp_BoolType(const char* constStr, SRT_Boolean* val);

/*-----------------------------------------------------------------------------
                                                FOR CALIBRATION
------------------------------------------------------------------------------*/

Err interp_CalibrationType(const char* constStr, SrtCalibrationType* val);

Err interp_inflation_indexed_bond_asset_swap(
    const char* constStr, SrtInflationIdxBondAssetSwapType* val);

#endif
