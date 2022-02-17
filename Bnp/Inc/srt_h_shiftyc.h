/* SRT_H_SHIFTYC.h */

/* ------------------------------------------------ */
/* ----- Bart Simpson      , Dec. 1994 ----- X-ENPC ----- */
/* ------------------------------------------------ */

#include "srt_h_all.h"

/* -------------------------------------------------------------------------- */

Err srt_f_parallelshift_yc(YC_Obj *yc, double shift_value);

/* -> Take yc as an input and shift it with the value 'shift_value
        The input is physically modified				 */

/* -------------------------------------------------------------------------- */

Err srt_f_update_yc(SrtUndPtr und, double shift_value);

/* -> Taking a market 'mkt and a value of shift 'shift_value      ,
        extract the yc from the market      , shift it      , put it back
        in the market and update everything that depends on it
      The market is physically modified                          */

/* -------------------------------------------------------------------------- */

YC_Obj *srt_f_dupYC_Obj(YC_Obj *old_YC_ptr);
