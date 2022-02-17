/* -------------------------------------------------------------------------

   FILE NAME	: swp_h_external_fct.h

   PURPOSE		: Defines the structure used to store all the pointer to
                  functions that have to be set before calling any SRT
                                  financial function.
                                  This include:
                                        - DiscFunc (df_curve_id/ccy      , start
   , end      , *df)
                                        - SpreadFunc (ref_rate_name      , start
   , end      , *spread)
                                        - VolFunc ( vol_curve_id      , start ,
   end , strike      , *vol )
                                        - ....

                                  Gives the functions to communicate with this
   structure
   ------------------------------------------------------------------------- */

#ifndef SWP_H_EXTERNAL_FCT_H
#define SWP_H_EXTERNAL_FCT_H

/* The DiscounFactor function type */

typedef Err (*DiscFuncType)(char *crv_id, double start, double end, double *df);

/* The SpreadFunction function type */

typedef Err (*SpreadFuncType)(char *rate_id, double start, double end,
                              double *spread);

/* The RateInformation function type */

typedef Err (*RateInfoFuncType)(char *rate_id, char *family_code, char *ccy,
                                char *tenor, char *basis, int *frequency,
                                char *calc_method, int *number_of_period);

/* The Volatility function type */

typedef Err (*VolFuncType)(char *vol_id, double start, double end,
                           double strike, double *vol, double *power);

/* The Correlation function type */

typedef Err (*CorrelFuncType)(char *vol_id, double start, double end,
                              double strike, double *vol);

/* The Volatility component function type */
/* ATMlog      , ATMNorm      , ATMBeta      , Alpha      , Beta      , Rho */
/* Log      , Norm for a given Strike */

typedef Err (*SABRVolFuncType)(char *vol_id, double start, double end,
                               double strike, double *vol, double *power,
                               int component);

/* The Historical Fixing function type */
typedef Err (*FixingFuncType)(long date, char *refratecode, double *fixing);

/* The Structure to store them all */

typedef struct _SrtExternalFunctions {
  DiscFuncType disc_func;

  SpreadFuncType spread_func;

  RateInfoFuncType rate_info_func;

  VolFuncType vol_func;

  CorrelFuncType correl_func;

  SABRVolFuncType sabrvol_func;

  FixingFuncType fixing_func;

} SrtExternalFunctions;

/* -----------------------------------------------------------------------------
   Resets to NULL the Official External Functions
   -----------------------------------------------------------------------------
 */
Err swp_f_ResetExternalFunctions();

/* -----------------------------------------------------------------------------
   Sets or Gets the Official Discount Function used in the pricing
   -----------------------------------------------------------------------------
 */
Err swp_f_SetDiscFunc(DiscFuncType disc_func);

Err swp_f_GetDiscFunc(DiscFuncType *disc_func);

/* -----------------------------------------------------------------------------
   Sets or Gets the Official Spread Function used in the pricing
   -----------------------------------------------------------------------------
 */
Err swp_f_SetSpreadFunc(SpreadFuncType spread_func);

Err swp_f_GetSpreadFunc(SpreadFuncType *spread_func);

/* -----------------------------------------------------------------------------
   Sets or Gets the Official Rate Information Function used in the pricing
   -----------------------------------------------------------------------------
 */
Err swp_f_SetRateInfoFunc(RateInfoFuncType rate_info_func);

Err swp_f_GetRateInfoFunc(RateInfoFuncType *rate_info_func);

/* -----------------------------------------------------------------------------
   Sets or Gets the Official Vol Interpolation Function used in the pricing
   -----------------------------------------------------------------------------
 */
Err swp_f_SetVolFunc(VolFuncType vol_func);

Err swp_f_GetVolFunc(VolFuncType *vol_func);

/* -----------------------------------------------------------------------------
   Sets or Gets the Official SABRVol Component Interpolation Function used in
   the pricing
   -----------------------------------------------------------------------------
 */
Err swp_f_SetSABRVolFunc(SABRVolFuncType sabrvol_func);

Err swp_f_GetSABRVolFunc(SABRVolFuncType *sabrvol_func);

/* -----------------------------------------------------------------------------
   Sets or Gets the Official Historical Fixing Function used in the pricing
   -----------------------------------------------------------------------------
 */
Err swp_f_SetFixingFunc(FixingFuncType fixing_func);

Err swp_f_GetFixingFunc(FixingFuncType *fixing_func);

#endif