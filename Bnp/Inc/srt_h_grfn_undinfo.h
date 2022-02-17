/* ==================================================================================

   FILENAME :  srt_h_grfn_undinfo.h

   PURPOSE :   the structure to store and the functions to work with the global
               underlyings info for a multi-underlyings evaluation in Grfn

   ==================================================================================
 */

#ifndef SRT_H_GRFN_UNDINFO_H
#define SRT_H_GRFN_UNDINFO_H

/* ---------------------------------------------------------------- */
/* Internal structure to store underlyings info                     */

typedef struct {
  char und_name[SRTBUFSZ];
  SRT_Boolean stochastic;
  char und_ccy[SRTBUFSZ];
} Und_Data;

/* ----------------------------------------------------------------------- */
/* All the underlyings informations for a multi-underlying Grfn evaluation */

typedef struct {
  int no_of_underlyings;
  int no_of_brownians;
  int numeraire_index; /* discounting underlying index */
                       /* it could be different from 0 (eg FX_STOCH_RATES) */
  SRT_Boolean jumping; /* whether jumping numeraire can be applied */
  Und_Data und_data[MAXUNDERLYING];
  SrtCorrLst *
      corr_ts; /* Points to the correlation matrix list used in the deal only */
  SRT_Boolean two_factor_model; /* If two factor model: do not correl */

  SRT_Boolean use_stochastic_vol; /* If stochastic vol  : do not correl */
  int use_corr_ts; /* Attempt to get rid of the correlation TS of sort ! */
  /* If set to 0      , we skip in MdlCom the corr TS checks... */

} SrtUndInfo;

/* Set all the SrtUndInfo flags from the list of undrlyings stored in the
 * SrtUndInfo */
Err analyze_und_info(SrtUndInfo *und_info);

/* Get the index in SrtUndInfo corresponding an underlying name */
Err get_index_from_und_info(SrtUndInfo *und_info, char *name, int *index);

/* Add some underlying which are not already in the SrtUndInfo list but are
   required for the computation (dom and for IR_UND for an FX_STOCH_RATES...) */
Err add_more_underlyings_to_und_info(SrtUndInfo *und_info);

#endif