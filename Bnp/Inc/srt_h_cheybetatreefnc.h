/* ============================================================================

   FILENAME:			srt_h_cheybetatreefnc.h

   DESCRIPTION:			utility functions for local vol Cheyette tree

   ============================================================================
 */
#ifndef srt_h_cheybetatreefnc_h
#define srt_h_cheybetatreefnc_h

typedef struct {
  long max_state_var_index;
  long min_state_var_index;
  long max_phi_index;
  long min_phi_index;
  double max_state_var;
  double min_state_var;
  double delta_state_var;
  double mid_phi;
  double max_phi;
  double min_phi;
  double delta_phi;
} SrtLocTreInf;

/* ---------------------------------------------------------------------------
 */
/* Find the local vol in a one factor model */
/* This function should replace at some stage the srt_f_chevaratsam function */
Err srt_f_one_fac_find_local_vol(SrtMdlType mdl_type, SrtStpPtr stp,
                                 SrtSample *cur_sample, int und_index,
                                 double *vol);

/*	Function:		srt_f_locvol_init_trinf
Object:			initialise the tree information: geomtry and limits */

Err srt_f_CheyBetaTree_init_trinf(
    SrtStpPtr stp,           /*	Step pointer */
    SrtGrfnParam *grfparam,  /*	Mdl discretisation parameters */
    SrtUndPtr und,           /*	Underlying */
    SrtLocTreInf *maxtrinf); /*	Will contain the smallest min_r_index
                                                                  and the
                                greates max_r_index (for allocation) */
/* ---------------------------------------------------------------------------
 */
/*	Function:		srt_f_locvol_make_r_grid
        Object:			make r grid  , once for all */

Err srt_f_CheyBetaTree_make_state_var_grid(SrtLocTreInf maxtrinf,
                                           double *r_grid);

/* ---------------------------------------------------------------------------
 */
/*	Function:		srt_f_locvol_make_phi_grid
        Object:			make phi grid  , for a given time step */

Err srt_f_CheyBetaTree_make_phi_grid(SrtLocTreInf *trinf, double *cur_phi_grid);

/* ---------------------------------------------------------------------------
 */
/*	Function:		srt_f_locvol_getmoments
        Object:			Compute moments and store in node
                                        Fills exp_sv and var_sv in node */

Err srt_f_locvol_getmoments(SrtTrinTreNdInf *node, SrtIRMTmInf *tminf);

/* ---------------------------------------------------------------------------
 */
/*	Function:		srt_f_locvol_recombine_in_r
        Object:			Compute recombination of the tree in r */

Err srt_f_CheyBetaTree_recombine_in_r(SrtTrinTreNdInf *node, double *r_grid,
                                      SrtLocTreInf *nxttrinf);

/* ---------------------------------------------------------------------------
 */
/*	Function:		srt_f_locvol_calc_prob
        Object:			Compute probabilities */

Err srt_f_CheyBetaTree_calc_prob(SrtTrinTreNdInf *node);

/* ---------------------------------------------------------------------------
 */
/*	Function:		srt_f_locvol_discount_assets
        Object:			Discount asset prices */

Err srt_f_CheyBetaTree_discount_assets(double *cur_assets,
                                       double ***next_assets,
                                       SrtTrinTreNdInf *node, long nassets,
                                       SrtLocTreInf *nxttrinf,
                                       double *next_r_grid, double next_phi,
                                       double *next_phi_grid);

#endif