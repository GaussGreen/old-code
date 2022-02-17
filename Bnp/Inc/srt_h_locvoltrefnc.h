#ifndef srt_h_locvoltrefnc_h
#define srt_h_locvoltrefnc_h

typedef struct {
  SrtSample cur_sam;
  SrtSample drift_sam;
  long son_index[3];
  double son_sv[3];
  double p[3];
  double df;
  double exp_sv;
  double var_sv;
  double time;
  double delta_t;
  SrtMdlType mdl_type;
  TermStruct *ts;
} SrtLocTreNdInf;

typedef struct {
  long max_r_index;
  long min_r_index;
  long max_phi_index;
  long min_phi_index;
  double max_r;
  double min_r;
  double delta_r;
  double mid_phi;
  double max_phi;
  double min_phi;
  double delta_phi;
} SrtLocTreInf;

Err srt_f_locvol_init_tminf(
    SrtStpPtr stp,  /*	Step pointer:
                                        the step list is supposed to have
                   been init */
    SrtUndPtr und); /*	Underlying pointer:
                                            contains the relevant info */

/*	Function:		srt_f_locvol_init_trinf
        Object:			init tree info in stp: allocate and populate */

Err srt_f_locvol_init_trinf(
    SrtStpPtr stp,           /*	Step pointer */
    SrtMdlPtr mdl,           /*	Mdl discretisation parameters */
    SrtUndPtr und,           /*	Underlying */
    SrtLocTreInf *maxtrinf); /*	Will contain the smallest min_r_index
                                                                  and the
                                greates max_r_index (for allocation) */

/*	Function:		srt_f_locvol_make_r_grid
        Object:			make r grid      , once for all */

Err srt_f_locvol_make_r_grid(SrtLocTreInf maxtrinf, double *r_grid);

/*	Function:		srt_f_locvol_make_phi_grid
        Object:			make phi grid      , for a given time step */

Err srt_f_locvol_make_phi_grid(SrtLocTreInf *trinf, double *cur_phi_grid);

/*	Function:		srt_f_locvol_getmoments
        Object:			Compute moments and store in node
                                        Fills exp_sv and var_sv in node */

Err srt_f_locvol_getmoments(SrtLocTreNdInf *node, SrtIRMTmInf *tminf);

/*	Function:		srt_f_locvol_recombine_in_r
        Object:			Compute recombination of the tree in r */

Err srt_f_locvol_recombine_in_r(SrtLocTreNdInf *node, double *r_grid,
                                SrtLocTreInf *trinf, SrtLocTreInf *nxttrinf);

/*	Function:		srt_f_locvol_calc_prob
        Object:			Compute probabilities */

Err srt_f_locvol_calc_prob(SrtLocTreNdInf *node);

/*	Function:		srt_f_locvol_discount_assets
        Object:			Discount asset prices */

Err srt_f_locvol_discount_assets(double *cur_assets, double ***next_assets,
                                 SrtLocTreNdInf node, long nassets,
                                 SrtLocTreInf *nxttrinf, double *next_r_grid,
                                 double next_phi, double *next_phi_grid);

#endif