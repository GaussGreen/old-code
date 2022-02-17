#ifndef SRT_H_IRMFCT_H
#define SRT_H_IRMFCT_H

/* ========================================================================

	IRM functions 

   ======================================================================== */



/*
 * set the Y_T_at_t_params needed to price a rate described
 * by a SrtRtFnc structure.  Will only allocate space if has not already
 * been allocated.
 */

SrtErr srt_f_rtinityp(SrtRtFnc rt, SrtUndPtr und, Ddate obsdate);

/*
 * eval the Y_T_at_t_params contained in a SrtRtFnc structure;
 * and set the dfs.  If df field is blank, will allocate it.
 */

SrtErr srt_f_rtevalyp(
		SrtRtFnc      rt, 
		SrtSample     *sam, 
		int           index, 
		SrtMdlDim     mdl_dim,
		SrtMdlType    mdl_type);



/* ======================================================================

	CHE specific functions 

   ====================================================================== */



/*  =====================================================================
    This function is defined in srt_f_chephiappx.c 
    ====================================================================== */	
void   srt_f_chelinrphi
  (
  SrtStpPtr top, 
  SrtStpPtr stp, 
  SrtSample *sam      /* answer returned here */
  ) 
;

/*  =====================================================================
    These two functions are defined in srt_f_chesammtx.c 
    ====================================================================== */	
void srt_f_chephilim
(
  SrtStpPtr top, 
  SrtStpPtr stp, 
  SrtIRMTmInf *tminf,
  SrtSample *sam,
  double *phimin,
  double *phimax
)
;

void srt_f_chesammatrix
(
  SrtStpPtr top,    /* ptr to top of lst */
  SrtStpPtr stp,    /* ptr to current step in lst */
  double *cur_r,    /* space for rates to be stored */ 
  double **cur_phi, /* space for phi's to be stored */
  SrtGrfnParam 	    *grfnparam/* info about model */
)
;


/* ======================================================================
    This function is defined in srt_f_vegatreche2dr.c 
   ====================================================================== */
Err srt_f_treche2dr
	(
	SrtUndPtr und, /* Underlying */
	SrtMdlPtr mdl, /* model parameters */
	SrtStpPtr stp, /* discretization of deal in time, wi/ events attached*/
	SrtDealDesc gd,/* deal descriptor structure */
	EvalEventFct evalcf, /* cashflow evaluator */
	double *price
	)
;

#endif
