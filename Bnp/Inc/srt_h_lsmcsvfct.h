/* ------------------------------------------------------------------------------
   AUTHOR: E. FOURNIE & al 

   DATE : Septembre 99

   FILENAME:    srt_h_cheysvtreefct.h

   PURPOSE:     

   MODIFICATION:

   ------------------------------------------------------------------------------ */

#ifndef SRT_H_CHEYSVTREEFCT_H
#define	SRT_H_CHEYSVTREEFCT_H

/* ------------------------------ STRUCTURES ------------------------------------ */
typedef struct 
{
/* state variables at time t- */
    double      x;
	double      phi;
	double      sigma;
	double      short_rate;
	double      df;
	double      numeraire;
} MCTreePointCheyBetaStochVol;

Err srt_f_alloc_lsmcsv_resources(int nreg, int node_dim);

Err srt_f_free_lsmcsv_resources(int nreg, int node_dim);

Err srt_f_lsmcsv_regr(      int nreg, 
							int node_dim, 
							int ntraj,
							int nt,                   
							MCTreePointCheyBetaStochVol  **grid,
							double           ***bwd_assets,
							double m1,
							double m2,
							double s11,
							double s12,
							double s21,
							double s22);

Err cheybetastochvol_simul(		  SrtGrfnParam    *grfnparam,
								  SrtStpPtr        step,
								  SrtUndInfo       *und_info,
								  MCTreePointCheyBetaStochVol     **grid,
								  long             nt);



#endif
/*--------------------------------- End of File -------------------------------------*/ 