/* ------------------------------------------------------------------------------
   AUTHOR: E. FOURNIE & al 

   DATE : MARS 98

   FILENAME:    srt_h_cheybeta2ftreefct.h

   PURPOSE:     New structure  of  random grid for the pink elephant 

   MODIFICATION:

   ------------------------------------------------------------------------------ */

#ifndef SRT_H_CHEYBETA2FTREEFCT_H
#define	SRT_H_CHEYBETA2FTREEFCT_H

#include "num_h_triangle.h"

/*------------------------------------------------------------------------------- */

#define    NB_NEIGHBOURS   5
#define    NB_BOUND_PTS	 4

/* ------------------------------ STRUCTURES ------------------------------------ */
typedef struct 
{
/* state variables at time t- */
    double      x1;
	double      x2;
	double 	    jump_x1;
	double      jump_x2;
/* non markovian variables time t-1 */
	double      phi1;
	double      phi2;
	double      phi1_2;

	double      fw1;
	double      fw2;
	double      short_rate;
/* discount factor from current step to the next one */	
	double      df;
	
} MCTreePoint;


/* ---------------------------------------------------------------------------------
   FUNCTION:     Err srt_f_cheybeta2d_make_grid()
   PURPOSE:     Build the regular grid based on simulation
   --------------------------------------------------------------------------------- */
Err  srt_f_cheybeta2ftree_make_grid( 
				SrtGrfnParam     *grfnparam,
				SrtStpPtr         step,
				SrtUndInfo       *und_info,
				MCTreePoint     **grid);

/* --------------------------------------------------------------------------------
   FUNCTION: Err srt_f_cheybeta2ftree_discount_to_node (...) 
   PURPOSE:  Backward induction for one given node : probabilities computed 
   -------------------------------------------------------------------------------- */
Err srt_f_cheybeta2ftree_discount_to_node ( 
						MCTreePoint       **grid,
						int                 index_t,
						int                 index_n,
						SrtStpPtr           step,
						double             *cur_asset,	       /* Current assets at node */
						SrtGrfnParam	   *grfnparam);


/* -------------------------------------------------------------------------------------------
   FUNCTION:  Err srt_f_cheybeta2ftree_compute_proba(...)
   PURPOSE:   compute the standards points and probas
   ------------------------------------------------------------------------------------------- */
Err srt_f_cheybeta2ftree_compute_proba(
					SrtStpPtr     step);

/* -------------------------------------------------------------------------------------------
   FUNCTION:  Err srt_f_alloc_computing_resources(...)
                      Err srt_f_free_computing_resources(...)
   PURPOSE:  alloc and free some data structures useful for the computation
   ------------------------------------------------------------------------------------------- */
Err srt_f_alloc_computing_resources(void);

Err srt_f_free_computing_resources(void);

/* -------------------------------------------------------------------------------------------
   FUNCTION:  void transfer_info_nodes(...)
   PURPOSE:  transfer nodes infos to triangulation structure
   ------------------------------------------------------------------------------------------- */
 void transfer_info_nodes(
	 MCTreePoint          *mesh,
	 double              **next_assets,
	 SrtTriangulationIO   *triang);


#endif
/*--------------------------------- End of File -------------------------------------*/ 