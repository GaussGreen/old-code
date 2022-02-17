/* ----------------------------------------------------------------------------------
   AUTHOR: E. FOURNIE & VE  

   DATE : JULY 98

   FILENAME:  srt_f_lgm2ftree.c

   PURPOSE:  main implementation of the 2 factor LGM tree in GRFN

   MODIFICATION:

   ----------------------------------------------------------------------------------*/
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_lgm2ftreefct.h"

#define SWAP(a,b) {void* temp=a; a=b; b=temp;}                            

#define ALLOC_MEM_ASSETS   {\
	cur_assets = f3tensor(- maxtrinf.max_index[0], maxtrinf.max_index[0],\
                                    - maxtrinf.max_index[1], maxtrinf.max_index[1],\
							        0, node_dim - 1);\
	next_assets = f3tensor(- maxtrinf.max_index[0], maxtrinf.max_index[0],\
                					  - maxtrinf.max_index[1], maxtrinf.max_index[1],\
									  0, node_dim - 1); }

#define FREE_MEM_ASSETS   {\
	free_f3tensor(cur_assets, - maxtrinf.max_index[0], maxtrinf.max_index[0],\
                						 - maxtrinf.max_index[1], maxtrinf.max_index[1],\
										 0, node_dim - 1);\
	free_f3tensor(next_assets, - maxtrinf.max_index[0], maxtrinf.max_index[0],\
                						   - maxtrinf.max_index[1], maxtrinf.max_index[1],\
										   0, node_dim - 1); }		
 
/* --------------------------------------------------------------------------------------------- */

Err srt_f_lgm2dtree(
					SrtUndPtr 			und,			   /* LGM2F underlying */
					SrtGrfnParam     *grfnparam,	 /* Discretisation parameters */
					SrtStpPtr 			stp,			    /* Pointer on step list */
					GrfnDeal	        *gd,				/* Deal description */
					EvalEventFct 	  evalcf,			 /* DF evaluation function */
					SrtIOStruct		    *iolist,		   /* Requests */
					SrtUndInfo 		    *und_info)		/* Underlying info (used in evalcf) */
{
Err 		   err = NULL;
int 			i, j; 
double		 cashflow;
double 		 zc_yield;
SrtStpPtr 	top, bot;
/* Information about the  tree */
SrtPentoTreeNodeInfo			 node;
SrtDiagTwoFacTreeInfo 		maxtrinf, *trinf;
/* Information needed for the model at a particular time*/
SrtIRMTmInf 		      *tminf;
/* Functions values */
double  					  ***next_assets, ***cur_assets;
/* Dimensions */
int							     node_dim = grfnparam->node_dim;
/* For basis calculus */
double						 y1, y2;

/*--------- END of declarations -------------*/
	top = stp = gototop (stp);
	bot = gotobot (stp);

	/* Allocates time info pointer at each step of the stp structure */
	if (err = srtstptminfalloc (stp, 1))   return err;

	/* Initialises the steps for an LGM2F model: allocates and populates time info */
	if (err = srt_f_irministp (stp, und, 0, und, und_info))   return err;
	
	/* Allocates and populates tree info */
	if (err = lgm2f_trelim (stp, &maxtrinf))   return err;

	/* Memory Allocation */
	ALLOC_MEM_ASSETS
  	if (!next_assets || !cur_assets)
	{
		FREE_MEM_ASSETS
		return serror ("allocation failure srt_f_lgm2dtree");
	}

	/* Backward induction ---- Back for the futures */
	for (stp = bot; stp != NULL; stp = stp->prev)
	{
		/* Get trinf and tminf */
    	trinf = (SrtDiagTwoFacTreeInfo*) (stp->trinf);
    	tminf = (SrtIRMTmInf*) (stp->tminf[0]);

		/* Stores PHIs independently of the node in LGM (depend only on time) */
		sam_get (node.cur_sam, 0, PHI1) = sam_get (tminf->fwd_sam, 0, PHI1);
		sam_get (node.cur_sam, 0, PHI2) =  sam_get (tminf->fwd_sam, 0, PHI2);
		sam_get (node.cur_sam, 0, CROSSPHI) = sam_get (tminf->fwd_sam, 0, CROSSPHI);

		/* Loop on XiYj  in dual basis */
		for (i = - trinf->max_index[0]; i <= trinf->max_index[0]; i++) 
			for (j = - trinf->max_index[1]; j <= trinf->max_index[1]; j++) 	 
		{
				/* Transfer statevars back into original basis */
				y1 = i * trinf->spacing[0];
				y2 = j * trinf->spacing[1];
				sam_get (node.cur_sam, 0, X1) = trinf->dual_basis[0][0] * y1 + trinf->dual_basis[0][1] * y2;
				sam_get (node.cur_sam, 0, X2) = trinf->dual_basis[1][0] * y1 + trinf->dual_basis[1][1] * y2;

				/* If it is not the last step, compute forward, connections and probas */
				if (stp->next)
					populate_lgm2f_tree_node (tminf, &node, stp->delta_t, stp->next->trinf);

				/* r1 = X1 + H11 + H12, r2 = X2 + H22 + H21,  r = F(0,t) + r1 + r2 */
      			sam_get (node.cur_sam, 0, X1) += tminf->rf.twof[0][0].H + tminf->rf.twof[0][1].H;
      			sam_get (node.cur_sam, 0, X2) += tminf->rf.twof[1][0].H + tminf->rf.twof[1][1].H;
      			sam_get (node.cur_sam, 0, SHORT_RATE) = sam_get (tminf->fwd_sam, 0, F_0_t)
														+ sam_get (node.cur_sam, 0, X1) + sam_get (node.cur_sam, 0, X2);

				/* If it is not the last step, discount to node */
      			if (stp->next)
      			{
					/* Compute discount factor for node */
					Y_T_at_t_compute(1, &node.cur_sam, &tminf->yp, &zc_yield, 0, TWO_FAC, LGM);
					node.df = exp(- zc_yield);

					/* Compute discounted asset prices at current node */
	        		lgm2f_tree_expectation(&node, next_assets, cur_assets[i][j], node_dim);
        		} 

				/* Evaluate cash flows in Grfn tableau */
				if (err = evalcf((GrfnEvent *) stp->e, &node.cur_sam, gd, (double*) cur_assets[i][j],
                				             (EvalEventDfsFct) srt_f_calc_grfn_event_dfs, und_info, &cashflow))
				{
					FREE_MEM_ASSETS 
					return err;
				}
		
		} /* End of loop on i, j */

		SWAP(cur_assets, next_assets);
	} /* End of backward induction */

	/* Stores the premium in the Input/Output list */
	err = srt_f_IOstructsetpremium((iolist), SRT_NO, next_assets[0][0][grfnparam->node_dim-1], "Tree_LGM_2f");
	
	/* Stores all the columns PV in the I/O list */
	err = srt_f_IOstructsetcolpvs((iolist), SRT_NO, (double*) next_assets[0][0], grfnparam->node_dim, ""	);
	
	/* Free the tensors previously allocated */
	FREE_MEM_ASSETS 

	return err;

} /* END Err srt_f_lgm2dtree(...) */


/* ----------------------------------------------------------------------------------*/

#undef FREE_MEM_ASSETS 
#undef ALLOC_MEM_ASSETS 
#undef SWAP


/* ========= END OF FILE =================================================== */