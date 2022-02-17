/* ----------------------------------------------------------------------------------------------

   FILENAME:  srt_f_twoundtree.c

   PURPOSE:  main implementation of the generic 2 underlying tree in GRFN

   AUTHOR : Eric Fournie

   MODIFICATION:
   -----------------------------------------------------------------------------------------------*/
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_stpcorata.h"
#include "srt_h_twoundtreefct.h"
                         
#define SWAP(a,b) {void* temp = a; a = b; b = temp;}                            

#define ALLOC_MEM_ASSETS   \
	cur_assets = f3tensor(- maxtrinf.max_index[0], maxtrinf.max_index[0], \
                                    - maxtrinf.max_index[1], maxtrinf.max_index[1], \
							        0, node_dim - 1);\
	next_assets = f3tensor(- maxtrinf.max_index[0], maxtrinf.max_index[0], \
                					  - maxtrinf.max_index[1], maxtrinf.max_index[1], \
									  0, node_dim - 1) 

#define FREE_MEM_ASSETS   \
	free_f3tensor(cur_assets, - maxtrinf.max_index[0], maxtrinf.max_index[0], \
                						 - maxtrinf.max_index[1], maxtrinf.max_index[1], \
										 0, node_dim - 1);\
	free_f3tensor(next_assets, - maxtrinf.max_index[0], maxtrinf.max_index[0], \
                						   - maxtrinf.max_index[1], maxtrinf.max_index[1], \
										   0, node_dim - 1) 		

/* ------------------------------------------------------------------------------------------------ */

Err srt_f_twoundtree( 
				SrtGrfnParam   *grfnparam, 
				SrtStpPtr 	stp,  
				GrfnDeal	   *gd,
				EvalEventFct 	evalcf,  
				SrtIOStruct	*iolist,  
				SrtUndInfo 	*und_info)
{
	/* ----- declare stuff ---- */
Err 		err = NULL;
SrtUndPtr 	und1,und2;
int 	     i, j; 
double	  cashflow;
SrtStpPtr 	top = NULL, bot = NULL;
/* Information about the  tree */
SrtPentoTreeNodeInfo			 node;
SrtDiagTwoFacTreeInfo 		maxtrinf, *trinf;

/* information needed for the model at a particular time for cases LGM-LGM / LGM-BS / BS-BS */
SrtIRMTmInf 	*tminf_irm1 = NULL, *tminf_irm2 = NULL;
SrtLogTmInf      *tminf_log1 = NULL, *tminf_log2 = NULL;

/* Functions values */
double  					  ***next_assets, ***cur_assets;

/* state variable values at current time step and previous time step */
double 	**cur_st_var = NULL, **prev_st_var = NULL;
int		   node_dim = grfnparam->node_dim;
             
SrtIOVal	*io_request_val_p = NULL;
SrtListAtom	*io_request = NULL; 

/* To reference to both underlyings */
SrtMdlType	mdl_type1, mdl_type2;

double        x1, x2, y1, y2;

/* --------------- end of declarations ------------------------------------ */


/* An easy check before messing up everything */
	if (und_info->no_of_underlyings != 2 )
	return serror ("Wrong number of und for two fac tree");

/* Stores pointers to first step and to last one */
	top = stp = gototop(stp);
	bot = gotobot(stp);

/* Allocate space for the time info to each step */
	err = srtstptminfalloc(top, und_info->no_of_underlyings); 
	if (err)
	return err;

/* Gets the models from the FIRST (== DOMESTIC) underlying and initialise the Time Info in the steps accordingly */
	und1 = lookup_und(und_info->und_data[0].und_name);
	err = get_underlying_mdltype(und1, &mdl_type1);
	if (err)
		return err;	
	if (mdl_type1 == LGM )       
	{ 
		err = srt_f_irministp(top, und1, 0, und1, und_info);
	}
	else if ((mdl_type1 == BLACK_SCHOLES)  || (mdl_type1 == EQ_STOCH_RATES))
	{ 
		err = srt_f_loginistp(top, und1, 0, und1, und_info);
	}
	if (err)
		return err;

/* Gets the models from the SECOND underlying and initialise the Time Info in the steps accordingly */
	und2 = lookup_und(und_info->und_data[1].und_name);
	err = get_underlying_mdltype(und2, &mdl_type2);
	if (err)
		return err;
	if (mdl_type2 == LGM )
	{ 
		err = srt_f_irministp(top,und2,1, und1, und_info);
	}
	else 
	if ((mdl_type2 == BLACK_SCHOLES) || (mdl_type2 == EQ_STOCH_RATES))
	{ 
		err = srt_f_loginistp(top,und2,1, und1, und_info);
	}
	if (err)
		return err;

/* Attaches the correlation to all the steps */
	err = srt_f_attach_correl_to_stp(top, und_info->corr_ts);
	if (err)
		return err;

/* Quanto adjustment for the second underlying if necessary */
	err = srt_f_steps_set_quanto_adjustment(top, und_info, 1, und2,und1);
	if (err)
	{
		return err;
	}

/* Builds the tree: Localisation and meshing of the state space */
	err = twound_trelim(top, &maxtrinf, mdl_type1, mdl_type2);
	if (err)
		return err;


/* Memory Allocation for the assets cubes */ 
	ALLOC_MEM_ASSETS;


/* Backward induction ---- Back for the futures */
	for (stp = bot; stp != NULL; stp = stp->prev) 
  	{
	/* Gets the Two Factor Tree Information for this time step */
		trinf = (SrtDiagTwoFacTreeInfo*) (stp->trinf);

	/* Initialisations which are independent of the level PHI and VAR */
	    if (mdl_type1 == LGM)
		{
			tminf_irm1 = (SrtIRMTmInf *) (stp->tminf[0]);
		    sam_get(node.cur_sam, 0, PHI) = sam_get(tminf_irm1->fwd_sam, 0, PHI);  
		}
	    else if ((mdl_type1 == BLACK_SCHOLES) || (mdl_type1 == EQ_STOCH_RATES))
		{
			tminf_log1 = (SrtLogTmInf *) (stp->tminf[0]);
		
		/* For non IRM, the discount factor can be set here for all nodes */
			if (stp->next)  
				node.df = ((SrtLogTmInf *) (stp->next->tminf[0]))->df / tminf_log1->df;		
		}
		
		if (mdl_type2 == LGM)
		{
			tminf_irm2 = (SrtIRMTmInf *) (stp->tminf[1]);
		    sam_get(node.cur_sam, 1, PHI) = sam_get(tminf_irm2->fwd_sam, 1, PHI);  
		
		/* Quanto correction (including the Delta t) */
			node.quanto_fwd = tminf_irm2->quanto_adjustment * tminf_irm2->ev.onef.sig * stp->delta_t;
		}
	    else if ((mdl_type2 == BLACK_SCHOLES) || (mdl_type2 == EQ_STOCH_RATES))
		{
			tminf_log2 = (SrtLogTmInf *) (stp->tminf[1]);
		
		/* Quanto correction (including the Delta t) */
			node.quanto_fwd = tminf_log2->quanto_adjustment * tminf_log2->int_sig_dt;
		}

	/* FIRST loop on the first state variable */
		for (i = - trinf->max_index[0]; i <= trinf->max_index[0]; i++) 
		{
		/* Second loop on the second state variable */
			for (j = - trinf->max_index[1]; j <= trinf->max_index[1]; j++) 	 
    		{
			/* Change in real coordinates */
				y1 = i * trinf->spacing[0];
				y2 = j * trinf->spacing[1];
				x1 = trinf->dual_basis[0][0] * y1 + trinf->dual_basis[0][1] * y2;
				x2 = trinf->dual_basis[1][0] * y1 + trinf->dual_basis[1][1] * y2;

			/* Sets the value of the FIRST state variable and rebuilds the underlying process value */    	
				if (mdl_type1 == LGM)	    		
				{
					sam_get(node.cur_sam, 0, STATEVAR) = x1;
					sam_get(node.cur_sam, 0, SHORT_RATE) = sam_get(node.cur_sam, 0, STATEVAR) 
										+ tminf_irm1->rf.onef.H + sam_get(tminf_irm1->fwd_sam, 0, F_0_t);
				}		
				else if ((mdl_type1 == BLACK_SCHOLES) || (mdl_type1 == EQ_STOCH_RATES))	
				{		
					sam_get(node.cur_sam, 0, STATEVAR) = x1;		
					sam_get(node.cur_sam, 0, SPOT) = tminf_log1->ajust_init_fwd_val * exp(x1);				
				}

			/* Sets the value of the SECOND state variable and rebuilds the underlying process value */    	
				if (mdl_type2 == LGM)		    
				{
                    sam_get(node.cur_sam, 1, STATEVAR) = x2;
					sam_get(node.cur_sam, 1, SHORT_RATE) = sam_get(node.cur_sam, 1, STATEVAR) 
												+ tminf_irm2->rf.onef.H	+ sam_get(tminf_irm2->fwd_sam, 1, F_0_t);
				}		    
				else if ((mdl_type2 == BLACK_SCHOLES) || (mdl_type2 == EQ_STOCH_RATES))		    
				{
					sam_get(node.cur_sam, 1, STATEVAR) = x2;
					sam_get(node.cur_sam, 1, SPOT) = tminf_log2->ajust_init_fwd_val * exp(x2);
				}

			/* If it is not the last step, populates the node and computes expectation */
				if (stp->next)
      		    {
				/* Compute forward, connections and probas */
					populate_twound_tree_node(mdl_type1, mdl_type2, stp, &node, stp->delta_t, stp->next->trinf);

      	      	/* Sets the discount factor in the node if an IRM (if not, set above) */
					if (mdl_type1 == LGM)
					{
					/* Adds and Removes the H function to the state var for the cal to Y_t_at_T... */
						node.df = exp(- sam_get(node.cur_sam, 0, SHORT_RATE) * stp->delta_t);
					}
					
				/* Compute discounted asset prices at current node */
	        		twound_tree_expectation(&node, next_assets, cur_assets[i][j], node_dim);
				}
				if (mdl_type1 == LGM)
				{
					sam_get(node.cur_sam, 0, STATEVAR) +=tminf_irm1->rf.onef.H;
				}
				if (mdl_type2 == LGM)
				{
					sam_get(node.cur_sam, 1, STATEVAR) +=tminf_irm2->rf.onef.H;
				}

			/* Evaluate cash flows in Grfn tableau */
				if (err = evalcf((GrfnEvent *) stp->e, 
						&node.cur_sam, 
						gd, 
						(double*) cur_assets[i][j],
						(EvalEventDfsFct) srt_f_calc_grfn_event_dfs, 
						und_info, 
						&cashflow))
				{
					FREE_MEM_ASSETS; 
					return err;
				}

			} /* END of loop on second underlying (j) */
		
		} /* End of loop on first underlying (i) */

		/* Swaps the assets: the new one becomes the old one */
		SWAP(cur_assets, next_assets);
	
	} /* END of backward induction */

/* Stores the premium in the Input/Output list */
	err = srt_f_IOstructsetpremium((iolist), SRT_NO, next_assets[0][0][grfnparam->node_dim-1], "Two Und Tree");
	
/* Stores all the columns PV in the I/O list */
	err = srt_f_IOstructsetcolpvs((iolist), SRT_NO, (double*) next_assets[0][0], grfnparam->node_dim, ""	);
     
/* Free the SrtStp list duplicated from stp */
	FREE_MEM_ASSETS;	

	return err;

} /* END Err srt_f_twoundtree(...) */

#undef FREE_MEM_ASSETS 
#undef ALLOC_MEM_ASSETS 
#undef SWAP

/* -------------------------------------------------------------------------- */





