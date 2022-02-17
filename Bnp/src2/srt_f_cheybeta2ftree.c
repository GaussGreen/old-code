/* -------------------------------------------------------------
	 
   FILENAME: srt_f_cheybeta2dtree.c

   PURPOSE: the a new implementation of pink elephant (or Malliavin stratification - copyright BFL) 
					for the Cheyette beta  2D in GRFN

  					First step:  generation of diffusion with MC optimized methods 

					second step: triangulation of random mesh O(nlog(n))

					third step: computation of probabilities solving Lagrangien problem for quadrature formula
									(exact formulas)

					fourth step:   pricing

  					fifth  step:  splitting up for non markovian part

  MODIFICATION:											  

 ----------------------------------------------------------------------- */

#include "srt_h_all.h"
#include "srt_h_cheybeta2ftree.h"
#include "srt_h_cheybeta2ftreefct.h"
#include "math.h"

/* --------------------------------------------------------------------------------- */

#define MAGICAL_NUMBER  3.40

#define SWAP(X,Y)  {\
	        void *_VPTR;\
	       _VPTR = X;\
	        X = Y;\
	        Y = _VPTR;}

#define FREE_CHEYBETA2F_MCTREE_MEMORY {\
	for(t = 0; t < num_time_pts; t++)\
		free(grid[t]);\
	free(grid);\
    free(in.pointlist);\
    free(in.pointattributelist);\
	free_dmatrix(next_assets, 0, grfnparams->num_MCarlo_paths - 1 + NB_BOUND_PTS, 0, grfnparams->node_dim - 1);\
	free_dmatrix(cur_assets,  0, grfnparams->num_MCarlo_paths - 1 + NB_BOUND_PTS, 0, grfnparams->node_dim - 1);\
}


/* -------------------------------------------------------------------------------
   FUNCTION: Err srt_f_cheybeta2dtree(...)
   PURPOSE:  The main functions to call when pricing backward with Cheyette Beta 2f
   ------------------------------------------------------------------------------- */

Err srt_f_cheybeta2dtree(
	 SrtUndPtr       und,                  /* underlying  pointer*/
	 SrtGrfnParam   *grfnparams,    /*grfn parameters*/
	 SrtStpPtr       stp,                   /* step pointer*/
	 GrfnDeal       *gd,                 /* GRFN deal description*/
	 EvalEventFct    evalcf,             /* cash-flow evaluation function*/
	 SrtIOStruct    *iolist,               /* list of requests*/ 
	 SrtUndInfo     *und_info           /* underlying info*/
	 )
{
Err                 err = NULL;
SrtStpPtr           top, bot;		  	               /*Step pointers*/
SrtSample           sam;		 	                     /* sample for evalcf */
SrtIRMTmInf        *tminf = NULL;	  	        /*Info for a particular time step*/
SrtListAtom        *io_request = NULL;	  	  /* For the storage of the premium ... */
SrtIOVal           *io_request_val_p = NULL;

long                t, i, num_time_pts, k;			 /*Indexes*/
double              cash_flow,  final_result;		  	/*Cash-flow*/
/* Arrays used to store the Grfn columns for all paths*/
double            **cur_assets = NULL, 
                  **next_assets = NULL;

/* A massive grid to store all paths */
MCTreePoint       **grid = NULL;

/* For the triangulation algorithm */
SrtTriangulationIO   in;



/* -------------------- STEP 1: INITIALISATION ---------------------------- */

	/* Set top (1st step), bot (last step) and stp (set to top) */
	top = stp = gototop(stp);
	bot = gotobot(stp);

	/* Create indexes for the time steps */
	num_time_pts = create_index(top) + 1;

	/* Need to create tminf structure in SrtStpPtr first;
		size = 1,only the domestic rate as underlying */
	if (err = srtstptminfalloc(stp,1))     return err;

	/* Allocates and populates time_info at each step for a CHEYBETA 2F model:*/
	err = srt_f_irministp(top, und, 0, und, und_info);
	if (err)	return err;


/* -------------- STEP 2 : PATHS GENERATION AND STORAGE ------------------- */

/* Sets the number of paths as necessarily 3.4 * N * ln N (with N time steps) */
	grfnparams->num_MCarlo_paths = (long) floor(MAGICAL_NUMBER * num_time_pts * log(num_time_pts));
	
/* Always use SPECTRUNC for the MC Scheme (better compact Brownian description) */
	grfnparams->sample_type = SPECTRUNC;	

/* Memory allocation for the grid */
	grid = (MCTreePoint**) srt_calloc(num_time_pts, sizeof(MCTreePoint*));
	for(t = 0 ; t < num_time_pts; t++)
		grid[t] = (MCTreePoint*) calloc(grfnparams->num_MCarlo_paths + NB_BOUND_PTS, sizeof(MCTreePoint));

	/* parameters for triangulation structure */
	in.numberofpoints = grfnparams->num_MCarlo_paths + NB_BOUND_PTS;
	in.pointlist = (double *) malloc(in.numberofpoints * 2 * sizeof(double));
	in.numberofpointattributes = grfnparams->node_dim;
	in.pointattributelist = (double *) malloc(in.numberofpoints * in.numberofpointattributes * sizeof(double));
	in.pointmarkerlist = (int *) NULL;
	in.numberofsegments = 0;
	in.numberofholes = 0;
	in.numberofregions = 0;

	/* Generate Monte Carlo paths, and store everything */
	err = srt_f_cheybeta2ftree_make_grid(grfnparams, top, und_info, grid);
	if (err) 
	{
		srt_free(grid);
		return err;
	}


/* -------------------- STEP 3 : BACKWARD INDUCTION  --------------------------- */

	/* Memory allocation for the asset grids ( npaths * Grfn dim ) */
	next_assets = dmatrix(0, grfnparams->num_MCarlo_paths - 1 + NB_BOUND_PTS, 0, grfnparams->node_dim -1);
	cur_assets = dmatrix(0, grfnparams->num_MCarlo_paths - 1 + NB_BOUND_PTS, 0, grfnparams->node_dim - 1); 
	if (!next_assets || !cur_assets || !grid)	
		return serror ("Big memory allocation failure in srt_f_cheybeta2dtree");

	
	/* allocates memory for rescaling, etc  */
	srt_f_alloc_computing_resources();

	/* LOOP ON TIME STEPS */
	t = num_time_pts - 1;
	for (stp = bot; stp != NULL; stp = stp->prev)
	{
		if (stp->next)
		{					
			/*--------- TRIANGULATION ---------- */
			/* compute transition probabilities */
			srt_f_cheybeta2ftree_compute_proba(stp);

			/* triangulation of forward position  */
			transfer_info_nodes(grid[t+1], next_assets, &in);
			triangulate_srt("pczQ", &in);

			/* LOOP FOR THE LINEAR PART */
			for ( i = 0; i < grfnparams->num_MCarlo_paths; i++)
			{
				/* Compute discounted asset prices at current node */
				err = srt_f_cheybeta2ftree_discount_to_node (grid, t, i, stp, cur_assets[i], grfnparams);
				if (err)   return err;
			} /* END loop i on  nodes */ 
		} 


	/* LOOP FOR THE CONSTRAIN PART */
		for ( i = 0; i < grfnparams->num_MCarlo_paths; i++)
		{
			/* Evaluate cash flows in Grfn tableau using passed function evalcf */
			sam_get(sam, 0, X1) = grid[t][i].x1 + grid[t][i].jump_x1;
			sam_get(sam, 0, X2) = grid[t][i].x2 + grid[t][i].jump_x2;
    		sam_get(sam, 0, SHORT_RATE) = grid[t][i].short_rate;
			sam_get(sam, 0, PHI1) = grid[t][i].phi1;
			sam_get(sam, 0, PHI2) = grid[t][i].phi2;
			sam_get(sam, 0, CROSSPHI) = grid[t][i].phi1_2;  

			err = evalcf((GrfnEvent *) stp->e, 
					&sam, 
					gd, 
					(double*) cur_assets[i], 
					(EvalEventDfsFct) srt_f_calc_grfn_event_dfs, 
					und_info, &cash_flow);
			if (err)
			{
				FREE_CHEYBETA2F_MCTREE_MEMORY;
				return err;
			}

		} /* END loop i on  nodes */ 

	/* Special treatment for boundaries condition */
		for ( i = grfnparams->num_MCarlo_paths; i < in.numberofpoints; i++)
			for (k = 0; k < grfnparams->node_dim; k++)
				cur_assets[i][k] = DBL_MAX; 

	/* If not at the last time step: df from t to t + dt */
		if (stp->next)
			triangledeinit();	
		
	/* Swaps the assets from cur to next */	
		SWAP (cur_assets, next_assets);
		t--;
	} /* END of loop on all time steps (backward induction) */

	/* free_memory discounting */
	srt_f_free_computing_resources();


/* -------------------- STEP 4 : PRICE COMPUTATION  --------------------------- */

/* The Grfn price is in the last column (all the paths give the same price) */
	final_result =  next_assets[0][grfnparams->node_dim - 1];
	

/* Stores the premium in the Input/Output list */
	err = srt_f_IOstructsetpremium(iolist,SRT_NO,final_result,"");
	if (err) 
	{
		FREE_CHEYBETA2F_MCTREE_MEMORY;
		return (err);
	}

/* Stores all the columns PV in the I/O list */
	err = srt_f_IOstructsetcolpvs(	
				iolist,
				SRT_NO,
				(double*) next_assets[0],
				grfnparams->node_dim,
				""	);
	if (err)
	{
		FREE_CHEYBETA2F_MCTREE_MEMORY;
		return err;
	}

/* Free the memory allocated for the arrays of assets, the mc grid
	    and the parameters for triangulation structure */
	FREE_CHEYBETA2F_MCTREE_MEMORY;
	
	return err;

} /* END of	Err srt_f_cheybeta2dtree(...) */

/*--------------------------------- End of File -------------------------------------*/