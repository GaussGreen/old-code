/******************************************************************************/
/*                                                                            */
/*      SYSTEM:         SRT     SORT, Fixed Income 2020 Addins                */
/*      SUB_SYSTEM:     SWT     Swap Tools                                    */
/*                                                                            */
/*      MODULE NAME:    SRT_F_VEGATRELGM.C                                    */
/*                                                                            */
/*      PURPOSE:        Functions to compute lgm tree                         */
/*                                                                            */
/*      AUTHORS:        G.Amblard, E.Auld, K.Chau, J.Malhi, A. Sahuguet       */
/*                                                                            */
/******************************************************************************/
/*                      Amendment History                                     */
/******************************************************************************/
/*                                                                            */
/*      AMENDED BY:     Jasbir S Malhi                                        */
/*                                                                            */
/*      DATE:           March 1995                                            */
/*                                                                            */
/*      REASON:         add greek functions                                   */
/*                                                                            */
/******************************************************************************/
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_lgmtreefct.h"
#include "srt_h_cheydynamics.h"

/* -------------------------------------------------------------------------- */

#define SWAP(X,Y){\
  	void *_VPTR;\
  	_VPTR 	= X;\
  	X 	= Y;\
  	Y 	= _VPTR;}

/* -------------------------------------------------------------------------- */

static SRT_Boolean is_within(
			Ddate 		today, 
			double 		left, 
			double 		date_or_time, 
			double 		right);


/* -------------------------------------------------------------------------- */
/* This function determines whether a request can be computed by the function
or not. Each computable request will be computed in srt_f_vegashifttrelgm1d. */
  

static SRT_Boolean srt_f_vegatrelgm1d_requestOK(
			int 		request_type);

/* -------------------------------------------------------------------------- */
/* This function defines the dimension of the tree using the parameters below.
However, when we deal with sigma shifts, the volatility at each date may be
different from one shifted price to another. This would lead to different
values of the u-parameter ( rate spacing between the nodes ). 
    see LGM_calc_delta_x() function
We want the tree to be frozen: the u is chosen constant for each shifted price.
The first time the function is called, it's for the premium - the real price,
without any shift. The value u computed then is stored and next calls to this
function will return this particular u instead of recomputing it.
The first time, *old_u = -1.						    */
 
static Err srt_f_vegashiftlgm1dtrelim
			(
  			SrtStpPtr 	stp,
  			SrtTreStpInf 	*maxtrinf,
  			double		*old_u);

/* -------------------------------------------------------------------------- */
/*This function computes the F, G, H and Psi for each step.		*/

static Err srt_f_vegashiftlgministp
			(
			SrtUndPtr   und,
			SrtUndInfo  *undinfo,
			SrtStpPtr 	stp,
			TermStruct 	*ts,
			double 		sigma_bucket_start, 
			double 		sigma_bucket_end,
			double 		sigma_shift,
			int 		sigma_shift_type,
			double 		tau_bucket_start, 
			double 		tau_bucket_end,
			double 		tau_shift,
			int 		tau_shift_type);

/* -------------------------------------------------------------------------- */

/* This function is the evolution of srt_f_trelgm1d in order to compute 
sigma and tau shifted prices with only 1 sweep of the tree  */

Err srt_f_vegashifttrelgm1d
	(	
	SrtUndPtr 	   und, 
	SrtGrfnParam   *grfnparam, 
	SrtStpPtr 	   stp, 
	GrfnDeal	   *gd,
	EvalEventFct   evalcf, 
	SrtIOStruct   *iolist,	
	SrtUndInfo    *und_info 
	)
{

/* ----- declare stuff ---- */
Err 		err = NULL;
int 		i; 
double 		zc_yield;
double		cashflow;
SrtStpPtr 	top,
		bot;

/* information about the  particular  node within the tree we are at */

SrtTrinTreNdInf node;

/* information about the dimensions of this tree */

SrtTreStpInf 	maxtrinf, *trinf;

/* information needed for the model at a particular time */

SrtIRMTmInf 	*nxttminf, *tminf;

/* val of assets being computed through the tree at cur and prev time step */
/* 3D arrays: r dimension, grfn_cols dimension, shift dimension */ 

double  	***prev_assets, 
		***cur_assets;

/* state variables at current time step and previous time step */

double 		*cur_x, 
		*prev_x;

/* amt of space for assets */

long 		sz;

/* current and constant value for rate spacing in the tree, for the prices
and all the schifted prices */

double 		old_u = -1;


/* number of extra calculation = 1 (for the price) +  
			 number of shifted values requested.
			The price can always be computed  */

int 		xtra_calc_num = 0;	

/* in order to know which shift corresponds to which calc.value */
                 
int 		xtra_calc_index; 

SrtStpPtr 	my_pointer;
SrtIOVal	*io_request_val_p;
SrtListAtom	*io_request; 


/* --------------- end of declarations --------------------------------------- */

	top = stp = gototop(stp);
	bot = gotobot(stp);



/*** Attach interesting und and mdl info to steps ***/
	
/* Initialise steps using the real SrtStpPtr that will be used in the tree
   to calculate premium  */

	if (err = srtstptminfalloc(stp,1))
	{
		return err;
	}
	if (err = srt_f_irministp(stp,und,0, und, und_info))
	{
		return err;
	}

	if (err =  srt_f_vegashiftlgm1dtrelim(stp, &maxtrinf, &old_u))
	{
 		return err;
	}
                                                                
/* Recreate SrtStpPtr for each greek request and add ptrs to io_request */

  	for (io_request = iolist->head;
			io_request != NULL;
			io_request = io_request->next)
	{
		io_request_val_p = (SrtIOVal*)(io_request->element->val.pval);

		if (srt_f_vegatrelgm1d_requestOK(io_request_val_p->type) == SRT_NO) 
	   {
	  	continue;
	   }
		if (io_request_val_p->type == IO_PREMIUM)
		{
		   ((SrtIOVal*)(io_request->element->val.pval))->pval 
					= (void*) stp;
		}  
		else
		{
                   if (err = srt_f_stpdup(stp, &my_pointer))
                   {
                        return err;
                   }
                   ((SrtIOVal*)(io_request->element->val.pval))->pval 
					= (void*) my_pointer;
                }
	}


/* For each requested greek calculation */


  	for (	io_request = iolist->head;
			io_request != NULL;
			io_request = io_request->next)
  	{

  	   io_request_val_p = (SrtIOVal*)(io_request->element->val.pval);

/* find out whether or not the requested value can be computed by this tree function */

  	   if (srt_f_vegatrelgm1d_requestOK(io_request_val_p->type) == SRT_NO) 
	   {
	  	continue;
	   }
	   else
	   {
		xtra_calc_num++;

	
		switch (io_request_val_p->type)
		{
		   case IO_PREMIUM:
			break;
		   case IO_SIGMA_SHIFT:
		   	my_pointer = (SrtStpPtr)(((SrtIOVal*)(io_request->element->val.pval))->pval); 

		   	if (err =  srt_f_vegashiftirministp(
			    		my_pointer, 
					und,
					und_info,
			    		io_request_val_p->bucket_start, 
					io_request_val_p->bucket_end,
					io_request_val_p->shift_value, 
					io_request_val_p->shift_type,
			    		-1, 
					-1,
					0, 
					SH_NONE  ))
		   	{
				return err;
		   	}
			break;
		   case IO_TAU_SHIFT:
		   	my_pointer = (SrtStpPtr)(((SrtIOVal*)(io_request->element->val.pval))->pval); 

		   	if (err =  srt_f_vegashiftirministp(
			    		my_pointer, 
					und,
						und_info,
			    		-1, 
					-1,
					0, 
					SH_NONE,
			    		io_request_val_p->bucket_start, 
					io_request_val_p->bucket_end,
					io_request_val_p->shift_value, 
					io_request_val_p->shift_type ))
		   	{
				return err;
		   	}
			break;
		  default:
			break;
		}/* END switch statement */ 

/* each request of the SrtIOStruct has now a field that points to a SrtStp structure.
For this structure, the parameters have been computed with the shift corresponding to the request */


	   }/* END srt_f_vegatrelgm1d_requestOK  == SRT_YES */ 
	   

 	} /* END of the for loop on io_request */

/* determine how much space we need to allocate,etc;
   attach tree geometry info to steps */

	

/* Memory Allocation -------------------------------------*/

  	sz = xtra_calc_num*(maxtrinf.max_x_index - maxtrinf.min_x_index+1)*
       		grfnparam->node_dim;

  	prev_assets = f3tensor(0,xtra_calc_num-1,
				maxtrinf.min_x_index,maxtrinf.max_x_index, 
      				0,grfnparam->node_dim-1);
  	cur_assets  = f3tensor(0,xtra_calc_num-1,
				maxtrinf.min_x_index,maxtrinf.max_x_index, 
      				0,grfnparam->node_dim-1);

  	prev_x = dvector(maxtrinf.min_x_index,maxtrinf.max_x_index);
  	cur_x = dvector(maxtrinf.min_x_index,maxtrinf.max_x_index);

  	if(!prev_assets || !cur_assets || !prev_x || !cur_x)
  	{
    		return serror("allocation failure srt_f_vegatrelgm.c");
  	}

  	memset(&cur_assets[0][maxtrinf.min_x_index][0],0,sz*sizeof(double));
/* Set to zero all the cell of the array */


/* Tree (going backwards in time) ---------------------------*/

  	for(stp = bot;stp != NULL; stp = stp->prev) 
  	{
    	   trinf = (SrtTreStpInf *)stp->trinf;
    	   tminf = (SrtIRMTmInf *)stp->tminf[0];

/** create grid of state variables **/
    	   cur_x[trinf->min_x_index]=trinf->xmin;
    	   for(i=trinf->min_x_index+1;i<=trinf->max_x_index;i++)
    	   {
      		cur_x[i] = cur_x[i-1] + trinf->u;
    	   }


/** For each x (i.e. for each node in the tree) **/

    	   for(i= trinf->min_x_index; i <= trinf->max_x_index; i++) 
    	   {
/* STATEVAR corresponds to r(t) - f(0,t)  and  X to SV - H(t) */
      		sam_get(node.cur_sam,0,STATEVAR) = cur_x[i] 
						+ tminf->rf.onef.H; 

/*    	
sam_get(tminf->fwd_sam,0,F_0_t) + tminf->H  +  ...
*/
	      	sam_get(node.cur_sam,0,SHORT_RATE) = 
	       				sam_get(node.cur_sam,0,STATEVAR)  
				      + sam_get(tminf->fwd_sam,0,F_0_t);

      		sam_get(node.cur_sam,0,PHI) = sam_get(tminf->fwd_sam,0,PHI);

     		xtra_calc_index = 0;

    		for (io_request = iolist->head;
	 			io_request != NULL;
	 			io_request = io_request->next) 
    		{
     		   io_request_val_p = 
				(SrtIOVal*)(io_request->element->val.pval);
    		   if (srt_f_vegatrelgm1d_requestOK(
					io_request_val_p->type) 
				== SRT_NO )
		   {
			continue;
		   } 

		   my_pointer = gotoindex( (SrtStpPtr)(io_request_val_p->pval), 
					stp->index);	
        	   tminf = (SrtIRMTmInf *)my_pointer->tminf[0];
		   sam_get(node.cur_sam,0,PHI) = sam_get(tminf->fwd_sam,0,PHI);

/** if not at the end: **/
      		   if(stp->next)
      		   {

        		nxttminf = (SrtIRMTmInf *)my_pointer->next->tminf[0];

/** compute discount factor for that stp **/ 

/*  Y_T_at_t_compute uses node->cur_sam...[STATEVAR] which correponds to
    r(t) - f(0,t)  that is X(t) + H(t) */
        		Y_T_at_t_compute( 1, 
                      			&node.cur_sam,
					&tminf->yp, 
                      			&zc_yield,
					0,
					ONE_FAC,
					LGM);

        		node.df = exp(- zc_yield );

/** calculate drift at this node **/

        		sam_get(node.drift_sam,0,PHI) 
					= sam_get(nxttminf->fwd_sam,0,PHI);
/*        		sam_get(node.drift_sam,0,SHORT_RATE) 
				= sam_get(node.cur_sam,0,SHORT_RATE) 
				+ sam_get(nxttminf->fwd_sam,0,F_0_t) 
				- sam_get(tminf->fwd_sam,0,F_0_t) 
				+ nxttminf->H - tminf->H 
				- tminf->lambda * my_pointer->delta_t 
					* (sam_get(node.cur_sam,0,SHORT_RATE) 
					- sam_get(tminf->fwd_sam,0,SHORT_RATE)
						- tminf->H);
*/
/*  The tree is built on X = r(t)-f(0,t)-H(t)
    Therefore, when drift_sam is used, STATEVAR corresponds to X  - 
    not to r-f(0,t) */
        		sam_get(node.drift_sam,0,STATEVAR) 
				= cur_x[i]* (1 - tminf->ev.onef.lambda 
							* my_pointer->delta_t);

/** calculate x index to center on in previous level **/

/*
        		srt_f_chetrerindex(my_pointer,prev_x,
						&node,mdl,CLOSESTINX);
*/
/* Modified node.mid_son_index  - Please note that this function uses
   only node->drift_sam...[STATEVAR] which is consistent with the previous
   assumption that the tree is build on X, not really on r(t)-f(0,t)*/
        		srt_f_trintrexindex(my_pointer,prev_x,
						&node,grfnparam,CLOSESTINX);

/** calculate variance at this node **/

        		node.var_at_sam = tminf->ev.onef.sig2 
							* my_pointer->delta_t;
/* CHECK if not better putting stdev_r */
/** calculate probabilities **/

/*
        		srt_f_chetreprob(my_pointer,prev_x,&node,mdl);
*/
        		srt_f_trintreprob(prev_x,&node, grfnparam);

/** calculate discounted prob weighted sum of previous cashflows **/

        		srt_f_trintredcntvec( my_pointer,
					      cur_assets[xtra_calc_index][i],
					      &node,
					      prev_assets[xtra_calc_index],
       					      grfnparam->node_dim);


      		   } /* END if(stp->next)  */

/** call cashflow function **/

      		   err = evalcf(
							(GrfnEvent *)my_pointer->e, 
							&node.cur_sam, 
							gd,
							(double *)cur_assets[xtra_calc_index][i],
							(EvalEventDfsFct)srt_f_calc_grfn_event_dfs,
							und_info,
							&cashflow);
			   if (err) return err;

    		   xtra_calc_index++;

    		} /* END of shift loop on io_request*/

	   }/* END of i loop from min_r_index to max_r_index */

    	   SWAP(cur_assets,prev_assets);
    	   SWAP(cur_x,prev_x);

  	}/* END of step loop from bot to top */

/* Stores all the requests in the IO list */
	xtra_calc_index = 0;
	for (io_request = iolist->head; io_request != NULL ; io_request = io_request->next)
	{
		io_request_val_p = (SrtIOVal*)(io_request->element->val.pval);
		if ( srt_f_vegatrelgm1d_requestOK(io_request_val_p->type) == SRT_NO ) 
		{
			continue;
		}
	
		strncpy(io_request_val_p->computation_origin, "Tree_LGM_1d", strlen("Tree_LGM_1d"));
		io_request_val_p->dval = prev_assets[xtra_calc_index][0][grfnparam->node_dim-1]; 

	/* free the SrtStp list duplicated from stp */
		if (io_request_val_p->type!=IO_PREMIUM) 
			free_list(io_request_val_p->pval);

	 	   xtra_calc_index++;
	} /* END of store all requests loop*/
                        
/* Stores all the columns PV in the I/O list */
	err = srt_f_IOstructsetcolpvs(	
				(iolist),
				SRT_NO,
				(double*) prev_assets[0][0],
				grfnparam->node_dim,
				""	);
	
/* Free memory */
	free_f3tensor(prev_assets,0,xtra_calc_num-1, 
			maxtrinf.min_x_index,maxtrinf.max_x_index, 
      			0,grfnparam->node_dim-1);
  	free_f3tensor(cur_assets, 0, xtra_calc_num-1, 
			maxtrinf.min_x_index,maxtrinf.max_x_index, 
      			0,grfnparam->node_dim-1);

  	free_dvector(prev_x, maxtrinf.min_x_index,maxtrinf.max_x_index);
  	free_dvector(cur_x, maxtrinf.min_x_index,maxtrinf.max_x_index);

/* Return a success or error message */
	return err;
}

/* Problemo:
probs are computed each time although sometimes it is useless.
*/




static SRT_Boolean is_within(
			Ddate 		today, 
			double 		left, 
			double 		date_or_time, 
			double 		right)
{
	double 		t;

/* if today = -1, the date is a real date -> no conversion */

	if (today == -1) 
	{
		t = date_or_time;
	}

/* if today != -1, the date is actually a time -> conversion into a date */

	else
	{
		t = (date_or_time / YEARS_IN_DAY) + today;
	} 
/* Anyway, t is now a date. */

	if (
		( (t>=left)  || (fabs(t-left) < EPSILON) )
		&& ( (t<=right) || (fabs(t-right) < EPSILON) )
   		) 
	{
		return SRT_YES;
	}
	else 
	{
		return SRT_NO;
	}
} 

/* -------------------------------------------------------------------------- */



static Err srt_f_vegashiftlgm1dtrelim
		(
  		SrtStpPtr 	stp,
  		SrtTreStpInf 	*maxtrinf,
  		double		*old_u	/* previous value of rate spacing u */	
)
{
  	Err 		err;
  	SrtStpPtr 	cur;
  	SrtTreStpInf 	*trinf;
  	SrtIRMTmInf 	*tminf;
  	double 		delta_x;
  
  	cur = gototop(stp);

 /* allocate space */

  	if (err = srtstpalloc(cur,sizeof(SrtTreStpInf),1,0))
  	{
		return err;
	}

  	if (*old_u == -1) 
	{ 
		LGM_calc_delta_x(cur,&delta_x, 0); 
		*old_u = delta_x; 
	}
  	else 
	{
		delta_x = *old_u;
	}

  	while(cur->next)
  	{
    		trinf = (SrtTreStpInf *)cur->trinf;
    		tminf = (SrtIRMTmInf *)cur->tminf[0];

    		trinf->min_x_index = -cur->index;
    		trinf->max_x_index =  cur->index;
    		trinf->u = delta_x;
    		trinf->xmin = - cur->index * delta_x ;   
/* trinf->rmin = sam_get(tminf->fwd_sam,0,F_0_t) + tminf->H  -  cur->index...;*/
    		cur = cur->next;
  	}
  	trinf = (SrtTreStpInf *)cur->trinf;
  	tminf = (SrtIRMTmInf *)cur->tminf[0];

  	trinf->min_x_index = -cur->index;
  	trinf->max_x_index =  cur->index;
  	trinf->u = delta_x;
    	trinf->xmin = - cur->index * delta_x ;   

/* Modify the call to TmInf */ 
/*trinf->rmin = sam_get(tminf->fwd_sam,0,F_0_t) + tminf->H  -  cur->index...; */

  	*maxtrinf = *trinf;
  	return NULL;
}

/* -------------------------------------------------------------------------- */

static Err srt_f_vegashiftlgministp
		(	
		SrtUndPtr    und,
		SrtUndInfo   *und_info,
		SrtStpPtr 	stp,
		TermStruct 	*ts,
		double 		sigma_bucket_start, 
		double 		sigma_bucket_end,
		double 		sigma_shift,
		int 		sigma_shift_type,
		double 		tau_bucket_start, 
		double 		tau_bucket_end,
		double 		tau_shift,
		int 		tau_shift_type  )
{   
  	SrtIRMTmInf *prvtminf;
	SrtListAtom	*ls;
	IrmTermStructVal	*tsval;
	Err err = NULL;

  	prvtminf = NULL;

/* According to the current shift (rate, sigma or tau)
we change -- temporarily -- the termstructure in order to 
compute G, H and Psi. */

	for (ls = ts->head; ls!= NULL; ls = ls->next)
	{
		tsval = (IrmTermStructVal*) (ls->element->val.pval);
		switch (sigma_shift_type)
		{
		  case SH_NONE:
			break;
		  case SH_ABSOLUTE:
			if ( is_within(	-1, sigma_bucket_start,
					tsval->date,
					sigma_bucket_end	) == SRT_YES )
				tsval->sig +=  sigma_shift;
			break;  
 	 	  case SH_RELATIVE:
			if ( is_within(	-1, sigma_bucket_start,
					tsval->date,
					sigma_bucket_end	) == SRT_YES )
				tsval->sig *= 1 + sigma_shift;
			break;
		  default:
			break;
		}/* end of switch */

		switch (tau_shift_type)
		{
		  case SH_NONE:
			break;
		  case SH_ABSOLUTE:
			if ( is_within(	-1, tau_bucket_start,
					tsval->date,
					tau_bucket_end	) == SRT_YES )
				tsval->tau +=  tau_shift;
			break;  
 	 	  case SH_RELATIVE:
			if ( is_within(	-1, tau_bucket_start,
					tsval->date,
					tau_bucket_end	) == SRT_YES )
				tsval->tau *= 1 + tau_shift;
			break;
		  default:
			break;
		}/* end of switch */

	} /* end of for ls=ts->head  loop */

/* ----------------------------------------------- */
/* only sigma and tau shifts have been implemented */
/* ----------------------------------------------- */

	srt_f_tsupdate(&ts);
/* the termstructure is updated since tau and sigma may have been shifted */


/* Now we need to update the F, G, H, Psi
   functions, according to the new values of sig and tau 
   (that might have been shifted...)
   This update is done using ts (that has just been updated...) 
*/
	if (err= srt_f_lgministp(stp,ts,0,ONE_FAC, und, und_info) )
		return(err);


/* We set the termstruct to its original value (no shift) */ 

	for (ls = ts->head; ls!= NULL; ls = ls->next)
	{
	   tsval = (IrmTermStructVal*) (ls->element->val.pval);
	   switch (sigma_shift_type)
	   {
		case SH_NONE:
			break;
		case SH_ABSOLUTE:
			if ( is_within( -1, sigma_bucket_start,
					tsval->date,
					sigma_bucket_end	) == SRT_YES )
				tsval->sig -=  sigma_shift;
			break;  
 	 	case SH_RELATIVE:
			if ( is_within(	-1, sigma_bucket_start,
					tsval->date,
					sigma_bucket_end	) == SRT_YES )
				tsval->sig /= 1 + sigma_shift;
			break;
		default:
			break;
	   }/* end of switch for sigma shift */

	   switch (tau_shift_type)
	   {
		case SH_NONE:
			break;
		case SH_ABSOLUTE:
			if ( is_within(	-1, tau_bucket_start,
					tsval->date,
					tau_bucket_end	) == SRT_YES )
				tsval->tau -=  tau_shift;
			break;  
 	 	case SH_RELATIVE:
			if ( is_within(	-1, tau_bucket_start,
					tsval->date,
					tau_bucket_end	) == SRT_YES )
				tsval->tau /= 1 + tau_shift;
			break;
		default:
			break;
	   }/* end of switch for tau shift */	
	} /* END of ls->head to tail loop*/
	
/* the termstructure is updated since tau and sigma may have been shifted */

	srt_f_tsupdate(&ts);

  	return NULL;       
}


/* -------------------------------------------------------------------------- */

/* initialises step list for tree */
Err srt_f_vegashiftirministp	(SrtStpPtr 	stp,					/* output */
								 SrtUndPtr 	und,
								 SrtUndInfo *und_info,              /* und info */
								 double 	sigma_bucket_start,		/* shift info */
								 double 	sigma_bucket_end,
								 double 	sigma_shift,
								 int 		sigma_shift_type,
								 double 	tau_bucket_start,
								 double 	tau_bucket_end,
								 double 	tau_shift,
								 int 		tau_shift_type)

{
  	SrtStpPtr 		top;
  	Err 			err		=	NULL;
  	SrtIRMTmInf 	*tminf, 
					*prvtminf;
  	TermStruct 		*ts;
  	Ddate 			today;
  	String 			yc_name,
					und_name;
	SrtMdlType      mdl_type;
	SrtCurvePtr     yldcrv;

	/* get info from und */
  	prvtminf =	NULL;
  	und_name =	get_underlying_name (und);
  	err = get_underlying_ts (und, &ts);
	if (err)
	{
		return err;
	}
	err = get_underlying_mdltype (und, &mdl_type);
	if (err)
	{
		return err;
	}
  	yc_name = get_ycname_from_irund(und) ;
  	yldcrv = lookup_curve(yc_name) ;
  	today = get_clcndate_from_yldcrv(yldcrv);
  	top = stp = gototop(stp);

	/* allocate space */
  	if (err = srtstptminfalloc (top, 1))
  	{
		return err ;
  	}
  	if (err=srtstpalloc (top, sizeof (SrtIRMTmInf), 0, 0))
  	{
		return err;
  	}

	/* populate srttimestps */
  	while (stp->next)
  	{
    	tminf = stp->tminf[0];
		tminf->df = swp_f_df (today, stp->ddate, yc_name);
    	tminf->ev.onef.sig = sqrt (find_sig2_interp (stp->time, stp->next->time, ts));
		
		/* shift */
		/* we perform the shift for sigma if we are in the right bucket */
		if (is_within (today, sigma_bucket_start, stp->next->time, sigma_bucket_end) == SRT_YES)
		{
        	switch (sigma_shift_type)
			{
				case SH_NONE:
					break;
				case SH_ABSOLUTE:
					tminf->ev.onef.sig +=  sigma_shift;
					break;  
  				case SH_RELATIVE:
					tminf->ev.onef.sig *= 1 + sigma_shift;
					break;
				default:
					break;
			}
		}

    	tminf->ev.onef.sig2 = tminf->ev.onef.sig * tminf->ev.onef.sig;
    	tminf->ev.onef.tau = find_tau (stp->next->time,ts);  

		/* we perform the shift for tau if we are in the right bucket */
		if (is_within (today, tau_bucket_start, stp->next->time, tau_bucket_end) == SRT_YES)
		{
			switch (tau_shift_type)
			{
				case SH_NONE:
					break;
				case SH_ABSOLUTE:
					tminf->ev.onef.tau +=  tau_shift;
					break;  
  				case SH_RELATIVE:
					tminf->ev.onef.tau *= 1 + tau_shift;
					break;
				default:
					break;
			}    
		}

    	tminf->ev.onef.lambda = 1.0 / tminf->ev.onef.tau;

		/* This is needed  only for cheyette */
    	if (mdl_type == CHEY)
  			{
				/* starting from step 2 */
				if (prvtminf)
    			{
					/* initialises the deterministic part of exp(r), exp(X) and exp(phi) */
      				srt_f_chedrfatsam (stp->prev, &(prvtminf->fwd_sam), &(tminf->fwd_sam), 0);
    			}
    			else
    			{
      		       	sam_get (tminf->fwd_sam, 0, PHI) = 0.0;
    			}
			}
			
			/* get IFR */
    		sam_get (tminf->fwd_sam, 0, F_0_t) = swp_f_zr (stp->ddate, stp->next->ddate, yc_name);

			/* get fwd disc to next date */
    		if (err = Y_T_at_t_param (stp->ddate, &stp->next->ddate, 1, und, &tminf->yp))
			{
				return err;
			}
		
    		prvtminf = tminf;
    		stp = stp->next;
  		}

	/* last stp */
  	tminf = stp->tminf[0];
  	tminf->df = swp_f_df(today, stp->ddate, yc_name);
  	tminf->ev.onef.sig = sqrt (find_sig2_interp(stp->time, stp->time+14.0, ts));
	if (is_within (today, sigma_bucket_start, stp->time, sigma_bucket_end) == SRT_YES)
	{
        switch (sigma_shift_type)
		{
			case SH_NONE:
				break;
			case SH_ABSOLUTE:
				tminf->ev.onef.sig +=  sigma_shift;
				break;  
  			case SH_RELATIVE:
				tminf->ev.onef.sig *= 1 + sigma_shift;
				break;
			default:
				break;
		}
	}
  	tminf->ev.onef.sig2 = tminf->ev.onef.sig * tminf->ev.onef.sig;
  	tminf->ev.onef.tau = find_tau(stp->time,ts);  
	if (is_within (today, tau_bucket_start, stp->time, tau_bucket_end) == SRT_YES)
	{
		switch (tau_shift_type)
		{
			case SH_NONE:
				break;
			case SH_ABSOLUTE:
				tminf->ev.onef.tau +=  tau_shift;
				break;  
  			case SH_RELATIVE:
				tminf->ev.onef.tau *= 1 + tau_shift;
				break;
			default:
				break;
		}
	}
  	tminf->ev.onef.lambda = 1.0/tminf->ev.onef.tau;
 	if (mdl_type == CHEY)
	{
  		if (prvtminf)
  		{
    		srt_f_chedrfatsam (stp->prev, &(prvtminf->fwd_sam), &(tminf->fwd_sam), 0);
  		}
  		else
	  	{
    		sam_get (tminf->fwd_sam, 0, PHI) = 0.0;
	  	}
	}
  	sam_get (tminf->fwd_sam,0,F_0_t) = swp_f_zr(stp->ddate, stp->ddate + 14.0, yc_name);

/* if this is the linear gauss markov model (which we regard as a
   subset of the cheyette model) then we can precompute phi and 
   several other variables, which we now do */
  	if( mdl_type == LGM)
  	{
    	if(err=srt_f_vegashiftlgministp		(und,
										und_info,
											top,ts,
											sigma_bucket_start,
											sigma_bucket_end,
											sigma_shift,
											sigma_shift_type,
											tau_bucket_start,
											tau_bucket_end,
											tau_shift,
											tau_shift_type))
		{
			return err;
		}
  	}
  	return err;
}

/* -------------------------------------------------------------------------- */
/* this function must answer the following question:
Can the function compute the request with type 'request_type ?
*/

static SRT_Boolean srt_f_vegatrelgm1d_requestOK(int request_type)

{
    	switch (request_type)
    	{
    	case IO_PREMIUM	:   return SRT_YES; break;
    	case IO_STDEV	:   return SRT_NO; break;
    	case IO_SIGMA_SHIFT :   return SRT_YES; break;
    	case IO_TAU_SHIFT	:   return SRT_YES; break;
    	case IO_RATE_SHIFT	:   return SRT_NO;  break;
    	default		:   return SRT_NO;  break;
    	}
}

