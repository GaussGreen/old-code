/* Black-Karasinski trinomial tree

   Assume the following dynamics:

   d(ln(r)) = (theta(t) - lambda(t)*ln(r))*dt + sigma(t)*dW
 ==      dr = (theta(t)-lambda(t)*ln(r)+sigma^2(t)/2)*r*dt + sigma(t)*r*dW


   Still use the name bdt... in most functions because of ease of 
   identification.

   K L Chau   Sep 1995
*/
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_lgmtreefct.h"

#define SWAP(X,Y){\
	void *_VPTR;\
	_VPTR = X;\
	X = Y;\
	Y = _VPTR;}

static Err bdttrelim
		(
  		SrtStpPtr 	stp,
  		SrtTreStpInf 	*maxtrinf
)
{
  	Err 		err;
  	SrtStpPtr 	cur;
  	SrtTreStpInf 	*trinf;
  	SrtIRMTmInf 	*tminf;
  	double 		delta_x;
  
  	cur = gototop(stp);

  	if (err = srtstpalloc(cur,sizeof(SrtTreStpInf),1,0))
		return err;

	LGM_calc_delta_x(cur,&delta_x, 0); 

  	while(cur->next)
  	{
    		trinf = (SrtTreStpInf *)cur->trinf;
    		tminf = (SrtIRMTmInf *)cur->tminf[0];

    		trinf->min_x_index = -cur->index;
    		trinf->max_x_index =  cur->index;
    		trinf->u = delta_x;
/* tree is centred on the forward */
    		trinf->xmin = - cur->index * delta_x 
				+ log(sam_get(tminf->fwd_sam,0,F_0_t));
    		cur = cur->next;
  	}

/* we can ignore the last step because it is a dummy step used to generate
   the final dt for forward induction.  Hence trinf and tminf are not
   required for this step because the tree nodes would not be used */

  	*maxtrinf = *trinf;
  	return NULL;
}

static void bdtdrfatsam
(
  SrtStpPtr stp, 
  SrtSample *cur_sam,
  SrtSample *drift_sam,
  double theta 
)
{
  double mu;
  SrtIRMTmInf *tminf;

  if(!stp->next)return;

  tminf = (SrtIRMTmInf *)stp->tminf[0];

  mu   =  theta - tminf->ev.onef.lambda * samptr_get(cur_sam,0,STATEVAR) ;
  
  samptr_get(drift_sam,0,STATEVAR) = samptr_get(cur_sam,0,STATEVAR) 
					+ mu*stp->delta_t ;
  return;
}

static void insert_extra_step(SrtStpPtr bot) 
{
  double dt,last_date,last_time ;

  dt = bot->prev->delta_t ;
  last_date = bot->ddate ;
  last_time = bot->time ;

  bot = insert_after(bot) ;
  bot->ddate = last_date+dt*DAYS_IN_YEAR ;
  bot->time  = last_time+dt ;
  bot->prev->delta_t = dt ;
  bot->prev->sqrt_delta_t = sqrt(dt) ;
}

static void gen_x
(
  SrtStpPtr stp,    /* ptr to current step in lst */
  double *cur_x    /* space for rates to be stored */ 
)
{
  int i;
  double delta_x;
  SrtTreStpInf *trinf;
  SrtIRMTmInf *tminf;

  trinf = stp->trinf;
  tminf = stp->tminf[0];

  if(trinf->max_x_index==trinf->min_x_index)
  {
    cur_x[trinf->min_x_index] = log(sam_get(tminf->fwd_sam,0,F_0_t));
    return;
  }

/** generate X=ln(r) **/
  delta_x = trinf->u ;
  cur_x[trinf->min_x_index] = trinf->xmin ;
  for(i=trinf->min_x_index+1 ; i<=trinf->max_x_index; i++)
    cur_x[i] = cur_x[i-1] + delta_x ;
  return;
}

static double calc_adj(SrtIRMTmInf *tminf, SrtTreStpInf *trinf,
		       double dt1,
		       double dt2, double *x, double *Q,
		       double df_at_2dt)
{
  double theta,sumtop,sumbot ;
  double a,b,lambda,tmp,r,dt1dt2 ;
  int i ;

  sumtop = sumbot = 0.0 ;
  lambda = tminf->ev.onef.lambda ;
  dt1dt2 = dt1*dt2 ;
  a = lambda*dt1dt2 ;
  b = tminf->ev.onef.sig2*dt1dt2/2 ;

  if (trinf->min_x_index==trinf->max_x_index)
  {
    r = exp(x[0]) ;
    tmp = exp(-r*(dt1+dt2)) ;
    sumtop = tmp*(1+(a*log(r)-b)*r) ;
    sumbot = tmp*dt1dt2*r ;
  }
  else
    for (i=trinf->min_x_index ; i<=trinf->max_x_index ; i++)
    {
      r = exp(x[i]) ;
      tmp = Q[i]*exp(-r*(dt1+dt2)) ;
      sumtop += tmp*(1+(a*log(r)-b)*r) ;
      sumbot += tmp*dt1dt2*r ;
    }

  theta = (sumtop - df_at_2dt)/sumbot ;
  return theta ;  
}

static void accumulate_Q(SrtTrinTreNdInf *node,double Q, double *next_Q)
{
  long index ;

  index = node->mid_son_index ;
  next_Q[index-1] += Q*node->p[0]*node->df ;
  next_Q[index]   += Q*node->p[1]*node->df ;
  next_Q[index+1] += Q*node->p[2]*node->df ;

  return ;
}

static void generate_prob(SrtStpPtr stp, SrtTrinTreNdInf *node, 
			   double varatsam, double adj, double cur_x,
			   double *next_x)
{
  double dt,r ;

  dt = stp->delta_t ;
/* set up state variables */
  sam_get(node->cur_sam,0,STATEVAR) = cur_x ; 
  sam_get(node->cur_sam,0,SHORT_RATE) = r = exp(cur_x) ;
  node->df = exp(-r*dt) ;
/* calculate drift and variance at node */
  bdtdrfatsam(stp,&(node->cur_sam),&(node->drift_sam),adj) ;
  node->var_at_sam = varatsam ;
/* establish connectivity to next layer of steps */
  srt_f_trintrexindex(stp,next_x,node,NULL,CLOSESTINX) ;
/* calculate probabilities and store results in node */
  srt_f_trintreprob(next_x,node,0) ;

  return ;
}

static void calc_prob_and_Q(SrtStpPtr stp, double *cur_x,
			    double *cur_Q, double adj, double *next_x, 
			    double *next_Q)
{
  SrtTrinTreNdInf node ;
  double dt,varatsam ;
  int i ;
  SrtTreStpInf *trinf,*nxttrinf ;
  SrtIRMTmInf  *tminf ;

  trinf = stp->trinf ;
  tminf = stp->tminf[0] ;
  nxttrinf = stp->next->trinf ;
  dt = stp->delta_t ;
  varatsam = tminf->ev.onef.sig2*dt ;

/* set every Q to zero at the next level */
  for (i=nxttrinf->min_x_index ; i<=nxttrinf->max_x_index ; i++)
    next_Q[i] = 0.0 ;
/* for each node, if Q is not zero, calc probs at the node */
  for (i=trinf->min_x_index ; i<=trinf->max_x_index ; i++)
  {
/* ignore nodes which are not connected at all */
    if (cur_Q[i]!=0)
    {
      generate_prob(stp,&node,varatsam,adj,cur_x[i],next_x) ;
/* accumulate Q for the 3 branches that are connected to the node */
      accumulate_Q(&node,cur_Q[i],next_Q) ;
    }
  }
  return ;
}

Err srt_f_trebdt1f
	(
	SrtUndPtr und, /* the underlying with its parameters */
	SrtGrfnParam   *grfnparam, /* model parameters */
	SrtStpPtr stp, /* discretization of deal in time, wi/ events attached*/
	GrfnDeal	   *gd,/* deal descriptor structure */
	EvalEventFct evalcf, /* cashflow evaluator */
	SrtIOStruct *iolist,
	SrtUndInfo *und_info
	)
{
  Err err=NULL;
  int i, stp_num, max_x_index=0; 
  SrtStpPtr top,bot,nxtstp;

  SrtTrinTreNdInf node;
  SrtTreStpInf maxtrinf,*trinf ;
  SrtIRMTmInf *tminf,*nxttminf;
  double *cur_x, *next_x, *cur_Q, *next_Q, *theta ;
  double **cur_assets,**next_assets ;
  long sz1,sz2;
  double dt1,dt2,df_at_2dt,varatsam,cashflow,price ;

  SrtIOVal *io_request_val_p ;
  SrtListAtom *io_request ;

  top = stp = gototop(stp);
  bot = gotobot(stp) ;

/* attach one more step at the end of step list, for use in forward induction
   because the discounting factor at that step is required */
  insert_extra_step(bot) ;

  err = srtstptminfalloc(stp,1) ;
  err = srt_f_irministp(stp,und,0, und, und_info);
  if(err)return err;
  err = bdttrelim(stp,&maxtrinf);
  if (err) return err;
  
/** Memory Allocation **/

  cur_assets	= dmatrix(maxtrinf.min_x_index,maxtrinf.max_x_index,
		 	  0,grfnparam->node_dim-1) ;
  next_assets	= dmatrix(maxtrinf.min_x_index,maxtrinf.max_x_index,
		 	  0,grfnparam->node_dim-1) ;
  cur_x		= dvector(maxtrinf.min_x_index,maxtrinf.max_x_index);
  next_x	= dvector(maxtrinf.min_x_index,maxtrinf.max_x_index);
  cur_Q		= dvector(maxtrinf.min_x_index,maxtrinf.max_x_index);
  next_Q	= dvector(maxtrinf.min_x_index,maxtrinf.max_x_index);
  theta		= dvector(0,maxtrinf.max_x_index) ;

  if(!cur_assets || !next_assets || !cur_x || !next_x || !cur_Q || !next_Q)
    return serror ("alloc failure in fwd tree");
           
  sz1 = maxtrinf.max_x_index-maxtrinf.min_x_index+1 ;
  sz2 = sz1*grfnparam->node_dim;
  memset(&cur_x[maxtrinf.min_x_index],   0,sz1*sizeof(double)) ;
  memset(&next_x[maxtrinf.min_x_index],  0,sz1*sizeof(double)) ;
  memset(&cur_Q[maxtrinf.min_x_index],   0,sz1*sizeof(double)) ;
  memset(&next_Q[maxtrinf.min_x_index],  0,sz1*sizeof(double)) ;
  memset(&cur_assets[maxtrinf.min_x_index][0],0,sz2*sizeof(double)) ;
  memset(&next_assets[maxtrinf.min_x_index][0],0,sz2*sizeof(double)) ;

/* Start at the top of the tree */

  trinf      = stp->trinf ;
  tminf      = stp->tminf[0] ;
  nxtstp     = stp->next ;
  nxttminf   = nxtstp->tminf[0] ;
  df_at_2dt  = ((SrtIRMTmInf *) nxtstp->next->tminf[0])->df ;
  cur_x[0]   = log(sam_get(tminf->fwd_sam,0,F_0_t)) ;
  cur_Q[0]   = 1.0 ;
  dt1        = stp->delta_t ;
  dt2        = nxtstp->delta_t ;

/* generate X at the next step */
  gen_x(nxtstp,next_x) ;

/* calculate adjustment */
  theta[0] = calc_adj(tminf,trinf,dt1,dt2,cur_x,cur_Q,df_at_2dt) ;

/* calculate probabilities at the three branches and get Arrow-Debreu prices */
  calc_prob_and_Q(stp,cur_x,cur_Q,theta[0],next_x,next_Q) ;

/* swap cur_x,Q with next_x,Q */
  SWAP(cur_x,next_x) ;
  SWAP(cur_Q,next_Q) ;
  memset(&next_x[maxtrinf.min_x_index],0,sz1*sizeof(double)) ;
  memset(&next_Q[maxtrinf.min_x_index],0,sz1*sizeof(double)) ;

/* move to next step and loop */

  for (stp = nxtstp,stp_num=1 ; stp->next->next !=NULL ; 
	stp = stp->next,stp_num++)
  {
    trinf      = stp->trinf ;
    tminf      = stp->tminf[0] ;
    nxtstp     = stp->next ;
    nxttminf   = nxtstp->tminf[0] ;
    df_at_2dt  = ((SrtIRMTmInf *) nxtstp->next->tminf[0])->df ;
    dt1        = stp->delta_t ;
    dt2        = nxtstp->delta_t ;

/* generate X */
    gen_x(nxtstp,next_x) ;

/* calculate adjustment */
    theta[stp_num] = calc_adj(tminf,trinf,dt1,dt2,cur_x,cur_Q,df_at_2dt) ;

/* calculate probabilities at the three branches and get Arrow-Debreu prices */
    calc_prob_and_Q(stp,cur_x,cur_Q,theta[stp_num],next_x,next_Q) ;

/* swap cur_r,phi,Q with next_x,phi,Q */
    SWAP(cur_x,next_x) ;
    SWAP(cur_Q,next_Q) ;
    memset(&next_x[maxtrinf.min_x_index],0,sz1*sizeof(double)) ;
    memset(&next_Q[maxtrinf.min_x_index],0,sz1*sizeof(double)) ;
  }

/* Evalulate cashflows at the bottom of the tree */
  trinf = stp->trinf ;
  tminf = stp->tminf[0] ;
  for (i=trinf->min_x_index ; i<= trinf->max_x_index ; i++)
  {
    sam_get(node.cur_sam,0,STATEVAR)   = cur_x[i] ; 
    sam_get(node.cur_sam,0,SHORT_RATE) = exp(cur_x[i]) ;
    err = evalcf(
				(GrfnEvent *)stp->e,
				&node.cur_sam,
				gd,
				(double *) cur_assets[i],
				(EvalEventDfsFct) srt_f_calc_grfn_event_dfs,
				und_info,
				&cashflow) ;
	if (err) return err;
  }
  SWAP(cur_x,next_x) ;
  SWAP(cur_assets,next_assets) ;

/* go backwards for all steps.  Since we start at the last but one step
   from maturity, we have to decrement the stp_num by 1 at start */

  for (stp=stp->prev, stp_num-- ; stp!=NULL ; stp=stp->prev, stp_num--)
  {
    trinf = stp->trinf ;
    tminf = stp->tminf[0] ;
    dt1 = stp->delta_t ;
/* generate all X in the same time step */
    gen_x(stp,cur_x) ;
/* In BK model, variance of X is only a function of t and not X, hence it
   needs to be calculated only once per step */
    varatsam = tminf->ev.onef.sig2*dt1 ;
/* for each X in the same time step */
    for (i=trinf->min_x_index ; i<= trinf->max_x_index ; i++)
    {
/* discount asset prices */
      generate_prob(stp,&node,varatsam,theta[stp_num],cur_x[i],next_x) ;
      srt_f_trintredcntvec(stp,cur_assets[i],&node,next_assets,
			   grfnparam->node_dim) ;
      err = evalcf(
					(GrfnEvent *)stp->e,
					&node.cur_sam,
					gd,
					(double *) cur_assets[i], 
					(EvalEventDfsFct) srt_f_calc_grfn_event_dfs,
					und_info,
					&cashflow);
	  if (err) return err;
    }
    SWAP(cur_assets,next_assets) ;
    SWAP(cur_x,next_x) ;
  }

   price = next_assets[0][grfnparam->node_dim-1] ;
   io_request = iolist->head ;
   io_request_val_p = (SrtIOVal*)(io_request->element->val.pval);
   strncpy(io_request_val_p->computation_origin,
	"Tre_BDT_1F", strlen("Tre_BDT_1F"));
   io_request_val_p->dval = price ;


/** free the world **/

  free_dmatrix(cur_assets,maxtrinf.min_x_index,maxtrinf.max_x_index,
		 	  0,grfnparam->node_dim-1) ;
  free_dmatrix(next_assets,maxtrinf.min_x_index,maxtrinf.max_x_index,
		 	  0,grfnparam->node_dim-1) ;
  free_dvector(cur_x,maxtrinf.min_x_index,maxtrinf.max_x_index);
  free_dvector(next_x,maxtrinf.min_x_index,maxtrinf.max_x_index);
  free_dvector(cur_Q,maxtrinf.min_x_index,maxtrinf.max_x_index);
  free_dvector(next_Q,maxtrinf.min_x_index,maxtrinf.max_x_index);
  free_dvector(theta,0,maxtrinf.max_x_index) ;

  return NULL;
}
