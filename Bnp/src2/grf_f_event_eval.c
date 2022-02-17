/******************************************************************************\

    FUNCTION        :   grfn_eval_event (evalcf)

    PURPOSE         :   Evaluate the cashflows for a GRFN event.


    DESCRIPTION     :   [1] Compute the discount factors. 
                        This means calling evaldf with sample and fwd_sample.
                        (Only if rate functions are used in the event)

                        [2] Evaluate the COMLL_PTR stored inside event. 
                        We use "reverse polish" notation. A stack in the form
                        of a static array is used to keep the result.
 
                        Results of each cell in the tableau will be saved in 
                        gd->cells; assignments in the GRFN language will result
                        in GRFN variables (stored inside gd) being modified, 
                        and pv_vect will be incremented (for tree computation).

    RETURNS         :   An error 
						Affect in the cash-flow input the value of the entire event,
						which will be the cashflow in the rightmost of the tableau.

                        Result = f(sample, fwd_sample, gd, [tree information])
                        
    PRAMETERS       :   GrfnEvent       *event
	                    Sample          *sample      
                        GrfnDeal        *gd       
                        Grfnpv_vectVec  pv_vect
                        DfEval          evaldf

    NOTE            :   If event == NULL or COMLL_PTR == NULL,
								cashflow = 0, return NULL (no error).

                        If pv_vect ==  NULL, use an internal array to take 
                        the place of pv_vect (For Monte-Carlo pv_vect can be 
                        NULL). COMLL_PTR inside event will still contain 
                        nodes of type COMLL_T_INCREMENT, because it is unknown
                        when an event is being constructed which numerical 
                        method will be used.

                        In order to avoid numerical overflow when computing 
                        swap rates when the value of the appropriate level 
                        payment is too small, it will approximate the swap rate
                        with the short rate.
  
    IMPORTANT       :   No checking for stack overflow (unlikely???)

    SPECIAL CASE    :   For Monte-Carlo last argument CAN be NULL.
                        => cash flow is returned.

                        [Optional] For Monte-Carlo computations pv_vect CAN be
                        substituted by a blank double array of size = Number of
                        cols in the tableau. Subsequently the contents of the
                        array will be incremented by the values of the
                        appropriatre cells. Consequently we could value several
                        deals at the same time (i.e. one for each column).
                        
                        For Tree last routine MUST be the vector being 
                        discounted through the tree.


    AND FINALLY     :   We always evaluate a COMLL from RIGHT to LEFT.
                                                             
                        On the stack f(x,y) will leave x above y.         

                                                             Stack
                                                              ---
                        So the stack will look like:        
                                                             | x |
                        The stack is NOT used for financial  | y |  
                        functions since their arguments are    .
                        known a-priori.

\******************************************************************************/


/******************************************************************************\
*                           Import Include Files                              *
\******************************************************************************/


#include "GRF_H_ALL.H>
#include "srt_h_sample.h"
#include "swp_h_cms.h"
#include "swp_h_cmsopt.h"
#include "math.h"
#include "opfnctns.h"
#include "srt_h_product_list.h"
#include "RainbowOpt.h"


/******************************************************************************\
*                           Macro Constants                                    *
\******************************************************************************/


#define STACK_SIZE  101                  
#define zero        0.0
#define one         1.0
#define TINY        1.0e-09
#define MAXDFS		1000

static double	pv_copy[STACK_SIZE];	/* Used for a copy of pv_vect vector  */
static double   *stack;					/* Pointer which describes the stack */


/* ----------------------------------------------------------------------------
                           Macro Functions  to Operate on the Stack
   ---------------------------------------------------------------------------- */


/* ------------ Take something off the stack  --------------- */
static double POP()
{
double x;
	
	x = *stack;
	stack--;
	
	return x;
}

/* ------------ Put something on the stack  ------------------ */
static void PUSH (double X)
{
	++stack;
	*stack = X;        
}


/* ----------------------------------------------------------------------------
               Macro Functions  When using or setting the PV of a column
   ---------------------------------------------------------------------------- */

#define PVREF(c1)   pv_copy[c1]         /* Gets a previously computed PV */
#define PVSET(c1,v) {curr_pv[c1]=(v);}  /* Set a new PV  */


/* -----------------------------------------------------------------------------
                           MAIN FUNCTION
   Loops throught the COMLL attached to the GrfnEvent to evaluate it, 
   using the stack to dump and collect values for operations...
   ---------------------------------------------------------------------------- */


Err grfn_eval_event_ammc  
    (   int               type_eval,        /* 0: fwd, 1: bwd, 2: nostatus */
		GrfnEvent        *event ,        /* Grfn Event as attached to a step   */ 
		SrtSample        *sample,        /* State variables                    */
        GrfnDeal         *gd,            /* GRFN Workspace                     */
        GrfnLatticeVec    pv_vect,        /* Discounted Vector in tree          */
		double           *fwd,
		double           *cur,
        EvalEventDfsFct   evaldf,         /* Function for calculation of dfs    */
        SrtUndInfo       *und_info,      /* Underlying Information             */
		double			 *cashflow )		/* CashFlow calculated via the comll  */

{
double	s[STACK_SIZE]; /* Stack to stock the different real results */

/*..Local variables..*/
double          rubbish;
double          coupon, clean, yield;
SwapDP          sdp;                /* Swap data parameters */
long            ipath;              /* Monte-Carlo path number */
static long     npath = 0;
Err             error;
COMLL_PTR       comll;
int             ix = 0;             /* Index of underlying curve */
int             i,ind, c1,c2, r1,r2;
int             call_put;
double          x,y,f,k,vol,delay,betavol,alpha,beta, rho;
double          spot, mat, dBarrier[2], dPar[16];
double			fwdx, fwdy, sigx, sigy, a, b, t;
int             nstep;
SrtDiffusionType difftype;
long            l;
double* tmp_row;
GrfnLatticeVec  curr_pv;                 /* State vector required for evaluation
                                       of a GRFN deal using a pv_vect.
                                       curr_pv[i] is the discounted expected
                                       value of ith level of pv_vect
                                     */
/* Dfs for COMLL_F_PAYOFF: */
double dfs_d[MAXDFS], dfs_f[MAXDFS];
double *dfs[2] = {dfs_d, dfs_f};
SrtProduct *product;

    
 /* For a NULL event or a NULL COMLL return 0 */
    if (event == NULL)
	{
		*cashflow = 0.0;
		return NULL;
	}

/* If this is a tree, make a copy of the pv_vect vector 
	to keep track of PV's before current event evaluations */
    if (pv_vect != NULL)
	{
        curr_pv = pv_vect;
        memcpy( pv_copy, pv_vect, gd->sswidth*sizeof(double) );
    }
    else
	{
	/* In this case, all the PV's will be computed on a local array INCREMENTED at each step */
		curr_pv = pv_copy;
    }
    
/* Initialise the stack to 0 (new event) */
	memset(s,0, STACK_SIZE * sizeof(double));
	
/* Make the static stack point to this STACK (for PUSH and POP) */
	stack = s;
  

/* Gets the comll from the evnt */
    comll = event->icom;
	if (comll == NULL)
	{
		*cashflow = 0.0;
		return NULL;
	}
  
 /* Calculate the discount factors (given that we are in the future)
    as seen from where we are in the future to the dates listed
    in this event */
	if (evaldf != NULL)
        evaldf(event, sample, und_info);
        
/* Loop through the command link list ... */
    while (comll != NULL) 
	{

		/* test if evaluation or skipping in ammc lsm */
		if ( (type_eval == 2) || ((type_eval == 0) && (comll->evalstatus == FWDCELL))
							|| ((type_eval == 1) && (comll->evalstatus == BWDCELL)) )
		{

			/* Look at the TYPE of the node and take appropriate action */
			switch(comll->type)
			{

				case COMLL_ASSIGN:
					
					/* Store value of stack in gd->cells - Note that
					   comll->gdcells points to gd->cells[i][j] */
					*((double *)comll->gdcells) = *stack;
					comll = comll->next;
					break;

				case COMLL_BRANCH:
					
					if ( POP() > 0.5 )      /* Values > 0.5 mean boolean TRUE */
						comll = comll->next;
					else
						comll = comll->jmp->next;
					break;
				case COMLL_CONSTF:
					x= POP();
					y = POP();
					PUSH(gd->cells[DTOL(y)][DTOL(x)]);
					comll = comll->next;
					break;
				case COMLL_VAR:         
				case COMLL_USERVAR:     
				case COMLL_CONST:
					
					x = *((double *)comll->gdcells); 
					PUSH(x);
					comll = comll->next;
					break;
      
				case COMLL_JUMP:
					
					comll = comll->jmp->next;
					break;

				case COMLL_POP:

					comll->dval = *stack;
					POP();
					comll = comll->next;
					break;
      
				case COMLL_RESET_STATE:

					memcpy(&gd->ssparam,&gd->hist_ssparam,sizeof(GrfnSSParam));
					comll = comll->next;
					break;

			 /* ----------------------------  Financial events ------------------------------------- */
 				case COMLL_F_PAYFROMPAST :
				/* Gets the KNOWN payment amount from the stack (paid in the future) */
					x = POP();
				/* Computes the PV with the Domestic Curve (through evaldf on top of function) */
					rubbish = x * event->df[0][comll->dfind[1]] ;
				/* Increment the PV of the past cashflows if we are in the past */
					if (gd->is_history) 
						gd->pv_of_past += rubbish;
				/* Keep track of this payment: store cashflow, pay date and event date in the GrfnDeal */
					gd->numpay += 1;
					gd->pay_dates = realloc(gd->pay_dates, gd->numpay * sizeof(long) );
					gd->pay_dates[gd->numpay-1] = DTOL(comll->dval);
					gd->event_pay_dates = realloc(gd->event_pay_dates, gd->numpay * sizeof(long) );
					gd->event_pay_dates[gd->numpay-1] = event->d;;
					gd->pay_amounts = realloc(gd->pay_amounts, gd->numpay * sizeof(double) );
					gd->pay_amounts[gd->numpay-1] = x;
				/* Puts the discounted amount back into the stack (fro proper event treatment) */	
					PUSH(rubbish);
					comll = comll->next;

				break;
 				case COMLL_F_PAYINPAST :
				/* Gets the KNOWN payment amount from the stack (paid in the past) */
					x = POP();
				/* Keep track of this payment: store cashflow, pay date and event date in the GrfnDeal */
					gd->numpay += 1;
					gd->pay_dates = realloc(gd->pay_dates, gd->numpay * sizeof(long) );
					gd->pay_dates[gd->numpay-1] = DTOL(comll->dval);
					gd->event_pay_dates = realloc(gd->event_pay_dates, gd->numpay * sizeof(long) );
					gd->event_pay_dates[gd->numpay-1] = event->d;
					gd->pay_amounts = realloc(gd->pay_amounts, gd->numpay * sizeof(double) );
					gd->pay_amounts[gd->numpay-1] = x;
				/* Puts a zero value in the stack (for proper event treatment) */	
					PUSH(0.0);
					comll = comll->next;

				break;
				case COMLL_F_PAY :
				case COMLL_F_PAYN :
				/* Gets the payment amount from the stack */
					x = POP();
				/* Computes the PV with the Domestic Curve (through evaldf on top of function) */
					rubbish = x * event->df[0][comll->dfind[1]] ;
				/* Puts the discounted amount back into the stack */	
					PUSH(rubbish);
					comll = comll->next;
				break;

				case COMLL_F_PVRANGE:
				case COMLL_F_LEVEL:     
            
					/* Get index of underlying */
					ix = comll->ivec[0];
                
					x = zero;

					for ( i=1; i < comll->dfindlen; i++ )
						x += event->df[ix][comll->dfind[i]] * comll->cvg[i-1];
                
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_F_PAYOFF:

				/* Get product and date from comll */
					memcpy(&product, comll->ptr, sizeof(SrtProduct*));
					l = (long)(comll->dval);

				/* Get domestic dfs */
					for (i=0; i < comll->dfindlen; i++)
						dfs[0][i] = event->df[0][comll->dfind[i]];

				/* Get foreign dfs if available */
					if ( comll->next && comll->next->type == COMLL_F_PAYOFF &&
						 comll->next->ptr == NULL )
					{
						comll = comll->next;
						ix = comll->ivec[0];
						for (i=0; i < comll->dfindlen; i++)
							dfs[1][i] = event->df[ix][comll->dfind[i]];
					}

				/* Evaluate payoff */

					error = product->Payoff(product, event->d, l, dfs, &x);
					if (error) return error;

					PUSH(x);
					comll = comll->next;
					break;

			    case COMLL_F_CMS:

				/* Get index of underlying */
					ix = comll->ivec[0];
					x = zero;

				/* Compute the level payment (x) */
					for ( i=1; i < comll->dfindlen; i++ ) 
					{
						x += event->df[ix][comll->dfind[i]] * comll->cvg[i-1];
					}

				/* Compute the floating leg value with all spreads (y): ini + fin + spreads */
					y = event->df[ix][comll->dfind[0]] - 
							event->df[ix][comll->dfind[comll->dfindlen-1]];
					for (i = 1; i < comll->dfindfloatlen; i++)
					{
						y += event->df[ix][comll->dfindfloat[i]] 
							* comll->spread[i-1] * comll->cvgfloat[i-1];
					}

				/* If rates are very high, x will be very small. In this case
					   replace the swap-rate with the short_rate, to avoid 
					   numerical overflow */

					if ( x < TINY )
					{
                		x = samptr_get(sample, ix, SHORT_RATE);
					}
					else
					{
					/* The value of the swap rate is floating divided by fixed leg  */
						x = y / x;
					}
				
				/* Adds the long spread to this rate */
					x += comll->long_spread;

				/* Get other arguments */
				/* they are store in sstore and dstore arrays */
				/* vol,delay, lognormal/normal in dstore */
					vol = comll->dstore[0];
					
					delay = comll->dstore[1];

					difftype = (SrtDiffusionType)(comll->dstore[2]);

				/* Compute the cms rate */
/*					error = swp_f_cmsratedp(
							(SwapDP *)comll->ptr, 
							event->d,
							x,
							vol,
							DTOL(delay), 
							DEFAULT_CMS_DELTA,
							MAX_CMS_SWAPS,
							difftype, 
							&x);	*/

					error = swp_f_cmsratedp_clsdfrm(
							(SwapDP *)comll->ptr, 
							event->d,
							x,
							vol,
							DTOL(delay), 
							difftype, 
							&x);
					
					if (error)
						return error;

				/* Puts the value into the stack and move on */
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_F_SWAP:
				
				/* Get index of underlying */
					ix = comll->ivec[0];
					x = zero;

				/* Compute the level payment (x) */
					for ( i=1; i < comll->dfindlen; i++ ) 
					{
						x += event->df[ix][comll->dfind[i]] * comll->cvg[i-1];
					}

				/* Compute the floating leg value with all spreads (y): ini + fin + spreads */
					y = event->df[ix][comll->dfind[0]] - 
							event->df[ix][comll->dfind[comll->dfindlen-1]];
					for (i = 1; i < comll->dfindfloatlen; i++)
					{
						y += event->df[ix][comll->dfindfloat[i]] 
							* comll->spread[i-1] * comll->cvgfloat[i-1];
					}

				/* If rates are very high, x will be very small. In this case
					   replace the swap-rate with the short_rate, to avoid 
					   numerical overflow */

					if ( x < TINY )
					{
                		x = samptr_get(sample, ix, SHORT_RATE);
					}
					else
					{
					/* The value of the swap rate is floating divided by fixed leg  */
						x = y / x;
					}
				
				/* Adds the long spread to this rate */
					x += comll->long_spread;

				/* Puts the value into the stack and move on */
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_F_DF:

					/* Get index of underlying */
					ix = comll->ivec[0];

					x = event->df[ix][comll->dfind[0]];
					if ( x < TINY )
					{
						x = zero;
					}
					else
					{
						x = event->df[ix][comll->dfind[1]] / x;
					}
                
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_F_FRA:

					/* Get index of underlying */
					ix = comll->ivec[0];

					x = event->df[ix][comll->dfind[1]];

					if ( x < TINY )
					{
						/*
						x = sample->und[ix].sv[0];
						*/
						x = samptr_get(sample, ix, SHORT_RATE);
					}
					else
					{
						x = event->df[ix][comll->dfind[0]]/x - one;
						x /= comll->cvg[0];
					}

					/* For FRA's, there is only one spread (one stored per period)
					   if it is not cash (otherwise comll->spread = NULL)
					*/
					if (comll->spread)
						x += comll->spread[0];
					
					PUSH(x);
					comll = comll->next;
					break;
    
				case COMLL_F_SWAPTION:

					/* Get index of underlying */

					ix = comll->ivec[0];

					x = zero;

					/* Get the level */
					for ( i=1; i < comll->dfindlen; i++ )
					{
						x += event->df[ix][comll->dfind[i]]*comll->cvg[i-1];
					}

					/* Get the swap rate */
					y = (event->df[ix][comll->dfind[0]] - 
						event->df[ix][comll->dfind[comll->dfindlen-1]])/x;
                
					call_put = comll->ivec[1];
                
					k = comll->dval;
                
					x *= DMAX( (double)call_put*(y-k) , 0 );
                
					PUSH(x);
					comll = comll->next;
					break;	

			/* ---------------------- Arithmetic ------------------------------ */
				case COMLL_A_PLUS:
                
					x  = POP();
					x += POP();
					PUSH(x);
					comll = comll->next;
					break;
      
				case COMLL_A_TIMES:

					x  = POP();
					x *= POP();
					PUSH(x);
					comll = comll->next;
					break;
      
				case COMLL_A_MINUS:
        
					x  = POP();
					x -= POP();
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_A_DIVIDE:
        
					x  = POP();
					x /= POP();
					PUSH(x);
					comll = comll->next;
					break;
      
				case COMLL_A_LOG:
        
					x = POP();
					PUSH(log(x));
					comll = comll->next;
					break;
      
				case COMLL_A_EXP:
                
					x = POP();
					x = exp(x);
					PUSH(x);
					comll = comll->next;
					break;
      
				case COMLL_A_FLOOR:
        
					x = POP();
					y = POP();
					if (y == 0)
					{
						x = 0;
					}
					else
					{
						x = x/y;
                		l = DTOL(x);
                		x = (double)l*y;
					}
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_A_CEILING:

					x = POP();
					y = POP();
					if (y == 0)
					{
						x = 0;
					}
					else
					{
						x = -x/y;
    					l = DTOL(x);
        				x = (double)l* -y;
					}
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_A_POW:

					x = POP();
					y = POP(); 
					x = pow(x,y);
					PUSH(x);
					comll = comll->next;
					break;
                   
				case COMLL_A_MAX:

					x = POP();
					y = POP();
					x = DMAX(x,y);
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_A_MIN:

					x = POP();
					y = POP();
					x = DMIN(x,y);
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_A_ABS:

					x = POP();
					x = fabs(x);
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_A_NORM:
        
					x = POP();
					x = norm(x);
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_A_BLKSCHLS:          /* Black-Scholes */

					f  = POP();                     /* FWD  */
					k  = POP();                     /* k    */
					x  = POP();                     /* vol  */
					y  = POP();                     /* time difference = mat */
					x *= sqrt(y);               	/* calc  */
					y  = (log(f/k) + 0.5*x*x) / x;  /* d1   */
					x  = y - x;                     /* d2   */

					if (POP() < 0.5)                /* call */
						x = f*norm(y)-k*norm(x);
					else                            /* put  */
						x = -f*norm(-y)+k*norm(-x);
                
					PUSH(x);
					comll = comll->next;
					break;  
				case COMLL_A_BSSABR:          /* Black-Scholes with SABR */

					f  = POP();                     /* FWD  */
					k  = POP();                     /* k    */
					betavol  = POP();               /* vol  */
					alpha = POP();                     
					beta = POP();
					rho = POP();
					y = POP();						/* time difference = mat */
					difftype = (SrtDiffusionType) POP();		/* 0 = Log, 1 = Norm, 2 = beta */
					call_put = (int) POP();

/*					if (error = srt_f_optbetastochvoltoblkvol(
						f,							
						k, 
						betavol,
						alpha,
						rho,
						y,  
						beta,
						&x)) return error;
*/

					if (  error = srt_f_optsarbvol( f, k, y, betavol, alpha, beta, rho, difftype, SRT_LOGNORMAL, &x )  )
						return error;
					if ( x <= 0.0 )
					{
						if ( call_put )					/* put */
							x = k > f ? k - f : 0.0;
						else
							x = f > k ? f - k : 0.0;
					}
					else
					{
						x *= sqrt(y);               	/* calc  */
						
						y  = (log(f/k) + 0.5*x*x) / x;  /* d1   */
						x  = y - x;                     /* d2   */

						if (call_put)                /* put */
							x = -f*norm(-y)+k*norm(-x);
						else                            /* call */
							x = f*norm(y)-k*norm(x);
					}


					PUSH(x);
					comll = comll->next;
					break;  


				case COMLL_A_BLKSCHLSNRM:          /* Black-Scholes Normal */

					f  = POP();                     /* FWD  */
					k  = POP();                     /* k    */
					x  = POP();                     /* vol  */
					y  = POP();                     /* time difference = mat */
					x *= sqrt(y);               	/* calc */
					y  = (f - k) / x;				/* d1   */
                
					if (POP() < 0.5)                /* call */
						x = x * gauss(y) + (f - k) * norm(y);
					else                            /* put  */
						x = x * gauss(y) - (f - k) * norm(-y);
                
					PUSH(x);
					comll = comll->next;
					break;  

				case COMLL_A_SPREADNUM:			/* Spread_numer */
                    
                    fwdx = POP();                     /* FWD X  */
					fwdy = POP();                     /* FWD Y  */
					sigx = POP();                     /* SigX   */
					sigy = POP();                     /* SigY   */
					rho  = POP();                     /* Rho   */
					a	 = POP();                     /* a   */
					b	 = POP();                     /* b   */
					k    = POP();                     /* K   */
					t	 = POP();                     /* t   */
					call_put = (int) POP();                 /* C/P  */
                    
					error = OptSpreadNew(
									fwdx,
									a,
									sigx,
									fwdy,
									-b,
									sigy,
									k,
									t, 
									rho,
									call_put,
									&x);

					if (error)
						return error;
					
					PUSH(x);
					comll = comll->next;
					break;
                

				case COMLL_A_NUMDAYS:			/* number of days in a corridor */
					
					f    = POP();
					spot = POP();
					dBarrier[0]  = POP();
					dBarrier[1] = POP();
					vol  = POP();
					delay = POP();
					mat  = POP();
					nstep= (int) POP();
					dPar[0] = POP();
					dPar[1] = POP();
					i = (int) POP();

					PUSH( prob_ExpectedDays(f, spot, dBarrier, vol, delay, mat, nstep, dPar, i) );
					comll = comll->next;
					break;  

				case COMLL_F_YIELD:             /* Bond yield */

					coupon = POP();
					clean  = POP();
					sdp = *((SwapDP *)comll->ptr);
 
					x = yield_fct(sdp, coupon, clean, 1.0, coupon); 

					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_F_CLEAN:             /* Bond clean price */

					coupon = POP();
					yield  = POP();
					sdp = *((SwapDP*)comll->ptr);
 
					x = clean_price_fct(sdp, coupon, yield, 1.0, coupon); 

					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_B_LT:        /* Less than */

					x = POP();
					y = POP();
					x = (x < y ? one  : zero);
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_B_GT:         /* Greater than */  

					x = POP();
					y = POP();
					x = (x > y ? one  : zero);
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_B_LE:        /* Less or equal */  

					x = POP();
					y = POP();
					if (fabs( y - x) <= DBL_EPSILON)
					{
						 x = one;
					}
					else
					{
						x = ( x <= y  ? one : zero );
					}
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_B_GE:           /* Greater or equal */  

					x = POP();
					y = POP();
					if (fabs(x - y) <= DBL_EPSILON)
					{
						 x = one;
					}
					else
					{
						x = ( x >= y  ? one : zero );
					}
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_B_EQ:            /* Equal */  

					x = POP();
					y = POP();
					if (fabs(x) <= DBL_EPSILON)
					{
						 x = ( ( fabs(y) <= DBL_EPSILON ) ? one  : zero );
					}
					else
					{
						x = ( ( fabs(y/x-1) <= DBL_EPSILON ) ? one : zero );
					}
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_B_AND:

					x = POP();
					y = POP();
					x = (y > .5 && x > .5 ? one  : zero);
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_B_OR:

					x = POP();
					y = POP();
					x = (y > .5 || x > .5 ? one  : zero);
					PUSH(x);
					comll = comll->next;
					break;

			/* ----------------------- Cell Functions ---------------------------------- */
				case COMLL_C_COLMAX:
				case COMLL_C_COLMIN:
				case COMLL_C_COLSUM:
				case COMLL_C_COLAVG:
				case COMLL_C_ROWMAX:
				case COMLL_C_ROWMIN:
				case COMLL_C_ROWSUM:
				case COMLL_C_ROWAVG:
				case COMLL_C_ROWMAXIDX: /* Cyril Godart Regis Benichou fecit per Seculum Secolorum */
				case COMLL_C_ROWMINIDX:

					c1 = comll->ivec[0];
					r1 = comll->ivec[1];
					r2 = comll->ivec[2];
					ind = 0;					/* used only for rowsort */
					grf_f_evlcllfct ( gd->cells, c1, r1, r2, ind, comll->type, &x); 
					PUSH(x);
					comll = comll->next;
					break;
        
				case COMLL_C_COLSORT:
				case COMLL_C_ROWSORT:

					c1 = comll->ivec[0];
					r1 = comll->ivec[1];
					r2 = comll->ivec[2];
					ind = comll->ivec[3];
					grf_f_evlcllfct ( gd->cells, c1, r1, r2, ind, comll->type, &x); 
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_C_ROWINTERP:
					r1  = comll->ivec[0];
					c1  = comll->ivec[1];
					c2  = comll->ivec[2];
					ind = comll->ivec[3];
				/* Get the value to interpolate from the stack */
					x = POP();
				/* Interpolate y=f(x) using the row as Y values and the auxiliary range as X values */
					grf_f_evlcllfct_interp(gd->aux[ind], gd->cells, r1, c1, c2, x, comll->type, &y); 
					PUSH(y);
					comll = comll->next;
					break;

			/* ------------------ State Variables ----------------------------- */
				case COMLL_S_SHORT_RATE:

					/* Get index of underlying */
					ix = comll->ivec[0];

					/*
					PUSH((sample->und[ix].sv[SHORT_RATE]));
					*/
					PUSH( samptr_get(sample, ix, SHORT_RATE) );
					comll = comll->next;

					break;

				case COMLL_S_PHI:

					/* Get index of underlying */
					ix = comll->ivec[0];
                
					/*
					PUSH((sample->und[ix].sv[PHI]));
					*/
					PUSH( samptr_get(sample, ix, PHI) );
					comll = comll->next;
					break;

				case COMLL_S_SIGMA:

					/* Get index of underlying */
					ix = comll->ivec[0];
                
					/*
					PUSH((sample->p[ix].s[SIGMA]));
					*/
					PUSH( samptr_get(sample, ix, SIGMA) );
					comll = comll->next;
					break;

				case COMLL_S_STATEVAR:

					/* Get index of underlying */
					ix = comll->ivec[0];
                
					/*
					PUSH((sample->p[ix].s[STATEVAR]));
					*/
					PUSH( samptr_get(sample, ix, STATEVAR) );
					comll = comll->next;
					break;

				case COMLL_S_NUMERAIRE:

					PUSH(sample->numeraire);
              
					comll = comll->next;
					break;

				case COMLL_S_MMACCOUNT:

					/* Get index of underlying */
					ix = comll->ivec[0];
                
					PUSH( samptr_get(sample, ix, BT) );
					comll = comll->next;
					break;

				case COMLL_S_STOCK:

					/* index of underlying = 0 */
					ix = 0;
                
					PUSH( samptr_get(sample, ix, SPOT) );
					comll = comll->next;
					break;

				case COMLL_S_UNDERLYING:

					/* Get index of underlying */
					ix = comll->ivec[0];
                
					PUSH( samptr_get(sample, ix, SPOT) );
					comll = comll->next;
					break;
			/* --------------------- Use or Increment of Column PV ------------------------------- */
				case COMLL_T_CUMREF:

					c1 = comll->ivec[0];
					x = fwd[c1];
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_T_CURREF:

					c1 = comll->ivec[0];
					x = cur[c1];
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_T_REF:

					c1 = comll->ivec[0];
					x = PVREF(c1);
					PUSH(x);
					comll = comll->next;
					break;

				// IMPLEMENTATION OF PVR IN GRFN
				case COMLL_T_RISKY_REF:
					return serror("No mention of PVR is allowed with this GRFN function !");
					break;

				case COMLL_T_SET:

					c1 = comll->ivec[0];
					PVSET(c1, *stack);
					comll = comll->next;
					break;
            
				case COMLL_T_INCREMENT:

					c1 = comll->ivec[0];
					PVSET(c1,PVREF(c1) + *stack);
					comll = comll->next;
					break;
				
				case COMLL_T_PVSOURCE_ASSIGN:
					
					return serror("No mention of SOURCE is allowed with this GRFN function !");
					break;

			/* -----------------------  Interpolation ---------------------------- */
				case COMLL_X_INTERP:

					r1 = comll->ivec[0];
					r2 = comll->ivec[1];
					x = POP();
					y = x;
					interp(gd->aux[r1],gd->aux[r2],
                		IMIN(gd->auxlen[r1],gd->auxlen[r2]),y,0,&x);
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_X_PVINTERP:
					ind = comll->ivec[0]; /* aux index */
//					r1  = comll->ivec[1]; /* row */
					c1  = comll->ivec[1]; /* first col */
					c2 = comll->ivec[2]; /* col skip */
				/* Get the value to interpolate from the stack */
					x = POP();
				/* Interpolate y=f(x) using the row as Y values and the auxiliary range as X values */
					
		/* Make a temporary copy of the current cell values of row # c  */
					ix = (int) gd->auxlen[ind];
					tmp_row = dvector(0,ix-1);
					for (i=0;i<ix;i++)
//				tmp_row[i] = CELL(c1+i*c2,r);
						tmp_row[i] = PVREF(c1+i*c2);
			
		/* Sort the values of this row (for cells from r1 to r2 ) */
		/*	qsort(tmp_row,c2-c1+1,sizeof(double),qcmp); */
					interp(gd->aux[ind], tmp_row, ix, x, 0, &y);
			
		/* Interpolate y=f(x) over the values of this row using the x_values input */
//					grf_f_evlcllfct_interp(, gd->cells, r1, c1, c2, x, comll->type, &y); 
					PUSH(y);
					comll = comll->next;
					break;


/*					c1 = comll->ivec[0];
					x = PVREF(c1);
					PUSH(x);
					comll = comll->next;
					break;

*/
				case COMLL_PRINT:

					if ( gd->outfileptr != NULL )
					{
						fprintf ( gd->outfileptr,
								  "%.6lf#%.6lf#%.6lf#\n",
								  samptr_get(sample, ix, 0),
								  samptr_get(sample, ix, 1),
								  *stack );
					}

					comll = comll->next;
					break;

			/* --------------------- Histogram ------------------------------------ */
				case COMLL_HIST:

					/* Get Monte-Carlo path number */
					ipath = sample->pathnum;

					/* Add data to histogram: This function creates the required
					   data structure (if necessary) and subsequently stores the
					   given value, x, in that structure (assuming POP has stored it
					   Note: comll->sval = Histogram label 
					*/

					if ((comll->prev != NULL) && (comll->prev->type == COMLL_POP)  )
						x = comll->prev->dval;
					else 
						return serror("the right syntax for hist is expr1|hist(\"name\")");

					error = srt_f_sendhistval(comll->sval, ipath, x);
					if (error)
						return error; 

					PUSH(x);
					comll = comll->next;
					break;
			/* -------------------- A Real Value : Push it in the Stack ----------- */
				case COMLL_REAL:
				default:

					PUSH(comll->dval);
					comll = comll->next;
					break;

			} /* END switch */
		}
		else 
		/* the cell cannot be evaluated in the current computing mode forward/backward */
		{
			comll = comll->nextcol;
		}

    } /* END while */

/* The Event Cash Flow is the last element of the Stack */    
 	if (stack != s)  *cashflow = POP();

/* Checks here that the stack is empty (otherwise something missing) */
	if (*stack != 0.0)
		smessage("Stack not empty in GrfnEvalEvent: %f",(double)*stack);

    return NULL;
}



/* -----------------------------------------------------------------------------
                           MODIFIED EVALUATION FOR AMMC
   ---------------------------------------------------------------------------- */


Err grfn_eval_event 
    (   GrfnEvent        *event ,        /* Grfn Event as attached to a step   */ 
		SrtSample        *sample,        /* State variables                    */
        GrfnDeal         *gd,            /* GRFN Workspace                     */
        GrfnLatticeVec    pv_vect,        /* Discounted Vector in tree          */
        EvalEventDfsFct   evaldf,         /* Function for calculation of dfs    */
        SrtUndInfo       *und_info,      /* Underlying Information             */
		double			 *cashflow )		/* CashFlow calculated via the comll  */

{

	/* standard call status nostatus */
	return(grfn_eval_event_ammc(2, event, sample, gd, pv_vect, 
					(double *) NULL, (double *) NULL, evaldf,und_info, cashflow));

}


/* ----------------------------------------------------------------------------
               Macro Functions  When using or setting the PV of a column (Credit)
   ---------------------------------------------------------------------------- */

#define PVRISKYREF(c1)			pvR_vect[c1]				/* Gets a previously computed PV */
#define PVSOURCESET(c1,v)		{pvSource_vect[c1]=(v);}	/* Set a new PVSource */
#define PVSOURCESTATUSSET(c1,v) {pvSource_status[c1]=(v);}	/* Set a new PVSource status*/


/* -----------------------------------------------------------------------------
                           MAIN FUNCTION FOR CREDIT DEALS
   ---------------------------------------------------------------------------- */

Err grfn_eval_event_credit  
    (   int               type_eval,        /* 0: fwd, 1: bwd, 2: nostatus */
		GrfnEvent        *event ,        /* Grfn Event as attached to a step   */ 
		SrtSample        *sample,        /* State variables                    */
        GrfnDeal         *gd,            /* GRFN Workspace                     */
        GrfnLatticeVec    pvNR_vect,        /* Discounted Vector in tree          */
		GrfnLatticeVec    pvR_vect,        /* Discounted Risky Vector in tree          */
		GrfnLatticeVec    pvSource_vect,        /* Discounted Source Vector in tree          */
		int				 *pvSource_status,        /* Discounted Vector in tree          */
		double           *fwd,
		double           *cur,
        EvalEventDfsFct   evaldf,         /* Function for calculation of dfs    */
        SrtUndInfo       *und_info,      /* Underlying Information             */
		double			 *cashflow )		/* CashFlow calculated via the comll  */


{
double	s[STACK_SIZE]; /* Stack to stock the different real results */

/*..Local variables..*/
double          rubbish;
double          coupon, clean, yield;
SwapDP          sdp;                /* Swap data parameters */
long            ipath;              /* Monte-Carlo path number */
static long     npath = 0;
Err             error;
COMLL_PTR       comll;
int             ix = 0;             /* Index of underlying curve */
int             i,ind, c1,c2, r1,r2;
int             call_put;
double          x,y,f,k,vol,delay,betavol,alpha,beta, rho;
double          spot, mat, dBarrier[2], dPar[16];
double			fwdx, fwdy, sigx, sigy, a, b, t;
int             nstep;
SrtDiffusionType difftype;
long            l;
GrfnLatticeVec  curr_pv;                 /* State vector required for evaluation
                                       of a GRFN deal using a pv_vect.
                                       curr_pv[i] is the discounted expected
                                       value of ith level of pv_vect
                                     */
/* Dfs for COMLL_F_PAYOFF: */
double dfs_d[MAXDFS], dfs_f[MAXDFS];
double *dfs[2] = {dfs_d, dfs_f};
SrtProduct *product;

    
 /* For a NULL event or a NULL COMLL return 0 */
    if (event == NULL)
	{
		*cashflow = 0.0;
		return NULL;
	}

/* If this is a tree, make a copy of the pv_vect vector 
	to keep track of PV's before current event evaluations */
    if (pvNR_vect != NULL)
	{
        curr_pv = pvNR_vect;
        memcpy( pv_copy, pvNR_vect, gd->sswidth*sizeof(double) );
    }
    else
	{
	/* In this case, all the PV's will be computed on a local array INCREMENTED at each step */
		curr_pv = pv_copy;
    }
    
/* Initialise the stack to 0 (new event) */
	memset(s,0, STACK_SIZE * sizeof(double));
	
/* Make the static stack point to this STACK (for PUSH and POP) */
	stack = s;
  

/* Gets the comll from the evnt */
    comll = event->icom;
	if (comll == NULL)
	{
		*cashflow = 0.0;
		return NULL;
	}
  
 /* Calculate the discount factors (given that we are in the future)
    as seen from where we are in the future to the dates listed
    in this event */
	if (evaldf != NULL)
        evaldf(event, sample, und_info);
        
/* Loop through the command link list ... */
    while (comll != NULL) 
	{

		/* test if evaluation or skipping in ammc lsm */
		if ( (type_eval == 2) || ((type_eval == 0) && (comll->evalstatus == FWDCELL))
							|| ((type_eval == 1) && (comll->evalstatus == BWDCELL)) )
		{

			/* Look at the TYPE of the node and take appropriate action */
			switch(comll->type)
			{

				case COMLL_ASSIGN:
					
					/* Store value of stack in gd->cells - Note that
					   comll->gdcells points to gd->cells[i][j] */
					*((double *)comll->gdcells) = *stack;
					comll = comll->next;
					break;

				case COMLL_BRANCH:
					
					if ( POP() > 0.5 )      /* Values > 0.5 mean boolean TRUE */
						comll = comll->next;
					else
						comll = comll->jmp->next;
					break;
				case COMLL_CONSTF:
					x= POP();
					y = POP();
					PUSH(gd->cells[DTOL(y)][DTOL(x)]);
					comll = comll->next;
					break;
				case COMLL_VAR:         
				case COMLL_USERVAR:     
				case COMLL_CONST:
					
					x = *((double *)comll->gdcells); 
					PUSH(x);
					comll = comll->next;
					break;
      
				case COMLL_JUMP:
					
					comll = comll->jmp->next;
					break;

				case COMLL_POP:

					comll->dval = *stack;
					POP();
					comll = comll->next;
					break;
      
				case COMLL_RESET_STATE:

					memcpy(&gd->ssparam,&gd->hist_ssparam,sizeof(GrfnSSParam));
					comll = comll->next;
					break;

			 /* ----------------------------  Financial events ------------------------------------- */
 				case COMLL_F_PAYFROMPAST :
				/* Gets the KNOWN payment amount from the stack (paid in the future) */
					x = POP();
				/* Computes the PV with the Domestic Curve (through evaldf on top of function) */
					rubbish = x * event->df[0][comll->dfind[1]] ;
				/* Increment the PV of the past cashflows if we are in the past */
					if (gd->is_history) 
						gd->pv_of_past += rubbish;
				/* Keep track of this payment: store cashflow, pay date and event date in the GrfnDeal */
					gd->numpay += 1;
					gd->pay_dates = realloc(gd->pay_dates, gd->numpay * sizeof(long) );
					gd->pay_dates[gd->numpay-1] = DTOL(comll->dval);
					gd->event_pay_dates = realloc(gd->event_pay_dates, gd->numpay * sizeof(long) );
					gd->event_pay_dates[gd->numpay-1] = event->d;;
					gd->pay_amounts = realloc(gd->pay_amounts, gd->numpay * sizeof(double) );
					gd->pay_amounts[gd->numpay-1] = x;
				/* Puts the discounted amount back into the stack (fro proper event treatment) */	
					PUSH(rubbish);
					comll = comll->next;

				break;
 				case COMLL_F_PAYINPAST :
				/* Gets the KNOWN payment amount from the stack (paid in the past) */
					x = POP();
				/* Keep track of this payment: store cashflow, pay date and event date in the GrfnDeal */
					gd->numpay += 1;
					gd->pay_dates = realloc(gd->pay_dates, gd->numpay * sizeof(long) );
					gd->pay_dates[gd->numpay-1] = DTOL(comll->dval);
					gd->event_pay_dates = realloc(gd->event_pay_dates, gd->numpay * sizeof(long) );
					gd->event_pay_dates[gd->numpay-1] = event->d;
					gd->pay_amounts = realloc(gd->pay_amounts, gd->numpay * sizeof(double) );
					gd->pay_amounts[gd->numpay-1] = x;
				/* Puts a zero value in the stack (for proper event treatment) */	
					PUSH(0.0);
					comll = comll->next;

				break;
				case COMLL_F_PAY :
				case COMLL_F_PAYN :
				/* Gets the payment amount from the stack */
					x = POP();
				/* Computes the PV with the Domestic Curve (through evaldf on top of function) */
					rubbish = x * event->df[0][comll->dfind[1]] ;
				/* Puts the discounted amount back into the stack */	
					PUSH(rubbish);
					comll = comll->next;
				break;

				case COMLL_F_PVRANGE:
				case COMLL_F_LEVEL:     
            
					/* Get index of underlying */
					ix = comll->ivec[0];
                
					x = zero;

					for ( i=1; i < comll->dfindlen; i++ )
						x += event->df[ix][comll->dfind[i]] * comll->cvg[i-1];
                
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_F_PAYOFF:

				/* Get product and date from comll */
					memcpy(&product, comll->ptr, sizeof(SrtProduct*));
					l = (long)(comll->dval);

				/* Get domestic dfs */
					for (i=0; i < comll->dfindlen; i++)
						dfs[0][i] = event->df[0][comll->dfind[i]];

				/* Get foreign dfs if available */
					if ( comll->next && comll->next->type == COMLL_F_PAYOFF &&
						 comll->next->ptr == NULL )
					{
						comll = comll->next;
						ix = comll->ivec[0];
						for (i=0; i < comll->dfindlen; i++)
							dfs[1][i] = event->df[ix][comll->dfind[i]];
					}

				/* Evaluate payoff */

					error = product->Payoff(product, event->d, l, dfs, &x);
					if (error) return error;

					PUSH(x);
					comll = comll->next;
					break;

			    case COMLL_F_CMS:

				/* Get index of underlying */
					ix = comll->ivec[0];
					x = zero;

				/* Compute the level payment (x) */
					for ( i=1; i < comll->dfindlen; i++ ) 
					{
						x += event->df[ix][comll->dfind[i]] * comll->cvg[i-1];
					}

				/* Compute the floating leg value with all spreads (y): ini + fin + spreads */
					y = event->df[ix][comll->dfind[0]] - 
							event->df[ix][comll->dfind[comll->dfindlen-1]];
					for (i = 1; i < comll->dfindfloatlen; i++)
					{
						y += event->df[ix][comll->dfindfloat[i]] 
							* comll->spread[i-1] * comll->cvgfloat[i-1];
					}

				/* If rates are very high, x will be very small. In this case
					   replace the swap-rate with the short_rate, to avoid 
					   numerical overflow */

					if ( x < TINY )
					{
                		x = samptr_get(sample, ix, SHORT_RATE);
					}
					else
					{
					/* The value of the swap rate is floating divided by fixed leg  */
						x = y / x;
					}
				
				/* Adds the long spread to this rate */
					x += comll->long_spread;

				/* Get other arguments */
				/* they are store in sstore and dstore arrays */
				/* vol,delay, lognormal/normal in dstore */
					vol = comll->dstore[0];
					
					delay = comll->dstore[1];

					difftype = (SrtDiffusionType)(comll->dstore[2]);

				/* Compute the cms rate */
/*					error = swp_f_cmsratedp(
							(SwapDP *)comll->ptr, 
							event->d,
							x,
							vol,
							DTOL(delay), 
							DEFAULT_CMS_DELTA,
							MAX_CMS_SWAPS,
							difftype, 
							&x);	*/

					error = swp_f_cmsratedp_clsdfrm(
							(SwapDP *)comll->ptr, 
							event->d,
							x,
							vol,
							DTOL(delay), 
							difftype, 
							&x);
					
					if (error)
						return error;

				/* Puts the value into the stack and move on */
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_F_SWAP:
				
				/* Get index of underlying */
					ix = comll->ivec[0];
					x = zero;

				/* Compute the level payment (x) */
					for ( i=1; i < comll->dfindlen; i++ ) 
					{
						x += event->df[ix][comll->dfind[i]] * comll->cvg[i-1];
					}

				/* Compute the floating leg value with all spreads (y): ini + fin + spreads */
					y = event->df[ix][comll->dfind[0]] - 
							event->df[ix][comll->dfind[comll->dfindlen-1]];
					for (i = 1; i < comll->dfindfloatlen; i++)
					{
						y += event->df[ix][comll->dfindfloat[i]] 
							* comll->spread[i-1] * comll->cvgfloat[i-1];
					}

				/* If rates are very high, x will be very small. In this case
					   replace the swap-rate with the short_rate, to avoid 
					   numerical overflow */

					if ( x < TINY )
					{
                		x = samptr_get(sample, ix, SHORT_RATE);
					}
					else
					{
					/* The value of the swap rate is floating divided by fixed leg  */
						x = y / x;
					}
				
				/* Adds the long spread to this rate */
					x += comll->long_spread;

				/* Puts the value into the stack and move on */
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_F_DF:

					/* Get index of underlying */
					ix = comll->ivec[0];

					x = event->df[ix][comll->dfind[0]];
					if ( x < TINY )
					{
						x = zero;
					}
					else
					{
						x = event->df[ix][comll->dfind[1]] / x;
					}
                
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_F_FRA:

					/* Get index of underlying */
					ix = comll->ivec[0];

					x = event->df[ix][comll->dfind[1]];

					if ( x < TINY )
					{
						/*
						x = sample->und[ix].sv[0];
						*/
						x = samptr_get(sample, ix, SHORT_RATE);
					}
					else
					{
						x = event->df[ix][comll->dfind[0]]/x - one;
						x /= comll->cvg[0];
					}

					/* For FRA's, there is only one spread (one stored per period)
					   if it is not cash (otherwise comll->spread = NULL)
					*/
					if (comll->spread)
						x += comll->spread[0];
					
					PUSH(x);
					comll = comll->next;
					break;
    
				case COMLL_F_SWAPTION:

					/* Get index of underlying */

					ix = comll->ivec[0];

					x = zero;

					/* Get the level */
					for ( i=1; i < comll->dfindlen; i++ )
					{
						x += event->df[ix][comll->dfind[i]]*comll->cvg[i-1];
					}

					/* Get the swap rate */
					y = (event->df[ix][comll->dfind[0]] - 
						event->df[ix][comll->dfind[comll->dfindlen-1]])/x;
                
					call_put = comll->ivec[1];
                
					k = comll->dval;
                
					x *= DMAX( (double)call_put*(y-k) , 0 );
                
					PUSH(x);
					comll = comll->next;
					break;	

			/* ---------------------- Arithmetic ------------------------------ */
				case COMLL_A_PLUS:
                
					x  = POP();
					x += POP();
					PUSH(x);
					comll = comll->next;
					break;
      
				case COMLL_A_TIMES:

					x  = POP();
					x *= POP();
					PUSH(x);
					comll = comll->next;
					break;
      
				case COMLL_A_MINUS:
        
					x  = POP();
					x -= POP();
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_A_DIVIDE:
        
					x  = POP();
					x /= POP();
					PUSH(x);
					comll = comll->next;
					break;
      
				case COMLL_A_LOG:
        
					x = POP();
					PUSH(log(x));
					comll = comll->next;
					break;
      
				case COMLL_A_EXP:
                
					x = POP();
					x = exp(x);
					PUSH(x);
					comll = comll->next;
					break;
      
				case COMLL_A_FLOOR:
        
					x = POP();
					y = POP();
					if (y == 0)
					{
						x = 0;
					}
					else
					{
						x = x/y;
                		l = DTOL(x);
                		x = (double)l*y;
					}
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_A_CEILING:

					x = POP();
					y = POP();
					if (y == 0)
					{
						x = 0;
					}
					else
					{
						x = -x/y;
    					l = DTOL(x);
        				x = (double)l* -y;
					}
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_A_POW:

					x = POP();
					y = POP(); 
					x = pow(x,y);
					PUSH(x);
					comll = comll->next;
					break;
                   
				case COMLL_A_MAX:

					x = POP();
					y = POP();
					x = DMAX(x,y);
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_A_MIN:

					x = POP();
					y = POP();
					x = DMIN(x,y);
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_A_ABS:

					x = POP();
					x = fabs(x);
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_A_NORM:
        
					x = POP();
					x = norm(x);
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_A_BLKSCHLS:          /* Black-Scholes */

					f  = POP();                     /* FWD  */
					k  = POP();                     /* k    */
					x  = POP();                     /* vol  */
					y  = POP();                     /* time difference = mat */
					x *= sqrt(y);               	/* calc  */
					y  = (log(f/k) + 0.5*x*x) / x;  /* d1   */
					x  = y - x;                     /* d2   */

					if (POP() < 0.5)                /* call */
						x = f*norm(y)-k*norm(x);
					else                            /* put  */
						x = -f*norm(-y)+k*norm(-x);
                
					PUSH(x);
					comll = comll->next;
					break;  
				case COMLL_A_BSSABR:          /* Black-Scholes with SABR */

					f  = POP();                     /* FWD  */
					k  = POP();                     /* k    */
					betavol  = POP();               /* vol  */
					alpha = POP();                     
					beta = POP();
					rho = POP();
					y = POP();						/* time difference = mat */
					difftype = (SrtDiffusionType) POP();		/* 0 = Log, 1 = Norm, 2 = beta */

/*					if (error = srt_f_optbetastochvoltoblkvol(
						f,							
						k, 
						betavol,
						alpha,
						rho,
						y,  
						beta,
						&x)) return error;
*/

					if (  error = srt_f_optsarbvol( f, k, y, betavol, alpha, beta, rho, difftype, SRT_LOGNORMAL, &x )  )
						return error;
						
					x *= sqrt(y);               	/* calc  */
					
					y  = (log(f/k) + 0.5*x*x) / x;  /* d1   */
					x  = y - x;                     /* d2   */

					if (POP() < 0.5)                /* call */
						x = f*norm(y)-k*norm(x);
					else                            /* put  */
						x = -f*norm(-y)+k*norm(-x);
                
					PUSH(x);
					comll = comll->next;
					break;  


				case COMLL_A_BLKSCHLSNRM:          /* Black-Scholes Normal */

					f  = POP();                     /* FWD  */
					k  = POP();                     /* k    */
					x  = POP();                     /* vol  */
					y  = POP();                     /* time difference = mat */
					x *= sqrt(y);               	/* calc */
					y  = (f - k) / x;				/* d1   */
                
					if (POP() < 0.5)                /* call */
						x = x * gauss(y) + (f - k) * norm(y);
					else                            /* put  */
						x = x * gauss(y) - (f - k) * norm(-y);
                
					PUSH(x);
					comll = comll->next;
					break;  
				
				case COMLL_A_SPREADNUM:			/* Spread_numer */
                    
                    fwdx = POP();                     /* FWD X  */
					fwdy = POP();                     /* FWD Y  */
					sigx = POP();                     /* SigX   */
					sigy = POP();                     /* SigY   */
					rho  = POP();                     /* Rho   */
					a	 = POP();                     /* a   */
					b	 = POP();                     /* b   */
					k    = POP();                     /* K   */
					t	 = POP();                     /* t   */
					call_put = (int) POP();                 /* C/P  */
                    
					error = OptSpreadNew(
									fwdx,
									a,
									sigx,
									fwdy,
									-b,
									sigy,
									k,
									t, 
									rho,
									call_put,
									&x);

					if (error)
						return error;
					
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_A_NUMDAYS:			/* number of days in a corridor */
					
					f    = POP();
					spot = POP();
					dBarrier[0]  = POP();
					dBarrier[1] = POP();
					vol  = POP();
					delay = POP();
					mat  = POP();
					nstep= (int) POP();
					dPar[0] = POP();
					dPar[1] = POP();
					i = (int) POP();

					PUSH( prob_ExpectedDays(f, spot, dBarrier, vol, delay, mat, nstep, dPar, i) );
					comll = comll->next;
					break;  

				case COMLL_F_YIELD:             /* Bond yield */

					coupon = POP();
					clean  = POP();
					sdp = *((SwapDP *)comll->ptr);
 
					x = yield_fct(sdp, coupon, clean, 1.0, coupon); 

					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_F_CLEAN:             /* Bond clean price */

					coupon = POP();
					yield  = POP();
					sdp = *((SwapDP*)comll->ptr);
 
					x = clean_price_fct(sdp, coupon, yield, 1.0, coupon); 

					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_B_LT:        /* Less than */

					x = POP();
					y = POP();
					x = (x < y ? one  : zero);
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_B_GT:         /* Greater than */  

					x = POP();
					y = POP();
					x = (x > y ? one  : zero);
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_B_LE:        /* Less or equal */  

					x = POP();
					y = POP();
					if (fabs( y - x) <= DBL_EPSILON)
					{
						 x = one;
					}
					else
					{
						x = ( x <= y  ? one : zero );
					}
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_B_GE:           /* Greater or equal */  

					x = POP();
					y = POP();
					if (fabs(x - y) <= DBL_EPSILON)
					{
						 x = one;
					}
					else
					{
						x = ( x >= y  ? one : zero );
					}
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_B_EQ:            /* Equal */  

					x = POP();
					y = POP();
					if (fabs(x) <= DBL_EPSILON)
					{
						 x = ( ( fabs(y) <= DBL_EPSILON ) ? one  : zero );
					}
					else
					{
						x = ( ( fabs(y/x-1) <= DBL_EPSILON ) ? one : zero );
					}
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_B_AND:

					x = POP();
					y = POP();
					x = (y > .5 && x > .5 ? one  : zero);
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_B_OR:

					x = POP();
					y = POP();
					x = (y > .5 || x > .5 ? one  : zero);
					PUSH(x);
					comll = comll->next;
					break;

			/* ----------------------- Cell Functions ---------------------------------- */
				case COMLL_C_COLMAX:
				case COMLL_C_COLMIN:
				case COMLL_C_COLSUM:
				case COMLL_C_COLAVG:
				case COMLL_C_ROWMAX:
				case COMLL_C_ROWMIN:
				case COMLL_C_ROWSUM:
				case COMLL_C_ROWAVG:
				case COMLL_C_ROWMAXIDX: /* Cyril Godart Regis Benichou fecit per Seculum Secolorum */
				case COMLL_C_ROWMINIDX:

					c1 = comll->ivec[0];
					r1 = comll->ivec[1];
					r2 = comll->ivec[2];
					ind = 0;					/* used only for rowsort */
					grf_f_evlcllfct ( gd->cells, c1, r1, r2, ind, comll->type, &x); 
					PUSH(x);
					comll = comll->next;
					break;
        
				case COMLL_C_COLSORT:
				case COMLL_C_ROWSORT:

					c1 = comll->ivec[0];
					r1 = comll->ivec[1];
					r2 = comll->ivec[2];
					ind = comll->ivec[3];
					grf_f_evlcllfct ( gd->cells, c1, r1, r2, ind, comll->type, &x); 
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_C_ROWINTERP:
					r1  = comll->ivec[0];
					c1  = comll->ivec[1];
					c2  = comll->ivec[2];
					ind = comll->ivec[3];
				/* Get the value to interpolate from the stack */
					x = POP();
				/* Interpolate y=f(x) using the row as Y values and the auxiliary range as X values */
					grf_f_evlcllfct_interp(gd->aux[ind], gd->cells, r1, c1, c2, x, comll->type, &y); 
					PUSH(y);
					comll = comll->next;
					break;

			/* ------------------ State Variables ----------------------------- */
				case COMLL_S_SHORT_RATE:

					/* Get index of underlying */
					ix = comll->ivec[0];

					/*
					PUSH((sample->und[ix].sv[SHORT_RATE]));
					*/
					PUSH( samptr_get(sample, ix, SHORT_RATE) );
					comll = comll->next;

					break;

				case COMLL_S_PHI:

					/* Get index of underlying */
					ix = comll->ivec[0];
                
					/*
					PUSH((sample->und[ix].sv[PHI]));
					*/
					PUSH( samptr_get(sample, ix, PHI) );
					comll = comll->next;
					break;

				case COMLL_S_SIGMA:

					/* Get index of underlying */
					ix = comll->ivec[0];
                
					/*
					PUSH((sample->p[ix].s[SIGMA]));
					*/
					PUSH( samptr_get(sample, ix, SIGMA) );
					comll = comll->next;
					break;

				case COMLL_S_STATEVAR:

					/* Get index of underlying */
					ix = comll->ivec[0];
                
					/*
					PUSH((sample->p[ix].s[STATEVAR]));
					*/
					PUSH( samptr_get(sample, ix, STATEVAR) );
					comll = comll->next;
					break;

				case COMLL_S_NUMERAIRE:

					PUSH(sample->numeraire);
              
					comll = comll->next;
					break;

				case COMLL_S_MMACCOUNT:

					/* Get index of underlying */
					ix = comll->ivec[0];
                
					PUSH( samptr_get(sample, ix, BT) );
					comll = comll->next;
					break;

				case COMLL_S_STOCK:

					/* index of underlying = 0 */
					ix = 0;
                
					PUSH( samptr_get(sample, ix, SPOT) );
					comll = comll->next;
					break;

				case COMLL_S_UNDERLYING:

					/* Get index of underlying */
					ix = comll->ivec[0];
                
					PUSH( samptr_get(sample, ix, SPOT) );
					comll = comll->next;
					break;
			/* --------------------- Use or Increment of Column PV ------------------------------- */
				case COMLL_T_CUMREF:

					c1 = comll->ivec[0];
					x = fwd[c1];
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_T_CURREF:

					c1 = comll->ivec[0];
					x = cur[c1];
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_T_REF:

					c1 = comll->ivec[0];
					x = PVREF(c1);
					PUSH(x);
					comll = comll->next;
					break;

				// IMPLEMENTATION OF PVR IN GRFN
				case COMLL_T_RISKY_REF:
				
					c1 = comll->ivec[0];
					x = PVRISKYREF(c1);
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_T_SET:

					c1 = comll->ivec[0];
					PVSET(c1, *stack);
					comll = comll->next;
					break;
            
				case COMLL_T_INCREMENT:

					c1 = comll->ivec[0];
					PVSET(c1,PVREF(c1) + *stack);
					comll = comll->next;
					break;
// IMPLEMENTATION OF PVR IN GRFN
				case COMLL_T_PVSOURCE_ASSIGN:

					c1 = comll->ivec[0];
					PVSOURCESTATUSSET(c1, POP());
					PVSOURCESET(c1, *stack);
					comll = comll->next;
					break;
				
			/* -----------------------  Interpolation ---------------------------- */
				case COMLL_X_INTERP:

					r1 = comll->ivec[0];
					r2 = comll->ivec[1];
					x = POP();
					y = x;
					interp(gd->aux[r1],gd->aux[r2],
                		IMIN(gd->auxlen[r1],gd->auxlen[r2]),y,0,&x);
					PUSH(x);
					comll = comll->next;
					break;

				case COMLL_X_PVINTERP:
					ind = comll->ivec[0];
//					r1  = comll->ivec[1];
					c1  = comll->ivec[2];
					c2 = comll->ivec[3];
				/* Get the value to interpolate from the stack */
					x = POP();
				/* Interpolate y=f(x) using the row as Y values and the auxiliary range as X values */
					y = gd->auxlen[ind];
					grf_f_evlcllfct_interp(gd->aux[ind], gd->cells, r1, c1, c2, x, comll->type, &y); 
					PUSH(y);
					comll = comll->next;
					break;


				case COMLL_PRINT:

					if ( gd->outfileptr != NULL )
					{
						fprintf ( gd->outfileptr,
								  "%.6lf#%.6lf#%.6lf#\n",
								  samptr_get(sample, ix, 0),
								  samptr_get(sample, ix, 1),
								  *stack );
					}

					comll = comll->next;
					break;

			/* --------------------- Histogram ------------------------------------ */
				case COMLL_HIST:

					/* Get Monte-Carlo path number */
					ipath = sample->pathnum;

					/* Add data to histogram: This function creates the required
					   data structure (if necessary) and subsequently stores the
					   given value, x, in that structure (assuming POP has stored it
					   Note: comll->sval = Histogram label 
					*/

					if ((comll->prev != NULL) && (comll->prev->type == COMLL_POP)  )
						x = comll->prev->dval;
					else 
						return serror("the right syntax for hist is expr1|hist(\"name\")");

					error = srt_f_sendhistval(comll->sval, ipath, x);
					if (error)
						return error; 

					PUSH(x);
					comll = comll->next;
					break;
			/* -------------------- A Real Value : Push it in the Stack ----------- */
				case COMLL_REAL:
				default:

					PUSH(comll->dval);
					comll = comll->next;
					break;

			} /* END switch */
		}
		else 
		/* the cell cannot be evaluated in the current computing mode forward/backward */
		{
			comll = comll->nextcol;
		}

    } /* END while */

/* The Event Cash Flow is the last element of the Stack */    
 	if (stack != s)  *cashflow = POP();

/* Checks here that the stack is empty (otherwise something missing) */
	if (*stack != 0.0)
		smessage("Stack not empty in GrfnEvalEvent: %f",(double)*stack);

// IMPLEMENTATION OF PVR IN GRFN
/* Assign the value of pvNR to pvR*/
	memcpy( pvR_vect, curr_pv, gd->sswidth*sizeof(double) );

    return NULL;
}



   
   
