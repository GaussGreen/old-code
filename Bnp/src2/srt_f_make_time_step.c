
/*****************************************************************************
   MODULE	: srt_F_MAKE_TIME_STEP.C
   FUNCTION	: srt_f_make_time_step
   AUTHOR	: E.Auld G.Amblard from function chell() by A.Litke 	
   DATE		: APRIL 94


<%%STA
  FUNCNAME        :srt_f_make_time_step 
  DESCRIPTION     : creates list of Steps to represent discretization of 
   underlying model over time, and populate them, 
   Steps are inserted for: 
   * times when cashflows occur 
   * times when changes occur in grfnparam parameters e.g. 
     term structure of volatility
   * increasing accuracy of discretization
     Discretization takes underlying mdl and numerical method into
     account. 

  MODIFIES	:*pointer_to_info  
  CALL          :

<%%END

   AMENDMENTS:
   Reference	:
   Author     	: K L Chau
   Date		: June 16, 1994
   Description	: ins_vol_dates no longer applicable; interpolate vol when
                  necessary
   Reference	:
   Author     	: E Auld
   Date		: July 13, 1994
   Description	: ins_vol_dates still an option, though in the long run will
		probably not be used; also changed insert vol dates not to
		insert last vol point (since vol is the same on either side
		of the last point, there is no information lost by not 
		inserting it)

   Reference	:
   Author     	: E Auld
   Date		: July 21, 1994
   Description	: undid previous change to insert vol dates, to prune
		sources of changed in value in deals currently priced
		by GRFN.  This remains a correct thing to do and will have
		to be redone when there is more time.

   Reference    :
   Author       : K L Chau
   Date         : August 16, 1994
   Description  : restructure program, so that the main routine is separated
                  into various subroutines, to enhance flexibility.
                  Some comments:
                  1) insert dates according to "sigma^2 * time".
                  2) code concerning the fields "center" and "dcenter" are
                     commented out in grfn_mdl_init_step.  If any problems arise
                     they should be put back in. 

   AMENDMENTS:
   Reference	:
   Author     	:E.AULD
   Date		:22 jul 94 (merged with this version 31 aug 94)
   Description	:removed adding of grfn events (to grfn_attach_events);
                 renamed file to srt_f_mdl_init_step.c, removed references
		 to grfn throughout file.


   AMENDMENTS:
   Reference	:
   Author     	:A. Sahuguet
   Date		:12 oct 94
   Description	:modified function in order to be used by the new 
		 TermStructure (2-linked list) 

   AMENDMENTS:
   Reference	:
   Author     	: K L Chau
   Date		: Dec 6, 1994
   Description	: use new srtmkt structure

   AMENDMENTS:
   Reference	:
   Author     	: O Van Eyseren
   Date		    : Aug 22, 1996
   Description	: renamed to srt_f_make_time_step

******************************************************************************/
#include "math.h"
#include "srt_h_all.h"

static int insert_ir_und_type_vol_dates(TermStruct *ts, SrtStpPtr step, double today);

/* --------------------------------------------------------------- */
/* PUT NODES ON FOR DATES OF GRFN DEAL */

static Err insert_eventdates(
		SrtStpPtr   step, 
		Date       *eventdates,
		long        numeventdates,
		double      today)
{
Err err = NULL;
long i;

/* If no event dates: return */
	if ( (!eventdates) || (numeventdates == 0) )
		return NULL;

/* We don't have to insert the first event date if it is today,
   because it should already have been inserted in the main loop */
	i = 0;

/* Insert all the event dates BEFORE and AFTER today (no need to reinsert today) */
	while ( i < numeventdates )
	{
	/* If date is before today, insert a step before */
		if ( eventdates[i] < today )
		{
			step = insert_before(step);
			step->prev->date = eventdates[i];
			step->prev->ddate = (double)step->prev->date;
			step->prev->time = (step->prev->ddate - today) * YEARS_IN_DAY;
		}
		else
	/* If date is after today, insert a step after */
		if ( eventdates[i] > today )
		{	
			step = insert_after(step);
			step->date = eventdates[i];
			step->ddate = (double)step->date;
			step->time = (step->ddate - today) * YEARS_IN_DAY;
		}

	/* Move on to next event date */
		i++;

	}
  
	return NULL;

}

/* -------------------------------------------------------------- */

SrtLst* get_next_sigma_date(SrtLst* ts)
{
SrtLst* ls;
	ls = ts;
	while ( (ls!=NULL)
		&& 
	( ((IrmTermStructVal*)ls->element->val.pval)->val_origin != SIGMA_DATE )
		&&
	( ((IrmTermStructVal*)ls->element->val.pval)->val_origin != BOTH_DATE )
	)
	ls = ls->next;

	return ls;
}

SrtLst* get_next_tau_date(SrtLst* ts)
{
SrtLst* ls;
	ls = ts;
	while ( (ls!=NULL)
		&& 
	( ((IrmTermStructVal*)ls->element->val.pval)->val_origin != TAU_DATE )
		&&
	( ((IrmTermStructVal*)ls->element->val.pval)->val_origin != BOTH_DATE )
	)
	ls = ls->next;

	return ls;
}

/* ----------------------------------------------------------------------- */

static Err insert_vol_dates(SrtUndInfo *und_info,
					 SrtStpPtr step, 
					 double today)
{
SrtUnderlyingType		und_type;
TermStruct		*ts;
SrtUndPtr		und;
Err err = NULL;
long i;
SrtLst	*ls;
Date   term_struct_date;
SrtMdlType	mdl_type;

	
	i = 0; /* Loop on the no of underlying */
	while(i < und_info->no_of_underlyings)
	{
		und = lookup_und(und_info->und_data[i].und_name);
		if(!und) return serror("can not get underlying ");

		err = get_underlying_mdltype(und,&mdl_type);
		if(err) return err;

		err = srt_f_interp_under(und->underl_lbl, &und_type);
		if(err) return err;

		err = get_underlying_ts(und,&ts);
		if(err) return err;

		if(und_type == INTEREST_RATE_UND)  
			insert_ir_und_type_vol_dates(ts,step,today);
		else if( (und_type == EQUITY_UND) && (mdl_type == EQ_STOCH_RATES_SRVGS) )
		{
			step = gototop(step);
			ls = ts->head;
			while(ls)
			{
				term_struct_date = (Date) ((EquityTermStructVal*)ls->element->val.pval)->date;

				while( (step->date < term_struct_date) && (step->next) ) step = step->next;

				if(step->date >= term_struct_date)
				{
					if(step->date == term_struct_date) ls = ls->next;
					else 
					{
						step = insert_before(step);
						step->prev->date = term_struct_date;
						step->prev->ddate = (double)step->prev->date;
						step->prev->time = (step->prev->ddate-today)*YEARS_IN_DAY;

						ls = ls->next;
					}
				}
				else ls = ls->next;
			}
		}

		i++;
	}

		return NULL;
}
static int insert_ir_und_type_vol_dates(TermStruct *ts, SrtStpPtr step, double today)
{

  SrtStpPtr top,bot;
  SrtLst* ls;
  IrmTermStructVal* tsval;

/* PUT DATES ON CURVE FOR TAU'S AND SIGMA'S WHERE FUNCTION CHANGES */
/* THIS MAKES SURE THAT TAU AND SIGMA ARE CONSTANT OVER A TIME STEP */ 

  top = gototop(step);
  bot = gotobot(step);
   
    step=top;
    ls = ts->head;
    ls = get_next_sigma_date(ls);
	
    while ( ( ls!= NULL ) && ( bot->date> (tsval = (IrmTermStructVal*)(ls->element->val.pval) )->date ) )
    {
      while (step->date < tsval->date)step=step->next; 
      if(step->date == tsval->date) 
	  {
		ls = ls->next; /* we take the next date */
		ls = get_next_sigma_date(ls); /* we check it's a sigma date */
	  }	
/* ONLY INSERT IF DATE NOT TOO CLOSE TO EVENT DATE. E.AULD */
/*      if(fabs((double)step->ddate-(double)sig_dates[i])<3.5)i++; */
      else
      {
        step=insert_before(step);
        step->prev->date = (Date)tsval->date;
        step->prev->ddate = (double)step->prev->date;
        step->prev->time = (step->prev->ddate-today)*YEARS_IN_DAY;
        ls = ls->next;
		ls = get_next_sigma_date(ls);
      }
    }

    step=top;
    ls = ts->head;
    ls = get_next_tau_date(ls);

    while( ( ls!= NULL ) && ( bot->date > (tsval = (IrmTermStructVal*)(ls->element->val.pval))->date ))
    {
      while (step->date < tsval->date)step=step->next; 
      if(step->date == tsval->date) 
		{
		ls = ls->next; /* we take the next date */
		ls = get_next_tau_date(ls); /* we check it's a tau date */
		}	
/* ONLY INSERT IF DATE NOT TOO CLOSE TO EVENT DATE. E.AULD */
/*      if(fabs((double)step->ddate-(double)sig_dates[i])<3.5)i++;*/
      else
      {
        step=insert_before(step);
        step->prev->date = (Date)tsval->date;
        step->prev->ddate = (double)step->prev->date;
        step->prev->time = (step->prev->ddate-today)*YEARS_IN_DAY;
        ls = ls->next;
	ls = get_next_tau_date(ls);
      }
    }

  return 0;
}

/* ----------------------------------------------------------------------- */

static Err insert_constrained_dates(
			SrtGrfnParam *grfnparam,
			SrtStpPtr step, 
			double today,
			SrtUndInfo *und_info)
{
Err			err = NULL;
long		i;
SrtUndPtr	und_info_und;
SrtMdlType	mdl_type;

/* Check if one of the underlying is of EQ_STOCH_RATES_SRVGS type */

	i = 0;
	while( (i < und_info->no_of_underlyings) && (grfnparam->insert_vol_dates == SRT_NO) )
	{
		und_info_und = lookup_und(und_info->und_data[i].und_name);
		err = get_underlying_mdltype(und_info_und, &mdl_type);
		if(err) return err;

		if(mdl_type == EQ_STOCH_RATES_SRVGS) grfnparam->insert_vol_dates = SRT_YES;

		i++;
	}

/* Check if it is needed to insert vol dates */ 

  if(
     (	
        !und_info->jumping ||
     	grfnparam->imp_type 	!= 	IMPMONTECARLO
     )  && (und_info->no_of_underlyings!=1)
	&& ((grfnparam->insert_vol_dates == SRT_YES) )          
    )
      if (err = insert_vol_dates(und_info, step, today))
	  return err;

  return NULL; 
}

/* ----------------------------------------------------------------------- */

/* -----------------------------------------------------------------------

   Different intermediate date insertion methods should be put here

   ------------------------------------------------------------------------ */

/* Inserts extra points (for discretisation) only for FUTURE dates (after today ) */
static int insert_extra_pts_in_future(
		SrtGrfnParam *grfnparam,
		SrtStpPtr     step,
		double        today,
		TermStruct   *ts)
{
/* PUT EXTRA DATES ON CURVE FOR MINIMUM DAYS PER TIME STEP  */

SrtStpPtr      top, bot;
long           i,n;
double         delta_t,max_time_per_slice;

    top = gototop(step);
    bot = gotobot(step); 	

/* Go to today  (SKIP EVENTS IN THE PAST) */
	step = top;
	while ((step->ddate < today ) && (step->next != NULL))
		step= step->next;

/* Insert dates according to time step division */
    max_time_per_slice=grfnparam->max_time_per_slice;
    delta_t= bot->time ;
	if(max_time_per_slice>(delta_t/(double)grfnparam->min_nodes_per_path))
      max_time_per_slice=(delta_t/(double)grfnparam->min_nodes_per_path);  

/* Insert all the extra steps */
    while(step->next!=NULL)
    {
		delta_t = step->next->time - step->time;
		n = DTOL(floor(delta_t/max_time_per_slice))-1;
		if(n>0)
			delta_t *= DAYS_IN_YEAR/(n+1);
		for(i=0;i<n;i++)
		{
			step = insert_after(step);
			step->date = step->prev->date + (long)(delta_t);
			step->ddate = step->prev->ddate + delta_t;
			step->time = (step->ddate-today)*YEARS_IN_DAY;
		}
		step = step->next;
	}
	return 0;	
}


/* ----------------------------------------------------------------------- */

/* -----------------------------------------------------------------------

   Set up the following parameters in each time step:
			delta_t, sqrt_delta_t 
   ----------------------------------------------------------------------- */

static Err populate_steps(SrtStpPtr step, double today)
{

/* Starts from the first time step */
    step = gototop(step);

/* Sets general parameters at each step  */
    while(step->next!=NULL)
    {
      step->delta_t = step->next->time-step->time;
      step->sqrt_delta_t = sqrt(step->delta_t);
      step = step->next;
    }

/* Return a success message */
  return NULL;
}

/* ----------------------------------------------------------------------- */

/* ------------------------------------------------------------------
	
    MAIN FUNCTION

   ------------------------------------------------------------------ */

Err srt_f_make_time_step
	(
	long            numeventdates,
	Date           *eventdates, 
	SrtUndPtr       und, 
	SrtStpPtr      *stpptr, 
	SrtGrfnParam   *grfnparam,
	SrtUndInfo     *und_info
	)
{                          
Err           err;
SrtStpPtr     top,step=NULL;
double        today ;
TermStruct   *ts ;

/* Get the underlying term structure */
	err = get_underlying_ts(und,&ts) ;

/*  Initialize the first node  as being TODAY */
	step = insert_after(step);
	step->ddate= today = get_today_from_underlying(und);
	step->time = (step->ddate-today)*YEARS_IN_DAY;
	step->date = DTOL(step->ddate);

/* Insert the event dates */
	err = insert_eventdates(
			step,
			eventdates,
			numeventdates,
			today);
	if (err)
		return err;

/* Insert constrained dates (vol dates etc. if needed) */
	err = insert_constrained_dates(
			grfnparam, 
			step, 
			today,
			und_info);

    if (err)
	   return ("Error in inserting constrained dates");


/* Insert intermediate dates for future dates */
/* Always insert extra points - they would be removed later if jumping
   numeraire in MC is used by rem_eventless_node in grfn_val_Grfndeal */
	if (insert_extra_pts_in_future(
			grfnparam,
			step,
			today,
			ts) > 1)
       return ("Error in inserting extra points"); 

/* Populate steps with the extra time information (delta t,...) */   
    err = populate_steps(
			step,
			today);
	if (err)
       return ("Error in populating steps");
    
/* Return the first time step */
	top = gototop(step);
    *stpptr = top;
    
	return(NULL);

} /* END Err srt_f_make_time_step(...) */

/* ----------------------------------------------------------------------------------- */