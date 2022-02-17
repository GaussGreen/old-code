/* =============================================================================
   
   FILENAME:    swp_f_swap_simple.c

   PURPOSE:     Provide a few function for simple swap calculations:
				- PV
				- Fwd Rate
				- IRR
				- Duration
				- Convexity
   ============================================================================= */

#include "swp_h_all.h"
#include "opfnctns.h"
#include "swp_h_swap_simple.h"


Err zcswap(
	SwapDP 		*sdp,
	double 		strike, 
	double 		ini, 
	double 		fin, 
	SwapOutput 	m, 
	String 	    ycname, 
	Date 		value_date, 
	double 		*answer)

{
SimpleSwap 	  s;
double 		  ldf;
double 		  pv;
double        frate;		
double        irr;
double        duration;
double        convexity;
double		  thirdmoment;
Date 		  today;
SrtErr		  err ;
SrtCurvePtr  crv;

	crv = lookup_curve(ycname);
	today = get_clcndate_from_yldcrv(crv);

	switch(m)
	{
		case COMPUTE_PV:
			ldf = swp_f_df((Ddate)today,(Ddate)value_date,ycname);
			s = make_SimpleSwap(sdp,strike,ini,fin,today,SWAP);
			pv = value_SimpleSwap(&s,crv)/ldf;
			free_inSimpleSwap(&s);
			*answer = pv;
			break;
		case COMPUTE_FWD_RATE:
			s = make_SimpleSwap(sdp,1.0,0.0,0.0,today,SWAP);
			frate = frate_SimpleSwap(&s,crv);
			free_inSimpleSwap(&s);
			*answer = frate;
			break;
		case COMPUTE_IRR:
			s   = make_SimpleSwap(sdp,strike,ini,fin,today,BOND);
			ldf = swp_f_df((Ddate)today,(Ddate)sdp->start,ycname);
			pv  = value_SimpleSwap(&s,crv)/ldf;
  			if(err=irr_range(
					s.dl.date, 
					s.cpn.d,
					s.cpn.len,
					s.dl.type,
					s.dl.prev,
					pv,
					s.sdp.basis_code,
					s.sdp.compd,
					&irr)) 
				return err;
			free_inSimpleSwap(&s);
			*answer = irr;
			break;
		case COMPUTE_DURATION:
			s   = make_SimpleSwap(sdp,strike,ini,fin,today,BOND);
			ldf = swp_f_df((Ddate)today,(Ddate)sdp->start,ycname);
			pv  = value_SimpleSwap(&s,crv)/ldf;
  			if(err=irr_range(
					s.dl.date, 
					s.cpn.d,
					s.cpn.len,
					s.dl.type,
					s.dl.prev,
					pv,
					s.sdp.basis_code,
					s.sdp.compd,
					&irr)) 
				return err;
  			if(err=duration_range(
					s.dl.date, 
					s.cpn.d,
					s.cpn.len,
					s.dl.type,
					s.dl.prev,
					s.sdp.basis_code,
					s.sdp.compd,
					irr,
					&duration))
				return err;
			free_inSimpleSwap(&s);
			*answer = duration;
			break;
		case COMPUTE_CONVEXITY:
			s   = make_SimpleSwap(sdp,strike,ini,fin,today,BOND);
			ldf = swp_f_df((Ddate)today,(Ddate)sdp->start,ycname);
			pv  = value_SimpleSwap(&s,crv)/ldf;
  			if (err=irr_range(
					s.dl.date, 
					s.cpn.d,
					s.cpn.len,
					s.dl.type,
					s.dl.prev,
					pv,
					s.sdp.basis_code,
					s.sdp.compd,
					&irr)) 
				return err;
  			if(err=convexity_range(
					s.dl.date, 
					s.cpn.d,
					s.cpn.len,
					s.dl.type,
					s.dl.prev,
					s.sdp.basis_code,
					s.sdp.compd,
					irr,
					&convexity))
				return err;
			free_inSimpleSwap(&s);
			*answer = convexity;
			break;
		case COMPUTE_THIRDMOMENT:
			s   = make_SimpleSwap(sdp,strike,ini,fin,today,BOND);
			ldf = swp_f_df((Ddate)today,(Ddate)sdp->start,ycname);
			pv  = value_SimpleSwap(&s,crv)/ldf;
  			if (err=irr_range(
					s.dl.date, 
					s.cpn.d,
					s.cpn.len,
					s.dl.type,
					s.dl.prev,
					pv,
					s.sdp.basis_code,
					s.sdp.compd,
					&irr)) 
				return err;
  			if(err=thirdmoment_range(
					s.dl.date, 
					s.cpn.d,
					s.cpn.len,
					s.dl.type,
					s.dl.prev,
					s.sdp.basis_code,
					s.sdp.compd,
					irr,
					&thirdmoment))
				return err;
			free_inSimpleSwap(&s);
			*answer = thirdmoment;
			break;
		case COMPUTE_MODIFIED_DURATION:
			s   = make_SimpleSwap(sdp,strike,ini,fin,today,BOND);
			ldf = swp_f_df((Ddate)today,(Ddate)sdp->start,ycname);
			pv  = value_SimpleSwap(&s,crv)/ldf;
  			if(err=irr_range(
					s.dl.date, 
					s.cpn.d,
					s.cpn.len,
					s.dl.type,
					s.dl.prev,
					pv,
					s.sdp.basis_code,
					s.sdp.compd,
					&irr)) 
				return err;
  			if(err=modified_duration_range(
					s.dl.date, 
					s.cpn.d,
					s.cpn.len,
					s.dl.type,
					s.dl.prev,
					s.sdp.basis_code,
					s.sdp.compd,
					irr,
					&duration))
				return err;
			free_inSimpleSwap(&s);
			*answer = duration;
			break;
		case COMPUTE_MODIFIED_CONVEXITY:
			s   = make_SimpleSwap(sdp,strike,ini,fin,today,BOND);
			ldf = swp_f_df((Ddate)today,(Ddate)sdp->start,ycname);
			pv  = value_SimpleSwap(&s,crv)/ldf;
  			if (err=irr_range(
					s.dl.date, 
					s.cpn.d,
					s.cpn.len, 
					s.dl.type,
					s.dl.prev,
					pv,
					s.sdp.basis_code,
					s.sdp.compd,
					&irr)) 
				return err;
  			if(err=modified_convexity_range(
					s.dl.date, 
					s.cpn.d,
					s.cpn.len, 
					s.dl.type,
					s.dl.prev,
					s.sdp.basis_code,
					s.sdp.compd,
					irr,
					&convexity))
				return err;
			free_inSimpleSwap(&s);
			*answer = convexity;
			break;
		case COMPUTE_MODIFIED_MATCHING_RATIO:
			s   = make_SimpleSwap(sdp,strike,ini,fin,today,BOND);
			ldf = swp_f_df((Ddate)today,(Ddate)sdp->start,ycname);
			pv  = value_SimpleSwap(&s,crv)/ldf;
  			if (err=irr_range(
					s.dl.date, 
					s.cpn.d,
					s.cpn.len, 
					s.dl.type,
					s.dl.prev,
					pv,
					s.sdp.basis_code,
					s.sdp.compd,
					&irr)) 
				return err;
  			if(err=modified_duration_range(
					s.dl.date, 
					s.cpn.d,
					s.cpn.len,
					s.dl.type,
					s.dl.prev,
					s.sdp.basis_code,
					s.sdp.compd,
					irr,
					&duration))
				return err;
			if(err=modified_convexity_range(
					s.dl.date, 
					s.cpn.d,
					s.cpn.len, 
					s.dl.type,
					s.dl.prev,
					s.sdp.basis_code,
					s.sdp.compd,
					irr,
					&convexity))
				return err;
			free_inSimpleSwap(&s);
			*answer = convexity/duration;
			break;
		default:
			return "Unknown argument to zcswap()";
	}

	
	return NULL;
}

/* --------------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------------
   This function is intended as the interface function for zcswap
   --------------------------------------------------------------------------------- */
Err swp_f_SimpleSwap(
		long        start,
		long        end_nfp,
		String      compStr,
		String      basisStr,
		double      coupon,
		double      initial,
		double      final,
		String      ycname,
		String      strMessage,
		long        value_date,
		double      *answer)
{
Err         err    = NULL;
SwapDP      swapdp;
SwapOutput  output;
SrtCurvePtr   yc_crv;
Date        clcn_date;
int         spot_lag;

/* Gets the swap yield curve to extract some date information */
	yc_crv = lookup_curve(ycname);
	if (!yc_crv)
		return serror("Could not find %s curve",ycname);

	clcn_date = get_clcndate_from_curve(yc_crv);
	spot_lag  = get_spotlag_from_curve(yc_crv);;

	err = swp_f_initSwapDP(	start, end_nfp, compStr, basisStr, &swapdp);
	if (err)
		return err;
	swapdp.spot_lag = spot_lag;

	err =swp_f_interp_swap_message( strMessage, &output);
	if (err)
		return err;

	
	if (!value_date)
		value_date = start;
	err = zcswap(&swapdp, coupon, initial, final, output, ycname, value_date, answer);
	if (err)
		return err;

	return NULL;
}

/* ------------------------------------------------------------------------------------- */


Err zcsens(String tenor, SrtCurvePtr yldcrv, double *answer)
{
	double 		dt,
			sens;
	double 		df;
	Date 		today;
	int 		nbrm=0;
	YCInstr 	instr;
	Err 		err;
	SrtCcyParam 	*ccy_param;


	ccy_param = get_ccyparam_from_yldcrv(yldcrv);   

	today = get_clcndate_from_yldcrv(yldcrv);

	if(err = string_to_YCInstr(tenor, &instr, NULL, 
			ccy_param, add_unit(today, ccy_param->spot_lag, SRT_BDAY, MODIFIED_SUCCEEDING)))
	{
		return err;
	}

	switch(instr.type){
	case SWAPINSTR:
		if(err = zcswap(&instr.dp.swapdates, 1.0, 0.0, 0.0,
				COMPUTE_PV, get_curve_name(yldcrv), today, &sens))
		{
			return err;
		}
		*answer = sens;
		break;
        default:
   		dt = coverage(instr.dp.futdates.start,
			instr.dp.futdates.end,
			instr.dp.futdates.basis_code);
		df = swp_f_df((Ddate)today, (Ddate)(instr.dp.futdates.end),get_curve_name(yldcrv));
		*answer = df * dt;
		break;
	}

	return NULL;
}


/* ------------------------------------------------------------------------------ */


SimpleSwap make_SimpleSwap(SwapDP *sdp, double strike, 
			double ini, double fin,
			Date today, StructType t)
{
	SimpleSwap cfs;
	Err			err=NULL;
	DateList start_date;
	DateList end_date;
	long i;


	cfs.sdp = *sdp;
	cfs.dl = SwapDP_to_DateList(sdp,MODIFIED_SUCCEEDING);
	
	start_date = new_DateList(cfs.dl.len - 1);
	end_date = new_DateList(cfs.dl.len - 1);
	for (i = 0; i< cfs.dl.len-1;i++)
	{
		start_date.date[i] = cfs.dl.date[i];
		end_date.date[i] = cfs.dl.date[i+1];
	}
	
	err = cpn_list(
			sdp->basis_code,
			sdp->compd,
			start_date,
			end_date,
			cfs.dl,
			strike,
			ini,
			fin,
			t,
			&(cfs.cpn));
	
	err = time_list(cfs.dl,today,&(cfs.t));
	
	swp_f_free_in_DateList(start_date);
	swp_f_free_in_DateList(end_date);

	return cfs;
}

/* ------------------------------------------------------------------------------ */


int free_inSimpleSwap(SimpleSwap *s)
{
	srt_free(s->dl.date);
	free_inDlist(&s->cpn);
	free_inDlist(&s->t);

	return( 0);
}
	
/* ------------------------------------------------------------------------------ */

double value_SimpleSwap(SimpleSwap *s, SrtCurvePtr m)
{
	Dlist df;
        int i;
	double pv = 0;
	Err		err = NULL;

	err= df_list(s->dl,get_curve_name(m),&df);
	for(i=0;i<s->dl.len;i++)	
		pv += df.d[i] * s->cpn.d[i];

	free_inDlist(&df);
	return pv;
}

/* ------------------------------------------------------------------------------ */

double frate_SimpleSwap(SimpleSwap *s, SrtCurvePtr m)
{
	Dlist df;
        int i;
	double pv = 0;
	Err		err = NULL;

	err = df_list(s->dl, get_curve_name(m),&df);
	for(i=0;i<s->dl.len;i++)	
		pv += df.d[i] * s->cpn.d[i];
	pv = (df.d[0] - df.d[df.len - 1])/pv;

	free_inDlist(&df);
	return pv;
}



/* --------------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------------
   This could be more useful when the PV is not quite the sum of discounted
   cashflows ( a bond, for instance) (we assume no initial exchange, and 1.0 as
   the final one)
   PV is the pv of the swap at the swap start date
   All the calculations are done on the CLEAN price (which is not bond math)
   ------------------------------------------------------------------------------ */
Err bond_modified_duration(
	SwapDP 		*sdp,
	double 		coupon, 
	double      pv,
	double 		*answer)
{
SimpleSwap 	  s;
double        irr;
double        duration;
Err           err = NULL;

	s   = make_SimpleSwap(sdp, coupon, 0.0, 1.0, sdp->start, BOND);
	if(err=irr_range(
			s.dl.date, 
			s.cpn.d,
			s.cpn.len,
			s.dl.type,
			s.dl.prev,
			pv,
			s.sdp.basis_code,
			s.sdp.compd,
			&irr)) 
		return err;
  	if(err=modified_duration_range(
			s.dl.date, 
			s.cpn.d,
			s.cpn.len,
			s.dl.type,
			s.dl.prev,
			s.sdp.basis_code,
			s.sdp.compd,
			irr,
			&duration))
		return err;
	free_inSimpleSwap(&s);
	
	*answer = duration;

	return NULL;
}

/* ------------------------------------------------------------------------------ */

Err bond_modified_convexity(
	SwapDP 		*sdp,
	double 		coupon, 
	double      pv,
	double 		*answer)
{
SimpleSwap 	  s;
double        irr;
double        convexity;
Err           err = NULL;

	s   = make_SimpleSwap(sdp, coupon, 0.0, 1.0, sdp->start, BOND);
	if(err=irr_range(
			s.dl.date, 
			s.cpn.d,
			s.cpn.len,
			s.dl.type,
			s.dl.prev,
			pv,
			s.sdp.basis_code,
			s.sdp.compd,
			&irr)) 
		return err;
  	if(err=modified_convexity_range(
			s.dl.date, 
			s.cpn.d,
			s.cpn.len, 
			s.dl.type,
			s.dl.prev,
			s.sdp.basis_code,
			s.sdp.compd,
			irr,
			&convexity))
		return err;
	free_inSimpleSwap(&s);
	
	*answer = convexity;

	return NULL;
}

/* ------------------------------------------------------------------------------ */





