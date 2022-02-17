/* ------------------------------------------------------------------------------ 
   FILENAME:      swp_f_bond_swaption.c

   PURPOSE:       Given a bond option, finds the swaption that best matches it 
   ------------------------------------------------------------------------------ */
#include "math.h"
#include "num_h_allhdr.h"
#include "swp_h_all.h"
#include "swp_h_bond_swaption.h"
#include "swp_h_swap_simple.h"
#include "opfnctns.h"

#define    FORWARD        1
#define    BACKWARD      -1
#define    NO_DIRECTION   0


/* We return the swaption maturity, strike and volatility, as well as vega ratio and 
hedge ratio.
If there is a problem to compute the volatility or the hedge ratio or the vega ratio,
we return just the swaption maturity and strike and put -1 for the volatility and 1 for
the hedge ratio and vega ratio */

Err swp_f_match_bond_swaption(
		long			 option_expiry,
		long             settle_bond_date,
		long             maturity_date,
		String           strBondComp,
		String           strBondBasis,
		double           coupon,  /* as 0.05 for 5% */
		double           fwd_clean_price,  /* as 1.00 for 100 */
		double           clean_strike,
		String           strCallPut,
		double           bond_option_value,
		String           yc_name,
		String           refRateCode,
		String           strSwapComp,
		String           strSwapBasis,
		String           strLogorNormal,

		long			 *swap_settle,
		long             *swap_maturity,
		double           *swap_coupon,
		double           *implied_vol,
		double           *vega_ratio,
		double           *hedge_ratio)
{
SwapDP            swap_dp;
SwapDP            bond_dp;
SrtCurvePtr       yc_crv = NULL;
Err               err = NULL;
double            acc_int;
double            fwd_dirty_price;
double            bond_duration;
double            bond_convexity;
double            bond_ratio;
double            bond_price_vol;
double            swap_duration;
double            swap_ratio;
double            swap_forward_level;
double            swaption_value;
double            swap_price_vol;
double            forward_df;
double            settle_bond_df;
double			  settle_swap_df;
double            option_maturity;
double            price_vol_ratio;
double            db, ds;
double            target;
long              new_direction, old_direction;
int               spot_lag;
Date              clcn_date;
Date              last_pay_date;
Date			  settle_swap_date;
SRT_Boolean           stop_flag;
SrtCallPutType    call_put;
SrtReceiverType   rec_pay;
SrtDiffusionType  log_or_norm;

double            real_swap_forward;
double            cash_swap_forward;
double            swaption_strike;

/* Initialisation of the returns */
	*hedge_ratio = 1;
	*vega_ratio = 1;
	*implied_vol = -1;

/* Gets the swap yield curve */
	yc_crv = lookup_curve(yc_name);
	if (!yc_crv)
		return serror("Could not find %s curve",yc_name);
	clcn_date = get_clcndate_from_yldcrv(yc_crv);
	spot_lag  = get_spotlag_from_curve(yc_crv);
	settle_swap_date = add_unit(option_expiry,spot_lag,SRT_BDAY, SUCCEEDING);
	*swap_settle = settle_swap_date;

/* Sets the SwapDP for the bond */
	err = swp_f_initSwapDP(	settle_bond_date, maturity_date, strBondComp, strBondBasis, &bond_dp);
	if (err)
		return err;

/* Compute the accrued interests to get the fwd dirty bond price */
	acc_int = acc_int_fct(bond_dp,coupon);
	fwd_dirty_price = fwd_clean_price + acc_int;

/* Compute the forward ratio convexity/duration for the forward bond */
	err = bond_modified_convexity( &bond_dp, coupon, fwd_dirty_price, &bond_convexity);
	if (err)
		return err;
	err = bond_modified_duration( &bond_dp, coupon, fwd_dirty_price, &bond_duration);
	if (err)
		return err;
	
	bond_ratio = bond_convexity / bond_duration;

	
/* Initialise a few flags before starting iterations on swap maturity */
	stop_flag = SRT_NO;
	old_direction = NO_DIRECTION;
	new_direction = NO_DIRECTION;

/* Loops on days until ratios are matched (Newton would be useless: days are integers) */
	*swap_maturity = maturity_date;
	while (stop_flag != SRT_YES)
	{
	
	/* Initialises a SwapDP for the swap with same maturity date as the bond */
		err = swp_f_initSwapDP(	settle_swap_date, *swap_maturity, strSwapComp, strSwapBasis, &swap_dp);
		if (err)
			return err;
		swap_dp.spot_lag = spot_lag;
	
	/* Sets the coupon of the swap so that forward values match */
		last_pay_date = bus_date_method(*swap_maturity, MODIFIED_SUCCEEDING);
		forward_df = swp_f_df(settle_swap_date, last_pay_date, yc_name);
		err = swp_f_Level_SwapDP( &swap_dp, yc_name, &swap_forward_level);
		swap_forward_level /= swp_f_df(clcn_date, settle_swap_date,yc_name);
		if (err)
			return err;
	
		*swap_coupon = ( fwd_clean_price /  clean_strike - forward_df ) / swap_forward_level;

	/* Compute the forward ratio convexity/duration for the forward swap */
		err = zcswap(&swap_dp, *swap_coupon, 0.0, 1.0, COMPUTE_MODIFIED_MATCHING_RATIO, 
			yc_name, settle_swap_date, &swap_ratio);
		if (err)
			return err;

	/* Determines which way to go: if ratio too high, increase maturity, else decrease it*/
		if (swap_ratio < bond_ratio)
			new_direction = FORWARD;
		else
		if (swap_ratio > bond_ratio)
			new_direction = BACKWARD;
		else
			stop_flag = SRT_YES;
		
		if (old_direction * new_direction == -1 )
		{
			stop_flag = SRT_YES;
		}
		else
		{
		/* Prepares the next iteration: moves by ONE day */
			if (new_direction == FORWARD)
				*swap_maturity = add_unit(*swap_maturity, 1, SRT_DAY, NO_BUSDAY_CONVENTION);
			else
			if (new_direction == BACKWARD)
				*swap_maturity = add_unit(*swap_maturity, -1, SRT_DAY, NO_BUSDAY_CONVENTION);
	
			old_direction = new_direction;
		}

	} /* END while loop on maturity date for the swap */

/* Computes the modified durations for the swap */	
	err = zcswap(&swap_dp, *swap_coupon, 0.0, 1.0, COMPUTE_MODIFIED_DURATION, 
		yc_name, settle_swap_date, &swap_duration);
	if (err)
		return err;

/* Checks the type of the bond option: sets the swaption type accordingly */
	err = interp_call_put (strCallPut, &call_put);
	if (call_put == SRT_CALL)
		rec_pay = SRT_RECEIVER;
	else
	if (call_put == SRT_PUT)
		rec_pay = SRT_PAYER;
	if (strLogorNormal[0]!= '\0')
		err = interp_diffusion_type(strLogorNormal, &log_or_norm);
	else
		log_or_norm = SRT_LOGNORMAL;

/* Gets the bond option implied PRICE volatility */
	option_maturity = (option_expiry - clcn_date) * YEARS_IN_DAY;
	settle_bond_df = swp_f_df(clcn_date, settle_bond_date, yc_name);
	target = bond_option_value  / settle_bond_df;
	err = srt_f_optimpvol( target, fwd_clean_price , clean_strike, option_maturity, 1.0 , 
						call_put,SRT_LOGNORMAL, &bond_price_vol);
	if (err)
		return err;
	
/* Deduces the swap PRICE volatility and the swaption price */	
	price_vol_ratio = bond_duration / swap_duration;
	swap_price_vol = bond_price_vol / price_vol_ratio;
	settle_swap_df = swp_f_df(clcn_date, settle_swap_date, yc_name);	
	swaption_value = settle_swap_df * srt_f_optblksch( fwd_clean_price/clean_strike, 1.0, swap_price_vol,
				option_maturity, 1.0, call_put, PREMIUM);

/* Gets the forwards with and without spreads (for strike correction) */
	err = swp_f_ForwardRate_SwapDP( &swap_dp, yc_name, refRateCode,&real_swap_forward);
	if (err)
		return err;
	err = swp_f_ForwardRate_SwapDP( &swap_dp, yc_name, "CASH", &cash_swap_forward);
	if (err)
		return err;

/* Make the strike correction to take the spreads into account */
	swaption_strike = *swap_coupon + real_swap_forward  - cash_swap_forward;
	*swap_coupon = swaption_strike;

/* Implies the swaption volatility that matches the swaption price*/
	err = swp_f_SwaptionImpliedVol_SwapDP( swaption_value, &swap_dp, swaption_strike,
				rec_pay, refRateCode, yc_name, log_or_norm, implied_vol);
	if (err)
		return err;

/* Computes the coefficient of hedging */
	db = bond_price_vol * bond_price_vol * option_maturity;
	db = (log(fwd_clean_price/clean_strike) + 0.5 * db) 
		/ sqrt(db);
	
	ds = bond_price_vol * bond_price_vol * option_maturity / (price_vol_ratio * price_vol_ratio);
	ds = (log(fwd_clean_price/clean_strike) + 0.5 * ds) 
		/ sqrt(ds);
 
	*hedge_ratio = bond_duration / swap_duration * clean_strike * norm(db) / norm(ds); 

/* Computes the vega ratio by shifting the bond price vol by 1% */
	bond_price_vol += 0.01;
	swap_price_vol = bond_price_vol / price_vol_ratio;
	swaption_value = settle_swap_df * srt_f_optblksch( fwd_clean_price/clean_strike, 1.0, swap_price_vol,
				option_maturity, 1.0, call_put, PREMIUM);
	
	err = swp_f_SwaptionImpliedVol_SwapDP( swaption_value, &swap_dp, *swap_coupon,
				rec_pay, refRateCode, yc_name, log_or_norm, vega_ratio);
	
	if (err)
		return err;
	*vega_ratio = (*vega_ratio - *implied_vol)/0.01;

	/* Return a success message */
	return NULL;
}


#undef    FORWARD 
#undef    BACKWARD
#undef    NO_DIRECTION