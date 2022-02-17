/* ===================================================================================
   FILENAME:      srt_h_lgmclsdfrm.h
   
   PURPOSE:       Compute swaption, bond option,caps/floors prices in LGM, via a 
                  closed form analytical solution.
   =================================================================================== */

#ifndef SRT_H_LGMCLSDFRM_H
#define SRT_H_LGMCLSDFRM_H


double srt_f_lgm_capfloor(
		TermStruct        *ts,
		double            *fixing_time,
		double            *period_time,   /* period_time[0] is the strike payment time */
		double            *df, 
		double            *payment,    /* This is cvg * ( cash_fwd + spread ) */
		int               num_caplets, 
		SrtReceiverType   rec_pay,
		SrtMdlDim         mdl_dim);

/* ---------------------------------------------------------------------------- */

double srt_f_newlgm_capfloor(
		TermStruct        *ts,
		double            *fixing_time,
		double            *period_time,   /* period_time[0] is the strike payment time */
		double            *df, 
		double            *payment,    /* This is cvg * ( cash_fwd + spread ) */
		int               num_caplets, 
		SrtReceiverType   rec_pay,
		SrtMdlDim         mdl_dim);

/* ---------------------------------------------------------------------------- */

double srt_f_lgm_coupon_bond_option(
			int               n,
			double            bond_strike,
			TermStruct        *ts,
			double            fixing_time,
			double            *coupon, 
			double            *pay_time, 
			double            *df, 
			SrtReceiverType   pay_rec,
			SrtMdlDim         mdl_dim);

/* ------------------------------------------------------------------------------ */	

double srt_f_newlgm_coupon_bond_option(
			int               n,
			double            bond_strike,
			TermStruct        *ts,
			double            fixing_time,
			double            *coupon, 
			double            *pay_time, 
			double            *df, 
			SrtReceiverType   pay_rec,
			SrtMdlDim         mdl_dim);

/* ------------------------------------------------------------------------------ */	

#endif

