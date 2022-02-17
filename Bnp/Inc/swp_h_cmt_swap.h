/* SWP_H_CMT_SWAP.H */

#ifndef SWP_H_CMT_SWAP_H
#define SWP_H_CMT_SWAP_H


Err swp_f_CMT_swap(
				long start,
				long end_or_nfp,
				String cmt_comp_str,
				String cmt_freq_str,
				double cmt_libor_spread,
				double strike,
				double initial_exchange,
				double final_exchange,
				String info_str,
				SrtCurvePtr cmt_crv,
				double *value);


#endif
