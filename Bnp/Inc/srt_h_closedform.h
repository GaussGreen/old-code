/* ===================================================================================
   FILENAME:      srt_h_closedform.h
   
   PURPOSE:       Compute swaption, bond option,caps/floors prices:
						- build the cash flow schedule
						- call the relevant pricing function (LGM or EtaBeta)
   =================================================================================== */

#ifndef SRT_H_CLOSEDFORM_H
#define SRT_H_CLOSEDFORM_H


Err srt_f_closed_form(    
		SrtUndPtr        und,
		SrtMdlType       mdl_type,
		SrtMdlDim        mdl_dim, 
		SwapDP           *sdp,
		double 		     strike, 
		double 		     bond_strike,
		SrtReceiverType  rec_pay,
		StructType 	     type,
		String           ref_rate_code,
		double 		     *answer);


#endif

