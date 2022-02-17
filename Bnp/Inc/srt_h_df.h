/* ===========================================================================
   
	 FILENAME:     srt_h_df.h

     PURPOSE:      The function that will have to be used everywhere to compute 
	               a discount factor for a given curve id (cash df)

   =========================================================================== */


#ifndef SRT_H_DF_H
#define SRT_H_DF_H


#ifndef SRT_DF_ERROR
#define SRT_DF_ERROR	DBL_MAX
#endif                             

/* --------------------------------------------------------------------------- 
	Note that curve_id can be either a full curve name or just a ccy 
	(the default cash curve attached to this ccy will then be picked)
   --------------------------------------------------------------------------- */

double 	srt_f_df( 
			Ddate       start, 
			Ddate       end,
			String      curve_id);


#endif