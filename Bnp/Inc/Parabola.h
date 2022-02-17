#ifndef PARABOLA_H
#define	PARABOLA_H

#include "utError.h"
#include "utTypes.h"

/*	Parabola functions	*/
/*	JB Oct04 */


Err conmpute_parabola_replication(
									char					*str, 
									SABR_VOL_TYPE			*val);

									/*	Setup parameters	*/
									Err op_sabr_set_param(	
									int						num_param, 
									char					**param_str, 
									char					**value_str,
									SABR_RISK_PARAM			*param);
