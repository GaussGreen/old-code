			/* SRT_H_VEGATRECHE2DR.h */

#include "srt_h_all.h"

#ifndef VEGATRECHE2DR_H
#define VEGATRECHE2DR_H

/* -------------------------------------------------------------------------- */

SrtStpPtr gotoindex( SrtStpPtr step, long index);

/* -> move along the step to fit index 'index */

/* -------------------------------------------------------------------------- */

void my_print_Node(FILE* out, long index, int r_index, int phi_index, 
						SrtTrinTreNdInf* node);

/* -> print a node */

/* -------------------------------------------------------------------------- */

/* double**** f4tensor(	long m_min, long m_max,
			long n_min, long n_max,
			long o_min, long o_max,
			long p_min, long p_max	);
*/
/* -> create a 4d array of double with the following dimensions:
    [m_min..m_max]
    [n_min..n_max]
    [o_min..o_max]
    [p_min..p_max]
In case of error, NULL is returned but every allocated and not used for the 
f4tensir has been properly freed 
*/

/* -------------------------------------------------------------------------- */
/*
void free_f4tensor( double**** t,
		    long m_min, long m_max,
		    long n_min, long n_max,
		    long o_min, long o_max,
		    long p_min, long p_max	);
*/
/* -> free memory for the 4d array t */ 

/* -------------------------------------------------------------------------- */

SRT_Boolean srt_f_treche_requestOK(int type);

/* -> return SRT_YES if the tree function using Cheyette model can compute
the request(s) of type 'type. */


/* -------------------------------------------------------------------------- */

Err srt_f_vegatreche2dr
	(
	SrtUndPtr und, /* market parameters */
	SrtGrfnParam	*grfnparam,	 /* model parameters */
	SrtStpPtr stp, /* discretization of deal in time, wi/ events attached*/
	GrfnDeal	   *gd,/* deal descriptor structure */
	EvalEventFct evalcf, /* cashflow evaluator */
	SrtIOStruct* iolist, /* list of price requests */
        SrtUndInfo *und_info
	);

/* -> fill,the request list with the shifted prices it can compute,
using a trinomial tree that fits the Cheyette model. */

#endif
