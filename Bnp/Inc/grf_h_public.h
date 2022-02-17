#ifndef GRF_H_PUBLIC_H
#define GRF_H_PUBLIC_H


/*===============================================================*/
/*===== public grfn types  ======================================*/
/*===============================================================*/


#include "grf_h_pubtypes.h"

/*===============================================================*/
/*===== printing and debugging  =================================*/
/*===============================================================*/


#include "grf_h_fltfl.h"


/*===============================================================*/
/*=============top level grfn function ==========================*/
/*===============================================================*/


Err grf_f_prgrdl(String		filename,
				 int		numeventdates,
				 Date		*eventdates,
				 long		nrows,
				 long		ncols,
				 GrfnCell	**sprdsht,
				 long		numgrng,
				 GrfnRng	*grng,
				 long		auxwidth,
				 long		*auxlen,
				 double		**aux);


/* sets *answer to be equal to d, prints out ref rate name,
 for testing historical call back function */
Err grf_f_testhistfct(Date d, char *ref_rate_name, double *answer);

#endif
