


























#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*  -----------------------------------------------------------------------
    IMSL Name:  C2INT/DC2INT (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 12, 1984

    Purpose:    Compute the cubic spline interpolant with the
                'not-a-knot' condition.

    Usage:      CALL C2INT (NDATA, XDATA, FDATA, BREAK, CSCOEF, IPVT)

    Arguments:
       NDATA  - Number of data points.  (Input)
                NDATA must be at least 2.
       XDATA  - Array of length NDATA containing the data point
                abscissas.  (Input)
                The data point abscissas must be distinct.
       FDATA  - Array of length NDATA containing the data point
                ordinates.  (Input)
       BREAK  - Array of length NDATA containing the breakpoints
                for the piecewise cubic representation.  (Output)
       CSCOEF - Matrix of size 4 by NDATA containing the local
                coefficients of the cubic pieces.  (Output)
       IPVT   - Work array of length NDATA.

    Remark:
       Note that column NDATA of CSCOEF is used as workspace.

    GAMS:       E1a

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_c2int(Mint *ndata, Mfloat xdata[], Mfloat fdata[], Mfloat break_[], Mfloat *cscoef, 
           Mint ipvt[])
#else
void imsl_c2int(ndata, xdata, fdata, break_, cscoef, ipvt)
	Mint            *ndata;
	Mfloat           xdata[], fdata[], break_[], *cscoef;
	Mint             ipvt[];
#endif
{
#define CSCOEF(I_,J_)	(cscoef+(I_)*(4)+(J_))
        Mint            _l0, _l1;
        Mfloat          _l2,_l3;


	imsl_e1psh("IMSL_C2INT ");
	/* CALL THE COMPUTATIONAL ROUTINE. */
         _l0 = 0;
         _l2 = F_ZERO;
         _l1 = 0;
         _l3 = F_ZERO;
	imsl_c2dec(ndata, xdata, fdata, &_l0, &_l2, &_l1, &_l3, break_, cscoef,
		   ipvt);

	imsl_e1pop("IMSL_C2INT ");
	return;
}				/* end of function */
