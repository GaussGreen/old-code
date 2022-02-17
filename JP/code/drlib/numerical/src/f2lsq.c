

/*   TO GET THIS PROGRAM TO COMPILE AND LINK IT WAS NECESSARY FOR
     ME TO PULL IN A COPY OF FLSFIT.C . (I NEEDED TO LINK WITH 
     IMSL_SROTM & IMSL_SROTMG)  ALL WARNINGS ISSUED BY Gnu-C
     APPEAR TO COME FROM FLSFIT.C !
*/

#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/* Structured by FOR_STRUCT, v0.2, on 08/24/90 at 15:02:32
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  F2LSQ/DF2LSQ (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 17, 1986

    Purpose:    Compute a least squares approximation with user-supplied
                functions.

    Usage:      CALL F2LSQ (F, INTCEP, NBASIS, NDATA, XDATA, FDATA, IWT,
                            WEIGHT, A, SSE, WK)

    Arguments:  (See FNLSQ)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_f2lsq(Mfloat (*f) (Mint, Mfloat), Mint *intcep, Mint *nbasis, Mint *ndata,
                    Mfloat xdata[], Mfloat fdata[], Mint *iwt,
	            Mfloat weight[], Mfloat a[], Mfloat *sse, Mfloat wk[])
#else
void imsl_f2lsq(f, intcep, nbasis, ndata, xdata, fdata, iwt,
        	   weight, a, sse, wk)
	Mfloat           (*f) ();
	Mint            *intcep, *nbasis, *ndata;
	Mfloat           xdata[], fdata[];
	Mint            *iwt;
	Mfloat           weight[], a[], *sse, wk[];
#endif
{
	Mint             _l0, _l1, i, id, idep, ido, idummy[1], ifrq, iind,
	                imax, imin, iprnt, ir, irank, istop, imsl_isub,
	                iwt2, ix, ixx, j, ldb, ldr, ldscpe, ldx, ncol,
	                nrmiss, nrow;
	Mfloat           dfe, scpe[1], tol;


	imsl_e1psh("imsl_f2lsq");
	scpe[0] = 0.0;

	/* Check INTCEP */
	if (*intcep < 0 || *intcep > 1) {
		imsl_e1sti(1, *intcep);

/*		imsl_ermes(5, 1, "The intercept option must be a zero or a one while INTCEP = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INTCEP_SHOULD_BE_ZERO);
	}
	/* Check NBASIS */
	if (*nbasis < 1) {
		imsl_e1sti(1, *nbasis);

/*		imsl_ermes(5, 2, "The number of basis functions must be at least one while NBASIS = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NBASIS_FCNS_TOO_SMALL);
	}
	/* Check NDATA */
	if (*ndata < 1) {
		imsl_e1sti(1, *ndata);

/*		imsl_ermes(5, 3, "The number of data points must be at least one while NDATA = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_AT_LEAST_1_PT);
	}
	/* Check IWT */
	if (*iwt < 0 || *iwt > 1) {
		imsl_e1sti(1, *iwt);

/*		imsl_ermes(5, 4, "The weighting option must be a zero or a one while IWT = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_BAD_IWT_OPTION);
	}
	if (imsl_n1rty(0) != 0)
		goto L_9000;

	nrow = 1;
	ncol = *nbasis + 1 + *iwt;
	ldx = 1;
	iind = -*nbasis;
	idep = -1;
	ifrq = 0;
	if (*iwt == 0) {
		iwt2 = 0;
	} else {
		iwt2 = *nbasis + 1;
	}
	imsl_isub = *intcep;
	tol = 100.0 * imsl_amach(4);
	ldb = *intcep + *nbasis;
	ldr = *intcep + *nbasis;
	ldscpe = 1;
	/* Partition workspace */
	ir = 1;
	ix = ir + imsl_ii_power(*nbasis + *intcep, 2);
	id = ix + *intcep + ncol;
	imin = id + *nbasis + *intcep;
	imax = imin + *nbasis + *intcep;

	for (i = 1; i <= *ndata; i++) {
		if (i == 1) {
			ido = 1;
		} else if (i == *ndata) {
			ido = 3;
		} else {
			ido = 2;
		}
		ixx = ix + *intcep;
		for (j = 1; j <= *nbasis; j++) {
			imsl_e1usr("ON");
			wk[ixx - 1] = (*f) (j, xdata[i - 1]);
			imsl_e1usr("OFF");
			ixx += 1;
		}
		if (*iwt == 1) {
			wk[ixx - 1] = weight[i - 1];
			ixx += 1;
		}
		wk[ixx - 1] = fdata[i - 1];
		/*
		 * Turn off printing and stopping of type 6 errors. First
		 * retreive current settings.
		 */
		imsl_e1pos(-6, &iprnt, &istop);
		/* Then set new values */
                _l0 = 0; _l1 = 0;
		imsl_e1pos(6, &_l0, &_l1);


		imsl_r2ivn(ido, nrow, ncol, &wk[ix + *intcep - 1], ldx, *intcep,
		iind, idummy, idep, idummy, ifrq, iwt2, imsl_isub, tol,
		a, ldb, &wk[ir - 1], ldr, &wk[id - 1], &irank, &dfe, scpe,
		ldscpe, &nrmiss, &wk[imin - 1], &wk[imax - 1], &wk[ix - 1]);
		/* Reset old values */
		imsl_e1pos(6, &iprnt, &istop);
		if (imsl_n1rty(0) == 4)
			goto L_30;
	}
	/* Clear warnings */
	if (imsl_n1rty(1) == 3 || imsl_n1rty(1) == 6)
		imsl_e1mes(0, 0, " ");
	if (irank < *intcep + *nbasis) {
		if (*intcep == 0) {

/*			imsl_ermes(3, 1, "Linear dependence of the basis functions was declared.  Appropriate elements of A are set to zero.");
*/
                imsl_ermes(IMSL_WARNING, IMSL_LINEAR_DEPENDENCE);
		} else {

/*			imsl_ermes(3, 2, "Linear dependence of the constant function (the intercept) and the basis functions was declared.  Appropriate elements of A are set to zero.");
*/
                imsl_ermes(IMSL_WARNING, IMSL_LINEAR_DEPENDENCE_CONST);
		}
	}
L_30:
	if (imsl_n1rty(1) == 4 && imsl_n1rcd(1) ==  IMSL_NONNEG_WEIGHT_REQUEST_1) {

/*		imsl_ermes(4, 1, "An element of the weight vector is not positive.  All elements must be positive.");
*/
                imsl_ermes(IMSL_FATAL, IMSL_NEGATIVE_WEIGHTS_2);
	}
	/* Exit section */
L_9000:
	imsl_e1pop("imsl_f2lsq");
	*sse = scpe[0];
	return;
}				/* end of function */
