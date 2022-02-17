#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*  -----------------------------------------------------------------------
    IMSL Name:  C2DEC/DC2DEC (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 1, 1984

    Purpose:    Compute the cubic spline interpolant with specified
                derivative endpoint conditions.

    Usage:      CALL C2DEC (NDATA, XDATA, FDATA, ILEFT, DLEFT, IRIGHT,
                            DRIGHT, BREAK, CSCOEF, IPVT)

    Arguments:  (See CSDEC)

    GAMS:       E1a

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_c2dec(Mint *ndata, Mfloat xdata[], Mfloat fdata[], Mint *ileft, Mfloat *dleft, Mint *iright,
     Mfloat *dright, Mfloat break_[], Mfloat *cscoef, Mint ipvt[])
#else
void imsl_c2dec(ndata, xdata, fdata, ileft, dleft, iright,
     dright, break_, cscoef, ipvt)
	Mint            *ndata;
	Mfloat           xdata[], fdata[];
	Mint            *ileft;
	Mfloat          *dleft;
	Mint            *iright;
	Mfloat          *dright, break_[], *cscoef;
	Mint             ipvt[];
#endif
{
#define CSCOEF(I_,J_)	(cscoef+(I_)*(4)+(J_))
	Mint             i, j, l, m;
	Mfloat           divdf1, divdf3, dtau, g;


	imsl_e1psh("IMSL_C2DEC");
	/* CHECK ARGUMENT NDATA */
	if (*ndata < 2) {
		imsl_e1sti(1, *ndata);

/*		imsl_ermes(5, 1, "The number of data points must be 2 or more while NDATA = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_AT_LEAST_2_PTS);
	}
	/* CHECK ARGUMENT ILEFT */
	if (*ileft < 0 || *ileft > 2) {
		imsl_e1sti(1, *ileft);

/*		imsl_ermes(5, 2, "The argument ILEFT = %(i1).  It must be 0, 1 or 2.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_WRONG_ILEFT_VALUE);
	}
	/* CHECK ARGUMENT IRIGHT */
	if (*iright < 0 || *iright > 2) {
		imsl_e1sti(1, *iright);

/*		imsl_ermes(5, 3, "The argument IRIGHT = %(i1).  It must be 0, 1 or 2.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_WRONG_IRIGHT_VALUE);
	}
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/* COPY AND SORT INPUT DATA */
	imsl_c1sor(*ndata, xdata, fdata, break_, cscoef, 4, ipvt);
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/*
	 * A TRIDIAGONAL LINEAR SYSTEM FOR THE UNKNOWN SLOPES S(I) OF F AT
	 * BREAK(I), I=1,...,N, IS GENERATED AND THEN SOLVED BY GAUSS
	 * ELIMINATION, WITH S(I) ENDING UP IN CSCOEF(2,I), ALL I.
	 * CSCOEF(3,.) AND CSCOEF(4,.) ARE USED INITIALLY FOR TEMPORARY
	 * STORAGE.
	 */
	l = *ndata - 1;
	/*
	 * COMPUTE FIRST DIFFERENCES OF TAU SEQUENCE AND STORE IN
	 * CSCOEF(3,.). ALSO COMPUTE FIRST DIVIDED DIFFERENCE OF DATA AND
	 * STORE IN CSCOEF(4,.)
	 */
	for (m = 2; m <= *ndata; m++) {
		*CSCOEF(m - 1, 2) = break_[m - 1] - break_[m - 2];
		*CSCOEF(m - 1, 3) = (*CSCOEF(m - 1, 0) - *CSCOEF(m - 2, 0)) /
			*CSCOEF(m - 1, 2);
	}
	/*
	 * CONSTRUCT FIRST EQUATION FROM THE BOUNDARY CONDITION, OF THE FORM
	 * CSCOEF(4,1)*S(1) + CSCOEF(3,1)*S(2) = DLEFT
	 */
	if (*ileft == 0) {
		if (*ndata == 2) {
			/*
			 * NO CONDITION AT LEFT END AND NDATA = 2
			 */
			*CSCOEF(0, 3) = F_ONE;
			*CSCOEF(0, 2) = F_ONE;
			*CSCOEF(0, 1) = F_TWO ** CSCOEF(1, 3);
		} else {
			/*
			 * NOT-A-KNOT CONDITION AT LEFT END AND NDATA .GT. 2
			 */
			*CSCOEF(0, 3) = *CSCOEF(2, 2);
			*CSCOEF(0, 2) = *CSCOEF(1, 2) + *CSCOEF(2, 2);
			*CSCOEF(0, 1) = ((*CSCOEF(1, 2) + F_TWO ** CSCOEF(0, 2)) ** CSCOEF(1, 3) *
					 *CSCOEF(2, 2) + imsl_fi_power(*CSCOEF(1, 2), 2) ** CSCOEF(2, 3)) / *CSCOEF(0, 2);
		}
	} else if (*ileft == 1) {
		/* SLOPE PRESCRIBED AT LEFT END */
		*CSCOEF(0, 3) = F_ONE;
		*CSCOEF(0, 2) = F_ZERO;
		*CSCOEF(0, 1) = *dleft;
	} else if (*ileft == 2) {
		/*
		 * SECOND DERIVATIVE PRESCRIBED AT LEFT END
		 */
		*CSCOEF(0, 3) = F_TWO;
		*CSCOEF(0, 2) = F_ONE;
		*CSCOEF(0, 1) = F_THREE ** CSCOEF(1, 3) - *CSCOEF(1, 2) / F_TWO ** dleft;
	}
	if (*ndata > 2) {
		/*
		 * IF THERE ARE INTERIOR KNOTS, GENERATE THE CORRESP.
		 * EQUATIONS AND CARRY OUT THE FORWARD PASS OF GAUSS
		 * ELIMINATION, AFTER WHICH THE M-TH EQUATION READS
		 * CSCOEF(4,M)*S(M) + CSCOEF(3,M)*S(M+1) = CSCOEF(2,M)
		 */
		for (m = 2; m <= l; m++) {
			g = -*CSCOEF(m, 2) / *CSCOEF(m - 2, 3);
			*CSCOEF(m - 1, 1) = g ** CSCOEF(m - 2, 1) + F_THREE * (*CSCOEF(m - 1, 2) *
			 *CSCOEF(m, 3) + *CSCOEF(m, 2) ** CSCOEF(m - 1, 3));
			*CSCOEF(m - 1, 3) = g ** CSCOEF(m - 2, 2) + F_TWO * (*CSCOEF(m - 1, 2) +
							     *CSCOEF(m, 2));
		}
		/*
		 * CONSTRUCT LAST EQUATION FROM THE SECOND BOUNDARY
		 * CONDITION, OF THE FORM (-G*CSCOEF(4,NDATA-1))*S(NDATA
		 * A-1)+CSCOEF(5
		 */
		if (*iright == 0) {
			if (*ndata == 3 && *ileft == 0) {
				/*
				 * NDATA=3 AND NOT-A-KNOT ALSO AT LEFT
				 */
				*CSCOEF(*ndata - 1, 1) = F_TWO ** CSCOEF(*ndata - 1, 3);
				*CSCOEF(*ndata - 1, 3) = F_ONE;
				g = -F_ONE / *CSCOEF(*ndata - 2, 3);
			}
		} else if (*iright == 1) {
			/*
			 * SLOPE IS PRESCRIBED AT RIGHT END, GO DIRECTLY TO
			 * BACKSUBSTITUTION, SINCE C ARRAY HAPPENS TO BE SET
			 * UP JUST RIGHT FOR IT AT THIS POINT
			 */
			*CSCOEF(*ndata - 1, 1) = *dright;
		} else if (*iright == 2) {
			/*
			 * SECOND DERIVATIVE PRESCRIBED AT RIGHT
			 */
			*CSCOEF(*ndata - 1, 1) = F_THREE ** CSCOEF(*ndata - 1, 3) + *CSCOEF(*ndata - 1, 2) /
				F_TWO ** dright;
			*CSCOEF(*ndata - 1, 3) = F_TWO;
			g = -F_ONE / *CSCOEF(*ndata - 2, 3);
		}
		if (*iright == 0 && !(*ndata == 3 && *ileft == 0)) {
			g = *CSCOEF(*ndata - 2, 2) + *CSCOEF(*ndata - 1, 2);
			*CSCOEF(*ndata - 1, 1) = ((*CSCOEF(*ndata - 1, 2) + F_TWO *
						   g) ** CSCOEF(*ndata - 1, 3) ** CSCOEF(*ndata - 2, 2) + imsl_fi_power(*CSCOEF(*ndata - 1, 2), 2) *
						  (*CSCOEF(*ndata - 2, 0) - *CSCOEF(*ndata - 3, 0)) / *CSCOEF(*ndata - 2, 2)) /
				g;
			g = -g / *CSCOEF(*ndata - 2, 3);
			*CSCOEF(*ndata - 1, 3) = *CSCOEF(*ndata - 2, 2);
		}
	} else {
		/* NDATA = 2 */
		if (*iright == 0) {
			if (*ileft > 0) {
				/*
				 * NDATA=2 AND NOT NOT-A-KNOT AT LEFT END
				 */
				*CSCOEF(*ndata - 1, 1) = F_TWO ** CSCOEF(*ndata - 1, 3);
				*CSCOEF(*ndata - 1, 3) = F_ONE;
			} else {
				/*
				 * NDATA=2 AND NOT-A-KNOT AT BOTH ENDPOINTS
				 */
				*CSCOEF(*ndata - 1, 1) = *CSCOEF(*ndata - 1, 3);
			}
		} else if (*iright == 1) {
			*CSCOEF(*ndata - 1, 1) = *dright;
		} else if (*iright == 2) {
			/*
			 * SECOND DERIVATIVE PRESCRIBED AT RIGHT
			 */
			*CSCOEF(*ndata - 1, 1) = F_THREE ** CSCOEF(*ndata - 1, 3) + *CSCOEF(*ndata - 1, 2) /
				F_TWO ** dright;
			*CSCOEF(*ndata - 1, 3) = F_TWO;
		}
		g = -F_ONE / *CSCOEF(*ndata - 2, 3);
	}

	if (*iright != 1 && ((*ndata > 2 || *iright != 0) || *ileft !=
			     0)) {
		/*
		 * COMPLETE FORWARD PASS OF GAUSS ELIMINATION
		 */
		*CSCOEF(*ndata - 1, 3) += g ** CSCOEF(*ndata - 2, 2);
		*CSCOEF(*ndata - 1, 1) = (g ** CSCOEF(*ndata - 2, 1) + *CSCOEF(*ndata - 1, 1)) /
			*CSCOEF(*ndata - 1, 3);
	}
	/* CARRY OUT BACK SUBSTITUTION */
	for (j = l; j >= 1; j--) {
		*CSCOEF(j - 1, 1) = (*CSCOEF(j - 1, 1) - *CSCOEF(j - 1, 2) ** CSCOEF(j, 1)) /
			*CSCOEF(j - 1, 3);
	}
	/*
	 * GENERATE CUBIC COEFFICIENTS IN EACH INTERVAL, I.E., THE DERIV.S AT
	 * ITS LEFT ENDPOINT, FROM VALUE AND SLOPE AT ITS ENDPOINTS.
	 */
	for (i = 2; i <= *ndata; i++) {
		dtau = *CSCOEF(i - 1, 2);
		divdf1 = (*CSCOEF(i - 1, 0) - *CSCOEF(i - 2, 0)) / dtau;
		divdf3 = *CSCOEF(i - 2, 1) + *CSCOEF(i - 1, 1) - F_TWO * divdf1;
		*CSCOEF(i - 2, 2) = F_TWO * (divdf1 - *CSCOEF(i - 2, 1) - divdf3) /
			dtau;
		*CSCOEF(i - 2, 3) = (divdf3 / dtau) * (F_SIX / dtau);
	}

L_9000:
	;
	imsl_e1pop("IMSL_C2DEC");
	return;
}				/* end of function */
