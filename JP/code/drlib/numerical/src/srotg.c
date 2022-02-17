#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/* -----------------------------------------------------------------------
    IMSL Name:  SROTG (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    August 9, 1986

    Purpose:    Construct a Givens plane rotation in single precision.

    Usage:      CALL SROTG (SA, SB, SC, SS)

    Arguments:
       SA     - First element of vector.  (Input/Output)
                On output, R = (+/-)SQRT(SA**2 + SB**2) overwrites SA.
       SB     - Second element of vector.  (Input/Output)
                On output, Z overwrites SB where Z is defined to be
                  SS        if ABS(SA) .GT. ABS(SB)
                  1.0/SC    if ABS(SB) .GE. ABS(SA) and SC .NE. 0.0
                  1.0       if SC .EQ. 0.0.
       SC     - Real scalar containing elements of the rotation matrix.
                (Output)
       SS     - Real scalar containing elements of the rotation matrix.
                (Output)

    Remark:
       SROTG constructs the Givens rotation
                           ( SC  SS )
                       G = (        ) , SC**2 + SS**2 = 1
                           (-SS  SC )
       which zeros the second element of (SA SB)(transpose).

    Keywords:   Level 1 BLAS; SROTG

    GAMS:       D1a8

    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Mathematical Support

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_srotg(Mfloat *sa, Mfloat *sb, Mfloat *sc, Mfloat *ss)
#else
void imsl_srotg(sa, sb, sc, ss)
	Mfloat          *sa, *sb, *sc, *ss;
#endif
{
	Mfloat           r, u, v;


	if (fabs(*sa) > fabs(*sb)) {
		/* HERE ABS(SA) .GT. ABS(SB) */
		u = *sa + *sa;
		v = *sb / u;
		/*
		 * NOTE THAT U AND R HAVE THE SIGN OF SA
		 */
		r = sqrt(.25 + imsl_fi_power(v, 2)) * u;
		/* NOTE THAT SC IS POSITIVE */
		*sc = *sa / r;
		*ss = v * (*sc + *sc);
		*sb = *ss;
		*sa = r;
	} else {
		/* HERE ABS(SA) .LE. ABS(SB) */
		if (*sb != F_ZERO) {
			u = *sb + *sb;
			v = *sa / u;
			/*
			 * NOTE THAT U AND R HAVE THE SIGN OF SB (R IS
			 * IMMEDIATELY STORED IN SA)
			 */
			*sa = sqrt(.25 + imsl_fi_power(v, 2)) * u;
			/* NOTE THAT SS IS POSITIVE */
			*ss = *sb / *sa;
			*sc = v * (*ss + *ss);
			if (*sc != F_ZERO) {
				*sb = F_ONE / *sc;
			} else {
				*sb = F_ONE;
			}
			/* HERE SA = SB = 0. */
		} else {
			*sc = F_ONE;
			*ss = F_ZERO;
			*sa = F_ZERO;
			*sb = F_ZERO;
		}
	}

	return;
}				/* end of function */
