#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*
  -----------------------------------------------------------------------
    IMSL Name:  C1DIM (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    August 16, 1985

    Purpose:    Call E1MES if any of the following inequalities are
                violated:  INT1 <= IARG1 <= IARG2 ; IARG2 > = 1.
                Increment the error code.

    Usage:      CALL C1DIM (INT1, IARG1, NMARG1, IARG2, NMARG2, NER)

    Arguments:
       INT1   - Lower bound for IARG1.  (Input)
                Usually INT1 = 1.  However, in routines that handle
                degenerate cases as not being an error, INT1 can be
                input as 0.
       IARG1  - Smaller argument.  (Input)
       NMARG1 - Character string that contains the name of the smaller
                argument.  (Input)
                Insert an asterisk (*) as the first character before
                the argument name if the check that INT1 <= IARG1 has
                already been performed.  See Remark.
       IARG2  - Larger argument.  (Input)
       NMARG2 - Character string that contains the name of the larger
                argument.  (Input)
       NER    - Error code.  (Input/Output)
                The input NER is the error code used by E1MES.  NER is
                incremented by 3 on output.

    Remark:
       C1DIM is useful for checking consistency of the logical and
       physical row dimension of a matrix.  The typical usage is
                CALL C1DIM (1, N, 'N', LDA, 'LDA', NER)
       where A is a matrix with N rows and LDA is the leading dimension
       of A specified in the calling program.  If an additional
       LD argument must also be checked that depends on N, use an
       asterisk in the second call to C1DIM.  This prevents two
       identical error messages indicating N must be greater or equal
       to INT from being issued.  For example a second call to check
       argument LDAINV is
                CALL C1DIM (1, N, '*N', LDAINV, 'LDAINV', NER)

    Chapter:    Reference Material - Error Handling (not documented)

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_c1dim(Mint int1, Mint iarg1, Mchar *nmarg1, Mint iarg2,
                Mchar *nmarg2, Mint *ner)
#else
void imsl_c1dim(int1, iarg1, nmarg1, iarg2, nmarg2, ner)
 	Mint            int1, iarg1;
	Mchar           *nmarg1;
	Mint            iarg2;
	Mchar           *nmarg2;
	Mint            *ner;
#endif
{


	if (*nmarg1 == '*'){
		*ner += 1;
		imsl_c1iarg(iarg2, nmarg2, 1, -1, ner);
		if (iarg2 >= 1) {

			imsl_c12ile(iarg1, (nmarg1 + 1), iarg2, nmarg2, ner);
		} else {
			*ner += 1;
		}
	} else {
		imsl_c1iarg(iarg1, nmarg1, int1, -1, ner);
		imsl_c1iarg(iarg2, nmarg2, 1, -1, ner);
		if (iarg2 >= 1) {

			imsl_c12ile(iarg1, nmarg1, iarg2, nmarg2, ner);
		} else {
			*ner += 1;
		}
	}
	return;
}				/* end of function */
