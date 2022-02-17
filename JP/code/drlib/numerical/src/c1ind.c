#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/* 
  -----------------------------------------------------------------------
    IMSL Name:  C1IND (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    August 16, 1985

    Purpose:    Call E1MES if either of the following is false:
                (1) INT1 .LE. IND, or (2) IND.LE.NVAR .AND. NVAR.GT.0.

    Usage:      CALL C1IND (INT1, IND, NMIND, NVAR, NMNVAR, NER)

    Arguments:
       INT1   - Lower bound for IND.  (Input)
                INT1 must be greater than or equal to zero.
       IND    - Index of variable.  (Input)
       NMIND  - Character string that contains the name of the index.
                (Input)
       NVAR   - Number of variables.  (Input)
       NMNVAR - Character string that contains the name of the number
                of variables.  (Input)
       NER    - Error code.  (Input/Output)
                The input NER is the error code used by E1MES.  NER is
                incremented by 1 on output.

    Remark:
       C1IND is useful for checking for consistency between an index
       and its range.  For example, in STAT routines IWT is an index
       which gives the column number in X where weights are found.  If
       IWT is zero no weights are given.  The usage of C1IND would be
                CALL C1IND (0, IWT, 'IWT', NVAR, 'NVAR', NER)

    Chapter:    STAT/LIBRARY Utilities (not documented)

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_c1ind(Mint int1, Mint ind, Mchar *nmind, Mint nvar, Mchar *nmnvar,
                Mint *ner)
#else
void imsl_c1ind(int1, ind, nmind, nvar, nmnvar, ner)
	Mint             int1, ind;
	Mchar           *nmind;
	Mint             nvar;
	Mchar           *nmnvar;
	Mint            *ner;
#endif
{


	imsl_c1iarg(ind, nmind, int1, -1, ner);
	if (nvar >= 1) {
		imsl_c12ile(ind, nmind, nvar, nmnvar, ner);
	} else {
		*ner += 1;
	}
	return;
}				/* end of function */
