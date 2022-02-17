#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/*Translated by FOR_C++, v0.1, on 06/08/90 at 16:03:25 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 06/08/90 at 16:03:21
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    Purpose:    Call E1MES if the inequality IARG1 <= IARG2 is violated.
                Increment the error code.

    Usage:      CALL C12ILE (IARG1, NMARG1, IARG2, NMARG2, NER)

    Arguments:
       IARG1  - Smaller argument.  (Input)
       NMARG1 - Character string that contains the name of the smaller
                argument.  (Input)
       IARG2  - Larger argument.  (Input)
       NMARG2 - Character string that contains the name of the larger
                argument.  (Input)
       NER    - Error code.  (Input/Output)
                The input NER is the error code used by E1MES.  NER is
                incremented by 1 on output.

    Copyright:  1984 by IMSL, Inc.  All rights reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied to
                this code.  No other warranty, expressed or implied, is
                applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_c12ile(Mint iarg1, Mchar *nmarg1, Mint iarg2, Mchar *nmarg2, Mint *ner)
#else
void imsl_c12ile(iarg1, nmarg1, iarg2, nmarg2, ner)
	Mint             iarg1;
	Mchar           *nmarg1;
	Mint             iarg2;
	Mchar           *nmarg2;
	Mint            *ner;
#endif
{


	if (iarg1 > iarg2) {
		imsl_e1sti(1, iarg1);
		imsl_e1sti(2, iarg2);
		imsl_e1stl(1, nmarg1);
		imsl_e1stl(2, nmarg2);

/*		imsl_ermes(5, *ner, "%(l1) = %(i1) and %(l2) = %(i2).  %(l1) must be less than or equal to %(l2).");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_INEQUALITY_VIOLATION);
	}
	*ner += 1;
	return;
}


/*-----------------------------------------------------------------------
    Purpose:    Call E1MES if an integer argument fails to fall in
                a specified interval.

    Usage:      CALL C1IARG (IARG, NMARG, INT1, INT2, NER)

    Arguments:
       IARG   - Argument to be checked.  (Input)
       NMARG  - Character string that contains the name of the
                argument to be checked.
       INT1   - Lower bound for IARG.  (Input)
       INT2   - Upper bound for IARG.  (Input)
                If INT2 is less than INT1, no upper bound is assumed.
       NER    - Error code.  (Input/Output)
                The input NER is the error code used by E1MES.  NER is
                incremented by 1 on output.
  -----------------------------------------------------------------------
 */
#ifdef ANSI
void imsl_c1iarg(Mint iarg, Mchar *nmarg, Mint int1, Mint int2, Mint *ner)
#else
void imsl_c1iarg(iarg, nmarg, int1, int2, ner)
	Mint             iarg;
	Mchar           *nmarg;
	Mint             int1, int2, *ner;
#endif
{


	if (int1 <= int2) {
		if (iarg < int1 || iarg > int2) {
			imsl_e1sti(1, iarg);
			imsl_e1sti(2, int1);
			imsl_e1sti(3, int2);
			imsl_e1stl(1, nmarg);

/*			imsl_ermes(5, *ner, "%(l1) = %(i1).  %(l1) must be greater than or equal to %(i2) and less than or equal to %(i3).");
*/
                        imsl_ermes(IMSL_TERMINAL,
			IMSL_INEQUALITY_VIOLATION_2);
		}
	} else {
		if (iarg < int1) {
			imsl_e1sti(1, iarg);
			imsl_e1sti(2, int1);
			imsl_e1stl(1, nmarg);

/*			imsl_ermes(5, *ner, "%(l1) = %(i1).  %(l1) must be greater than or equal to %(i2).");
*/
                        imsl_ermes(IMSL_TERMINAL,
			IMSL_INEQUALITY_VIOLATION_3);
		}
	}
	*ner += 1;
	return;
}
 
