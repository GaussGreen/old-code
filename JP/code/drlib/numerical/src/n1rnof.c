#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/* Structured by FOR_STRUCT, v0.2, on 01/12/90 at 18:09:54
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  N1RNOF

    Computer:   FORC/SINGLE

    Revised:    August 2, 1985

    Purpose:    Set or retrieve the error checksum number-of-failures
                flag.

    Usage:      N1RNOF (IOPT)

    Arguments:
       IOPT   - Integer specifying the desired option.  (Input)
                If IOPT=1 the number-of-failures flag is increased by 1.
                If IOPT=2 the number-of-failures flag value is returned
                          in N1RNOF and the flag is set to zero.
       N1RNOF - Integer function. (Output)  The number-of-failures flag
                value is returned in N1RNOF.

    Copyright:  1985 by IMSL, Inc.  All rights reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */

#ifdef ANSI
Mint imsl_n1rnof(Mint iopt)
#else
Mint imsl_n1rnof(iopt)
    Mint	     iopt;
#endif
{
    Mint   n1rnof_v;
    static Mint	nof = 0;

    if (iopt == 1) {
	/* Increment no.-of-failures flag */
	nof++;
	n1rnof_v = nof;
    } else if (iopt == 2) {
	/* Retrieve no.-of-failures flag */
	n1rnof_v = nof;
	/* Clear no.-of-failures flag */
	nof = 0;
    } else
/*	(8, 1, "The argument passed to N1RNOF must be 1 or 2.");
*/
        imsl_ermes(8, IMSL_BAD_ARGUMENT_TO_N1RNOF);
    return (n1rnof_v);
}
