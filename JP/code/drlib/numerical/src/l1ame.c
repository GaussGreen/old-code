#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/* Structured by FOR_STRUCT, v0.2, on 07/30/90 at 11:42:40
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  L1AME (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    August 3, 1987

    Purpose:    Compares CA to CB regardless of case.

    Usage:      L1AME (CA, CB)

    Arguments:
       CA     - Character*1 scalar.  (Input)
       CB     - Character*1 scalar.  (Input)
                On entry, CA and CB specify characters to be compared.

    GAMS:       R1

    Chapters:   MATH/LIBRARY Reference Material

    Copyright:  1987 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
/* CA_S and CB_S are not used, but leave the calling sequence intact. */
#ifdef ANSI
Mint imsl_l1ame(Mchar *ca, Mint ca_s, Mchar *cb, Mint cb_s)
#else
Mint imsl_l1ame(ca, ca_s, cb, cb_s)
	Mchar           *ca;
	Mint            ca_s;
	Mchar           *cb;
	Mint            cb_s;
#endif
{
	/*
	 * Test if the characters are equal
	 */

	return (tolower(*ca) == tolower(*cb));
}
