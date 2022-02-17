#ifndef CX_ALIBCONV_H
#define CX_ALIBCONV_H

#include "cxutils.h"

#include <alib/cdate.h>

/*f
 * Converts a string representation ("2Y", "3M", "1S", etc) to a TDateInterval
 */
int CxStringToDateInterval
(char          *input,      /* (I) String w/ 1A, 3M, 4D, etc */
 TDateInterval *interval);  /* (O) Value read from file */

/*f
 * Formats a TDateInterval. Can be called twice from the same print
 * statement, but not more.
 */
char* CxFormatDateInterval
(TDateInterval interval);  /* (I) */

/*
** calendar coercion uses the string version
*/
CxTCalendar* CxCalendarFromAlibHolidayName(char*);

/*f
***************************************************************************
** Converts a TCurve into an ALIB custom curve.
**
** The resulting TCurve will have fBaseDate = baseDate.
**
** This TCurve will have the correct interpolation if you use it using
** standard ALIB functions, but will not necessarily have good values in
** the dates and rates parameters.
**
** This curve will only be good for a single session since it does not
** have serialisation built in at present.
***************************************************************************
*/
TCurve* CxZeroCurveToAlibCustomCurve (TCurve *in);

#endif
