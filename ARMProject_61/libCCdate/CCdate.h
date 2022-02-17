/* -------------------------------------------------------------------------

   File: CCdate.h
   Path: /home/julia/projets/dev/Common/libdate/SCCS/s.CCdate.h
   Description: bib de date
   Created: 96/10/06
   Author: Jacques WERNERT
   Modified: 00/03/20 15:24:01
   Last maintained by: Jacques WERNERT
   Revision: 1.17

   -------------------------------------------------------------------------

   Note:

   ------------------------------------------------------------------------- */

#ifndef _CCdate_h
# define _CCdate_h

#include "CCcommon.h"
SCCS_ID (CCdate_h_SccsId, "@(#)CCdate.h	1.17, modified 00/03/20");

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#ifdef Solaris
#include <pthread.h>
#endif

#define DAT_now                          (time (NULL))
#define DAT_mins_since_midnight_to_gmt(mins, today) DAT_secs_since_midnight_to_gmt(60*(mins), today)
#define DAT_gmt_to_mins_since_midnight(time_t) (DAT_DAT_secs_since_midnight_to_gmt(time_t) / 60)
#define DAT_secs_value(time)             (time % 60)
#define DAT_mins_from_secs(secs)         (secs / 60)

#define DAT_error                        (time_t)-1

#define DAT_STRING_SIZE                  (26)
#define DAT_FSTRING_SIZE                 (21)
#define DAT_FSTRING_FORMAT               "%h %d %Y %I:%M %p"

CCEXTERN_FUNCTION (time_t DAT_gmt_to_secs_since_midnight,
		 (time_t timestamp));
CCEXTERN_FUNCTION (time_t DAT_secs_since_midnight_to_gmt,
		 (time_t secs, time_t today));

CCEXTERN_FUNCTION (char *DAT_gmt_to_string, (time_t timestamp));
CCEXTERN_FUNCTION (char *DAT_gmt_to_string_r, (time_t timestamp));
CCEXTERN_FUNCTION (char *DAT_gmt_to_mystring, 
		 (time_t timestamp, char *buf, int bufsize));
CCEXTERN_FUNCTION (char *DAT_gmt_to_syb_string, (time_t timestamp));
CCEXTERN_FUNCTION (char *DAT_gmt_to_syb_string_r, (time_t timestamp));
CCEXTERN_FUNCTION (char *DAT_gmt_to_XLstring, (time_t timestamp));

CCEXTERN_FUNCTION (time_t DAT_local_to_gmt, (time_t local));
CCEXTERN_FUNCTION (time_t DAT_gmt_to_local, (time_t gmt));

CCEXTERN_FUNCTION (time_t DAT_ssdate_to_gmt, (double ssdate));
CCEXTERN_FUNCTION (double DAT_gmt_to_ssdate, (time_t timestamp));

CCEXTERN_FUNCTION (int DAT_is_week_end, (time_t timestamp));

CCEXTERN_FUNCTION (time_t DAT_struct_to_gmt, (int year, int month, int day, 
					    int hour, int minute, int second));
CCEXTERN_FUNCTION (int DAT_gmt_to_struct, (time_t gmt, 
					 int *year, int *month, int *day, 
					 int *hour, int *minute, int *second));

CCEXTERN_FUNCTION (time_t DAT_fstring_to_gmt, (char *source, char *format));
CCEXTERN_FUNCTION (char *DAT_gmt_to_fstring, (time_t gmt, char *format));

CCEXTERN_FUNCTION (time_t DAT_gmt_to_secs_until_midnight, (void));

CCEXTERN_FUNCTION (time_t DAT_adgdate_to_gmt, (int adgdate));
CCEXTERN_FUNCTION (int DAT_gmt_to_adgdate, (time_t gmt));

#define DAT_ymd_to_ssdate(y, m, d) DAT_gmt_to_ssdate (DAT_struct_to_gmt ((y), (m), (d), 0, 0, 0))

#define DAT_gmt_to_XLdate(gmt)   	DAT_gmt_to_ssdate (gmt)
#define DAT_gmt_to_AXdate(gmt)   	DAT_gmt_to_ssdate (gmt)
#define DAT_gmt_to_WZdate(gmt)   	DAT_gmt_to_ssdate (gmt)
#define DAT_XLdate_to_gmt(excel) 	DAT_ssdate_to_gmt (excel)
#define DAT_AXdate_to_gmt(applix) 	DAT_ssdate_to_gmt (applix)
#define DAT_WZdate_to_gmt(wingz) 	DAT_ssdate_to_gmt (wingz)

CCEXTERN_FUNCTION (char *DAT_gmt_to_fstring_r, (time_t gmt, char *format));
CCEXTERN_FUNCTION (time_t DAT_fstring_to_gmt_r, (char *source, char *format));

CCEXTERN_FUNCTION (long DAT_DMYtoJulian, (int d, int m, int y, double* Julian));
CCEXTERN_FUNCTION (long DAT_JulianToDMY, (double Julian, int* d, int* m, int* y));

CCEXTERN_FUNCTION (long DAT_ssdate_to_struct, (double ssdate, int* y, int* m, int* d));
CCEXTERN_FUNCTION (double DAT_struct_to_ssdate, (int y, int m, int d));

CCEXTERN_FUNCTION (double XLDateToJulian, (double d ) );
CCEXTERN_FUNCTION (double JulianToXLDate, (double d ) );
CCEXTERN_FUNCTION (int isValidXLDate, (double d ) );


/*
 * anciens noms
 */
#define DAT_secs_since_midnight2gmt(secs, today) DAT_secs_since_midnight_to_gmt(secs, today)

#endif _CCdate_h

/* EOF CCdate.h */

