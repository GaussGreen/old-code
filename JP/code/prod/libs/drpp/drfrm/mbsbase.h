/***************************************************************************
 *	SCCS Keyword Information
 *	------------------------
 *	Module name	:  mbsbase.h
 *	Company Name	:  JP Morgan Securities Inc.
 *	Author		:  Robert Lenk
 *			   Davis (Shuenn-Tyan) Lee
 *			   Derivatives Research
 *	Version		:  1.64
 *	Extracted	:  5/28/97 at 16:58:25
 *	Last Updated	:  5/28/97 at 16:58:22
 ***************************************************************************
 *      Base library for DR/MBS analytics (err handling, etc.)
 *      Libraries required:  GTO analytics lib
 *      Compiler directives: DEBUG_MBS
 *
 *      Copyright 1996 J.P. Morgan & Co. Incorporated. All rights reserved.
 ***************************************************************************/
#ifndef __mbsbase_h
#define __mbsbase_h

/* Define all the usual plain C hdrs here,
   to save having to do this in each module;
   exclude: malloc, ... */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <time.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Include declarations from other headers in mbs base lib */
#include "mbsutils.h"

/* Use std GTO lib macros */
#include "cgeneral.h"
#include "macros.h"
#include "cmemory.h"
#include "cerror.h"
#include "smooth.h"


/*****************************************************************************
 *  Misc constants/macros
 *****************************************************************************/
/* Std app label for error msgs */
#define DEF_MBSBASE_LIBNAME  "MBS Analytics "

/* std string lengths */
#define	SHORTSTRLEN	32	/* for brief strings */
#define	LINESTRLEN	132	/* reasonable max length for single
				   output line */
#define	MAXPATHLEN	256	/* generous max len for pathspec */
#define	ERRSTRLEN	MAXPATHLEN+2*LINESTRLEN	/* large enough for pathspec +
						   2 lines */
#define	GLOB_ERRSTRLEN	6*LINESTRLEN		/* size of global warn/error
						   string */
/* max length of single line in perl data file */
#define   MAX_PERL_LEN         16384

#define   HOLIDAY_FILE   "holiday.dat"   /* name of MBS analytics holiday file */
#define   DRDATA_SIG_FILE "drdata.sig"   /* file checked in $DRDATA */

/* Hard-coded PC DRDATA location (when DRDATA not defined) */
#define PC_DRDATADIR            "p:\\dr_data"

/* Default PC home dir */
#define   DOS_HOME_DIR    "c:"

/* Path delim */
#define   MBS_PATH_DELIM     '\\'     /* backslash for DOS paths */
#define FIMS_DATADIR_SUFFIX   '\\'

/* GTO makes it hard to access these strings, so keep local copies */
#define LOC_GTO_ACT_365_STR       "Actual/365"
#define LOC_GTO_ACT_365F_STR      "Actual/365F"
#define LOC_GTO_ACT_360_STR       "Actual/360"
#define LOC_GTO_B30_360_STR       "30/360"
#define LOC_GTO_B30E_360_STR      "30E/360"
#define LOC_GTO_ACT_365FJ_STR     "Actual/365FJ"
#define LOC_GTO_B30E_360I_STR     "30E/360I"
#define LOC_GTO_B30EP_360_STR     "30E+/360"
#define LOC_GTO_B30_360F_STR      "30/360F"


/*****************************************************************************
 *  Common inline functions
 *****************************************************************************/
#ifndef ROUND
  #define ROUND(x) ((long)((long) floor(x) + (((x)-floor(x) >= 0.5) ? 1. : 0.)))
#endif

#ifndef SQR
  #define SQR(x)      ((x)*(x))
#endif

#ifndef HIGHINT
  #define HIGHINT(x)  ((long) ceil((double) (x)))
#endif

#ifndef LOWINT
  #define LOWINT(x)   ((long) floor((double) (x)))
#endif

#ifndef MIN
  #define MIN(x,y)   (((x) < (y)) ? (x) : (y))
#endif

#ifndef MAX
  #define MAX(x,y)   (((x) > (y)) ? (x) : (y))
#endif

#ifndef COLLAR
  #define COLLAR(x,xmin,xmax) MIN(MAX((x),(xmin)),(xmax))
#endif

/* Macro to perform MBS/ARM cpn reset: rounding, life/sticky collar, etc 
 * Note: call this macro as function, but do not encase w/XONFAIL */
#ifndef MBS_CPN_RESET
#define MBS_CPN_RESET(xNewCpn,xPrevCpn,xCapSprd,xFlrSprd,              \
                      xLifeCap,xLifeFlr,xRnd,xModNewCpn)  {            \
    xModNewCpn = MAX(MIN(((fabs((xRnd))<1.e-7) ? (xNewCpn) : (xRnd)*   \
                 (floor((xNewCpn)/(xRnd)) +                            \
                  ((((xNewCpn)/(xRnd))-floor((xNewCpn)/(xRnd)) >= 0.5) ? \
                   1. : 0.))),                                         \
                MIN((xLifeCap),(xPrevCpn)+(xCapSprd))),                \
            MAX((xLifeFlr),MAX(0.,(xPrevCpn)+(xFlrSprd))));            \
}
#endif
/* Same as above, but using smoothed min/max functions */
#ifndef MBS_CPN_SMOOTHRESET
#define MBS_CPN_SMOOTHRESET(xNewCpn,xPrevCpn,xCapSprd,xFlrSprd,        \
                            xLifeCap,xLifeFlr,xRnd,xMaxDiffs,xModNewCpn) { \
    double xTmpCpn;                                                    \
    if(GtoMinSmooth(((fabs((xRnd))<1.e-7) ? (xNewCpn) : (xRnd) *       \
                     (floor((xNewCpn)/(xRnd)) +                        \
                      ((((xNewCpn)/(xRnd))-                            \
                        floor((xNewCpn)/(xRnd)) >= 0.5) ? 1. : 0.))),  \
                    xMaxDiffs, MIN((xLifeCap),(xPrevCpn)+(xCapSprd)),  \
                    &xTmpCpn) ISNT SUCCESS ||                          \
       GtoMaxSmooth(xTmpCpn, xMaxDiffs,                                \
                    MAX((xLifeFlr),MAX(0.,(xPrevCpn)+(xFlrSprd))),     \
                    &xModNewCpn) ISNT SUCCESS) {                       \
        MbsMesg(MERRTRACE,FALSE,routine);                              \
        status = FAILURE;  goto done;  }                               \
}
#endif



/*****************************************************************************
 *  Error-handling stuff
 ****************************************************************************/
/* Status codes for MbsMesg() */
#ifndef MINFO
  #define   MINFO    0
  #define   MWARN    1
  #define   MERROR   2
  #define   MERRTRACE   3
#endif

/* 
 * Error-handling macros: 
 */


/* F_INIT:
 * Macro of needed err-handling stuff at start of function;
 * Note that default status is now "FAILURE" (GTO convention)
 */
#ifndef F_INIT
#define F_INIT(x)               \
   int status = FAILURE;        \
   static char outmesg[256];    \
   static char routine[] = x
#endif


/* F_END:
 * Macro of (minimum) err-handling stuff for end of function;
 * Note that (GTO convention) if code has reached "done"
 * without errors, set status to success
 */
#ifndef F_END
#define F_END                   \
     status = SUCCESS;          \
   done:                        \
     return status
#endif


/* XONTRUE: 
 * if error condition (x) true, 
 * prints error msg and jumps to end of curr func 
 */
#ifndef XONTRUE
#define XONTRUE(x,y) { 	                    \
	if ((x)) {                          \
	  MbsMesg(MERROR, TRUE, y);  	    \
	  MbsMesg(MERRTRACE,FALSE,routine); \
          status = FAILURE;	            \
	  goto done;  	                    \
	}			            \
}
#endif


/* XONFAIL: 
 * Detects error from func return (x),
 * and jumps to end of current function 
 */
#ifndef XONFAIL
#define XONFAIL(x) {                     \
	if ((x) != SUCCESS) {              \
	  MbsMesg(MERRTRACE,FALSE,routine);   \
	  status = FAILURE;	           \
	  goto done;                       \
      }                                  \
}
#endif



/* XONFAILMSG: 
 * like XONFAIL, but also prints err msg 
 */
#ifndef XONFAILMSG
#define XONFAILMSG(x,y) {                \
	if ((x) != SUCCESS) {              \
	  MbsMesg(MERROR, TRUE, y);  	   \
	  MbsMesg(MERRTRACE,FALSE,routine);   \
	  status = FAILURE;	           \
	  goto done;                       \
      }                                  \
}
#endif


/* XONRANGE:
 * checks if value exceeds valid (inclusive) range; 
 * generates error if not
 */
#ifndef XONRANGE
#define XONRANGE(x,xlow,xhigh,xname) {           \
       if( (double) (x) < (double) (xlow) ||     \
           (double) (x) > (double) (xhigh) )     \
       {                                         \
           sprintf(outmesg,                      \
                   "\nBad value for %s(%lf): must be in range %lf to %lf", \
               xname,(double) (x), (double) (xlow),(double) (xhigh));  \
           XONTRUE(TRUE,outmesg);                \
       }                                         \
}
#endif



/****************************************************************************
 *  Structure declarations
 ****************************************************************************/

/* Structure for elements of behavior linked list */
typedef struct {
    char behavior_label[LINESTRLEN];        /* label itself */
    void *prev;				    /* ptr to prev item */
    void *next;				    /* ptr to next item */
} BEHAVIOR_STRUC;


/****************************************************************************
 *  STATIC vars
 ****************************************************************************/
/* pointer to start of behavior linked list */
static BEHAVIOR_STRUC *p_behavior_start = NULL; 


/*****************************************************************************
 *   Public functions
 ****************************************************************************/

/****************************************************************************
 *  MbsBaseInit
 *  Newer (8/96) version of init function for mbsbase library
 *  Calling code should call this lib when starting up
 *  Notes:
 *   - can turn on/off msgs to stdout & log file independently
 *   - uses our code to send msgs to stdout
 *   - uses GTO funcs to send msgs to log file
 *   (and therefore this func will turn on/off GTO error logging)
 ***************************************************************************/
void MbsBaseInit
   (char     *appLabel,	        /* (I) brief label (for output msgs) */
    TBoolean  enableStdOutMsgs, /* (I) TRUE=enable output msgs to stdout */
    TBoolean  enableLogMsgs,    /* (I) TRUE=enable output msgs to logfile */
    char     *logFileName);     /* (I) Log file name (if logging enabled);
                                 * if NULL or empty, current/default
                                 * GTO log file name will be used;
                                 * SPECIAL CASE: if name begins with '@',
                                 * will use subsequent label together with
                                 * unique run id to make log file name */


/*****************************************************************************
 *  MbsMesg
 *  Newer version (8/96) of mesg function
 *  Now uses GTO msg functions exclusively
 *****************************************************************************/
void MbsMesg
   (long      msgStatus,         /* (I) MINFO=info mesg(log), MWARN=warning
				   (stderr), MERROR=error(stderr) */
    TBoolean timeStamp,         /* (I) false=no timestamp, true=also print
				   timestamp */
    char    *mesg);             /* (I) msg to print */

/****************************************************************************
 *  mbs_init
 *  Older version of mbsbase init function, for backward compatibility
 ***************************************************************************/
void mbs_init(char *app_label,	/* (I) brief label for app/run (for output
				   msgs) */
	      long enable_mesgs,	/* (I) true enables mbs_mesg() output to
				   stdout/stderr, false prevents any
				   stdout/stderr output			*/
	      char *logfile);   /* (I) name for log file (NULL for none,
				   "" for default) */

/*****************************************************************************
 *  mbs_mesg
 *  Older function to output message (either normal or error msg)
 *  Included for backward compatibility
 *****************************************************************************/
void mbs_mesg(long msgStatus,	/* (I) MINFO=info mesg(log), MWARN=warning
				   (stderr), MERROR=error(stderr) */
	      long timestamp,	/* (I) false=no timestamp, true=also print
				   timestamp */
	      char *mesg);      /* (I) msg to print */


/*****************************************************************************
 *  MbsDbg
 *  Outputs string only in debug mode
 *  (Simplified, special-purpose version of MbsMesg)
 *****************************************************************************/
void MbsDbg(char *msg);

/*****************************************************************************
 *  MbsGetLastErrMesg
 *  Copies most recent err mesg to output string
 *****************************************************************************/
int MbsGetLastErrMesg(long maxlen,            /* (I) max # chars to copy */
                      char *last_errmesg);	/* (O) mesg copied to this
						   string */

/*****************************************************************************
 *  MbsSetPrimaryDataDir
 *  Sets/overrides path to Primary data directory 
 *****************************************************************************/
int MbsSetPrimaryDataDir
   (char *dataDir);         /* (I) name of directory */

/*****************************************************************************
 *  MbsSetDrDataDir
 *  Sets/overrides path to DRDATA data directory 
 *****************************************************************************/
int MbsSetDrDataDir
   (char *dataDir);         /* (I) name of directory */

/*****************************************************************************
 *  MbsSetFimsDataDir
 *  Sets/overrides path to FIMS $DATA data directory 
 *****************************************************************************/
int MbsSetFimsDataDir
   (char *dataDir);         /* (I) name of FIMS $DATA directory */

/***************************************************************
 * MbsGetRunIdString
 * Returns pointer to run id string for this job
 ***************************************************************/
int EXPORT MbsGetRunIdString
   (char  **runIdStr);            /* (O) returns ptr to run id string */



#ifdef __cplusplus
}
#endif

#endif
