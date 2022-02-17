/******************************************************************************
 * Module:	Q3
 * Submodule:
 * File: aliblog.h	
 * Function:	
 * Author:	Interest Rates DR
 * Revision:	$Header$
 *****************************************************************************/
#ifndef	ALIBLOG_H
#define	ALIBLOG_H


# include "cgeneral.h"
# include "cerror.h"
# include "cmemory.h"
# include <stdio.h>              /* FILE  */
# include "macros.h"
# include "cdate.h"	         /* TDate */
# include "dateconv.h"           /* TMonthDayYear */

/*----------------------------------------------------------------------------
 * Variable type definitions.
 */

/* variable type specification */
typedef	long	DVType;

#define	Q3_NULL_T		((DVType) 0)
						
/* basic LIL types */
#define	Q3_DOUBLE_L		((DVType) 80)
#define	Q3_LONG_L		((DVType) 81)
#define	Q3_TDATE_L		((DVType) 82)
#define	Q3_TDATEINTERVAL_L	((DVType) 83)
#define	Q3_CHAR_BLOCK_L	        ((DVType) 84)	
#define	Q3_CHAR_L		((DVType) 85)

/* derived LIL types */
#define	Q3_PERCENT_L		((DVType) 91)	/* Q3_DOUBLE_L, but %  */
#define	Q3_BPOINT_L		((DVType) 92)	/* Q3_DOUBLE_L, but in bp  */

/* unwrap LIL arguments */				       
#define	CASTSTR(ptr, nidx)	(&(((char*)ptr)[WRAP_STR_IDX(nidx+1)]))
#define	CASTTYPE(ptr, nidx, type)	(((type *) ptr)[nidx+1])

/*----------------------------------------------------------------------------
 * loglil.c
 */

int	Q3LilVectLogging         (FILE *fp, 
				  DVType varType, 
				  void *lptr,
				  char  *argName);

int	Q3LilVectLoggingFile     (char *fnam,		
				  char *mode,	      
				  char *funcName, 
				  ...);

int	Q3LilVectLoggingFp       (FILE *fp,
				  char *funcName,  
				  ...);

int     Q3LilVectSize            (DVType varType,
				  void *lptr,
				  char *argName);

int     Q3FPrintf                (FILE *fp,
				  char *fmt,
				  ...);

char*   Q3TDatePrint             (char *string,
				  TDate aDate);

#endif



