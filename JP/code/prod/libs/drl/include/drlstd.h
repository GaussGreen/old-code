/************************************************************************
 * Module:	DRL - DR C Utilities
 * Submodule:	
 * File:	drlstd.h
 * Function:	Standard definitions
 *		MUST BE INCLUDED BEFORE ANY SYSTEM HEADER FILE
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#ifndef	_drlstd_H
#define	_drlstd_H

/*
 * C Library Compatibility
 */
#if defined(DRL_CLIB)

#if !defined(DBL_EPSILON)
# include <float.h>	/* otherwise DBL_EPSILON redefined by DRL_CLIB */
# include <limits.h>	/* otherwise MAX/MIN redefined */
# undef	DBL_EPSILON	/* USE ALIB DBL_EPSILON !!! */
#endif

# include "exponent.h"	/* For GtoExp */
# include "cgeneral.h"
# include "cerror.h"
# include "cmemory.h"
# include "macros.h"
# include "ldate.h"
# include "bastypes.h"

typedef	TDate	DDate;
typedef	TDateInterval	DInterval;    

# include "drlerr.h"	/* Macros DrlErrMsg, etc. */

typedef long		DDayCount;
#define DRL_ACT_ACT	GTO_ACT_ACT
#define DRL_ACT_365F	GTO_ACT_365F
#define DRL_ACT_360	GTO_ACT_360
#define DRL_B30_360	GTO_B30_360
#define DRL_B30E_360	GTO_B30E_360

/*
 * The following flags can be defined for compilation:
 * DRL_CLIB for linking to the C Analytics Library
 * UNIX to compile under UNIX.
 * _WINDLL to compile under Windows.
 * _WIN32 to compile under Windows NT.
 */

# if defined (UNIX) 
#  undef exp		/* do not use GtoExp */
# endif


#else /*!DRL_CLIB*/

/*---------------------------------------------------------------------*
 * Stand alone library compilation
 *---------------------------------------------------------------------*/

# include <float.h>
# include <stddef.h>
# include <stdlib.h>

/*
 * Standard macros used
 */

# ifndef TRUE
#  define TRUE 1
# endif
# ifndef FALSE
#  define FALSE 0
# endif

# ifndef SUCCESS
#  define SUCCESS 0
# endif

# ifndef FAILURE
#  define FAILURE -1
# endif

# ifndef MAX
#  define MAX(a,b) ((a) > (b) ? (a) : (b))
# endif

# ifndef MIN
#  define MIN(a,b) ((a) < (b) ? (a) : (b))
# endif


#define IS_ALMOST_ZERO(x) ((x) < DBL_EPSILON && (x) > -DBL_EPSILON)
#define ARE_ALMOST_EQUAL(x,y) IS_ALMOST_ZERO((x)-(y))
#define IS_BETWEEN(x,a,b) ((a) < (b) ? \
                           ((x) >= (a) && (x) <= (b)) : \
                           ((x) >= (b) && (x) <= (a)))

#define IS_ZERO_TO_TOL(x,tol) ((x) < (tol) && (x) > -1.*(tol) )
#define ARE_EQUAL_TO_TOL(x,y,tol) IS_ZERO_TO_TOL((x)-(y),(tol))

#ifndef PROGRAM_BUG
#define PROGRAM_BUG()   DrlErrMsg("Program bug:%s line %d\n",__FILE__,__LINE__)
#endif

#ifndef NOT_IMPLEMENTED
#define NOT_IMPLEMENTED	DrlErrMsg("Not implemnted:%s line %d\n",__FILE__,__LINE__)
#endif

#define	MALLOC	malloc
#define	FREE	free



#ifndef NEW
#define NEW(t)              (t *) MALLOC(sizeof(t))
#define NEW_ARRAY(t,n)      (t *) MALLOC(sizeof(t)*(n))
#define FREE_ARRAY(ptr)     FREE(ptr)
#endif

/*
 * Standard types
 */

typedef	int	DBoolean;

/*t-------------------------------------------------------------
 * Type for day count convention.
 * 
 * <br>The following macros are defined:<br>
 * <PRE>
 * #define DRL_ACT_ACT           'A'
 * #define DRL_ACT_365F          '5'
 * #define DRL_ACT_360           '0'
 * #define DRL_B30_360           '3'
 * </PRE>
 */
typedef char		DDayCount;
/*e*/


#define DRL_ACT_ACT		'A'		/* Actual/365 */
#define DRL_ACT_365F		'5'		/* Actual/365 Fixed */
#define DRL_ACT_360		'0'		/* Actual/360 */
#define DRL_B30_360		'3'		/* 30/360 */
#define DRL_B30E_360		'E'		/* 30E/360 */
#define DRL_ACT_365FJ		'J'		/*  */
#define DRL_B30E_360I		'I'
#define DRL_B30_360_FIXED	'D'		/* For bond coupon */
#define DRL_B30EP_360		'P'		/* 30E+/360 */
#define DRL_ACT_ACT_FRF		'F'		/* For French TAM  */


/*#define	DRL_TDATE_MDY*/
#ifndef	DRL_TDATE_MDY
/*t-------------------------------------------------------------
 * Basic date type.
 * <br>
 * Uses SRM3 date long format YYYYMMDD.
 */
typedef	long	DDate;		/* YYYYMMDD format */
/*e*/

#else
typedef	struct	{ int	m; int	d; int	y; } DDate;
#endif


/*t-------------------------------------------------------------
 * Definition of a date interval : mimic ALIB.
 */
typedef struct
{
    int prd;        /* number of periods from offset date*/
    char prd_typ;   /* type of periods                   */
                    /* D - day; M - month; W - week;     */
                    /* Q - 3 months; S - 6 months;       */
                    /* A - 12 months; I - IMM period     */
                    /* F - max of next IMM or +3 month   */
    int flag;       /* 0 - offset is value date
                       -1 - offset is the previous date in the date array
                       x - any other number is index into array of intervals.
                           the date at that location is an offset */
} DInterval;    
/*e*/

#endif	/*!DRL_CLIB*/



/*----------------------------------------------------------------------
 * Variable type definitions.
 */

			/* variable type specification */
typedef	long	DVType;
						/* Data types */
#define	DRL_NULL_T		((DVType) 0)
						/* LIL types */
#define	DRL_DOUBLE_L		((DVType) 80)
#define	DRL_FLOAT_L		((DVType) 80)	/* same as DRL_DOUBLE_L */
#define	DRL_LONG_L		((DVType) 81)
#define	DRL_INT_L		((DVType) 81)	/* same as DRL_LONG_L */
#define	DRL_TDATE_L		((DVType) 82)
#define	DRL_TDATEINTERVAL_L	((DVType) 83)
#define	DRL_CHAR_BLOCK_L	((DVType) 84)	
#define	DRL_CHAR_L		((DVType) 85)
						/* Derived LIL types */
#define	DRL_PERCENT_L		((DVType) 91)	/* DRL_FLOAT_L, but %  */
#define	DRL_BPOINT_L		((DVType) 92)	/* DRL_FLOAT_L, but in bp  */
						/* C types */
#define	DRL_POINTER_T		((DVType) 02)

#define	DRL_DOUBLE_T		((DVType) 10)
#define	DRL_FLOAT_T		((DVType) 11)
#define	DRL_INT_T		((DVType) 12)
#define	DRL_LONG_T		((DVType) 13)
#define	DRL_CHAR_T		((DVType) 14)	/* single char */
#define	DRL_STRING_T		((DVType) 20)	/* char pointer (ie char *p) */
#define	DRL_CHAR_ARRAY_T	((DVType) 21)	/* char array (ie char p[..] */
#define	DRL_TDATE_T		((DVType) 30)
#define	DRL_TDATEINTERVAL_T	((DVType) 31)
#define	DRL_TDAYCOUNT_T		((DVType) 32)
					/* Types build form basic types: */
#define	DRL_PERCENT_T		((DVType) 33)	/* double multilpied by 100 in I/Os */
#define	DRL_CUR_T		((DVType) 34)	/* currency format I/Os */
#define	DRL_CURK_T		((DVType) 35)	/* currency format in K I/Os */
#define	DRL_CURM_T		((DVType) 36)	/* currency format in M I/Os */
#define	DRL_BOOLEAN_T		((DVType) 37)	/* boolean (TRUE,FALSE) as int */

					/* Objects types */
#define	DRL_CVAR_T		((DVType) 129)
#define	DRL_CARRAY_T		((DVType) 130)
#define	DRL_CMATRIX_T		((DVType) 131)
#define	DRL_CVECTOR_T		((DVType) 134)
#define	DRL_CDVECTOR_T		((DVType) 132)
#define	DRL_CDMATRIX_T		((DVType) 133)
					/* complex data structures */
					/* LIL types */
#define	DRL_LILVAR_L		((DVType) 257)
#define	DRL_LILARRAY_L		((DVType) 258)
#define	DRL_LILVECT_L		((DVType) 259)
#define	DRL_LILMATR_L		((DVType) 260)
#define	DRL_LILVECTARRAY_L	((DVType) 261)
#define	DRL_LILVECTRANGE_L	((DVType) 262)
#define	DRL_LILMATR2VECT_L	((DVType) 263)

/*e*/





/*----------------------------------------------------------------------
 * The following macros are recommended for wrappers for compatibility
 * with future LIL,EIL, etc. releases
 */

#define DRL_ARGSIZE(arg)	((unsigned) (arg)[0])
#define DRL_IS_SCALAR(arg)	(DRL_ARGSIZE(arg) == 1)
#define DRL_ISNT_SCALAR(arg)	(DRL_ARGSIZE(arg) > 1)
#define DRL_IS_VECTOR(arg)	(DRL_ARGSIZE(arg) >= 1)


/*
 * Use this to access the ith string in a Wrapper string array.
 */

#ifndef DRL_WRAP_STR_BYTES
# define DRL_WRAP_STR_BYTES	128
#endif
#ifndef DRL_WRAP_STR_IDX
# define DRL_WRAP_STR_IDX(idx) (((idx)-1)*DRL_WRAP_STR_BYTES + 1)
#endif


/*
 * For an array in LIL: goes by lines first
 *	-----------> 		----------------->
 *	|     Y			0		NY-1
 *	|			----------------->
 *	| X			NY+1		2*NY-1
 *	|			  	....
 *	|
 *	V
 * Element of (X=i, Y=j) is located in the
 * wrapper array at (j + NY*i).
 */
#define	DRL_WRAP_MATR_IDX(iIdxMax, jIdxMax, iIdx, jIdx) \
		((jIdx)+(iIdx)*(jIdxMax)+1)


/*
 * Useful macros to check arg len. Assume that
 * char *routine and label done are deined.
 */
#define	DRL_WRAP_CHECK_VECTOR(arg)	\
		{if (DRL_ARGSIZE(arg) < 1) {DrlErrMsg(\
		"%s: argument `%s' is not a vector (len=%d)\n",\
		routine, #arg, DRL_ARGSIZE(arg)); goto done;}}

#define	DRL_WRAP_CHECK_SCALAR(arg)	\
		{if (DRL_ARGSIZE(arg) != 1) {DrlErrMsg(\
		"%s: argument `%s' is not a scalar\n",\
		routine, #arg); goto done;}}


#define	DRL_WRAP_CHECK_VECTOR_LEN(arg, len)	\
		{if (DRL_ARGSIZE(arg) != (len)) {\
		DrlErrMsg("%s: argument `%s' is a vector "\
		"of length %d (expected %d)\n",routine, \
		#arg, DRL_ARGSIZE(arg), (len)); goto done;}}


#define	ASSERT_OR_DONE(cond) 	{if (!(cond)) {DrlErrMsg(\
				"%s: assertion `%s' failed.\n",\
				routine, #cond); goto done;}}

#define	IF_FAILED_DONE(statement)	\
				{if ((statement) != SUCCESS) {goto done;};}


#endif	/* _drlstd_H */
