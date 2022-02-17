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
#if defined(CLIB)

#if !defined(DBL_EPSILON)
# include <float.h>	/* otherwise DBL_EPSILON redefined by CLIB */
# include <limits.h>	/* otherwise MAX/MIN redefined */
#undef	DBL_EPSILON	/* USE ALIB DBL_EPSILON !!! */
#endif

# include "exponent.h"	/* For GtoExp */
# include "cgeneral.h"
# include "cerror.h"
# include "cmemory.h"
# include "macros.h"
# include "ldate.h"	/* has TDate */

/*t-@CDOC(idxn=TDayCount,catn=structdef)
 * Type for day count convention.
 */
typedef long		TDayCount;
/*e*/

/*
 * The following flags can be defined for compilation:
 * CLIB for linking to the C Analytics Library
 * UNIX to compile under UNIX.
 * _WINDLL to compile under Windows.
 * _WIN32 to compile under Windows NT.
 */


/*
 * The macro  DLL_EXPORT defined for platform compatibility.
 */
# if defined(_WINDLL) && (!defined(_WIN32) && !defined(WIN32))
#  define DLL_EXPORT(type)	type __export
# elif defined(_WIN32) || defined(WIN32)
#  define DLL_EXPORT(type)	__declspec(dllexport) type
# else
#  define DLL_EXPORT(type)	type
# endif
/*e*/


/* misc */
# ifdef	CLIB7X
#  define	GtoErrMsgSet(x)	(!(x))
# endif

# if defined (UNIX) 
#  undef exp		/* do not use GtoExp */
# endif


#else /*!CLIB*/

# include <float.h>
# include <stddef.h>

# if defined(_WINDLL) && (!defined(_WIN32) && !defined(WIN32))
#  define DLL_EXPORT(type)	type __export
# elif defined(_WIN32) || defined(WIN32)
#  define DLL_EXPORT(type)	__declspec(dllexport) type
# else
#  define DLL_EXPORT(type)	type
# endif


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


#define IS_ALMOST_ZERO(x) ((x) < DBL_EPSILON && (x) > -DBL_EPSILON)
#define ARE_ALMOST_EQUAL(x,y) IS_ALMOST_ZERO((x)-(y))
#define IS_BETWEEN(x,a,b) ((a) < (b) ? \
                           ((x) >= (a) && (x) <= (b)) : \
                           ((x) >= (b) && (x) <= (a)))

#define IS_ZERO_TO_TOL(x,tol) ((x) < (tol) && (x) > -1.*(tol) )
#define ARE_EQUAL_TO_TOL(x,y,tol) IS_ZERO_TO_TOL((x)-(y),(tol))

#ifndef PROGRAM_BUG
#define PROGRAM_BUG()   GtoErrMsg("Program bug:%s line %d\n",__FILE__,__LINE__)
#endif

#ifndef NEW
#define NEW(t)              (t *) MALLOC(sizeof(t))
#define NEW_ARRAY(t,n)      (t *) MALLOC(sizeof(t)*(n))
#define FREE_ARRAY(ptr)     FREE(ptr)
#endif


typedef	long	TDate;

/* Definition of a date interval */
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
} TDateInterval;    

#define	TDATEINTERVAL_NIL	{0, '\0', 0}


/* TDateList is a list of dates.  */
typedef struct
{
    int    fNumItems;
    TDate *fArray;
} TDateList;

/* routine name subsitution */
# define GtoErrMsg	PrintErrMsg
# define GtoErrMsgSet	PrintErrMsgSet
DLL_EXPORT(int)	PrintErrMsg(char *fmt,...);
DLL_EXPORT(int)	PrintErrMsgSet(int errMsgStatus);

#endif	/*!CLIB*/



/*----------------------------------------------------------------------
 * Variable Types for LIL interface
 */

/*t-@CDOC(idxn="LIL_Types_Definitions",catn=comm)
 * These type definitions are used in LIL (Lotus Interface Layer)
 * wrapper routines used for addins.
 */

typedef	double		FloatL;
typedef	TDate		TDateL;
typedef	double		TDateIntervalL;
typedef	long int	IntL;
typedef	char		CharBlockL;
typedef	char		CharL;

/*e*/

/*t-@CDOC(idxn=Varianble_Types_Definitions,catn=comm)
 * Variable type definitions.
 */

			/* variable type specification */
typedef	long	TVType;
						/* Data types */
#define	DRL_NULL_T		((TVType) 0)
						/* LIL types */
#define	DRL_DOUBLE_L		((TVType) 80)
#define	DRL_FLOAT_L		((TVType) 80)	/* same as DRL_DOUBLE_L */
#define	DRL_LONG_L		((TVType) 81)
#define	DRL_INT_L		((TVType) 81)	/* same as DRL_LONG_L */
#define	DRL_TDATE_L		((TVType) 82)
#define	DRL_TDATEINTERVAL_L	((TVType) 83)
#define	DRL_CHAR_BLOCK_L	((TVType) 84)	
#define	DRL_CHAR_L		((TVType) 85)
						/* Derived LIL types */
#define	DRL_PERCENT_L		((TVType) 91)	/* DRL_FLOAT_L, but %  */
#define	DRL_BPOINT_L		((TVType) 92)	/* DRL_FLOAT_L, but in bp  */
						/* C types */
#define	DRL_POINTER_T		((TVType) 02)

#define	DRL_DOUBLE_T		((TVType) 10)
#define	DRL_FLOAT_T		((TVType) 11)
#define	DRL_INT_T		((TVType) 12)
#define	DRL_LONG_T		((TVType) 13)
#define	DRL_CHAR_T		((TVType) 14)	/* single char */
#define	DRL_STRING_T		((TVType) 20)	/* char pointer (ie char *p) */
#define	DRL_CHAR_ARRAY_T	((TVType) 21)	/* char array (ie char p[..] */
#define	DRL_TDATE_T		((TVType) 30)
#define	DRL_TDATEINTERVAL_T	((TVType) 31)
#define	DRL_TDAYCOUNT_T		((TVType) 32)
					/* Types build form basic types: */
#define	DRL_PERCENT_T		((TVType) 33)	/* double multilpied by 100 in I/Os */
#define	DRL_CUR_T		((TVType) 34)	/* currency format I/Os */
#define	DRL_CURK_T		((TVType) 35)	/* currency format in K I/Os */
#define	DRL_CURM_T		((TVType) 36)	/* currency format in M I/Os */
#define	DRL_BOOLEAN_T		((TVType) 37)	/* boolean (TRUE,FALSE) as int */

					/* Objects types */
#define	DRL_CVAR_T		((TVType) 129)
#define	DRL_CARRAY_T		((TVType) 130)
#define	DRL_CMATRIX_T		((TVType) 131)
#define	DRL_CVECTOR_T		((TVType) 134)
#define	DRL_CDVECTOR_T		((TVType) 132)
#define	DRL_CDMATRIX_T		((TVType) 133)
					/* complex data structures */
					/* LIL types */
#define	DRL_LILVAR_L		((TVType) 257)
#define	DRL_LILARRAY_L		((TVType) 258)
#define	DRL_LILVECT_L		((TVType) 259)
#define	DRL_LILMATR_L		((TVType) 260)
#define	DRL_LILVECTARRAY_L	((TVType) 261)
#define	DRL_LILVECTRANGE_L	((TVType) 262)
#define	DRL_LILMATR2VECT_L	((TVType) 263)

/*e*/


/*----------------------------------------------------------------------
 * The following macros are recommended for wrappers for compatibility
 * with future LIL,EIL, etc. releases
 */
#ifndef	ARGSIZE	/* Defined in macros.h */
#define ARGSIZE(arg) ((unsigned) (arg)[0])
#define IS_SCALAR(arg) (ARGSIZE(arg) == 1)
#define ISNT_SCALAR(arg) (ARGSIZE(arg) > 1)
#define IS_VECTOR(arg) (ARGSIZE(arg) >= 1)
#endif	/* ARGSIZE */

/*
 * Use this to access the ith string in a Wrapper string array.
 */
#ifdef	_SKIP
#ifndef WRAP_STR_BYTES
# define WRAP_STR_BYTES	128
#endif
#ifndef WRAP_STR_IDX
# define WRAP_STR_IDX(idx) ((idx-1)*WRAP_STR_BYTES + 1)
#endif
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
#define	WRAP_MATR_IDX(iIdxMax, jIdxMax, iIdx, jIdx) \
		((jIdx)+(iIdx)*(jIdxMax)+1)



/*
 * Useful macros to check arg len. Assume that
 * char *routine and label done are deined.
 */
#define	WRAP_CHECK_SCALAR(arg)	{if (ARGSIZE(arg) != 1) {GtoErrMsg(\
				"%s: argument `%s' is not a scalar\n",\
				routine, #arg); goto done;}}

#define	WRAP_CHECK_VECTOR(arg)	{if (ARGSIZE(arg) < 1) {GtoErrMsg(\
				"%s: argument `%s' is not a vector (len=%d)\n",\
				routine, #arg, ARGSIZE(arg)); goto done;}}

#define	WRAP_CHECK_VECTOR_LEN(arg, len)	{if (ARGSIZE(arg) != (len)) {\
				GtoErrMsg("%s: argument `%s' is a vector "\
				"of length %d (expected %d)\n",routine, \
				#arg, ARGSIZE(arg), (len)); goto done;}}


#define	ASSERT_OR_DONE(cond) 	{if (!(cond)) {GtoErrMsg(\
				"%s: assertion `%s' failed.\n",\
				routine, #cond); goto done;}}

#define	IF_FAILED_DONE(statement)	\
				{if ((statement) != SUCCESS) {goto done;};}


/* Input volatility type: percentage or basis point vol */
typedef enum { 
	NORMVOL,
	LOGVOL
} KVolType;


#endif	/* _drlstd_H */
