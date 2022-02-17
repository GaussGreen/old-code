/*********************************************************************************
 * CRXUTILIO.H 
 * crx io utils
 *
 ********************************************************************************/

#ifndef __CRXUTILIO_H__
#define __CRXUTILIO_H__

#ifdef __cplusplus
extern "C" {
#endif

#include <common/include/drmacros.h>
#include "t_curve.h"
#include "l_date.h"
#include "crxerror.h"
#include "crxutil.h"
    
#include <alib/cerror.h>
    
typedef	long		TDayCount;
    
#define	READ_DATA(type,ptr,str)	\
    { if (CrxFScanVType(fp, type, (void*) ptr) != SUCCESS) \
	  { DR_Error("%s: can't read %s.\n", routine, str); \
								 goto RETURN;}}

#define NEW_ARRAY(t,n)      (t *) MALLOC(sizeof(t)*(n))
#define FREE_ARRAY(ptr)     FREE(ptr)

/*
 * Useful macros to check arg len. Assume that
 * char *routine and label done are deined.
 */
/* Defined in ALIB
#ifndef ARGSIZE
#define ARGSIZE(arg) ((unsigned) (arg)[0])
#define IS_SCALAR(arg) (ARGSIZE(arg) == 1)
#define ISNT_SCALAR(arg) (ARGSIZE(arg) > 1)
#define IS_VECTOR(arg) (ARGSIZE(arg) >= 1)
#endif
*/

#ifndef WRAP_CHECK_SCALAR
#define WRAP_CHECK_SCALAR(arg)  {if (ARGSIZE(arg) != 1) {DR_Error(\
                "%s: argument `%s' is not a scalar\n",\
                routine, #arg); goto RETURN;}}

#define WRAP_CHECK_VECTOR(arg)  {if (ARGSIZE(arg) < 1) {DR_Error(\
                "%s: argument `%s' is not a vector (len=%d)\n",\
                routine, #arg, ARGSIZE(arg)); goto RETURN;}}

#define WRAP_CHECK_VECTOR_LEN(arg, len) {if (ARGSIZE(arg) != (len)) {\
                DR_Error("%s: argument `%s' is a vector "\
                "of length %d (expected %d)\n",routine, \
                #arg, ARGSIZE(arg), (len)); goto RETURN;}}


#define ASSERT_OR_RETURN(cond)    {if (!(cond)) {DR_Error(\
                "%s: assertion `%s' failed.\n",\
                routine, #cond); goto RETURN;}}

#define IF_FAILED_RETURN(statement)   \
                {if ((statement) != SUCCESS) {goto RETURN;};}

#endif


/* variable type specification */
typedef	long	TVType;
						/* Data types */
#define	DRL_NULL_T		    ((TVType) 0)
						/* LIL types */
#define	DRL_DOUBLE_L		((TVType) 80)
#define	DRL_FLOAT_L		    ((TVType) 80)	/* same as DRL_DOUBLE_L */
#define	DRL_LONG_L		    ((TVType) 81)
#define	DRL_INT_L		    ((TVType) 81)	/* same as DRL_LONG_L */
#define	DRL_TDATE_L		    ((TVType) 82)
#define	DRL_TDATEINTERVAL_L	((TVType) 83)
#define	DRL_CHAR_BLOCK_L	((TVType) 84)	
#define	DRL_CHAR_L		    ((TVType) 85)
						/* Derived LIL types */
#define	DRL_PERCENT_L	    ((TVType) 91)	/* DRL_FLOAT_L, but %  */
#define	DRL_BPOINT_L	    ((TVType) 92)	/* DRL_FLOAT_L, but in bp  */
						/* C types */
#define	DRL_POINTER_T	    ((TVType) 02)

#define	DRL_DOUBLE_T	    ((TVType) 10)
#define	DRL_FLOAT_T		    ((TVType) 11)
#define	DRL_INT_T		    ((TVType) 12)
#define	DRL_LONG_T		    ((TVType) 13)
#define	DRL_CHAR_T		    ((TVType) 14)	/* single char */
#define	DRL_STRING_T		((TVType) 20)	/* char pointer (ie char *p) */
#define	DRL_CHAR_ARRAY_T	((TVType) 21)	/* char array (ie char p[..] */
#define	DRL_TDATE_T		    ((TVType) 30)
#define	DRL_TDATEINTERVAL_T	((TVType) 31)
#define	DRL_TDAYCOUNT_T		((TVType) 32)
					/* Types build form basic types: */
#define	DRL_PERCENT_T		((TVType) 33)	/* double multilpied by 100 in I/Os */
#define	DRL_CUR_T		    ((TVType) 34)	/* currency format I/Os */
#define	DRL_CURK_T		    ((TVType) 35)	/* currency format in K I/Os */
#define	DRL_CURM_T		    ((TVType) 36)	/* currency format in M I/Os */
#define	DRL_BOOLEAN_T		((TVType) 37)	/* boolean (TRUE,FALSE) as int */

					/* Objects types */
#define	DRL_CVAR_T		    ((TVType) 129)
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

int		CrxFScanVType(FILE *fp, TVType type, void *valPtr);
int		CrxFPrintVType(FILE *fp, TVType type, void *valPtr);

int     CrxVTypeScan(char *str, TVType type, void *p);
char*   CrxVTypePrint(char *str, TVType type, void *p);

int	CrxLilVectArrayFpRead(
	FILE *fp,                             /* (I) file pointer */
	int numItems,                         /* (I) num elements in each column */
	int numVect,                          /* (I) num of columns */
	TVType *varType,                      /* (I) variable types [0..numVect-1] */
	void ***lptr);                        /* (I) array of & of pointers [0..numVect-1] */

int CrxLilVectArrayFpReadV(
	FILE *fp,                             /* (I) file pointer */
	int numItems,                         /* (I) num elements to be read */
	/* TVType varType, void *lptr,
	 * ...
	 * TVType varType, void *lptr,
	 * DRL_NULL_T (last argument MUST be DRL_NULL_T)
	 */
	...);

/*f----------------------------------------------------------------------
 * Returns the name of a type "type".
 */
char* CrxVTypeName(TVType type);

/*f----------------------------------------------------------------------
 * Returns SUCCESS if type "type" is a valid C type.
 */
int CrxVTypeCheckValidCType(TVType type);

/*f----------------------------------------------------------------------
 * Returns 0 if type <i> type</i> is a valid LIL type.
 */
int CrxVTypeCheckValidLType(TVType type);
    
/*f----------------------------------------------------------------------
 * Returns the size of <i> type</i>.
 */
size_t CrxVTypeSizeof(TVType type);

/**----------------------------------------------------------------------
 * Returns the element with offset "offset" in an array
 * "p" of type "type".
 */
void* CrxVTypeOffsetVect(void *p, int offset, TVType type);



/**-------------------------------------------------------------
 * Date routines : scan date in YYYYMDD format.
 *
 * Scan the string "string" for a date. Recognized format
 * is YYYYMMDD. Puts the
 * result in "aDate". Returns 0 if success.
 */

int CrxTDateScanYMD(char *string, TDate *aDate);


/*f-------------------------------------------------------------
 * I/O: read next line (skip comments).
 */
char* CrxFGetLine(char *s, int n, FILE *fp, int *line);


/*f---------------------------------------------------------------------
 * Scans a long value in a string <i> s</i> by trying to match
 * the string with different values. On successful return, the
 * value is put in <i> value</i>.
 * <br> <b> Example:</b>
 * \begin{verbatim}
 * char  *buf;
 * long  value;
 * ...
 * errCode = CrxStrLongValueScan(buf, "stub type" &value,
 *              "B", GTO_BOND_STUB,
 *              "S", GTO_SIMPLE_STUB,
 *              NULL);
 * ...
 * \end{verbatim}
 */
int CrxStrLongValueScan( 
        char *s,                /* (I) input string (unchanged) */
        char *errName,          /* (I) variable name for debugging */
        long *value,            /* (O) value */
        /* char *varName1, long varValue1,
         * ...
         * char *varNameN, long varValueN,
         * NULL)                        LAST ARGUMENT MUST BE NULL
         */ ...);


/**---------------------------------------------------------
 *  I/O: advance to specific token (skip comments).
 *  Returns SUCCESS/FAILURE.
 */
int CrxFAdvanceToToken(FILE *fp, char *token);



/**---------------------------------------------------------
 */
int CrxTCurveFileRead(TCurve **that, char *fnam);



/*f-------------------------------------------------------------
 * Sorts a vector <i> p</i> with <i> numItem</i> elements
 * of type <i> type</i> in ascending order.
 * If the flag <i> removeMultiple</i> is TRUE, removes multiple
 * instances of a same element, and <i> numItems</i>
 * may be changed on exit.
 */


int
CrxVTypeVectSort(
    void *p,             /* (B) vector */
    int *numItems,       /* (B) num elem in vector */
    TVType type,         /* (I) vector type */
    int removeMultiple); /* (I) TRUE=remove multiple items */


/* simple set of DR wrapper utility functions for which we changed the
   prefix and compiled locally */

typedef int (*CrxTStringConverter)(const char*, void*);

/* none of the StringTo... functions reports errors in cases of failure */
int CrxStringToSize (const char* str, size_t *val);
int CrxStringToInt (const char* str, int *val);
int CrxStringToLong (const char* str, long *val);
int CrxStringToDouble (const char* str, double *val);
int CrxStringToDate (const char* str, TDate *val);
int CrxStringToBoolean (const char* str, TBoolean *val);
int CrxStringToString (const char* src, char **val);

int CrxReadSize (FILE*, size_t*);
int CrxReadInt (FILE*, int*);
int CrxReadLong (FILE*, long*);
int CrxReadDouble (FILE*, double*);
int CrxReadDate (FILE*, TDate*);
int CrxReadBoolean (FILE*, TBoolean*);
int CrxReadString (FILE*, char**);
int CrxReadArrays (FILE*, size_t, char*, ...);
int CrxReadGenericData (FILE*, const char*, CrxTStringConverter, void*);
char* CrxReadToken (FILE*, size_t, char*);

#define CrxDateToYYYYMMDD CrxTDate2DrDate
    
/* end of extern "C" scope */
#ifdef __cplusplus
}
#endif

#endif
