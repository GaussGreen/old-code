/****************************************************************************/
/*      Standard input output for yield and volatility curves.              */
/****************************************************************************/
/*      STDINPUT.c                                                          */
/****************************************************************************/


#include <ctype.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>

#include <irx/macros.h>
#include <irx/date.h>
#include <irx/convert.h>

#include "irx/dateutils.h"
#include "irx/irxutilio.h"


/*f--------------------------------------------------------------
 * I/O: read variable type (skip comments).
 *
 * <br><br>
 * Scans a variable type value of type <i> type</i>
 * in the file <i> fp</i> and puts the result in <i> valPtr</i>.
 * Returns SUCCESS/FAILURE.
 */
int irxFScanVType(FILE *fp, IrxTVType type, void *valPtr)
{
static        char        buf[1024];
        if (irxFScanString(fp, buf) != SUCCESS)
                return(FAILURE);

        if (irxVTypeScan(buf, type, valPtr) != SUCCESS)
                return(FAILURE);

        return(SUCCESS);
}




/*f----------------------------------------------------------------------
 * Allocates a vector of type "type" of length "size".
 */
void* irxVTypeVectAlloc(int size, IrxTVType type)
{
static        char        routine[] = "irxVTypeVectAlloc";
        int         i;
        void        *p = NULL;
        size_t      typeSize;

        switch (type) {
        case IRX_DOUBLE_T:
        case IRX_PERCENT_T:
        case IRX_BPOINT_T:
                p = (void*) NEW_ARRAY(double, size+1);
                if (p) ((double*)p)[0] = (double) size;
                return(p);
        case IRX_INT_T:
                p = (void*) NEW_ARRAY(int, size+1);
                if (p) ((int*)p)[0] = (int) size;
                return(p);
        case IRX_LONG_T:
                p = (void*) NEW_ARRAY(long, size+1);
                if (p) ((long*)p)[0] = (long) size;
                return(p);
        case IRX_TDATE_T:
                p = (void*) NEW_ARRAY(IrxTDate, size+1);
                if (p) ((IrxTDate*)p)[0] = (IrxTDate)size;
                return(p);

        case IRX_CHAR_BLOCK_T:
                p = (void*) NEW_ARRAY(char, 128*size+1);
                if (p) ((char*)p)[0] = (char)((unsigned char) size);
                return(p);

                /**
                 ** Other
                 **/
        case IRX_STRING_T:
                /* for IRX_STRING_T, allocate char string */
                p = (void*) NEW_ARRAY(char*, size);
                if (p == NULL) {
                        irxError("%s: malloc IRX_STRING_T "
                                "length %d failed.\n", routine, (int) (size));
                        return(NULL);
                }
                    for (i=0; i<=size-1; i++) {
                        ((char**)p)[i] = NEW_ARRAY(char,  WRAP_STR_BYTES);
                        if (((char**)p)[i] == NULL) {
                                irxError("%s: malloc failed.\n", routine);
                                return(NULL);
                        }
                }
                break;


        default:
                /* check type size */
                if (irxVTypeCheckValidCType(type) != 0) {
                        irxError("%s: bad C type.\n", routine);
                        return(NULL);
                }

                 /* get type size */
                if ((typeSize = irxVTypeSizeof(type)) <= 0) {
                        irxError("%s: bad type size.\n", routine);
                        return(NULL);
                }

                /* allocate memory */
                if ((p = (void*) irxMemAlloc((size_t) (typeSize*size))) == NULL) {
                        irxError("%s: malloc length %d failed.\n", routine, 
                                        (int) (typeSize*size));
                        return(NULL);
                }
                break;
        }

        return(p);
}


/*f----------------------------------------------------------------------
 * Frees an array of type "type" of length "size".
 */

int irxVTypeVectFree(void* p, int size, IrxTVType type)
{
        int        i;

        /* nothing to do */
        if (p == NULL) return(0);

        switch (type) {
                /**
                 ** LIL counted arrays
                 **/
        case IRX_DOUBLE_T:
        case IRX_PERCENT_T:
        case IRX_BPOINT_T:
                FREE((double*) p);
                return(SUCCESS);
        case IRX_INT_T:
                FREE((int*) p);
                return(SUCCESS);
        case IRX_LONG_T:
                FREE((long*) p);
                return(SUCCESS);
        case IRX_TDATE_T:
                FREE((IrxTDate*) p);
                return(SUCCESS);
        case IRX_CHAR_BLOCK_T:
                FREE((char*) p);
                return(SUCCESS);
        case IRX_STRING_T:
                    for (i=0; i<=size-1; i++) {
                        if (((char**)p)[i] != NULL)
                                FREE((void*) ((char**)p)[i]);
                    }
                FREE((void*) p);
                break;
        default:
                if (irxVTypeCheckValidCType(type) != 0) {
                        irxError("irxVTypeVectFree: bad C type\n");
                        return(1);
                }

                FREE((void*) p);

                return(0);
        }
        return(SUCCESS);
}

/*f--------------------------------------------------------------
 * I/O: read char string (skip comments).
 *
 * <br><br>
 * Scans a string value in the file pointer <i> fp</i> and puts
 * the result in <i> s</i>
 * (a string being a block of characters either enclosed
 * in double brackets or not containing any space).
 * Skips all end of lines after '\#' is encountered.
 * Returns SUCCESS/FAILURE.
 */

int irxFScanString(FILE *fp, char *s)
{
        char       *q;
        int        c;

        q = s;

        /* skip over '#' commented lines */
        while (((c = getc(fp)) == '#') || (isspace(c))) {
            if (c == '#') {
                while (((c = getc(fp)) != EOF) && (c != '\n'));
                if (c == EOF) return(FAILURE);
            }
        }

        if (c == EOF)
                return(FAILURE);

        if (c != '"') {
                /*return (fscanf(fp, "%s", s) != 1);*/
                *q++ = c;
                while (((c = getc(fp)) != EOF) && (!isspace(c)))
                        *q++ = c;
                *q = '\0';
        } else {
                while (((c = getc(fp)) != EOF) && (c != '"'))
                        *q++ = c;
                *q = '\0';
        }

        return(SUCCESS);
}



/*f--------------------------------------------------------------
 * I/O: print variable type.
 *
 * <br><br>
 * Prints a variable type value <i> valPtr</i> of type <i> type</i>
 * in the file <i> fp</i>.
 * Returns SUCCESS/FAILURE.
 */

int irxFPrintVType(FILE *fp, IrxTVType type, void *valPtr)
{
static        char        buf[256];

        strcpy(buf, irxVTypePrint(NULL, type, valPtr));

        if (fp == NULL) {
                irxError(buf);
        } else {
                fputs(buf, fp);
        }

        return(SUCCESS);
}


/*f----------------------------------------------------------------------
 * Returns the element with offset "offset" in an array
 * "p" of type "type".
 */

static void* irxVTypeOffsetVect(void *p, int offset, IrxTVType type)
{

        switch (type) {
        case IRX_POINTER_T:
                return(((void**)p) + offset);
        case IRX_DOUBLE_T:
        case IRX_PERCENT_T:
        case IRX_CUR_T:
        case IRX_CURK_T:
        case IRX_CURM_T:
                return(((double*)p) + offset);
        case IRX_FLOAT_T:
                return(((float*)p) + offset);
        case IRX_INT_T:
        case IRX_BOOLEAN_T:
                return(((int*)p) + offset);
        case IRX_LONG_T:
                return(((long*)p) + offset);
        case IRX_CHAR_T:
                return(((char*)p) + offset);
        case IRX_CHAR_BLOCK_T:
        case IRX_STRING_T:
                return(((char**)p) + offset);
        case IRX_TDATE_T:
                return(((IrxTDate*)p) + offset);
        case IRX_TDATEINTERVAL_T:
                return(((IrxTDateInterval*)p) + offset);
        case IRX_TDAYCOUNT_T:
                return(((IrxTDayCountConv*)p) + offset);

        default:
                irxError("OffsetType: bad type %ld.\n", (long)type);
                return(NULL);
        }
}


/*f----------------------------------------------------------------------
 * Reads a file <i> fp</i> an array of <i> numVect</i> columns
 * each containing <i> numItems</i> elements,
 * of the form
 * \begin{verbatim}
 *   <elem #1        of vect #1>  ...  <elem #1        of vect #numVect> 
 *   ...
 *   <elem #numItems of vect #1>  ...  <elem #numItems of vect #numVect> 
 *   end
 * \end{verbatim}
 */

int irxVectArrayFpRead(
        FILE *fp,                             /* (I) file pointer */
        int numItems,                         /* (I) num elements in each column */
        int numVect,                          /* (I) num of columns */
        IrxTVType *varType,                      /* (I) array of variable types [0..numVect-1] */
        void ***lptr)                         /* (I) array of void* pointers [0..numVect-1] */
{
static  char       routine[] = "irxVectArrayFpRead";
        int        status = FAILURE;

        int        idxV, idx;
        void        *vptr;


        /* reset to NULL */
        for (idxV=0; idxV<=numVect-1; idxV++)
                *(lptr[idxV]) = (void*) NULL;

        /* malloc */
        for (idxV=0; idxV<=numVect-1; idxV++) {
                *(lptr[idxV]) = irxVTypeVectAlloc(numItems, varType[idxV]);
                if (*(lptr[idxV]) == NULL) goto RETURN;
        }



        /* read elements */
        for (idx=0; idx<=numItems-1; idx++) {
            for (idxV=0; idxV<=numVect-1; idxV++) {

                /* access element #idx of array #idxV */
                vptr = irxVTypeOffsetVect(*(lptr[idxV]), idx, varType[idxV]);
                if (vptr == NULL) {
                        goto RETURN;
                }

                /* read element */
                if (irxFScanVType(fp, varType[idxV], vptr) != SUCCESS) {
                    irxError("%s: can't read element #%d/%d "
                        "in column # %d/%d "
                        "(expecting type %s).\n", routine,
                        idx+1, numItems, idxV+1, numVect,
                        irxVTypeName(varType[idxV]));
                    goto RETURN;
                }
            }
        }


        /* made it through */
        status = SUCCESS;
RETURN:
        if (status != 0) {
                for (idxV=0; idxV<=numVect-1; idxV++)
                    irxVTypeVectFree(*(lptr[idxV]), numItems, varType[idxV]);
                irxError("%s: failed.\n", routine);
        }
        return(status);
}


/*f----------------------------------------------------------------------
 * A version of <i> irxVectArrayFpRead</i> that allows the call
 * with variable number of arguments (one does not need to
 * set up arrays as in <i> irxVectArrayFpRead</i>. \\
 * <b> Example:</b> to read 
 * \begin{verbatim}
 * # dates and notional
 *    01-Jan-1998    100.0
 *    01-Jul-1998     80.0
 *    01-Jan-1999     60.0
 *    01-Jul-1999     40.0
 *    01-Jan-1999     20.0
 * \end{verbatim}
 * call the routine as in 
 * \begin{verbatim}
 *     IrxTDate      *dates = NULL;
 *     double        *notionals = NULL;
 *     ...
 *     status = irxVectArrayFpReadV(
 *                    stdin,
 *                    5,
 *                    IRX_TDATE_T,  (void*) &dates,
 *                    IRX_DOUBLE_T, (void*) &notionals,
 *                    IRX_NULL_T);
 *     ...
 * \end{verbatim}
 */

int irxVectArrayFpReadV(
        FILE *fp,                             /* (I) file pointer */
        int numItems,                         /* (I) num elements to be read */
        /* IrxTVType varType, void *lptr,
         * ...
         * IrxTVType varType, void *lptr,
         * IRX_NULL_T (last argument MUST be IRX_NULL_T)
         */
        ...)
{
static        char      routine[] = "irxVectArrayFpReadV";
        int             status = FAILURE;

        va_list         ap;
#undef  MAX_ARGS
#define MAX_ARGS        32
        IrxTVType       varType[MAX_ARGS];
        void            **lptr[MAX_ARGS];
        int             numVect = 0;


        /* get arguments */
        va_start(ap, numItems);
        while ((varType[numVect] = (IrxTVType) va_arg(ap, IrxTVType))
                        != IRX_NULL_T) {
                lptr[numVect] = (void**) va_arg(ap, void*);
                numVect++;
                if (numVect >= MAX_ARGS) {
                        irxError("%s: too many arguments (max %d).\n",
                                routine, MAX_ARGS);
                        goto RETURN;
                }
        }
        va_end(ap);

        /* call routine */
        if (irxVectArrayFpRead(
                fp,
                numItems,
                numVect,
                varType,
                lptr) != SUCCESS)
                        goto RETURN;

        /* made it through */
        status = SUCCESS;
RETURN:
        if (status != 0) {
                irxError("%s: failed.\n", routine);
        }
        return(status);
#undef  MAX_ARGS
}

/*f----------------------------------------------------------------------
 * Returns the name of a type "type".
 */
char* irxVTypeName(IrxTVType type)
{
        switch (type) {
        case IRX_DOUBLE_T:               return("double");
        case IRX_FLOAT_T:                return("float");
        case IRX_INT_T:                  return("int");
        case IRX_LONG_T:                 return("long");
        case IRX_CHAR_T:                 return("char");

        case IRX_STRING_T:               return("char*");
        case IRX_CHAR_ARRAY_T:           return("char[]");
        case IRX_TDATE_T:                return("IrxTDate");
        case IRX_TDATEINTERVAL_T:        return("IrxTDateInterval");
        case IRX_TDAYCOUNT_T:            return("IrxTDayCountConv");

        case IRX_PERCENT_T:              return("percent");
        case IRX_CUR_T:                  return("currency");
        case IRX_CURK_T:                 return("currencyK");
        case IRX_CURM_T:                 return("currencyM");
        case IRX_BOOLEAN_T:              return("TBoolean");

        case IRX_CHAR_BLOCK_T:           return("char_block");

        case IRX_BPOINT_T:               return("bpoint");


        default:                return("ERROR");
        }
}



/*f----------------------------------------------------------------------
 * Returns SUCCESS if type "type" is a valid C type.
 */
int irxVTypeCheckValidCType(IrxTVType type)
{
        switch (type) {
        /* C types */
        case IRX_POINTER_T:
        case IRX_DOUBLE_T:
        case IRX_FLOAT_T:
        case IRX_INT_T:
        case IRX_LONG_T:
        case IRX_CHAR_T:
        case IRX_STRING_T:
        case IRX_TDATE_T:
        case IRX_TDATEINTERVAL_T:
        case IRX_TDAYCOUNT_T:
        case IRX_PERCENT_T:
        case IRX_CUR_T:
        case IRX_CURK_T:
        case IRX_CURM_T:
        case IRX_BOOLEAN_T:
        case IRX_CHAR_BLOCK_T:
        case IRX_BPOINT_T:
                return(SUCCESS);
        default:
                return(FAILURE);
        }
}



/*f----------------------------------------------------------------------
 * Returns the size of <i> type</i>.
 */
size_t irxVTypeSizeof(IrxTVType type)
{
static        char        routine[] = "irxVTypeSizeof";
        switch (type) {
        /* C types */
        case IRX_DOUBLE_T:
        case IRX_PERCENT_T:
        case IRX_CUR_T:
        case IRX_CURK_T:
        case IRX_CURM_T:
                return sizeof(double);
        case IRX_FLOAT_T:
                return sizeof(float);
        case IRX_INT_T:
        case IRX_BOOLEAN_T:
                return sizeof(int);
        case IRX_LONG_T:
                return sizeof(long);
        case IRX_CHAR_T:
                return sizeof(char);
        case IRX_STRING_T:
        case IRX_CHAR_BLOCK_T:
                return sizeof(char*);
        case IRX_TDATE_T:
                return sizeof(IrxTDate);
        case IRX_TDATEINTERVAL_T:
                return sizeof(IrxTDateInterval);
        case IRX_TDAYCOUNT_T:
                return sizeof(IrxTDayCountConv);
        default:
                irxError("%s: bad type %ld.\n", routine, (long) type);
                return((size_t) 0);
        }
}


/*f-------------------------------------------------------------
 * Scans a double value in a currency format
 * $\pm XXX,XXX,XXX.XX$ in the string <i> s</i>
 * and puts the result in <i> x</i>. 
 * If the string containts the percent sign, the output is 
 * divided by 100.
 * Returns 0 on success.
 */

int irxCurScan(char *s, double *x)
{
#undef        DELETECHAR
#define        DELETECHAR(s,p)                {char *q; for(q=p; *(q+1)!= '\0'; q++) \
                                *q=*(q+1); *q = '\0';}
        char        buf[256],
                *p;
        int        isPercent = FALSE;

        strcpy(buf, s);
        if (strchr(buf, '%')) isPercent = TRUE;

        for(p = buf; *(p+1) != 0; p++) {
                if (*p == ',') DELETECHAR(s, p);
        }

        if (sscanf(buf, "%lf", x) != 1) return (FAILURE);

        if (isPercent) *x *= 1e-2;
        return (SUCCESS);
}



/*f-------------------------------------------------------------
 * Prints a double value <i> x</i> in a currency format
 * $\pm XXX,XXX,XXX.XX$ in the string <i> s</i>.
 * If <i> s</i> is NULL, uses a static buffer and returns it.
 */

char* irxCurPrint(char *s, double x, int ndigits)
{
#undef        INSERTCHAR
#define       INSERTCHAR(s,p,c)  {char *q; for(q=s+(int)strlen(s);q>=p;q--)\
                                 *(q+1)=*q; *p = c;}
static        char          fmt[32];
              char          *p ;
              int           i;

/* buffers for printing */
#undef        BUF_IDX
#define       BUF_IDX       32
static        char          tmpBuf[BUF_IDX][64];
static        int           tmpBufIdx=0;

#undef        BUF_NEXT
#define       BUF_NEXT(str) str = (str != NULL ? str : tmpBuf[tmpBufIdx++]);\
                            tmpBufIdx = (tmpBufIdx > BUF_IDX-1 ? 0 : tmpBufIdx)


        /*s = (s == NULL ? buf : s);*/
        BUF_NEXT(s);

        if (ndigits <= 0) {
                sprintf(s, "%.0f", x);

                /* goto point */
                p = s + strlen(s);

        } else {
                sprintf(fmt, "%%.%df", ndigits);
                sprintf(s, fmt, x);

                /* goto point */
                p = s;
                while (*p != '.') p++;
        }


        if (p == s) goto RETURN;

        for(;;) {
            for (i=0; i<=2; i++) {
                p--;
                if ((p == s) || !(isdigit(*p))) goto RETURN;
            }
            if (!isdigit(*(p-1))) goto RETURN;

            if (isdigit(*(p-1))) INSERTCHAR(s, p, ',');

                /*if (((p-1) == s) || !(isdigit(*(p-1)))) goto RETURN;
                INSERTCHAR(s, p, ',');*/
        }


RETURN:
        return(s);

#undef        INSERTCHAR
#undef        BUF_IDX
#undef        BUF_NEXT
}





/*f-------------------------------------------------------------
 * This routine scans for strings
 * in a character string pointer by <i> p</i>.
 * The read strings can be enclosed in quotes (') or double quotes ("),
 * which are removed.
 * On exit, the result string is put in <i> arg</i> and returned (if <i> arg</i> 
 * is NULL, a static copy is returned).
 * The pointer <i> p</i> is advanced to the remaining unread portion 
 * the scanned string. It should not be changed between
 * successive calls to the routine.
 * For example,
 * \begin{verbatim}
 *          char    *p, *q, *string = "1.23 'hello world'";
 *          ...
 *          / * p points stringto be scanned * /
 *          p = &string[0];
 *          q = irxStrScanString(&p, NULL);        / * q is "1.23" * /
 *          / * p now points to the unread portion of string  * /
 *          ...
 *          q = irxStrScanString(&p, NULL);        / * q is "hello world" * /
 *          ...
 * \end{verbatim}
 */

char* irxStrScanString(char **p, char *arg)
{
const        char        EOS = '\0';                /* end of string */
static        char        arg1[256];
        char        *s, *q, c;

        arg = (arg != NULL ? arg : arg1);
        q = arg;
        s = *p;

        while (((c = *s++) != EOS) && (isspace(c)));

        if (c == EOS) return(NULL);


        if (c == '"') {
                while (((c = *s++) != EOS) && (c != '"'))
                        *q++ = c;
                if (c == EOS) return(NULL);
                *q = '\0';
        } else if (c == '\'') {
                while (((c = *s++) != EOS) && (c != '\''))
                        *q++ = c;
                if (c == EOS) return(NULL);
                *q = '\0';
        } else {
                *q++ = c;
                while (((c = *s++) != EOS) && (!isspace(c)))
                        *q++ = c;
                *q = '\0';
        }
        *p = s;
        return(arg);
}


/*f----------------------------------------------------------------------
 * Scans in a formatted string <i> str</i> a type <i> type</i> and
 * puts the result in <i> p</i>.\\
 * {\it Warning: for LIL types, <i> p</i> has to point to the 
 * right element of the array, that one accesses using 
 * <i> irxVTypeOffsetVect</i>, even if it is the first element.}
 * Returns 0 if scan successful.
 */

int irxVTypeScan(char *str, IrxTVType type, void *p)
{
static        char        routine[] = "irxVTypeScan";
        int     status    = FAILURE;

        char        *buf2, buf[256];

        switch (type) {
        case IRX_DOUBLE_T:
                if (irxCurScan(str, (double*)p) != SUCCESS) goto RETURN;
                break;
        case IRX_PERCENT_T:
                if (sscanf(str, "%lf", (double*)p) != 1) goto RETURN;
                *((double*)p) *= 1e-2;
                break;
        case IRX_FLOAT_T:
                if (sscanf(str, "%f",  (float*)p) != 1) goto RETURN;
                break;
        case IRX_INT_T:
                if (sscanf(str, "%d",  (int*)p) != 1) goto RETURN;
                break;
        case IRX_LONG_T:
                if (sscanf(str, "%ld",  (long*)p) != 1) goto RETURN;
                break;
        case IRX_CHAR_T:
                if (sscanf(str, "%s",  buf) != 1) goto RETURN;
                if (sscanf(buf, "%c",  (char*)p) != 1) goto RETURN;
                break;
        case IRX_STRING_T:
                if (sscanf(str, "%s",  ((char**)p)[0]) != 1) goto RETURN;
                break;
        case IRX_CHAR_ARRAY_T:
                strncpy(buf, str, sizeof(buf));
                buf2 = &buf[0];
                if (irxStrScanString(&buf2, ((char*)p)) == NULL)
                        goto RETURN;
                break;
        case IRX_TDATE_T:
                if (sscanf(str, "%s",  buf) != 1) goto RETURN;
                if (irxStringToDate(buf, (IrxTDate*)p) != 0) goto RETURN;
                break;
        case IRX_TDATEINTERVAL_T:
                if (sscanf(str, "%s",  buf) != 1) goto RETURN;
                if (irxStringToDateInterval(buf, (IrxTDateInterval*)p)
                                != SUCCESS) goto RETURN;
                break;
        case IRX_TDAYCOUNT_T:
                if (sscanf(str, "%s",  buf) != 1) goto RETURN;
                if (irxDayCountConvFromString(buf, (IrxTDayCountConv*)p) 
                                != SUCCESS) goto RETURN;
                break;

        case IRX_CUR_T:
                if (irxCurScan(str, (double*)p) != SUCCESS) goto RETURN;
                break;
        case IRX_CURK_T:
                if (irxCurScan(str, (double*)p) != SUCCESS) goto RETURN;
                *((double*)p) *= 1e3;
                break;
        case IRX_CURM_T:
                if (irxCurScan(str, (double*)p) != SUCCESS) goto RETURN;
                *((double*)p) *= 1e6;
                break;
        case IRX_BOOLEAN_T:
                /* T,F or integer value */
                if (sscanf(str, "%s",  buf) != 1) goto RETURN;
                if (toupper(buf[0]) == 'T') {
                        *((int*) p) = 1;
                } else if (toupper(buf[0]) == 'F') {
                        *((int*) p) = 0;
                } else {
                        if (sscanf(buf, "%d",  (int*)p) != 1)
                                goto RETURN;
                }
                break;
        case IRX_CHAR_BLOCK_T:
                if (sscanf(str, "%s",  ((char*)p)) != 1) goto RETURN;
                break;
        case IRX_BPOINT_T:
                if (sscanf(str, "%lf", (double*)p) != 1) goto RETURN;
                *((double*)p) *= 1e-4;
                break;
        default:
                irxError("%s: bad type %ld.\n", routine, (long) type);
                return(2);
        }

        status = SUCCESS; 

RETURN:
        if (status != SUCCESS)
        {
            irxError("%s: can't scan `%s' (type %s)\n",
                    routine, str, irxVTypeName(type));
        }

        return status;
}


/*f----------------------------------------------------------------------
 * Prints in a formatted string <i> str</i> a pointer <i> </i>p
 * of type <i> type</i>.  Returns <i> str</i>.\\
 * {\it Warning: for LIL types, <i> p</i> has to point to the 
 * right element of the array, that one accesses using 
 * <i> irxVTypeOffsetVect</i>, even if it is the first element.}
 */

char* irxVTypePrint(char *str, IrxTVType type, void *p)
{
static        char        routine[] = "irxVTypePrint";
static        char        buf[64];

        str = (str == (char*)NULL ? buf : str);

        if (p == NULL) {
                sprintf(str, "NULL");
                return(str);
        }

        switch (type) {
        case IRX_POINTER_T:
                sprintf(str, "%p", *((void**) p));
                break;
        case IRX_DOUBLE_T:
                sprintf(str, "%12.8f", *((double*) p));
                break;
        case IRX_FLOAT_T:
                sprintf(str, "%f",  *((float*) p));
                break;
        case IRX_INT_T:
                sprintf(str, "%d",  *((int*) p));
                break;
        case IRX_LONG_T:
                sprintf(str, "%ld",  *((long*) p));
                break;
        case IRX_CHAR_T:
                sprintf(str, "%c",  *((char*) p));
                break;
        case IRX_STRING_T:
                sprintf(str, "%s",  ((char**)p)[0]);
                break;
        case IRX_CHAR_ARRAY_T:
                sprintf(str, "%s",  ((char*)p));
                break;
        case IRX_TDATE_T:
                sprintf(str, "%10s", irxFormatDate(*((IrxTDate*) p)));
                break;
        case IRX_TDATEINTERVAL_T:
                sprintf(str, "%6s",
                        irxDateIntervalToString(((IrxTDateInterval*) p)));
                break;
        case IRX_TDAYCOUNT_T:
                sprintf(str, "%10s",
                        irxDayCountConvToString(*((IrxTDayCountConv*) p)));
                break;

        case IRX_PERCENT_T:
                sprintf(str, "%1.8f", *((double*) p) * 1e2);
                break;
        case IRX_CUR_T:
                irxCurPrint(str, *((double*) p), 0);
                break;
        case IRX_CURK_T:
                irxCurPrint(str, *((double*) p)*1e-3, 0);
                break;
        case IRX_CURM_T:
                irxCurPrint(str, *((double*) p)*1e-6, 0);
                break;
        case IRX_BOOLEAN_T:
                sprintf(str, "%s", (*((int*) p) != 0 ? "T" : "F"));
                break;

        case IRX_BPOINT_T:
                sprintf(str, "%12.8f", *((double*) p) * 1e4);
                break;

        default:
                irxError("%s: bad type %ld.\n", routine, (long)type);
                strcpy(str, "ERR");
                break;
        }

        return(str);
}


/*f-------------------------------------------------------------
 * I/O: read next line (skip comments).
 *
 * <br><br>
 * Acts like the standard ANSI routine <i> fgets</i>, except
 * that is discards all lines starting by '\#' (commented)
 * and counts number of lines: the counter <i> line</i>
 * is incremented every time a line is read (useful for
 * error reporting).
 */

static char*
irxFGetLine(char *s, int n, FILE *fp, int *line)
{
        int     acceptFlag = 0;
        char    *p;

        do {
                /* read next line */
                if (fgets(s, n, fp) == NULL) return(NULL);
                if (line != NULL) (*line)++;

                /* check for empty line (comment # ;) */
                p = s;
                while ((isspace(*p)) && (*p != '\0')) p++;
                if ((*p != '#') && (*p != ';') && (*p != '\0'))
                        acceptFlag = 1;
        } while (acceptFlag == 0);
        return(s);
}

