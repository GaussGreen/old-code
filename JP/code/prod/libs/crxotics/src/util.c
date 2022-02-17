/*********************************************************************************
 *
 * Misc utils
 *
 ********************************************************************************/
#include <string.h>

#include "transition.h"
#include "stdarg.h"
#include "cerror.h"
#include "srm123.h"

static int one = 1;

#define ERRBUFFSIZE    1024L        /* size of error message */
static char gErrMsg[ERRBUFFSIZE];


#if(defined(WIN32)||defined(WINNT)||defined(_WIN32))
#define vsnprintf   _vsnprintf
#endif



static DRDiagnosticMsgCB ErrorCallbackFnp = NULL;
static DRDiagnosticMsgCB WarningCallbackFnp = NULL;
static DRDiagnosticMsgCB TraceCallbackFnp = NULL;

/*********************************************************************************
 * set the error callback function as porposed by Jerry.
 *
 ********************************************************************************/
void CRXErrorSetCallBack(DRDiagnosticMsgCB errorCallbackFnp)
{
        /* Initialize error msg */
        strcpy (gErrMsg, "");
        ErrorCallbackFnp = errorCallbackFnp;
}

/*********************************************************************************
 * place holder for other additional functions proposed by Jerry.
 *
 ********************************************************************************/
void CRXWarningSetCallBack(DRDiagnosticMsgCB f)
{
        WarningCallbackFnp=f;
}
void CRXTraceSetCallBack(DRDiagnosticMsgCB f)
{
        TraceCallbackFnp=f;
}



/*********************************************************************************
 *   Print an error message.
 *
 ********************************************************************************/
void CRXError(const char* format, ...)
{
        int handled = 0;
        va_list args;

        va_start(args, format);
        if (ErrorCallbackFnp != NULL)
        {
                handled = (*ErrorCallbackFnp)(format, args);
                (*ErrorCallbackFnp)("\n", NULL);
        }

        if (!handled)
        {
                vfprintf(stdout,format,args);
                fprintf(stdout,"\n");
                fflush(stdout);
        }
        va_end(args);
}

/*********************************************************************************
 *       Store the error message in a string buffer
 *
 ********************************************************************************/
int CRXErrorCallback(const char* format, va_list ap)
{

    char    newErr[ERRBUFFSIZE];
    long    oldErrSize, newErrSizeMax;

    /* Existing error msg size */
    oldErrSize = strlen(gErrMsg);

    /* Max size of new error msg */
    newErrSizeMax = ERRBUFFSIZE - oldErrSize;

    /* Error bufffer full (including \0), no action */
    if (newErrSizeMax <= 1) 
        return TRUE;



    /* Print new error msg in the buffer, up to the full buffer including \0 
     */
    vsnprintf(newErr, newErrSizeMax, format, ap);

    /* Append to the existing error msg */
    strcat(gErrMsg, newErr);

    return TRUE;

}


/*********************************************************************************
 *       Retrieve the error message in a string buffer
 *
 ********************************************************************************/
const char* CRXErrorRetrieve()
{
    /* In case of truncation, separated by new line */
    if (gErrMsg[strlen(gErrMsg)-1] != '\n')
        gErrMsg[strlen(gErrMsg)-1] = '\n';

    return &gErrMsg[0];
}

static int  ExtendTCurve(T_CURVE *tc, long fromDate, long toDate);


/*********************************************************************************
 * Date routines : convert from YYYYMMDD long format.
 *                              
 * Converts a DR date type (encoded as YYYYMMDD in a long)
 * to a TDate. Returns -1 if failed.
 *
 ********************************************************************************/
TDate DrDate2TDate(long drDate)
{
    static  char routine[] = "DrDate2TDate";
    TDate       tDate;
    TMonthDayYear mdy;
    
    mdy.year  =  drDate/10000;            /* Truncation forced                  */
    mdy.month = (drDate - 10000*mdy.year)/100;
    mdy.day   = (drDate - 10000*mdy.year - 100*mdy.month);

    if (GtoMDYToDate(&mdy, &tDate) != SUCCESS) {
        CRXError("%s: failed.\n", routine);
        return(-1L);
    }
    return(tDate);
}   

/*********************************************************************************
 *      Given a DateList and its size
 *      sort its contents in 1st: date ascending order, 
 *                           2nd: value ascending order. (if not NULL)
 *      Returns SUCCESS or FAILURE
 ********************************************************************************/
int     CRXSortDateList(
    int      NbDates,                     /* (I) Num of dates                   */
    TDate   *DateList,                    /* (I/O) Date list                    */
    long    *SuppValue)                   /* (I) sort value                     */
{
    static char routine[] = "CRXSortDateList";
    
    int     status = FAILURE;
    TDate   tmpDate;
    long    tmpValue;  /* temp storage */
    int     i, j;
    long    *Val = NULL;

#undef  IS_GREATER
#define IS_GREATER(dateA, valA, dateB, valB)          \
                ((dateA > dateB) || ((dateA == dateB) && (valA > valB)))


    if (DateList == NULL) goto RETURN;
    
    if (SuppValue == NULL) Val = DateList;  /* trick to sort dates only */
                    else   Val = SuppValue;

    for (j = 1 ; j < NbDates; j++)
    {
        tmpDate = DateList[j];
        tmpValue = Val[j];

        i = j-1;
        while ((i >= 0) && 
               IS_GREATER(DateList[i],Val[i],tmpDate,tmpValue))
        {
            DateList[i+1] = DateList[i];
            Val[i+1]      = Val[i];
            i--;
        }

        DateList[i+1] = tmpDate;
        Val[i+1] = tmpValue;

    }  /* for j */

    status = SUCCESS;
RETURN:
    if(status == FAILURE)
    {
        CRXError("%s: Failed.\n", routine);
    }

    return status;
#undef IS_GREATER
}

/*********************************************************************************
 *      Merge and sort two date lists, redundant entries are removed
 *
 ********************************************************************************/
int CRXSortMergeDateList(
	int     *NbMerge,                     /* (O) Num of points in merged list   */
	TDate   *DLMerge,                     /* (O) Merged & sorted date list      */
	int      Nb1,                         /* (I) Num of points in 1st list      */
	TDate   *DL1,                         /* (I) 1st date list                  */
	int      Nb2,                         /* (I) Num of points in 2nd list      */
	TDate   *DL2)                         /* (I) 2nd date list                  */
{
    static char routine[] = "CRXSortMergeDateList";

	int    status = FAILURE;    
	int    i, idx;
	TDate  DLLocal[MAX_NB];


	/* merge */
	for(i = 0;i < Nb1;i++)
	{
		DLLocal[i] = DL1[i];
	}
	for(i = 0;i < Nb2;i++)
	{
		DLLocal[Nb1+i] = DL2[i];
	}
	*NbMerge = Nb1 + Nb2;
	
	/* sort */
	if(SortDateList(*NbMerge,
					DLLocal,
					NULL) == FAILURE)
	{
        goto RETURN;
	}

	/* remove redundant entries */

	DLMerge[0] = DLLocal[0];
	idx = 1;
	for(i = 1;i < Nb1 + Nb2;i++)
	{
		if(DLLocal[i] != DLLocal[i-1])
		{
			DLMerge[idx] = DLLocal[i];
			idx++;
		}
	}
	*NbMerge = idx;
	
	status = SUCCESS;
RETURN:
    if(status == FAILURE)
    {
        CRXError("%s: Failed.\n", routine);
    }

	return status;
}

/*********************************************************************************
 * Date routines : convert to YYYYMMDD long format.
 *                              
 * Converts a date in TDate type to DR date (encoded as YYYYMMDD in a long)
 *
 ********************************************************************************/
long TDate2DrDate(TDate tDate)
{
    static  char routine[] = "TDate2DrDate";
    long        drDate;
    TMonthDayYear   mdy;
    
    if (GtoDateToMDY(tDate, &mdy) != SUCCESS) {
        return(-1L);
    }
    
    drDate  = (long) mdy.year * 10000L;
    drDate += (long) mdy.month * 100L;
    drDate += (long) mdy.day;
    
    return(drDate);
}

/*********************************************************************************
 * Convert DR DCC to Alib DCC
 *
 ********************************************************************************/
int AlibDCC2DrDCC(
    char *dcc,                            /* (O) Dr DCC                         */
    long AlibDCC)                         /* (I) Alib DCC                       */
{
    static char routine[] = "AlibDCC2DrDCC";
    
    int      status = FAILURE;

    switch(AlibDCC)
    {
    case GTO_B30_360:
        *dcc = '3';
        break;
    case GTO_ACT_365F:
        *dcc = '5';
        break;
    case GTO_ACT_360:
        *dcc = '0';
        break;
    default:
        CRXError("%s: Unknown DCC.\n", routine);
        goto RETURN;
    }

    status = SUCCESS;
RETURN:
    if(status == FAILURE)
    {
        CRXError("%s: Failed.\n", routine);
    }

    return status;
}

/*********************************************************************************
 * Convert DR DCC to Alib DCC, return double
 *
 ********************************************************************************/
double AlibDCC2DrDCCDouble(long AlibDCC)                         /* (I) Alib DCC                       */
{
    static char routine[] = "AlibDCC2DrDCC";
    
    switch(AlibDCC)
    {
    case GTO_B30_360:
        return 3;
        break;
    case GTO_ACT_365F:
        return 5;
        break;
    case GTO_ACT_360:
        return 0;
        break;
    default:
        return -1;
    }
}

/*********************************************************************************
 * Convert Freq (char) to Alib Basis (double)
 *
 ********************************************************************************/
char   Basis2Freq(double Basis)
{

    switch((int)(Basis)){
    case 1:
        return 'A';
    case 2:
        return 'S';
    case 4:
        return 'Q';
    case 12:
        return 'M';
    default:
        return 'N';
    }
}

/*********************************************************************************
 *      Calculates a deterministic zero bond price
 *      Returns -999.99 if failure
 *
 *      Supports
 *      (ZeroInterpTypeFlag = 0) Linear Zero Cpn:
 *      -- linear zero cpn interp
 *      -- flat zero cpn extrapolation on both sides
 *
 *      (ZeroInterpTypeFlag = 1) Flat Fwd:
 *      -- flat fwd interp and extrapolation
 *
 ********************************************************************************/
static int  ZeroInterpTypeFlag = 1;       /* 0=Linear Zero Cpn; 1=Flat Fwd      */

static double   CRXZeroPrice(
    long     MatDate,                     /* (I) Mat date of the zero       */
    long     ValueDate,                   /* (I) Value date                 */
    int      NbZero,                      /* (I) Number of zeros in curve   */
    long    *ZeroDates,                   /* (I) maturity dates in zero crv */
    double  *ZeroRates)                   /* (I) Zero rates                 */
{
    static char routine[] = "CRXZeroPrice";
    
    double Price = -999.99;
    double t,  ZR, 
           t1, ZR1, Z1, 
           t2, ZR2, Z2;
    int    idx;

    /* basic checks */
    if ((ZeroDates == NULL) || (ZeroRates == NULL)) goto RETURN;
    if (NbZero <= 0) goto RETURN;
    if (MatDate == ValueDate) return(1.0);

    t   = Daysact(ValueDate, MatDate)/365.0;
    idx = GetDLOffset(NbZero, ZeroDates, MatDate, CbkHIGHER);

    if ((idx == 0) || (NbZero == 1)) /* MatDate <= 1st ZeroDate or crv has only 1 pt */
    {
        t1  = 0;
        ZR1 = ZeroRates[0];

        t2  = Daysact(ValueDate, ZeroDates[0])/365.0;
        ZR2 = ZeroRates[0];
    }
    else if (idx<0) /* i.e. all zero dates are < MatDate, flat fwd extrap */
    {
        switch (ZeroInterpTypeFlag)
        {
        case 0: /* Linear Zero Cpn */
            t1  = Daysact(ValueDate, ZeroDates[NbZero-1])/365.0;
            ZR1 = ZeroRates[NbZero-1];

            t2  = t;
            ZR2 = ZeroRates[NbZero-1];
            break;
        
        case 1: /* Flat Fwd */    
            t1  = Daysact(ValueDate, ZeroDates[NbZero-2])/365.0;
            ZR1 = ZeroRates[NbZero-2];

            t2  = Daysact(ValueDate, ZeroDates[NbZero-1])/365.0;
            ZR2 = ZeroRates[NbZero-1];
            break;        
        
        default:
            
            goto RETURN;
        }
    }
    else
    {
        t1  = Daysact(ValueDate, ZeroDates[idx-1])/365.0;
        ZR1 = ZeroRates[idx-1];

        t2  = Daysact(ValueDate, ZeroDates[idx])/365.0;
        ZR2 = ZeroRates[idx];
    }
    
    switch (ZeroInterpTypeFlag)
    {
    case 0: /* Linear Zero Cpn */
        if (linterp(t, &ZR,
                    t1,  t2,
                    ZR1, ZR2) == FAILURE) goto RETURN;
                    
        Price = pow(1.0 + ZR, -t);                    
        break;
        
    case 1: /* Flat Fwd */    
        if (IS_EQUAL(t1, t2)) goto RETURN;
        Z1  = pow(1.0 + ZR1, -t1); 
        Z2  = pow(1.0 + ZR2, -t2);
        Price = Z1 * pow(Z2/Z1, (t-t1)/(t2-t1));        
        break;
        
    default:
        goto RETURN;
    }

RETURN:

    if (Price < 0.0)
    {
        CRXError("%s: failed.", routine);
    }
    return (Price);

} 

/*********************************************************************************
 * Convert Dr T_CURVE to Alib TCurve
 *
 ********************************************************************************/
int TCurve2DrCurve(
    T_CURVE  *DrCurve,                    /* (O) Dr curve                       */
    TDate     today,                      /* (I) today                          */
    TCurve   *AlibCurve)                  /* (I) Alib curve                     */
{
    static char routine[] = "TCurve2DrCurve";
    
    int     status = FAILURE;
    int     i;
    
    if(GtoEncapsulatedCurveSetRates(AlibCurve) == FAILURE)
    {
        CRXError("Encapsulation failed.\n");
    }


    if(today == AlibCurve->fBaseDate)
    {
        DrCurve->Today = today;
    }

    DrCurve->Today = TDate2DrDate(AlibCurve->fBaseDate);
    DrCurve->ValueDate = DrCurve->Today;
    
    DrCurve->SpotDays = 0;
    
    if(AlibDCC2DrDCC(&DrCurve->CurveDCC,AlibCurve->fDayCountConv) == FAILURE)
    {
        goto RETURN;
    }
    
    DrCurve->CurveFreq = Basis2Freq(AlibCurve->fBasis);

    DrCurve->NbZero = AlibCurve->fNumItems-1;

    for(i = 1;i < DrCurve->NbZero+1;i++)
    {
        DrCurve->ZeroDate[i-1] = TDate2DrDate(AlibCurve->fArray[i].fDate);

        if(DrCurve->ZeroDate[i-1] < 0)
        {
            goto RETURN;
        }
        
        DrCurve->Zero[i-1]     = AlibCurve->fArray[i].fRate;
    }


    DrCurve->InterpType = SRM3_FLATFWD_INTERP; /* Warning! */

    if(ExtendTCurve(DrCurve,TDate2DrDate(today) ,0) != SUCCESS)
    {
        goto RETURN;
    }
    
    DrCurve->ValueDate = TDate2DrDate(today);

    status = SUCCESS;
RETURN:
    if(status == FAILURE)
    {
        CRXError("%s: Failed.\n", routine);
    }
    return status;
}

/*********************************************************************************
 * flat extend zero curve between fromDate and toDate dates
 * code modified from function InitYieldCurve() in irdiffuse.c
 *
 ********************************************************************************/
static int  ExtendTCurve(
    T_CURVE   *tc,                        /* (O) extended tcurve                */
    long      fromDate,                   /* (I) from date                      */
    long      toDate)                     /* (I) to date                        */
{
    int      i, nbZero = tc->NbZero;
    double   T2Vrate   = tc->Zero[0];
    double   stmZero;
    double   AZero, ARate;
    T_CURVE  tmpTc;
    
    /* no need to extend if value date equals today */
    if (tc->ValueDate == fromDate )
        return SUCCESS;
    
    if ( tc->ValueDate < fromDate )
    {
        CRXError("ExtendTCurve. Curve value date %ld < from date %ld",
                 tc->ValueDate, fromDate);
        return FAILURE;
    }
    
    stmZero   = pow((1.0 + T2Vrate), -Daysact(fromDate,tc->ValueDate)/365.0);
    
    /* keep a copy of the old zero curve */
    tmpTc = *tc;
    
    /* adjust zero rates to be from from date */
    for (i=0; i<nbZero; i++)
    {
        AZero = CRXZeroPrice(tmpTc.ZeroDate[i],
                             tmpTc.ValueDate,
                             tmpTc.NbZero,
                             tmpTc.ZeroDate,
                             tmpTc.Zero);

        if (AZero <= TINY)
        {
            CRXError("ExtendTCurve. AZero < TINY for date %ld",
                     tc->ZeroDate[i]);
            return FAILURE;
        }

        AZero *= stmZero;   /* convert it to a zero from fromDate */

        ARate = pow(AZero, -365.0/Daysact(fromDate, tc->ZeroDate[i])) - 1.0;

        tc->Zero[i] = ARate;
    }

    /* set today and value date to fromDate */
    tc->Today = tc->ValueDate = fromDate;

    if (toDate != 0 && tc->ZeroDate[tc->NbZero-1] < toDate)
    {
        tc->NbZero += 1;
        tc->ZeroDate[tc->NbZero-1] = toDate;
        tc->Zero[tc->NbZero-1] = tc->Zero[tc->NbZero-2];
    }
    return SUCCESS;
}


/*********************************************************************************
 * n-th eigenvector for Upper-Triangular matrix, given eigenval
 * which is (n,n) element of the matrix
 *
 ********************************************************************************/
int MatNthEigVecUT(
    double        *evec,                  /* (O) ptr to eigenvec                */
    int            idx,                   /* (I) index of the evec to calc      */
    double         eval,                  /* (I) eigenval specified             */
    Mat           *B,                     /* (I) matrix                         */
    double        *pEvec)                 /* (I) array to eigenvecs             */
{
    int                status  = FAILURE;
    char               routine[] = "MatNthEigVecUT";
    int                i, j;
    int                n = B->size;
    double             *ptr = evec;
    double             sum;
    int                istart;
    
    /* initialize */
    PtrSet(ptr,n,0);
    
    if(idx >= n){
        CRXError("%s: Mat idx too large.\n",routine);
        goto RETURN;
    }

    ptr[idx] = 1.0;
    if(idx == 0)
    {
        status = SUCCESS;
        goto RETURN;
    }
            
    for(i = idx - 1;i >= 0;i--)
    {
        istart = i * n;
            
        if(fabs(B->ptr[istart+i] - eval) < SMALL){
            PtrCopy(ptr,&pEvec[istart],n);
            status = SUCCESS;
            goto RETURN;
            
        } else {
            sum = 0.0;
            for(j = i+1;j <= idx;j++)
            {
                sum += B->ptr[istart + j] * ptr[j];
            }

            ptr[i] = -sum / (B->ptr[istart+i] - eval);
        }
    }
    
    /* normalize */
    sum = 0.0;
    for(i = 0;i < n;i++)
    {
        sum += ptr[i] * ptr[i];
    }
    sum = sqrt(sum);
    for(i = 0;i < n;i++)
    {
        ptr[i] /= sum;
    }
    
    status = SUCCESS;
RETURN:
    if(status == FAILURE)
    {
        CRXError("%s: Failed.\n", routine);
    }

    return status;
}

/*********************************************************************************
 * set size for mat
 *
 ********************************************************************************/
int MatSetSize(
    Mat    *A,                            /* (I/O) output matrix                */
    int    size)                          /* (I) size                           */
{
    A->size = size;
    A->size2 = size * size;

    return SUCCESS;
}

/*********************************************************************************
 * A = B * C
 *
 ********************************************************************************/
int MatMult(
    Mat    *A,                            /* (O) output matrix                  */
    Mat    *B,                            /* (I) input matrix                   */
    Mat    *C)                            /* (I) input matrix                   */
{
    static char routine[] = "MatMult";

    int     status = FAILURE;
    int     i, j, k;
    Mat     mat;
    double  dtmp;

    MatSetSize(&mat, B->size);
    
    if(B->size != C->size){
        CRXError("%s: Mat size not matched.\n", routine);
        goto RETURN;
    }

    for(i = 0;i < mat.size;i++)
    {
        for(j = 0;j < mat.size;j++)
        {
            dtmp = 0.0;
            for(k = 0;k < B->size;k++)
            {
                dtmp += B->ptr[i * C->size + k] * C->ptr[k * C->size + j];
            }
            mat.ptr[i * C->size + j] = dtmp;
        }
    }


    for(i = 0;i < mat.size2;i++)
    {
        A->ptr[i] = mat.ptr[i];
    }
    
    A->size = B->size;
    A->size2 = B->size2;    

    status = SUCCESS;
RETURN:
    if(status == FAILURE)
    {
        CRXError("%s: Failed.\n", routine);
    }
    return status;
}

/*********************************************************************************
 * scale an array
 *
 ********************************************************************************/
int PtrScale(
    double    *ptr,                       /* (I/O) array of doubles             */
    int        size,                      /* (I) size of ptr                    */
    double     scale)                     /* (I) scale factor                   */
{
    BLAS_DECL(dscal)(&size,&scale,ptr,&one);
    return SUCCESS;
}


/*********************************************************************************
 * copy one ptr to another
 *
 ********************************************************************************/
int PtrCopy(
    double    *to,                        /* (O) destination                    */
    double    *from,                      /* (I) source                         */
    int        size)                      /* (I) size of array                  */
{
    BLAS_DECL(dcopy)(&size,from,&one,to,&one);
    
    return SUCCESS;
}

/*********************************************************************************
 * A = B * C
 * Upper Triangle matrix multiplication only 
 *
 ********************************************************************************/
int MatMultUT(
    Mat    *A,                            /* (O) product                        */
    Mat    *B,                            /* (I) matrix B                       */
    Mat    *C)                            /* (I) matrix C                       */
{
    static char routine[] = "MatMultUT";
    
    int     status = FAILURE;
    Mat     mat;
    char    side = 'r';
    char    uplo = 'l';
    char    transa = 'n';
    char    diag = 'n';
    int     m, n, lda, ldb;
    double  alpha = 1.0;

    
    MatSetSize(&mat, B->size);
    
    m = C->size;
    n = C->size;
    lda = m;
    ldb = m;

    if(B->size != C->size){
        CRXError("%s: Mat size not matched.\n", routine);
        goto RETURN;
    }

    PtrCopy(mat.ptr,C->ptr,mat.size2);
    
    BLAS_DECL(dtrmm)(&side,&uplo,&transa,&diag,&m,&n,&alpha,B->ptr,&lda,mat.ptr,&ldb);
    
    PtrCopy(A->ptr,mat.ptr,mat.size2);
    A->size = B->size;
	A->size2 = B->size2;

    status = SUCCESS;
RETURN:
    if(status == FAILURE)
    {
        CRXError("%s: Failed.\n", routine);
    }
    return status;
}
            
/*********************************************************************************
 * X = B * Y
 * X: Nx1, B: NxN, Y: Nx1
 *
 ********************************************************************************/
int MatMultVec(
    double     *X,                        /* (O) vector Y                       */
    Mat        *B,                        /* (I) mat                            */
    double     *Y)                        /* (I) vector                         */
{
    double ptr[MAX_NB];
    int    i, j;

    for(i = 0;i < B->size;i++)
    {
        ptr[i] = Y[i];
    }

    
    for(i = 0;i < B->size;i++)
    {
        X[i] = 0.0;
        for(j = 0;j < B->size;j++)
        {
            X[i] += B->ptr[i * B->size + j] * ptr[j];
        }
    }

    return SUCCESS;
}

/*********************************************************************************
 * X = B * Y
 * X: Nx1, B: NxN, Y: Nx1
 *
 ********************************************************************************/
int MatMultVecUT(
    double    *X,                         /* (O) vector Y                       */
    Mat       *B,                         /* (I) mat                            */
    double    *Y)                         /* (I) vector                         */
{
    char    uplo = 'l';
    char    trans = 't';
    char    diag = 'n';
    int     one = 1;
    double  tmpY[MAX_NB];

    PtrCopy(tmpY,Y,B->size);
    
    BLAS_DECL(dtrmv)(&uplo,&trans,&diag,&(B->size),B->ptr,&(B->size),tmpY,&one);

    PtrCopy(X,tmpY,B->size);

    return SUCCESS;
}

/*********************************************************************************
 * X = Y * B
 * X: 1xN, B: NxN, Y: 1xN
 *
 ********************************************************************************/
int VecMultMat(
    double    *X,                         /* (O) vector Y                       */
    double    *Y,                         /* (I) vector                         */
    Mat       *B)                         /* (I) mat                            */
{
    double ptr[MAX_NB];
    int    i, j;

    for(i = 0;i < B->size;i++)
    {
        ptr[i] = Y[i];
    }
        
    for(i = 0;i < B->size;i++)
    {
        X[i] = 0.0;
        for(j = 0;j < B->size;j++)
        {
            X[i] += ptr[j] * B->ptr[j * B->size + i];
        }
    }

    return SUCCESS;
}

/*********************************************************************************
 * set an arry to a certain value
 *
 ********************************************************************************/
int PtrSet(
    double       *ptr,                    /* (O) array                          */
    int           size,                   /* (I) size of array                  */
    double        val)                    /* (I) value to be set to             */
{
    int i;
    for(i = 0;i < size;i++)
    {
        ptr[i] = val;
    }
    
    return SUCCESS;
}

/*********************************************************************************
 *     Computes exp(t*H), the matrix exponential of a general matrix in
 *     full, using the irreducible rational Pade approximation to the 
 *     exponential function exp(x) = r(x) = (+/-)( I + 2*(q(x)/p(x)) ),
 *     combined with scaling-and-squaring.
 *     adapted from expokit.
 *
 ********************************************************************************/
int MatExp(
    double       *expH,                   /* (O) output exp(tH)                 */
    int           ideg,                   /* (I) degree of diagonal Pade.       */
    double        t,                      /* (I) time-scale                     */
    double       *H,                      /* (I) matrix                         */
    int           m)                      /* (I) order of H                     */
{
    static char     routine[] = "MatExp";
    
    int      status = FAILURE;
    int      mm,i,j,k,ih2,ip,iq,iused,ifree,iodd,icoef,iput,iget;
    double   hnorm,scale,scale2,cp,cq;
    double   *wsp;
    int      *ipiv;
    char     cn = 'N';
    double   alpha, beta;
    int      iflag, ns;
    
    
    wsp   = (double *)malloc((4*m*m + ideg + 1)*sizeof(double));
    ipiv  = (int *)   malloc(m * sizeof(int));

    /* initialise pointers ... */
    mm = m*m;
    icoef = 0;
    ih2 = icoef + (ideg+1);
    ip  = ih2 + mm;
    iq  = ip + mm;
    ifree = iq + mm;

    /* scaling: seek ns such that ||t*H/2^ns|| < 1/2; 
       and set scale = t/2^ns ... */

    for(i = 0;i < m;i++)
    {
        wsp[i] = 0.0;
        for(j = 0;j < m;j++)
        {
            wsp[i] += fabs(H[i * m + j]);
        }
    }
    
    hnorm = 0.0;

    for(i = 0;i < m;i++)
    {
        hnorm = MAX(hnorm, wsp[i]);
    }
    
    hnorm = fabs( t*hnorm );
    ns = 0;
    if(hnorm > SMALL){    
        ns = MAX( 0,(int)(log(hnorm)/log(2.0))+2 );
    }
    scale = t / pow(2,ns);
    scale2 = scale*scale;
    
    /*---  H2 = scale2*H*H ... */
    beta = 0.0;
    BLAS_DECL(dgemm)(&cn,&cn,&m,&m,&m,&scale2,H,&m,H,&m,&beta,&wsp[ih2],&m);
    
    /*---  compute Pade coefficients ... */
    i = ideg+1;
    j = 2*ideg+1;
    wsp[icoef] = 1.0;
    for(k = 1;k <= ideg;k++)
    {
        wsp[icoef + k] = wsp[icoef+k-1] * (double)(i-k)/(double)(k*(j-k));
    }

    /*---  initialize p (numerator) and q (denominator) ... */
    cp = wsp[icoef+ideg-1];
    cq = wsp[icoef+ideg];

    for(j = 0;j < m;j++)
    {
        for(i = 0;i < m;i++)
        {
            wsp[ip + j * m + i] = 0.0;
            wsp[iq + j * m + i] = 0.0;
        }
        
        wsp[ip + j * (m+1)] = cp;
        wsp[iq + j * (m+1)] = cq;
    }

    /*---  Apply Horner rule ... */
    iodd = 1;
    k = ideg - 1;
    do {
        iused = iodd*iq + (1-iodd)*ip;

        alpha = 1.0; beta = 0.0;
        BLAS_DECL(dgemm) (&cn,&cn,&m,&m,&m,&alpha,&wsp[ih2],&m,&wsp[iused],&m,&beta,&wsp[ifree],&m);

        for(j = 0;j < m;j++)
        {
            wsp[ifree+j*(m+1)] += wsp[icoef+k-1];
        }

        ip = (1-iodd)*ifree + iodd*ip;
        iq = iodd*ifree + (1-iodd)*iq;
        ifree = iused;
        iodd = 1-iodd;
        k = k-1;
    } while(k > 0);
    
    /*---  Obtain (+/-)(I + 2*(p\q)) ... */
    if(iodd == 1)
    {
        beta = 0.0;
        BLAS_DECL(dgemm)(&cn,&cn,&m,&m,&m,&scale,H,&m,&wsp[iq],&m,&beta,&wsp[ifree],&m);
        iq = ifree;
    } else {
        beta = 0.0;
        BLAS_DECL(dgemm) (&cn,&cn,&m,&m,&m,&scale,H,&m,&wsp[ip],&m,&beta,&wsp[ifree],&m);
        ip = ifree;
    }

    alpha = -1.0;
    BLAS_DECL(daxpy)(&mm,&alpha,&wsp[ip],&one,&wsp[iq],&one);
    BLAS_DECL(dgesv)(&m,&m,&wsp[iq],&m,ipiv,&wsp[ip],&m,&iflag);
    if(iflag != 0){
        CRXError("%s: failed.\n",routine);
        goto RETURN;
    }

    alpha = 2.0;
    BLAS_DECL(dscal)(&mm,&alpha,&wsp[ip],&one);

    for(j = 0;j < m;j++)
    {
        wsp[ip + j *(m+1)] += 1.0;
    }
    iput = ip;
    if(ns == 0 && iodd == 1)
    {
        alpha = -1.0;
        BLAS_DECL(dscal)(&mm,&alpha,&wsp[ip],&one);
    } else {
        /*--   squaring : exp(t*H) = (exp(t*H))^(2^ns) ...*/
        iodd = 1;
        for(k = 0;k < ns;k++)
        {
            iget = iodd*ip + (1-iodd)*iq;
            iput = (1-iodd)*ip + iodd*iq;
            
            alpha = 1.0; beta = 0.0;
            BLAS_DECL(dgemm)(&cn,&cn,&m,&m,&m,&alpha,&wsp[iget],&m,&wsp[iget],&m,&beta,&wsp[iput],&m);
            iodd = 1-iodd;
        }
    }

    PtrCopy(expH,&wsp[iput],m*m);

    status = SUCCESS;
RETURN:
    if(ipiv != NULL) free(ipiv);
    if(wsp  != NULL) free(wsp);

    if(status == FAILURE)
    {
        CRXError("%s: Failed.\n", routine);
    }

    return status;
}

/*********************************************************************************
 *      Smooth version of the max(a,b)
 *      This is the based on the integral of the DrSmoothStep function
 *      Returns b                   if a <= b - step
 *              'smooth curve'      if b - step < a < b + step
 *              a                   if a >= b + step
 *
 ********************************************************************************/
double  CRXSmoothMax(
    double    a,                          /* (I) Argument 1                     */
    double    b,                          /* (I) Argument 2                     */
    double    Step)                       /* (I) smooth in [-step, step]        */
{
    double smoothValue, x, xSquare;

    Step = ABS(Step);

    /* do the limiting case when step = 0.0 */
    if (Step < TINY)
    {
        smoothValue = (a < b) ? b : a;
        return (smoothValue);
    }

    /* do the smooth case */
    x = (a - b)/Step;

    if (x <= -1.0)
    {
        smoothValue = b;
    }
    else if (x >= 1.0)
    {
        smoothValue = a;
    }
    else
    {
        xSquare = x * x;
        smoothValue = ((xSquare-5.) * xSquare + 15.) * xSquare * 0.03125 + 
                      0.5 * x + 0.15625;
	    smoothValue = Step * smoothValue + b;
    } 
    return (smoothValue);
}



/* Coefficients for inverse cumulative normal */
#define   a0     2.50662823884
#define   a1   -18.61500062529
#define   a2    41.39119773534
#define   a3   -25.44106049637
#define   b0    -8.47351093090
#define   b1    23.08336743743
#define   b2   -21.06224101826
#define   b3     3.13082909833
#define   c0     0.3374754822726147
#define   c1     0.9761690190917186
#define   c2     0.1607979714918209
#define   c3     0.0276438810333863
#define   c4     0.0038405729373609
#define   c5     0.0003951896511919
#define   c6     0.0000321767881768
#define   c7     0.0000002888167364
#define   c8     0.0000003960315187

/* coefficients for double precision cumulative normal -- courtesy of ALIB */

#define SQRT2   1.414213562373095049     /* sqrt(2) */
#define SQRTPI  1.772453850905516027     /* sqrt(pi) */

/* Coefficients in expression of erf(x) for -0.46875<=x<=0.46875 */
#define P10 3209.377589138469472562    /* Numerator */
#define P11 377.4852376853020208137
#define P12 113.8641541510501556495
#define P13 3.161123743870565596947
#define P14 0.1857777061846031526730
#define Q10 2844.236833439170622273   /* Denominator */
#define Q11 1282.616526077372275645
#define Q12 244.0246379344441733056
#define Q13 23.60129095234412093499
#define Q14 1.0

/* Coefficients in expression of erfc(x) for 0.46875<=x<=4.0 */
#define P20 1230.33935479799725272  /* Numerator */
#define P21 2051.07837782607146532
#define P22 1712.04761263407058314
#define P23 881.952221241769090411
#define P24 298.635138197400131132
#define P25 66.1191906371416294775
#define P26 8.88314979438837594118
#define P27 0.564188496988670089180
#define P28 2.15311535474403846343e-8
#define Q20 1230.33935480374942043  /* Denominator */
#define Q21 3439.36767414372163696
#define Q22 4362.61909014324715820
#define Q23 3290.79923573345962678
#define Q24 1621.38957456669018874
#define Q25 537.181101862009857509
#define Q26 117.693950891312499305
#define Q27 15.7449261107098347253
#define Q28 1.0

/* Coefficients in expression of erfc(x) for x>= 4.0 */
#define P30 -6.58749161529837803157E-4    /* Numerator */
#define P31 -1.60837851487422766278E-2
#define P32 -1.25781726111229246204E-1
#define P33 -3.60344899949804439429E-1
#define P34 -3.05326634961232344035E-1
#define P35 -1.63153871373020978498E-2
#define Q30  2.33520497626869185443E-3    /* Denominator */
#define Q31  6.05183413124413191178E-2
#define Q32  5.27905102951428412248E-1
#define Q33  1.87295284992346047209
#define Q34  2.56852019228982242072
#define Q35  1.0

/*f----------------------------------------------------------------------------
 * Cumulative error function weighted by exponential factor. Computes
 * exp(b+a^2/2)erfc(x-a)
 *
 * Alib comment:  The routine has a relative accuracy no worse than 1.0E-14, 
 * where relative accuracy is defined as (computed - truth)/truth, and truth
 * comes from a continued fraction calculation.  This is essentially 
 * accurate to the next to last decimal digit of machine accuracy on the Sun.
 */
double ExpCErrFcn (double a, double b, double x)
{
    double numerator;            /* numerator of polynomial in expression */
    double denominator;          /* denominator of polynomial in expression */
    double y,y2;                 /* y = abs(x)/sqrt(2), y2 = y*y */
    double erfc;                 /* return value */

    double W = b + 0.5 * a * a;

    y  = fabs(x-a) / SQRT2;
    y2 = y * y;

    if (y < 0.46875) 
    {
        numerator   = P10 + y2*(P11 + y2*(P12 + y2*(P13 +y2*P14)));
        denominator = Q10 + y2*(Q11 + y2*(Q12 + y2*(Q13 +y2*Q14)));
        erfc = exp(W) * (1 - y * numerator / denominator);
    }
    else if (y < 4.0) 
    {
        numerator   = P20 + y*(P21 + y*(P22 + y*(P23 +
                                                 y*(P24 + y*(P25 + y*(P26 + y*(P27 + y*P28)))))));
        denominator = Q20 + y*(Q21 + y*(Q22 + y*(Q23 +
                                                 y*(Q24 + y*(Q25 + y*(Q26 + y*(Q27 + y*Q28)))))));
        erfc = exp(-y2 + W) * numerator / denominator;
    }
    else /* (y > 4.0) */ 
    {
        double z2 = 1/y2; 
        numerator   = P30 + z2*(P31 + z2*(P32 + z2*(P33 + z2*(P34 +z2*P35))));
        denominator = Q30 + z2*(Q31 + z2*(Q32 + z2*(Q33 + z2*(Q34 +z2*Q35))));
        erfc = (exp(-y2 + W)/y) * 
            (1.0 / SQRTPI + numerator / (denominator * y2)); 
    }
    
    return erfc;
} /*ExpCerrFcn */

/*f----------------------------------------------------------------------------
 * Double precision cumulative normal function
 */
double NormCum (double x)
{
    return (x>0.0) ? 1. - 0.5 * ExpCErrFcn(0,0,x) : 0.5 * ExpCErrFcn(0,0,x);
}


/*f----------------------------------------------------------------------------
 * Inverse cumulative normal distribution
 *
 * Based on Risk Magazine
 */
double NormCumInv (double prob)
{
    double  t;
    double  x;
    double  r;
   
    t = (prob < 0.5) ? (1. - prob) : prob;
   
    x = t - 0.5;
    if (fabs (x) < 0.42) 
    {
        r = x * x;
        r = x * (((a3 * r + a2) * r + a1) * r + a0) 
            / ((((b3 * r + b2) * r + b1) * r + b0) * r + 1.);
        return (prob < 0.5) ? -r : r;
    } 
    else 
    {
        r = t;
        if (x > 0.) r = 1. - t;
        r = log (- log (r));
        r = c0 + r * (c1 + r * (c2 + r * (c3 + r * (c4 + r * (c5 + r 
                                                              * (c6 + r * (c7 + r * c8)))))));
        if (x < 0.) r = -r;
        return (prob < 0.5) ? -r : r;
    }
} /* NormCumInv */


/*********************************************************************************
 * Returns the value ln [gamma(xx ]for xx > 0 .
 *
 ********************************************************************************/
static double gammln(double xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
                          24.01409824083091,-1.231739572450155,
                          0.1208650973866179e-2,-0.5395239384953e-5};
    int j;
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

/*********************************************************************************
 * Used by betai Evaluates continued fraction for incom lete beta function by
 * modified Lentz
 *
 ********************************************************************************/
#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
static double betacf(double a, double b, double x)
{
    int m,m2;
    double aa,c,d,del,h,qab,qam,qap;
    
    qab=a+b;
    qap=a+1.0;
    qam=a-1.0;
    c=1.0; 
    d=1.0-qab*x/qap;
    if (fabs(d) < FPMIN) d=FPMIN;
    d=1.0/d;
    h=d;
    for (m=1;m<=MAXIT;m++) {
        m2=2*m;
        aa=m*(b-m)*x/((qam+m2)*(a+m2));
        d=1.0+aa*d; 
        if (fabs(d) < FPMIN) d=FPMIN;
        c=1.0+aa/c;
        if (fabs(c) < FPMIN) c=FPMIN;
        d=1.0/d;
        h *= d*c;
        aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
        d=1.0+aa*d; 
        if (fabs(d) < FPMIN) d=FPMIN;
        c=1.0+aa/c;
        if (fabs(c) < FPMIN) c=FPMIN;
        d=1.0/d;
        del=d*c;
        h *= del;
        if (fabs(del-1.0) < EPS) break; 
    }
    if (m > MAXIT) CRXError("a or b too big, or MAXIT too small in betacf");
    return h;
}

/*********************************************************************************
 * Returns the value of the beta function B(z, w)
 *
 ********************************************************************************/
double Beta(double z, double w)
{
    return exp(gammln(z)+gammln(w)-gammln(z+w));
}

/*********************************************************************************
 * Returns the value of the beta dist function B(z, w) pdf
 *
 ********************************************************************************/
double BetaPDF(double x, double z, double w)
{
    return exp( (z - 1.0) * log(x) + (w - 1.0) * log(1.0 - x)
                - gammln(z) - gammln(w) + gammln(z+w));
}

/*********************************************************************************
 * Returns the incom lete beta function I x (a, b).
 *
 ********************************************************************************/
double Betai(double x, double a, double b)
{
    double bt;

    if (x < 0.0 || x > 1.0)
        CRXError("Bad x in routine betai");
    if (x == 0.0 || x == 1.0)
        bt=0.0;
    
    else                                  /* Factors in front of the continued fraction. */
        bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));

    
    if (x < (a+1.0)/(a+b+2.0))            /* Use continued fraction directly. */
        return bt*betacf(a,b,x)/a;
    else
        return 1.0-bt*betacf(b,a,1.0-x)/b;
}
