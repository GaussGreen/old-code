/*********************************************************************************
 * CRXUTIL.C
 * crx utils
 *
 ********************************************************************************/

#include <string.h>
#include <math.h>

#include "crxutil.h"
#include "crcrv.h"



/*f-------------------------------------------------------------
 * Creates and returns a TDate corresponding to
 * month "m", day "d" and year "y".
 */

TDate
CrxTDateMake(int m, int d, int y)
{
        TMonthDayYear   mdy;
        TDate           theDate;

        mdy.month = (long) m;
        mdy.day = (long) d;
        mdy.year = (long) y;

        if (GtoMDYToDate(&mdy, &theDate) != 0) {
                return 0;
        } else {
                return theDate;
        }

}



/********************************************************************************
 * Date routines : convert from YYYYMMDD long format.
 *                              
 * Converts a DR date type (encoded as YYYYMMDD in a long)
 * to a TDate. Returns -1 if failed.
 ********************************************************************************/
TDate CrxDrDate2TDate(long drDate)
{
    static  char routine[] = "CrxDrDate2TDate";
    TDate       tDate;
    TMonthDayYear mdy;
    
    mdy.year  =  drDate/10000;                /* Truncation forced      */
    mdy.month = (drDate - 10000*mdy.year)/100;
    mdy.day   = (drDate - 10000*mdy.year - 100*mdy.month);

    if (GtoMDYToDate(&mdy, &tDate) != SUCCESS) {
        DR_Error("%s: failed.\n", routine);
        return(-1L);
    }
    return(tDate);
}   


/********************************************************************************
 * Date routines : convert to YYYYMMDD long format.
 *                              
 * Converts a date in TDate type to DR date (encoded as YYYYMMDD in a long)
 * Returns -1 if failed.
 ********************************************************************************/

long CrxTDate2DrDate(TDate tDate)
{
    static  char routine[] = "CrxTDate2DrDate";
    long        drDate;
    TMonthDayYear   mdy;
    
    if (GtoDateToMDY(tDate, &mdy) != SUCCESS) {
        DR_Error("%s: failed.\n",routine);
        return(-1L);
    }
    
    drDate  = (long) mdy.year * 10000L;
    drDate += (long) mdy.month * 100L;
    drDate += (long) mdy.day;
    
    return(drDate);
}

/********************************************************************************
 * Convert Freq (char) to Alib Basis (double)
 ********************************************************************************/
double   CrxFreq2Basis(char Freq)
{
    double      basis;

    switch(Freq){
    case 'A':
        basis = 1;
        break;
    case 'S':
        basis = 2;
        break;
    case 'Q':
        basis = 4;
    case 'M':
        basis = 12;
        break;
    default:
        DR_Error("Error converting freq to basis.\n");
        return -1;
    }
    
    return basis;
}


/********************************************************************************
 * Convert Freq (char) to Alib Basis (double)
 ********************************************************************************/
char   CrxBasis2Freq(double Basis)
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

/********************************************************************************
 * Convert DR DCC to Alib DCC
 ********************************************************************************/
int AlibDCC2DrDCC(char *dcc, long AlibDCC)
{
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
        goto RETURN;
    }

    status = SUCCESS;
RETURN:
    return status;
}

/********************************************************************************
 * Convert DR DCC to Alib DCC
 ********************************************************************************/
long CrxDr2AlibDCC(char DrDCC)
{
    long AlibDCC = FAILURE;
    switch(DrDCC)
    {
    case '3':
        AlibDCC = GTO_B30_360;
        break;
    case '5':
        AlibDCC = GTO_ACT_365F;
        break;
    case '0':
        AlibDCC = GTO_ACT_360;
        break;
    default:
        DR_Error("Error converting DCC.\n");
    }

    return AlibDCC;
}

/********************************************************************************
 * Convert Dr T_CURVE to Alib TCurve
 ********************************************************************************/
int CrxDrCurve2TCurve(TCurve  **AlibCurve,
                      T_CURVE * DrCurve)
{
    int         status = FAILURE;
    int         i;
    TDate       TBaseDate;
    TDate       TZeroDate[MAX_ZRATE];
    TCurve      *tcurve = NULL;

    TBaseDate = CrxDrDate2TDate(DrCurve->ValueDate);

    for(i = 0;i < DrCurve->NbZero;i++){
        TZeroDate[i] = CrxDrDate2TDate(DrCurve->ZeroDate[i]);
    }
    
    tcurve = GtoMakeTCurve(TBaseDate,
                           TZeroDate,
                           DrCurve->Zero,
                           DrCurve->NbZero,
                           CrxFreq2Basis(DrCurve->CurveFreq),
                           CrxDr2AlibDCC(DrCurve->CurveDCC));

    if(!tcurve){
        DR_Error("Conversion to Alib TCurve failed.\n");
        goto RETURN;
    }

    
    *AlibCurve = tcurve;

    status = SUCCESS;
 RETURN:
    return status;
}



/*****  CrxSortDateList  ******************************************/
/*
 *      Given a DateList and its size
 *      sort its contents in 1st: date ascending order, 
 *                           2nd: value ascending order. (if not NULL)
 *      Returns SUCCESS or FAILURE
 */

int     CrxSortDateList(int    NbDates, 
                        long  *DateList,
                        long  *SuppValue)
{
    TDate   tmpDate;
    long   tmpValue;  /* temp storage */
    int    i, j;
    long  *Val = NULL;

#undef  IS_GREATER
#define IS_GREATER(dateA, valA, dateB, valB)          \
                ((dateA > dateB) || ((dateA == dateB) && (valA > valB)))


    if (DateList == NULL) return (FAILURE);
    
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

    return (SUCCESS);

#undef IS_GREATER

} /* CrxSortDateList */




int CrxSortMergeDateList(
	int   *NbMerge,				         /* Num of points in merged list   */
	long  *DLMerge,					 /* Merged & sorted date list      */
	int   Nb1,							 /* Num of points in 1st list      */
	long  *DL1,						 /* 1st date list                  */
	int   Nb2,							 /* Num of points in 2nd list      */
	long  *DL2)						 /* 2nd date list                  */
{
	int   i, idx;
	TDate  DLLocal[MAXPOINTS];
	int   status = FAILURE;

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
	if(CrxSortDateList(*NbMerge,
					DLLocal,
					NULL) == FAILURE)
	{
		goto FREE_MEM_AND_RETURN;
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
 FREE_MEM_AND_RETURN:
	return status;
}



int CrxTDateLinearInterp1d(
    TDate *xa,  /* (I) array of X(i) (i=0,..,m-1) */
    double *ya, /* (I) array of Y(i) (i=0,..,m-1) */
    int m,      /* (I) arrays size */
    TDate x,    /* (I) point to intepolate */
    double *y)  /* (O) interpolated value */
{
static  char    routine[] = "CrxLinearInterp1d";

    int i;
    double  u;

    if (m < 1) {
        DR_Error("%s: not enough points (%d).\n",
            routine, m);
        return(1);
    }

    if (m == 1) {
        *y = ya[0];
        return(0);
    }

    if (x <= xa[0]) {
        *y = ya[0];
        return(0);
    }
    if (x >= xa[m-1]) {
        *y = ya[m-1];
        return(0);
    }


    for (i=0; (xa[i] <= x) && (i<=m-1) ; i++);
    i--;
    if (i > m-1) i--;

    u = ((double)x - (double)xa[i]) / ((double)xa[i+1] - (double)xa[i]) ;

    *y =    (1.0-u) * ya[i  ] +
        u   * ya[i+1] ;

    return(SUCCESS);
}

