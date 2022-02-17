/****************************************************************************/
/* Purpose */
/****************************************************************************/

/****************************************************************************/
/* EXTERNAL ROUTINES USED */
/****************************************************************************/
/*

*/
/****************************************************************************/
/* Headers */
/****************************************************************************/
#include "SRT_H_ALL.H>
#include "opfnctns.H>
#include "srt_h_lgmtypes.h"
#include "srt_h_lgmprotos.h"
#include "srt_h_lgmUSprotos.h"
#include "math.h"

/* Add a maximum number of periods = 100 Year of maturity */
static int LGMUTIL_MAXNB_PERIOD	= 210;

/****************************************************************************/
/* Prototypes of private functions */
/****************************************************************************/
/* Gets sigma(t) = sigma[j] for t[j]<t<t[j+1] for j=0,..., ntau-1
from zeta[j] = zeta(t) at t=tau[j] for j=0,..., numS 
and kap[i] = kappa(t) for tkap[i]<t<t[i+1] */
static void GetSigFromZeta(double g0, long ntau, double *sig,
						   Date *tau, double *zeta,
						   long numK, Date *t, double *kap);

/* compute integral of exponential and its derivative */
static double expint(double arg);
static double expintder(double arg);

/* Gets kappa(t) = kappa[i] for t[i]<t<t[i+1] for i=0,..., n-1
from G[i] = G(t) at t=t[i] for i=0,..., n */
static void GetKapFromG(double b, long n, Date *t, double *G, double *kap);

/* Use Newton to solve (exp(theta)-1)/theta = ratio */
static double FindTheta(double ratio);

/* compute the time-average value of kappa(t) */
static double computeavgkap(long n, Date *t, double *kap);

/* find largest deviation from the average */
static long GetWorstDif(long n, double *Arr, double avg);

static double GetDeltaLogb(long n, double *Arr, Date *tArr, double avg, long iw);

/****************************************************************************/
/**** Comments ****/

/*******************************************/
/******* Utilities for Deal structures *****/
/*******************************************/
/* Create an array of (exer date, end date, fixed rate) for
conveniently passing exer boundary info */ 
/* Returns pointer to structure; on failure returns NULL */
LGMSwptnsPtr LGMCreateSwptns(long nSwptns)
{	LGMSwptns *ptr = NULL;
	long n;
	
	n = nSwptns;
	if (n<1)
		n=1;
	ptr = (LGMSwptnsPtr) srt_calloc(1, sizeof(LGMSwptns));
	if (ptr==NULL)
		return(ptr);

	ptr->n = nSwptns;
	ptr->tEx = (Date*) srt_calloc(n, sizeof(Date));
	ptr->tEnd = (Date*) srt_calloc(n, sizeof(Date));
	ptr->Rfix = (double*) srt_calloc(n, sizeof(double));

	if (ptr->tEx==NULL ||ptr->tEnd==NULL || ptr->Rfix==NULL)
	{	LGMFreeSwptns(&ptr);
		return(NULL);
	}
	return(ptr);
}

/*******************************************/
/* Free a LGMSwptn structure */
void LGMFreeSwptns(LGMSwptns **ptr)
{	if ((*ptr) == NULL) return;

	if( (*ptr)->tEx != NULL) srt_free((*ptr)->tEx);
	if( (*ptr)->tEnd != NULL) srt_free((*ptr)->tEnd);
	if( (*ptr)->Rfix != NULL) srt_free((*ptr)->Rfix);

	(*ptr)->tEx = NULL;
	(*ptr)->tEnd = NULL;
	(*ptr)->Rfix = NULL;

	srt_free((*ptr));
	(*ptr) = NULL;
	return;
}

/*******************************************/
/* Free contents of SrtLgmRefSwptnData structure */
void LGMFreeRefSwptnData(SrtLgmRefSwptnData *ptr)
{	if (ptr==NULL) return;

	if(ptr->refShortSwptnArr!=NULL) srt_free(ptr->refShortSwptnArr);
	if(ptr->refLongSwptnArr!=NULL) srt_free(ptr->refLongSwptnArr);
	if(ptr->refCapSwptnArr!=NULL) srt_free(ptr->refCapSwptnArr);
	if(ptr->refFixSwptnArr!=NULL) srt_free(ptr->refFixSwptnArr);

	ptr->refShortSwptnArr = NULL;
	ptr->refLongSwptnArr = NULL;
	ptr->refCapSwptnArr = NULL;
	ptr->refFixSwptnArr = NULL;

	ptr->NrefShortSwptn = 0;	
	ptr->NrefLongSwptn = 0;
	ptr->NrefCapSwptn = 0;
	ptr->NrefFixSwptn = 0;
	ptr->isSemi = 0;	

	return;
}

/*******************************************/
/* Create a General European stucture with nPay payments */ 
/* Returns pointer to structure; on failure returns NULL */
SrtGenEurPtr LGMCreateGenEur(long nPay)
{	SrtGenEur *ptr;
	long np;

	np = nPay;					/* minimum array size */
	if (np<1)
		np=1;
	ptr = (SrtGenEur*) srt_calloc(1, sizeof(SrtGenEur));
	if (ptr==NULL)
		return (ptr);

	ptr->nPay = nPay;
	ptr->Payment = (double*) srt_calloc(np, sizeof(double));
	ptr->tPay = (Date*) srt_calloc(np, sizeof(Date));
	if (ptr->Payment==NULL || ptr->tPay==NULL)
	{	LGMFreeGenEur(&ptr);
		return (NULL);
	}
	return (ptr);
}

/*******************************************/
/* Frees the general European structure pointed to by *ptr */
void LGMFreeGenEur(SrtGenEur **ptr)
{	if ((*ptr)==NULL) return;

	if ((*ptr)->Payment !=NULL) srt_free((*ptr)->Payment);
	if ((*ptr)->tPay !=NULL) srt_free((*ptr)->tPay);

	(*ptr)->Payment = NULL;
	(*ptr)->tPay = NULL;

	srt_free((*ptr));
	*ptr = NULL;

	return;
}

/*******************************************/
/* Check validity of general European deal */
LGMErr LGMValidGenEur(SrtGenEur *ptr, Date tfirst, long *nExleft)
{	long i, nPay;
	
/* Check exercise */
	if (ptr->tEx > tfirst+11000)
		return ("No valid exercise date");

	*nExleft=1;
	if (ptr->tEx < tfirst)
		*nExleft = 0;				/* exercise date has passed */

/* Check pay dates */
	nPay = ptr->nPay;
	if (nPay<1)
		return ("Deal has no payments");

	if (ptr->tPay[0] < ptr->tEx)
		return("Pay date before exercise");

	for (i=1; i<nPay; i++)
	{	if (ptr->tPay[i-1] > ptr->tPay[i])
			return("Pay dates out of order");
	}	
	return (NULL);
}

/*******************************************/
/* Create a simple MidAtlantic stucture with nEx exercise dates and nPay payments */ 
/* Returns pointer to structure; on failure returns NULL */
SrtSimMidAtPtr LGMCreateSimMidAt(long nEx, long nPay)
{	SrtSimMidAt *ptr;
	long ne, np;

	ne=nEx;
	if (ne<1)
		ne=1;		/* minimum array size */
	np=nPay;
	if (np<1)
		np=1;		/* minimum array size */

	ptr = (SrtSimMidAt*) srt_calloc(1, sizeof(SrtSimMidAt));
	if (ptr==NULL)
		return (NULL);

	ptr->nPay = nPay;
	ptr->Payment = (double*) srt_calloc(np, sizeof(double));
	ptr->tPay = (Date*) srt_calloc(np, sizeof(Date));
	ptr->nEx = nEx;
	ptr->tEx = (Date*) srt_calloc(ne, sizeof(Date));
	ptr->tStart = (Date*) srt_calloc(ne, sizeof(Date));
	ptr->Strike = (double*) srt_calloc(ne, sizeof(double));
	ptr->FirstPay = (long*) srt_calloc(ne, sizeof(long));
	ptr->RedFirstPay = (double*) srt_calloc(ne, sizeof(double));
	
	if (ptr->Payment==NULL || ptr->tPay==NULL     ||
		ptr->tEx==NULL	   || ptr->tStart==NULL	  ||
		ptr->Strike==NULL  || ptr->FirstPay==NULL ||
		ptr->RedFirstPay==NULL )
	{	LGMFreeSimMidAt(&ptr);
		return (NULL);
	}
	return (ptr);
}

/*******************************************/
/* Frees simple MidAtlantic structure */
void LGMFreeSimMidAt(SrtSimMidAt **ptr)
{	if ((*ptr)==NULL) return;

	if ((*ptr)->Payment != NULL) srt_free((*ptr)->Payment);
	if ((*ptr)->tPay != NULL) srt_free((*ptr)->tPay);
	if ((*ptr)->tEx != NULL) srt_free((*ptr)->tEx);
	if ((*ptr)->tStart != NULL) srt_free((*ptr)->tStart);
	if ((*ptr)->Strike != NULL) srt_free((*ptr)->Strike);
	if ((*ptr)->FirstPay != NULL) srt_free((*ptr)->FirstPay);
	if ((*ptr)->RedFirstPay != NULL) srt_free((*ptr)->RedFirstPay);

	(*ptr)->Payment = NULL;
	(*ptr)->tPay = NULL;
	(*ptr)->tEx = NULL;
	(*ptr)->tStart = NULL;
	(*ptr)->Strike = NULL;
	(*ptr)->FirstPay = NULL;
	(*ptr)->RedFirstPay = NULL;



	srt_free((*ptr));
	*ptr = NULL;
	return;
}

/*******************************************/
/* Check validity of Simple MidAtlantic deal */
LGMErr LGMValidSimMidAt(SrtSimMidAt *ptr, Date tfirst, long *nExleft)
{	long i, j, nPay, nEx;

	nPay = ptr->nPay;
	nEx = ptr->nEx;

	if (nPay<1 || nEx<1)
		return ("no pay date or no exer date");

/* Check pay dates */
	for (i=1; i<nPay; i++)
	{	if (ptr->tPay[i] < ptr->tPay[i-1])
			return ("Pay dates out of order");
	}

/* Check exercise dates */
	for (j=1; j<nEx; j++)
	{	if (ptr->tEx[j] <= ptr->tEx[j-1])
			return("exer dates out of order");
	}

/* Check settlement dates, first paydates, and strikes */
	for (j=0; j<nEx; j++)
	{	i = ptr->FirstPay[j];
		if (i<0 || i>=nPay)
			return("no first pay date upon exercise");

		if (ptr->tEx[j] > ptr->tStart[j] ||
				ptr->tStart[j] > ptr->tPay[i]  )
			return ("settlement date too early or late");

		if (ptr->Strike[j] < 0)
			return ("strike is negative");
	}

/* Find first exercise after tfirst */
	if (ptr->tEx[nEx-1] < tfirst)
	{	ptr->FirstExer = nEx;
		*nExleft=0;
	}
	else
	{	for (j=0; j<nEx && ptr->tEx[j] < tfirst; j++)		
			;
		ptr->FirstExer = j;
		*nExleft = nEx - j;
	}
	return (NULL);
}

/*******************************************/
/* Create a general MidAtlantic stucture with nEx exercise 
	dates and nPay[j] payments at exercise j, j=0,...,nEx-1 
/* Returns pointer to structure; on failure returns NULL */
SrtGenMidAtPtr LGMCreateGenMidAt(long nEx, long *nPayArr)
{	SrtGenMidAt *ptr;
	long j, ne, np;

	ptr = (SrtGenMidAt*) srt_calloc(1, sizeof(SrtGenMidAt));
	if (ptr==NULL)
		return (NULL);

	ne=nEx;
	if (ne<1)
		ne=1;

	ptr->nEx = nEx;
   	ptr->tEx = (Date*) srt_calloc(ne, sizeof(Date));
	ptr->nPay = (long*) srt_calloc(ne, sizeof(long));
	ptr->Payment = (double**) srt_calloc(ne, sizeof(double*));
	ptr->tPay = (Date**) srt_calloc(ne, sizeof(Date*));

	if (ptr->tEx==NULL || ptr->nPay==NULL ||
		ptr->Payment==NULL || ptr->tPay==NULL )
	{	srt_free(ptr->tEx); srt_free(ptr->nPay);
		srt_free(ptr->Payment); srt_free(ptr->tPay);
		srt_free(ptr);
		return (NULL);
	}

	for (j=0; j<ne; j++)
	{	if (nEx<1)
			np=1;
		else
		{	np=nPayArr[j];
			if (np<1)
				np=1;
		}
		ptr->nPay[j] = nPayArr[j];
		ptr->Payment[j]	= (double*) srt_calloc(np,sizeof(double));
		ptr->tPay[j] = (Date*) srt_calloc(np, sizeof(Date));
	}	
	for (j=0; j<ne; j++)
	{	if (ptr->Payment[j]==NULL || ptr->tPay[j]==NULL)
		{	LGMFreeGenMidAt(&ptr);
			return (NULL);
		}
	}
	return (ptr);
}

/************************************************************/
/* Frees the General MidAtlantic structure */
void LGMFreeGenMidAt(SrtGenMidAt **ptr)
{	long j, ne;

	if ((*ptr)==NULL) return;

	ne = (*ptr)->nEx;
	if (ne<1) ne=1;

	for (j=0; j<ne; j++)
	{	if ( (*ptr)->Payment!=NULL && (*ptr)->Payment[j]!=NULL)
		{	srt_free((*ptr)->Payment[j]);
			(*ptr)->Payment[j] = NULL;
		}
		if ( (*ptr)->tPay!=NULL && (*ptr)->tPay[j]!=NULL)
		{	srt_free((*ptr)->tPay[j]);
			(*ptr)->tPay[j] = NULL;
		}
	}
	if ((*ptr)->tEx) srt_free((*ptr)->tEx);
	if ((*ptr)->nPay) srt_free((*ptr)->nPay);
	if ((*ptr)->Payment) srt_free((*ptr)->Payment);
	if ((*ptr)->tPay) srt_free((*ptr)->tPay);

	(*ptr)->tEx = NULL;
	(*ptr)->nPay = NULL;
	(*ptr)->Payment = NULL;
	(*ptr)->tPay = NULL;

	srt_free((*ptr));
	*ptr=NULL;
	return;
}

/************************************************************/
/* Check validity of General MidAtlantic deal */
LGMErr LGMValidGenMidAt(SrtGenMidAt *ptr, Date tfirst, long *nExleft)
{	long i, j, nPay, nEx;

	nEx = ptr->nEx;
	if (nEx<1)
		return ("no exer date");

/* Check order of exercise dates */
	for (j=1; j<nEx; j++)
	{	if (ptr->tEx[j] <= ptr->tEx[j-1])
			return("exer dates out of order");
	}

/* For each exercise date, check the order of the pay dates */
	for (j=0; j<nEx; j++)
	{	nPay = ptr->nPay[j];
		if (nPay < 1)
			return("no pay dates for an exercise");

		if (ptr->tPay[j][0] < ptr->tEx[j])
			return("Pay date before exercise date");

		for (i=1; i<nPay; i++)
		{	if (ptr->tPay[j][i] < ptr->tPay[j][i-1])
				return ("Pay dates out of order");
		}
	}
	
/* Find first exercise */
	if (ptr->tEx[nEx-1]<tfirst)
	{	ptr->FirstExer = nEx;
		*nExleft=0;
	}
	else
	{	for (j=0; j<nEx && ptr->tEx[j]<tfirst; j++)
			;
		ptr->FirstExer = j;
		*nExleft = nEx - j;
	}
	return (NULL);
}

/************************************************************/
/* Create a simple American stucture with nfix fixed leg periods
and nflt payments */ 
/* Returns pointer to structure; on failure returns NULL */
SrtSimAmerPtr LGMCreateSimAmer(long nfix, long nflt)
{	SrtSimAmer *ptr;
	long nfixAlloc, nfltAlloc;

	ptr = (SrtSimAmer*) srt_calloc(1, sizeof(SrtSimAmer));
	if (ptr==NULL)
		return (NULL);

	ptr->nfix = nfix;
	ptr->nflt = nflt;
	ptr->EarlyFlagFix = 1;
	ptr->EarlyFlagFlt = 1;
	ptr->ResetFlt = 0;

	nfixAlloc= max(1, nfix);	/* minimum array size */
	nfltAlloc= max(1, nflt);	/* minimum array size */

	ptr->tfixStart = (Date*) srt_calloc(nfixAlloc, sizeof(Date));
	ptr->tfixEnd = (Date*) srt_calloc(nfixAlloc, sizeof(Date));
	ptr->tfixPay = (Date*) srt_calloc(nfixAlloc, sizeof(Date));
	ptr->fixCoupon = (double*) srt_calloc(nfixAlloc, sizeof(double));
	ptr->ExtraPrem = (double*) srt_calloc(nfixAlloc, sizeof(double));

	ptr->tfltFixing = (Date*) srt_calloc(nfltAlloc, sizeof(Date));
	ptr->tfltPay = (Date*) srt_calloc(nfltAlloc, sizeof(Date));
	
	if (ptr->tfixStart==NULL || ptr->tfixEnd==NULL		||
		ptr->tfixPay==NULL	 || ptr->fixCoupon==NULL	||
		ptr->ExtraPrem==NULL || ptr->tfltFixing==NULL	||
		ptr->tfltPay==NULL )
	{	LGMFreeSimAmer(&ptr);
		return (NULL);
	}
	return (ptr);
}

/************************************************************/
/* Frees simple American structure */
void LGMFreeSimAmer(SrtSimAmer **ptr)
{	if ((*ptr)==NULL) return;

	if ((*ptr)->tfixStart) srt_free((*ptr)->tfixStart);
	if ((*ptr)->tfixEnd) srt_free((*ptr)->tfixEnd);
	if ((*ptr)->tfixPay) srt_free((*ptr)->tfixPay);
	if ((*ptr)->fixCoupon) srt_free((*ptr)->fixCoupon);
	if ((*ptr)->ExtraPrem) srt_free((*ptr)->ExtraPrem);
	if ((*ptr)->tfltFixing) srt_free((*ptr)->tfltFixing);
	if ((*ptr)->tfltPay) srt_free((*ptr)->tfltPay);

	(*ptr)->tfixStart = NULL;
	(*ptr)->tfixEnd = NULL;
	(*ptr)->tfixPay = NULL;
	(*ptr)->fixCoupon = NULL;
	(*ptr)->ExtraPrem = NULL;
	(*ptr)->tfltFixing = NULL;
	(*ptr)->tfltPay = NULL;

	srt_free((*ptr));
	*ptr = NULL;
	return;
}

/************************************************************/
/* Check validity of simple American deal */
LGMErr LGMValidSimAmer(SrtSimAmer *ptr, Date tfirst)
{	long i, j, nfix, nflt;
	Date tfirstExer, tfirstStart;

	nfix = ptr->nfix;
	nflt = ptr->nflt;

	if (nfix<1 || nfix>601 || nflt<1 || nflt>601)
		return ("too few or too many dates");

	tfirstExer = max(tfirst, ptr->tFirstExer);
	if (tfirstExer >= ptr->tfltFixing[nflt-1])
		return("last floating rate has already been fixed");

	if (ptr->CalBusLag!=0) 
		ptr->CalBusLag=1;
	if (ptr->CalBusLag==1)
		tfirstStart = add_unit(tfirstExer, ptr->lagExerSettle, SRT_BDAY, ptr->convSettle);
	else
		tfirstStart = add_unit(tfirstExer, ptr->lagExerSettle, SRT_DAY, ptr->convSettle);
	if (tfirstStart >= ptr->tfixEnd[nfix-1] || tfirstStart >= ptr->tfltPay[nflt-1])
		return("no underlying left");

	for (i=0; i<(nfix-1)
			&& ptr->tfixStart[i] < ptr->tfixEnd[i]
			&& ptr->tfixStart[i] < ptr->tfixStart[i+1]
			&& ptr->tfixEnd[i] < ptr->tfixEnd[i+1]
			&& ptr->tfixPay[i] <= ptr->tfixPay[i+1]
			&& ptr->tfixStart[i] < ptr->tfixPay[i]; i++)
        ;
    if (ptr->tfixStart[nfix-1] >= ptr->tfixEnd[nfix-1] ||
		ptr->tfixStart[nfix-1] >= ptr->tfixPay[nfix-1])
        i = nfix-2;
	
    for (j=0; j<(nflt-1)
			&& ptr->tfltPay[j]<ptr->tfltPay[j+1]
			&& ptr->tfltFixing[j]<ptr->tfltFixing[j+1]
			&& ptr->tfltFixing[j]<ptr->tfltPay[j]; j++)
        ;

	if (ptr->tfltFixing[nflt-1] >= ptr->tfltPay[nflt-1])
		j=-1;

    if (i<(nfix-1) || j<(nflt-1))
		return("dates are wrong");

	return (NULL);
}

/*******************************************/
/* Allocate a Bermudan inverse floater with nEx exercise dates and n coupon periods */ 
/* Returns pointer to structure; on failure returns NULL */
SrtCallInvFltPtr LGMCreateCallInvFlt(long nEx, long n)
{	
	SrtCallInvFlt *ptr;
	long nExAlloc, nAlloc;

	nExAlloc = max (1, nEx);
	nAlloc = max (0, n);

	ptr = (SrtCallInvFlt*) srt_calloc (1, sizeof (SrtCallInvFlt));
	if (!ptr) return NULL;

	ptr->nEx = nEx;
	ptr->tEx = (Date*) srt_calloc (nExAlloc, sizeof(Date));
	ptr->tSet = (Date*) srt_calloc (nExAlloc, sizeof(Date));
	ptr->iSet = (long*) srt_calloc (nExAlloc, sizeof(long));
	ptr->strike = (double*) srt_calloc (nExAlloc, sizeof(double));

	ptr->nCpn = n;
	ptr->tCpnStart = (Date*) srt_calloc (nAlloc, sizeof(Date));
	ptr->tCpnPay = (Date*) srt_calloc (nAlloc, sizeof(Date));
	ptr->a = (double*) srt_calloc (nAlloc, sizeof(double));
	ptr->gear = (double*) srt_calloc (nAlloc, sizeof(double));
	ptr->cvg = (double*) srt_calloc (nAlloc, sizeof(double));
	ptr->lcvg = (double*) srt_calloc (nAlloc, sizeof(double));
	ptr->cap_str = (double*) srt_calloc (nAlloc, sizeof(double));

	if (!ptr->tEx || !ptr->iSet || !ptr->strike 
		|| !ptr->tCpnStart || !ptr->tCpnPay
		|| !ptr->a || !ptr->gear || !ptr->cvg || !ptr->lcvg || !ptr->cap_str)
	{	
		LGMFreeCallInvFlt (&ptr);
		return NULL;
	}

	ptr->fwdVolCpn = 0;
	ptr->fwdVolEx = 0;
	ptr->fwdVolMat = NULL;

	return ptr; 
}

/* Frees a Bermudan inverse floater structure */
void LGMFreeCallInvFlt (SrtCallInvFlt **ptr)
{	
	if ((*ptr) == NULL) return;

	if ((*ptr)->tEx != NULL) srt_free ((*ptr)->tEx);
	if ((*ptr)->tSet != NULL) srt_free ((*ptr)->tEx);
	if ((*ptr)->iSet != NULL) srt_free ((*ptr)->iSet);
	if ((*ptr)->strike != NULL) srt_free ((*ptr)->strike);
	if ((*ptr)->tCpnStart != NULL) srt_free ((*ptr)->tCpnStart);
	if ((*ptr)->tCpnPay != NULL) srt_free ((*ptr)->tCpnPay);
	if ((*ptr)->a != NULL) srt_free ((*ptr)->a);
	if ((*ptr)->gear != NULL) srt_free ((*ptr)->gear);
	if ((*ptr)->cvg != NULL) srt_free ((*ptr)->cvg);
	if ((*ptr)->lcvg != NULL) srt_free ((*ptr)->lcvg);
	if ((*ptr)->cap_str != NULL) srt_free ((*ptr)->cap_str);
	if ((*ptr)->fwdVolMat != NULL) free_dmatrix ((*ptr)->fwdVolMat,0,(*ptr)->fwdVolCpn,0,(*ptr)->fwdVolEx);

	(*ptr)->tEx = NULL;
	(*ptr)->tSet = NULL;
	(*ptr)->iSet = NULL;
	(*ptr)->strike = NULL;
	(*ptr)->tCpnStart = NULL;
	(*ptr)->tCpnPay = NULL;
	(*ptr)->a = NULL;
	(*ptr)->gear = NULL;
	(*ptr)->cvg = NULL;
	(*ptr)->lcvg = NULL;
	(*ptr)->cap_str = NULL;
	(*ptr)->fwdVolMat = NULL;

	srt_free ((*ptr));
	*ptr = NULL;
	return;
}

/*******************************************/
/* Check validity of a Bermudan inverse floater */
LGMErr LGMValidCallInvFlt(SrtCallInvFlt *ptr, Date tfirst, long *EffnEx)
{	long i, j, n, nEx;

	n = ptr->nCpn;
	nEx = ptr->nEx;
	if (n<1 || nEx<1)
		return ("no cpn periods or no exer dates");

/* Check coupons */
	for (i=1; i<n; i++)
	{	if (ptr->tCpnStart[i] <= ptr->tCpnStart[i-1] || ptr->tCpnPay[i] <= ptr->tCpnPay[i-1])
			return ("Coupon dates out of order");
		if (ptr->a[i]<=0 || ptr->a[i-1]<=0 || ptr->gear[i]<0 || ptr->gear[i-1]<0)
			return ("handle is wrong");
		if (ptr->cap_str[i]<=0 || ptr->cap_str[i-1]<=0)
			return ("cap strike is incorrect");
	}

/* Check exercises */
	for (j=1; j<nEx; j++)
	{	if (ptr->tEx[j] <= ptr->tEx[j-1])
			return("exer dates out of order");
	}

/* find first exercise date on or after tfirst */
	for (j=0; j<nEx && ptr->tEx[j]<tfirst; j++);
	ptr->FirstEx = j;
	*EffnEx = ptr->nEx - j;
	return (NULL);
}


/*******************************************
********************************************
********************************************/
/* Allocate a Callable Time Swap with nEx exercise dates and n coupon periods */ 
/* Returns pointer to structure; on failure returns NULL */
SrtCallTimeSwapPtr LGMCreateCallableTimeSwap(long nEx, long n, long observation_freq)
{	
	SrtCallTimeSwap *ptr;
	long nExAlloc, nAlloc,i,j;

	nExAlloc = max (1, nEx);
	nAlloc = max (0, n);

	ptr = (SrtCallTimeSwap*) srt_calloc (1, sizeof (SrtCallTimeSwap));
	if (!ptr) return NULL;

	ptr->observation_freq = observation_freq;

	ptr->nEx = nEx;
	ptr->tEx = (Date*) srt_calloc (nExAlloc, sizeof(Date));
	ptr->iSet = (long*) srt_calloc (nExAlloc, sizeof(long));
	ptr->tSet = (Date*) srt_calloc (nExAlloc, sizeof(Date));
	ptr->strike = (double*) srt_calloc (nExAlloc, sizeof(double));

	ptr->nCpn = n;
	ptr->tCpnStart = (Date*) srt_calloc (nAlloc, sizeof(Date));
	ptr->tCpnPay = (Date*) srt_calloc (nAlloc, sizeof(Date));
	ptr->tCpn = (double*) srt_calloc (nAlloc, sizeof(double));
	ptr->gear = (double*) srt_calloc (nAlloc, sizeof(double));
	ptr->tCvgCpn = (double*) srt_calloc (nAlloc, sizeof(double));
	ptr->tFundingPayment = (double*) srt_calloc(nAlloc,sizeof(double));
	
	ptr->barriers = (double**) srt_calloc(nAlloc, sizeof(double*));
	for (i=0;i<nAlloc;i++)
		(ptr->barriers)[i] = (double*) srt_calloc(2,sizeof(double));

	ptr->tForwardTS = (double**) srt_calloc(2,sizeof(double*));
	for (i=0;i<=1;i++)
		(ptr->tForwardTS)[i] = (double*) srt_calloc(nExAlloc,sizeof(double));

	ptr->observationdays = (Date***) srt_calloc(nAlloc, sizeof(Date**));
	for (i=0;i<nAlloc;i++)
	{
		(ptr->observationdays)[i] = (Date**) srt_calloc(observation_freq,sizeof(Date*));
		for (j=0;j<observation_freq;j++)
		{
			(ptr->observationdays)[i][j] = (Date*) srt_calloc(2,sizeof(Date));
		}
	}

	/* Just for debug 
	ptr->reval_times = (double*) dvector(0,19);
	*/

	ptr->ratiodays_for_subperiod = (double**) srt_calloc(nAlloc, sizeof(double*));
	for (i=0;i<nAlloc;i++)
		(ptr->ratiodays_for_subperiod)[i] = (double*) srt_calloc(observation_freq,sizeof(double));


	if (!ptr->tEx || !ptr->iSet ||  !ptr->tCpnStart || !ptr->tCpnPay
		|| !ptr->tCpn || !ptr->gear || !ptr->tCvgCpn || !ptr->barriers  
		|| !ptr->observationdays || !ptr->ratiodays_for_subperiod)
	{	
		LGMFreeCallableTimeSwap (&ptr);
		return NULL;
	}

	ptr->fwdVolCpn = 0;
	ptr->fwdVolEx = 0;
	ptr->fwdVolMat = NULL;

	return ptr; 
}

/* Frees a Callable time swap structure */
void LGMFreeCallableTimeSwap (SrtCallTimeSwap **ptr)
{
int i,j;
	if ((*ptr) == NULL) return;

	if ((*ptr)->tEx != NULL) srt_free ((*ptr)->tEx);
	if ((*ptr)->iSet != NULL) srt_free ((*ptr)->iSet);
	if ((*ptr)->tSet != NULL) srt_free ((*ptr)->tSet);
	if ((*ptr)->strike != NULL) srt_free ((*ptr)->strike);
	if ((*ptr)->tCpnStart != NULL) srt_free ((*ptr)->tCpnStart);
	if ((*ptr)->tCpnPay != NULL) srt_free ((*ptr)->tCpnPay);
	if ((*ptr)->gear != NULL) srt_free ((*ptr)->gear);
	if ((*ptr)->tCvgCpn != NULL) srt_free ((*ptr)->tCvgCpn);
	if ((*ptr)->tCpn != NULL) srt_free ((*ptr)->tCpn);
	if ((*ptr)->tFundingPayment != NULL) srt_free((*ptr)->tFundingPayment);
	for (i=0;i<(*ptr)->nCpn;i++)
	{
		if (((*ptr)->barriers)[i] != NULL) srt_free(((*ptr)->barriers)[i]);
	}
	if ((*ptr)->barriers != NULL) srt_free((*ptr)->barriers);
	for (i=0;i<=1;i++)
	{
		if (((*ptr)->tForwardTS)[i] != NULL) srt_free(((*ptr)->tForwardTS)[i]);
	}
	if ((*ptr)->tForwardTS != NULL) srt_free((*ptr)->tForwardTS);
	for (i=0;i<(*ptr)->nCpn;i++)
	{
		for (j=0;j<(*ptr)->observation_freq;j++)
		{
			if (((*ptr)->observationdays)[i][j] != NULL) srt_free(((*ptr)->observationdays)[i][j]);
		}	
		if (((*ptr)->observationdays)[i] != NULL) srt_free(((*ptr)->observationdays)[i]);
	}
	if ((*ptr)->observationdays != NULL) srt_free((*ptr)->observationdays);
	for (i=0;i<(*ptr)->nCpn;i++)
	{
		if (((*ptr)->ratiodays_for_subperiod)[i] != NULL) srt_free(((*ptr)->ratiodays_for_subperiod)[i]);
	}
	if ((*ptr)->ratiodays_for_subperiod != NULL) srt_free((*ptr)->ratiodays_for_subperiod);
	if ((*ptr)->fwdVolMat != NULL) free_dmatrix ((*ptr)->fwdVolMat,0,(*ptr)->fwdVolCpn,0,(*ptr)->fwdVolEx);
	/* Just for debug 
	if ((*ptr)->reval_times != NULL) free_dvector((*ptr)->reval_times,0,19);
	*/


	(*ptr)->tEx = NULL;
	(*ptr)->iSet = NULL;
	(*ptr)->tSet = NULL;
	(*ptr)->strike = NULL;
	(*ptr)->tCpnStart = NULL;
	(*ptr)->tCpnPay = NULL;
	(*ptr)->gear = NULL;
	(*ptr)->tCvgCpn = NULL;
	(*ptr)->tCpn = NULL;
	(*ptr)->tFundingPayment = NULL;
	(*ptr)->tForwardTS = NULL;
	(*ptr)->barriers = NULL;
	(*ptr)->observationdays = NULL;
	(*ptr)->ratiodays_for_subperiod = NULL;
	(*ptr)->fwdVolMat = NULL;
	/* Just for debug 
	(*ptr)->reval_times = NULL;
	*/
	srt_free ((*ptr));
	*ptr = NULL;
	return;
}


/*******************************************/
/* Check validity of a Callable Time Swap */
LGMErr LGMValidCallTimeSwap(SrtCallTimeSwap *ptr, Date tfirst, long *EffnEx)
{	long i, j, n, nEx;

	n = ptr->nCpn;
	nEx = ptr->nEx;
	if (n<1 || nEx<1)
		return ("no cpn periods or no exer dates");

/* Check coupons */
	for (i=1; i<n; i++)
	{	if (ptr->tCpnStart[i] <= ptr->tCpnStart[i-1] || ptr->tCpnPay[i] <= ptr->tCpnPay[i-1])
			return ("Coupon dates out of order");
	/* Negative gearings useful for callable inverse floaters on range accrual */
	/*
	if (ptr->gear[i-1]<0)
			return ("handle is wrong");
			*/
	}

/* Check exercises */
	for (j=1; j<nEx; j++)
	{	if (ptr->tEx[j] <= ptr->tEx[j-1])
			return("exer dates out of order");
	}

/* find first exercise date on or after tfirst */
	for (j=0; j<nEx && ptr->tEx[j]<tfirst; j++);
	ptr->FirstEx = j;
	*EffnEx = ptr->nEx - j;
	return (NULL);
}

/*******************************************/
/* Create a Bermudan cap floater with nEx exercise dates
and n coupon periods */ 
/* Returns pointer to structure; on failure returns NULL */
SrtCallCapFltPtr LGMCreateCallCapFlt(long nEx, long n)
{	SrtCallCapFlt *ptr;
	long nExAlloc, nAlloc;

	nExAlloc = max(1,nEx);
	nAlloc = max(0,n);

	ptr = (SrtCallCapFlt*) srt_calloc(1, sizeof(SrtCallCapFlt));
	if (ptr==NULL) return (NULL);

	ptr->nCpn = n;
	ptr->tCpn = (Date*) srt_calloc(nAlloc+1, sizeof(Date));
	ptr->amax = (double*) srt_calloc(nAlloc+1, sizeof(double));
	ptr->amin = (double*) srt_calloc(nAlloc+1, sizeof(double));
	ptr->marg = (double*) srt_calloc(nAlloc+1, sizeof(double));

	ptr->nEx = nEx;
	ptr->tEx = (Date*) srt_calloc(nExAlloc, sizeof(Date));
	ptr->iSet = (long*) srt_calloc(nExAlloc, sizeof(long));
	ptr->strike = (double*) srt_calloc(nExAlloc, sizeof(double));
	
	if (ptr->tCpn==NULL || ptr->amax==NULL || ptr->amin==NULL || ptr->marg==NULL ||
		ptr->tEx==NULL	|| ptr->iSet==NULL || ptr->strike==NULL)
	{	LGMFreeCallCapFlt(&ptr);
		return (NULL);
	}
	return (ptr);
}

/*******************************************/
/* Frees a Bermudan cap floater structure */
void LGMFreeCallCapFlt(SrtCallCapFlt **ptr)
{	if ((*ptr)==NULL) return;

	if ((*ptr)->tCpn != NULL) srt_free((*ptr)->tCpn);
	if ((*ptr)->amax != NULL) srt_free((*ptr)->amax);
	if ((*ptr)->amin != NULL) srt_free((*ptr)->amin);
	if ((*ptr)->marg != NULL) srt_free((*ptr)->marg);
	if ((*ptr)->tEx != NULL) srt_free((*ptr)->tEx);
	if ((*ptr)->iSet != NULL) srt_free((*ptr)->iSet);
	if ((*ptr)->strike != NULL) srt_free((*ptr)->strike);

	(*ptr)->tCpn = NULL;
	(*ptr)->amax = NULL;
	(*ptr)->amin = NULL;
	(*ptr)->marg = NULL;
	(*ptr)->tEx = NULL;
	(*ptr)->iSet = NULL;
	(*ptr)->strike = NULL;

	srt_free((*ptr));
	*ptr = NULL;
	return;
}

/*******************************************/
/* Check validity of a Bermudan cap floater */
/* set the iSet[] array */
LGMErr LGMValidCallCapFlt(SrtCallCapFlt *ptr, Date tfirst, long *EffnEx)
{	long i, j, n, nEx;

	n = ptr->nCpn;
	nEx = ptr->nEx;
	if (n<1 || nEx<1)
		return ("no cpn periods or no exer dates");

	if (ptr->cpnBasis<0 || ptr->cpnBasis>=LASTBASISCODE ||
		ptr->aBasis<0 ||   ptr->aBasis>=LASTBASISCODE )
		return ("unknown basis");

/* Check coupons */
	for (i=1; i<n; i++)
	{	if (ptr->tCpn[i] <= ptr->tCpn[i-1])
			return ("Coupon dates out of order");

		if (ptr->marg[i]<-1.0 || ptr->marg[i]>=1.0)
			return ("margin is incorrect");
	}

	if (ptr->lvg<=0.0 || ptr->lvg>=100.0)
		return ("bad leverage");

/* Check exercises */
	for (j=1; j<nEx; j++)
	{	if (ptr->tEx[j] <= ptr->tEx[j-1])
			return("exer dates out of order");
	}

/* set the iSet[j] array (settlement upon exercise) */
	i=0;
	for (j=0; j<nEx; j++)
	{	for (; i<n && ptr->tCpn[i] < ptr->tEx[j]; i++)
			;
		if (i<n)
			ptr->iSet[j] = i;
		else
			return ("bad settlement");
		if (ptr->strike[j]<=0.0 || ptr->strike[j]>=10.0)
			return ("bad strike");
		
	}

/* find first exercise date on or after tfirst */
	for (j=0; j<nEx && ptr->tEx[j]<tfirst; j++)
		;
	ptr->FirstEx = j;
	*EffnEx = ptr->nEx - j;
	return (NULL);
}

/************* Utilities for term structures ****************/
/* ALL term structures should be verified by LGMVerifyLGM_TS or
LGMVerifySigKapTS before use, since on error the rouintes that
get Zeta, G, Sigma, and Kappa will crash */ 
/************************************************************/
/* Get zeta at thedate from term structure */
/* On error, crashes */
double LGMZetaFromTS(Date thedate, Date tNow, LGM_TS *tsPtr)
{	double dt1, dt2, slope, zeta1, zeta2;
	long j, lastj;

	if (thedate<=tNow)		/* zeta = 0 if today or before */
	{	zeta2 = 0.0;
		return (zeta2);
	}
	
/* interpolate to find zeta(tNow) */
	lastj = tsPtr->numZ - 1;
	for (j=1; j < lastj && tsPtr->zdate[j] < tNow; j++)
		;

	dt1 = (double) (tsPtr->zdate[j] - tsPtr->zdate[j-1]);
	slope = (tsPtr->zeta[j] - tsPtr->zeta[j-1])/dt1;
	dt2 = (double)(tNow - tsPtr->zdate[j-1]);
	zeta1 = tsPtr->zeta[j-1] + slope*dt2;

/* interpolate to find zeta(thedate) */
	for (j=1; j < lastj && tsPtr->zdate[j] < thedate; j++)
		;

	dt1 = (double) (tsPtr->zdate[j] - tsPtr->zdate[j-1]);
	slope = (tsPtr->zeta[j] - tsPtr->zeta[j-1])/dt1;
	dt2 = (double)(thedate - tsPtr->zdate[j-1]); 
	zeta2 = tsPtr->zeta[j-1] + slope*dt2;

	if (zeta2<=zeta1)
		return (0.);
	else
		return (zeta2 - zeta1);
}

/******************************************************/
/* Get G at thedate from term structure */
/* On error, crashes */
double LGMGFromTS(Date thedate, LGM_TS *tsPtr)
{	long i, lasti;
	double dt, dt1, Gatdate;
	double slope;

/* interpolate to find G(thedate) */
	lasti = tsPtr->numG - 1;
	for (i=1; i<lasti && tsPtr->Gdate[i] < thedate; i++)
		;

	dt1 = (double)(tsPtr->Gdate[i] - tsPtr->Gdate[i-1]);
	slope = (tsPtr->G[i] - tsPtr->G[i-1])/dt1;
	dt = (double)(thedate - tsPtr->Gdate[i-1]); 
	Gatdate = tsPtr->G[i-1] + slope*dt;

	return (Gatdate);
}


/******************************************************/
/* Get sigma at thedate from term structure */
/* On error, crashes */
double LGMSigFromTS(Date thedate, SigKapTS *tsPtr)
{	double Sigatdate;
	long j, lastj;

/* find Sigma(thedate) from piecewise constant curve */
	lastj = tsPtr->numS - 1;
	for (j=0; j<lastj && tsPtr->sdate[j] < thedate; j++)
		;
	Sigatdate = tsPtr->sig[j];
	return (Sigatdate);
 }

/******************************************************/
/* Get kappa at the date from term structure */
/* On error, crashes */
double LGMKapFromTS(Date thedate, SigKapTS *tsPtr)
{	double Kapatdate;
	long j, lastj;

/* find Kappa(thedate) from piecewise constant curve */
	lastj = tsPtr->numK - 1;
	for (j=0; j<lastj && tsPtr->kdate[j] < thedate; j++)
		;
	Kapatdate = tsPtr->kap[j];
	return (Kapatdate);
 }

/******************************************************/
/* Routine to verify zeta-`G term structure is sensible */
/* Also re-scales and cleans up the zeta-G term structure */
LGMErr LGMVerifyLGM_TS(LGM_TS *tsPtr)
{	long j, i, m, n;
	double slope, dt, glast;
	double *G, *zeta;

	if (tsPtr==NULL)
		return ("no term structure");

	m= tsPtr->numZ;
	n= tsPtr->numG;
	zeta = tsPtr->zeta;
	G = tsPtr->G;

	if (m<2 || n<2)
		return ("too few dates in term struc");

/* check dates */	
	for (j=1; j<m; j++)
	{	if (tsPtr->zdate[j] <= tsPtr->zdate[j-1])
			return ("zeta dates out of order");
	}

	for (i=1; i<n; i++)
	{	if (tsPtr->Gdate[i] <= tsPtr->Gdate[i-1])
			return ("G dates out of order");
	}

/* check zeta values, and  change to feasible set if needed	*/
	for (j=1; j<m; j++)
	{	if (zeta[j] < zeta[j-1])
		{
			smessage("calibration failed for diag. swaption : %d \n", j); 
			/* zeta(t) must be increasing */
			zeta[j] = zeta[j-1];
		}
	}

/* check G values, and change to feasible set if needed */
	if (G[0] < G[n-1])
	{	for (i=0; i<n; i++)
			G[i] = -G[i];
	}
	for (i=1; i<n; i++)
	{	if (G[i] > G[i-1])			/* G(t) must be decreasing */
			G[i] = G[i-1];
	}

/* Re-scale term structure */
	glast = G[n-1];
	dt = (double)(tsPtr->Gdate[n-1] - tsPtr->Gdate[0]);
	slope = (G[0] - glast)*365.0/dt;
	if (fabs(slope) < 1.e-8)
		return("G values too small");

	for (i=0; i<n; i++)
		G[i] = (G[i] - glast)/slope;

	for (j=0; j<m; j++)
		zeta[j] = zeta[j] * slope * slope;

	return (NULL);
}

LGMErr LGMVerifyLGM_TS_For_DiagSwpts(LGMCalSet *CSPtr,LGM_TS *tsPtr)
{
	long j, i, m, n,most_exp_index, n_swpts;
	double slope, dt, glast;
	double *G, *zeta;

	if (tsPtr==NULL)
		return ("no term structure");

	m = tsPtr->numZ;
	n = tsPtr->numG;
	zeta = tsPtr->zeta;
	G = tsPtr->G;

	if (m<2 || n<2)
		return ("too few dates in term struc");

/* check dates */	
	for (j=1; j<m; j++)
	{	if (tsPtr->zdate[j] <= tsPtr->zdate[j-1])
			return ("zeta dates out of order");
	}

	for (i=1; i<n; i++)
	{	if (tsPtr->Gdate[i] <= tsPtr->Gdate[i-1])
			return ("G dates out of order");
	}


	/* find the most exp swpts index */
	most_exp_index = 1;
	n_swpts = CSPtr->nEx;
	
	for(i = 2; i <= n_swpts; i++)
	{
		if(CSPtr->Vlong[i] > CSPtr->Vlong[most_exp_index] ) most_exp_index = i;
	}

/* check zeta values, and  change to feasible set if needed	*/
	for (j = (most_exp_index+1); j<m; j++)
	{	if (zeta[j] < zeta[j-1])
		{
			smessage("calibration failed for diag. swaption : %d \n", j);
			/* zeta(t) must be increasing */
			zeta[j] = zeta[j-1];
		}
	}

	for(j = most_exp_index; j > 1 ; j--)
	{
		if(zeta[j] < zeta[j-1])
		{
			smessage("calibration failed for diag. swaption : %d \n", j);
			zeta[j-1] = zeta[j]; 
		}
	}

/* check G values, and change to feasible set if needed */
	if (G[0] < G[n-1])
	{	for (i=0; i<n; i++)
			G[i] = -G[i];
	}
	for (i=1; i<n; i++)
	{	if (G[i] > G[i-1])			/* G(t) must be decreasing */
		G[i] = G[i-1];
	}

/* Re-scale term structure */
	glast = G[n-1];
	dt = (double)(tsPtr->Gdate[n-1] - tsPtr->Gdate[0]);
	slope = (G[0] - glast)*365.0/dt;
	if (fabs(slope) < 1.e-8)
		return("G values too small");

	for (i=0; i<n; i++)
		G[i] = (G[i] - glast)/slope;

	for (j=0; j<m; j++)
		zeta[j] = zeta[j] * slope * slope;

	return (NULL);
}



/*****************************************************************/
/* Routine to verify that sigma-kappa term structure is sensible */
LGMErr LGMVerifySigKapTS(SigKapTS *tsPtr)
{	long j, i, m, n;

	if (tsPtr==NULL)
		return ("no term structure");
	
	m = tsPtr->numS;
	n = tsPtr->numK;

	if (m<1 || n<1)
		return ("too few dates in term structure");

/* check dates */	
	for (j=1; j<m; j++)
	{	if (tsPtr->sdate[j] <= tsPtr->sdate[j-1])
			return ("sigma dates out of order");
	}

	for (i=1; i<n; i++)
	{	if (tsPtr->kdate[i] <= tsPtr->kdate[i-1])
			return ("kappa dates out of order");
	}

/* check sigma values */
	for (j=0; j<m; j++)
	{	if (tsPtr->sig[j]<0)
			tsPtr->sig[j] = -(tsPtr->sig[j]);
	}
	return (NULL);
}

/******************************************************/
/*  Create LGM zeta-G term structure */
/* On failure, returns NULL */

LGM_TSPtr LGMCreateLGM_TS(long numZ, long NumG)
{	LGM_TS *ptr;




	if (numZ<2 || numZ>600 || NumG<2 || NumG>600)
		return (NULL);

	ptr = (LGM_TS*) srt_calloc(1, sizeof(LGM_TS));
	if (ptr==NULL)
		return (NULL);

	ptr->numZ = numZ;
	ptr->numG = NumG;
	ptr->zdate = (Date*) srt_calloc(numZ, sizeof(Date));
	ptr->zeta = (double*) srt_calloc(numZ, sizeof(double));
	ptr->Gdate = (Date*) srt_calloc(NumG, sizeof(Date));
	ptr->G = (double*) srt_calloc(NumG, sizeof(double));
	
	if (ptr->zdate==NULL || ptr->zeta==NULL ||
			ptr->Gdate==NULL || ptr->G==NULL )

	{
		LGMFreeLGM_TS(&ptr);
	}

	return (ptr);
}

LGM_TSPtr LGMCreateLGM2F_TS(long numZ, long numG)
{	
	LGM_TS *ptr;
	long i;

	if (numZ<2 || numZ>600 || numG<2 || numG>600)
		return (NULL);

	ptr = (LGM_TS*) srt_calloc(1, sizeof(LGM_TS));
	if (ptr==NULL)
		return (NULL);

	ptr->numZ = numZ;
	ptr->numG = numG;
	ptr->zdate = (Date*) srt_calloc(numZ, sizeof(Date));
	ptr->Gdate = (Date*) srt_calloc(numG, sizeof(Date));

	/* TWO FACTOR CASE */
	ptr->Zeta1 = (double*) srt_calloc(numZ, sizeof(double));
	ptr->Zeta2 = (double*) srt_calloc(numZ, sizeof(double));
	ptr->Zeta12 = (double*) srt_calloc(numZ, sizeof(double));
	ptr->H1 = (double*) srt_calloc(numG, sizeof(double));
	ptr->H2 = (double*) srt_calloc(numG, sizeof(double));


	ptr->Zeta1Scenari = (double**) srt_calloc(2*numZ+1, sizeof(double*));
	ptr->Zeta2Scenari = (double**) srt_calloc(2*numZ+1, sizeof(double*));
	ptr->Zeta12Scenari = (double**) srt_calloc(2*numZ+1, sizeof(double*));
	ptr->H1Scenari = (double**) srt_calloc(2*numZ+1, sizeof(double*));
	ptr->H2Scenari = (double**) srt_calloc(2*numZ+1, sizeof(double*));
	for ( i = 0; i <= 2*numZ; i++)
	{
		ptr->Zeta1Scenari[i] = (double*) srt_calloc(numZ, sizeof(double));
		ptr->Zeta2Scenari[i] = (double*) srt_calloc(numZ, sizeof(double));
		ptr->Zeta12Scenari[i] = (double*) srt_calloc(numZ, sizeof(double));
	}
	for ( i = 0; i <= 2*numZ; i++)
	{
		ptr->H1Scenari[i] = (double*) srt_calloc(numG, sizeof(double));
		ptr->H2Scenari[i] = (double*) srt_calloc(numG, sizeof(double));
	}

	if (ptr->zdate==NULL ||	ptr->Gdate==NULL || ptr->H1==NULL || ptr->H2==NULL
		|| ptr->Zeta1Scenari==NULL || ptr->Zeta2Scenari==NULL || ptr->Zeta12Scenari==NULL
		|| ptr->H1Scenari==NULL || ptr->H2Scenari==NULL)
	{
		LGMFreeLGM2F_TS(&ptr);
	}

	return (ptr);
}

/******************************************************/
/* free LGM zeta G term structure */
void LGMFreeLGM_TS(LGM_TS **ptr)
{	if ((*ptr)==NULL) return;

	if ((*ptr)->zdate) srt_free((*ptr)->zdate);
	if ((*ptr)->zeta) srt_free((*ptr)->zeta);
	if ((*ptr)->Gdate) srt_free((*ptr)->Gdate);
	if ((*ptr)->G) srt_free((*ptr)->G);

	(*ptr)->zdate = NULL;
	(*ptr)->zeta = NULL;
	(*ptr)->Gdate = NULL;
	(*ptr)->G = NULL;

	srt_free((*ptr));
	*ptr = NULL;
	return;
}

void LGMFreeLGM2F_TS(LGM_TS **ptr)
{	
	int numZ, numG, i;

	if ((*ptr)==NULL) return;

	if ((*ptr)->Gdate) {srt_free((*ptr)->Gdate);(*ptr)->Gdate = NULL;}
	if ((*ptr)->zdate) {srt_free((*ptr)->zdate); (*ptr)->zdate = NULL;}

	if ((*ptr)->Zeta1) {srt_free((*ptr)->Zeta1);(*ptr)->Zeta1 = NULL;}
	if ((*ptr)->Zeta2) {srt_free((*ptr)->Zeta2);(*ptr)->Zeta2 = NULL;}
	if ((*ptr)->Zeta12) {srt_free((*ptr)->Zeta12);(*ptr)->Zeta12 = NULL;}
	if ((*ptr)->H1) {srt_free((*ptr)->H1);(*ptr)->H1 = NULL;}
	if ((*ptr)->H2) {srt_free((*ptr)->H2);(*ptr)->H2 = NULL;}

	numZ = (*ptr)->numZ;
	numG = (*ptr)->numG;
	for ( i = 0; i <= 2*numZ; i++)
	{
		srt_free((*ptr)->Zeta1Scenari[i]);
		srt_free((*ptr)->Zeta2Scenari[i]);
		srt_free((*ptr)->Zeta12Scenari[i]);
		srt_free((*ptr)->H1Scenari[i]);
		srt_free((*ptr)->H2Scenari[i]);
	}
	srt_free((*ptr)->Zeta1Scenari);(*ptr)->Zeta1Scenari = NULL;
	srt_free((*ptr)->Zeta2Scenari);(*ptr)->Zeta2Scenari = NULL;
	srt_free((*ptr)->Zeta12Scenari);(*ptr)->Zeta12Scenari = NULL;
	srt_free((*ptr)->H1Scenari);(*ptr)->H1Scenari = NULL;
	srt_free((*ptr)->H2Scenari);(*ptr)->H2Scenari = NULL;

	srt_free((*ptr));
	*ptr = NULL;
	return;
}

/******************************************************/
/* Create SigmaKappa term structure */
/* On failure, returns NULL */
SigKapTSPtr LGMCreateSigKapTS(long NumS, long NumK)
{	SigKapTS *ptr;

	if (NumS<1 || NumS>600 || NumK<1 || NumK>600)
		return (NULL);

	ptr = (SigKapTS*) srt_calloc(1, sizeof(SigKapTS));
	if (ptr==NULL)
		return (NULL);

	ptr->numS = NumS;
	ptr->numK = NumK;
	ptr->sdate = (Date*) srt_calloc(NumS, sizeof(Date));
	ptr->sig = (double*) srt_calloc(NumS, sizeof(double));
	ptr->kdate = (Date*) srt_calloc(NumK, sizeof(Date));
	ptr->kap = (double*) srt_calloc(NumK, sizeof(double));
	
	if (ptr->sdate==NULL || ptr->sig==NULL ||
			ptr->kdate==NULL || ptr->kap==NULL )
	{
			LGMFreeSigKapTS(&ptr);
	}

	return (ptr);
}

/* free sigma kappa term struture */
void LGMFreeSigKapTS(SigKapTS **ptr)
{	if ((*ptr)==NULL) return;

	if ((*ptr)->sdate) srt_free((*ptr)->sdate);
	if ((*ptr)->sig) srt_free((*ptr)->sig);
	if ((*ptr)->kdate) srt_free((*ptr)->kdate);
	if ((*ptr)->kap) srt_free((*ptr)->kap);
	
	(*ptr)->sdate = NULL;
	(*ptr)->sig = NULL;
	(*ptr)->kdate = NULL;
	(*ptr)->kap = NULL;

	srt_free((*ptr));
	*ptr = NULL;
	return;
}

/******************************************************************/
/* Routines to convert between zetaG term structure and
  sigma kappa term structures*/
/******************************************************************/
/* This routine creates a LGM zeta G term structure, and fills it
in with values equivalent to the sigma kappa term structure.
On error, it returns a NULL */
LGM_TSPtr LGMConvertSigKaptoZG(Date tNow, SigKapTS *SKtsPtr)
{	LGMErr error;
	LGM_TS* LGMtsPtr;
	Date   tprev;
	long i, j, i1, j1, numG, numZ, numK, numS;
	double dt, gprev, cumul, kap;

	error = NULL;

/* check out sigma-kappa ts */
	error = LGMVerifySigKapTS(SKtsPtr);
	if (error!=NULL)
		return(NULL);

	numK = SKtsPtr->numK;
	numS = SKtsPtr->numS;

/* find all dates after today (tNow) */
	for (j1=0; j1 < numS && SKtsPtr->sdate[j1] <= tNow; j1++)
		;
	for (i1=0; i1 < numK && SKtsPtr->kdate[i1] <= tNow; i1++)
		;

	numG = numK + 1 - i1;
	numZ = numS + 1 - j1;

	LGMtsPtr = LGMCreateLGM_TS(numZ, numG);
	if (LGMtsPtr==NULL)
		return(NULL);

/* fill in G dates */
	LGMtsPtr->Gdate[0] = tNow;
	for (i=i1; i < numK; i++)
		LGMtsPtr->Gdate[i+1-i1] = SKtsPtr->kdate[i];

/* construct G(t) */	
	LGMtsPtr->G[0] = 0.0;
	gprev = 1.0;
	for (i=1; i<numG; i++)
	{	dt = (LGMtsPtr->Gdate[i] - LGMtsPtr->Gdate[i-1])/365.0;
		kap = SKtsPtr->kap[i+i1-1];
		LGMtsPtr->G[i] = LGMtsPtr->G[i-1] - gprev*dt*expint(-kap*dt);
		gprev = gprev*exp(-kap*dt);
	}

/* fill zeta dates */
	LGMtsPtr->zdate[0] = tNow;
	for (j=j1; j < numS; j++)
		LGMtsPtr->zdate[j+1-j1] = SKtsPtr->sdate[j];

/* construct zeta(t) */
	LGMtsPtr->zeta[0] = 0.0;

	gprev = 1.0;
	tprev = tNow;
	i=i1;
	for (j=1; j <= numZ; j++)
	{	cumul = 0.;
		for ( ; i < numK && SKtsPtr->kdate[i] <= LGMtsPtr->zdate[j]; i++)
		{	dt = (double)(SKtsPtr->kdate[i] - tprev);
			kap = SKtsPtr->kap[i];
			cumul = cumul + dt*expint(2.0*dt*kap)/(gprev*gprev);
			gprev = gprev * exp(-dt*kap);
			tprev = SKtsPtr->kdate[i];
		}
		if (i < numK)
			kap = SKtsPtr->kap[i];
		else
			kap = SKtsPtr->kap[numK-1];
		dt = (double)(LGMtsPtr->zdate[j] - tprev);
		cumul = cumul + dt*expint(2.0*dt*kap)/(gprev*gprev);
		gprev = gprev * exp(-dt*kap);
		tprev = LGMtsPtr->zdate[j];
		LGMtsPtr->zeta[j] = LGMtsPtr->zeta[j-1] + 
							SKtsPtr->sig[j-j1-1] * SKtsPtr->sig[j-j1-1] * cumul;
	}

/* validate new ts, and return */
	error = LGMVerifyLGM_TS(LGMtsPtr);
	if (error!=NULL)
		LGMFreeLGM_TS(&LGMtsPtr);
	return (LGMtsPtr);
}

/******************************************************************/
/* This routine creates a sigma kappa term structure, and fills it
in with values equivalent to the zeta G term structure.
Of all the sigma-kappa term structures that match the LGM zeta and
G values, it chooses the least oscillatory
On error, it returns a NULL */

/*	Modified AS -still doesn't work in all cases, so we keep track of
	the original Autocal TS, in order to do tests directly on it	*/

static double integ (double dt, double arg)
{
	if (fabs (arg) > 1.0e-06)
	{
		return (exp (arg * dt) - 1.0) / arg;
	}

	return dt;
}

static double dg_func (double dt, double lam, double alpha, double beta)
{
	return alpha * beta * integ (dt, lam);
}

#define MAXITER 50

static double find_lambda (double dt, double dg, double alpha, double beta)
{
	double a[3] = {0.05, 0.10, 0.15}, b[3];
	double nstop;
	int i;

	for (i=0; i<2; i++)
	{
		b[i] = dg_func (dt, a[i], alpha, beta);
		if (fabs (b[i] - dg) < 1.0e-06)
		{
			return a[i];
		}
	}

	i = 0;
	nstop = 0.0;
	while	(nstop < 1.0 && i < MAXITER)
	{
		b[2] = dg_func (dt, a[2], alpha, beta);
		i++;
		newton (dg, 2.0, a, b, &nstop);
	}

	return a[2];
}

static double find_alpha (double dt12, double dt23, double dg12, double dg23, double *lambda, double *beta)
{
	double a[3] = {0.90, 1.00, 1.10}, b[3];
	double nstop;
	double l12, l23;
	int i;

	for (i=0; i<2; i++)
	{
		l23 = find_lambda (dt23, dg23, a[i], 1.0);
		l12 = find_lambda (dt12, dg12, a[i], exp (l23 * dt23));
		b[i] = l23 - l12;
		if (fabs (b[i]) < 1.0e-06)
		{
			return a[i];
		}
	}

	i = 0;
	nstop = 0.0;
	while	(nstop < 1.0 && i < MAXITER)
	{
		l23 = find_lambda (dt23, dg23, a[2], 1.0);
		l12 = find_lambda (dt12, dg12, a[2], exp (l23 * dt23));
		b[2] = l23 - l12;
		i++;
		newton (0.0, 2.0, a, b, &nstop);
	}

	*lambda = l23;
	*beta = exp (l23 * dt23);
	return a[2];
}

static double interpolate_lambda (long t1, long t2, long *ti, double *l, int n)
{
	int i = 0;

	if (n < 2)
	{
		return l[0];
	}

	while (i < n-1 && ti[i] <= t1)
	{
		i++;
	}

	if (i == n - 1 || ti[i] >= t2)
	{
		return l[i];
	}

	return	((double) (ti[i] - t1) / (t2 - t1)) * l[i] 
		+	((double) (t2 - ti[i]) / (t2 - t1)) * l[i+1];
}

static double find_sigma (double dt, double dz, double lambda, double alpha, double beta)
{
	return sqrt (alpha * alpha * beta * beta * dz / integ (dt, - 2 * lambda));
}

SigKapTSPtr LGMConvertZGtoSigKap (Date tNow, LGM_TS *input)
{	
	SigKapTSPtr	output;
	Err			error;
	
	long		i;
	double		alpha, lambda, beta;
	
	int			numG;
	long		*Gdate;
	double		*G;

	numG = input->numG;
	Gdate = input->Gdate;
	G = input->G;

	/*	GET THE EQUIVALENT TAU TERM STRUCTURE	*/

	/*	Number of Gs has to be at least two	*/

	if (numG < 2)
	{
		return NULL;
	}

	/*	Case 1: there is only two Gs, thus a single mean-reversion that can be chosen
		arbitrarily, so we choose 0	*/

	if (numG == 2)
	{
		output = LGMCreateSigKapTS (input->numZ-1, 1);
		if (!output)
		{	
			return NULL;
		}
	
		output->kdate[0] = Gdate[1];
		output->kap[0] = 0.00;

		alpha = (G[0] - G[1]) / ((Gdate[1] - Gdate[0]) * YEARS_IN_DAY);
	}

	else

	/*	Case 2: there is more than two Gs	*/
	{
		output = LGMCreateSigKapTS (input->numZ-1, numG-1);
		if (!output)
		{	
			return NULL;
		}

		/*	Set alpha so that the two last lambdas match	*/
		alpha = find_alpha (	(Gdate[numG-2] - Gdate[numG-3]) * YEARS_IN_DAY,
								(Gdate[numG-1] - Gdate[numG-2]) * YEARS_IN_DAY,
								G[numG-3] - G[numG-2],
								G[numG-2] - G[numG-1],
								&lambda,
								&beta);

		/*	Find the corresponding lambdas	*/
		beta = 1.0;
		for (i = numG - 1; i >= 1; i--)
		{
			output->kdate[i-1] = Gdate[i];
			output->kap[i-1] = find_lambda (	(Gdate[i] - Gdate[i-1]) * YEARS_IN_DAY,
												G[i-1] - G[i],
												alpha,
												beta);
			beta *= exp (output->kap[i-1] * (Gdate[i] - Gdate[i-1]) * YEARS_IN_DAY);
		}
	}

	/*	GET THE EQUIVALENT SIGMA TERM STRUCTURE	*/

	lambda = interpolate_lambda (	input->zdate[input->numZ-1],
									Gdate[numG-1],
									output->kdate,
									output->kap,
									output->numK);

	beta = exp (lambda * (Gdate[numG-1] - input->zdate[input->numZ-1]) * YEARS_IN_DAY);

	for (i = input->numZ-1; i >= 1; i--)
	{
		lambda = interpolate_lambda (	input->zdate[i],
										input->zdate[i-1],
										output->kdate,
										output->kap,
										output->numK);

		output->sdate[i-1] = input->zdate[i];
		output->sig[i-1] = find_sigma (	(input->zdate[i] - input->zdate[i-1]) * YEARS_IN_DAY,
										input->zeta[i] - input->zeta[i-1],
										lambda,
										alpha,
										beta);

		beta *= exp (lambda * (input->zdate[i] - input->zdate[i-1]) * YEARS_IN_DAY);
	}

	error = LGMVerifySigKapTS (output);
	if (error)
	{
		LGMFreeSigKapTS (&output);
		return NULL;
	}

	return output;
}

#undef MAXITER 

/******************************************************************/
/* Gets sigma(t) = sigma[j] for tau[j]<t<tau[j+1] for j=0,..., ntau-1 */
/* from zeta[j] = zeta(t) at t=tau[j] for j=0,..., ntau */
/* and from kap[i] = kappa(t) for t[i]<t<t[i+1] for i=0,...,numK */
/* WARNING: zeta must be scaled so that littleg(t[0]) = g0 */
static void GetSigFromZeta(double g0, long ntau, double *sig,
						   Date *tau, double *zeta,
						   long numK, Date *t, double *kap)
{	Date tprev;
	long i, j;
	double dt, gprev, cumul;

	tprev = t[0];
	gprev = g0;
	i=1;

	for (j=0; j<=ntau; j++)
	{	cumul = 0.0;
		for ( ; i<=numK && t[i]<=tau[j]; i++)
		{	dt = ((double)(t[i]-tprev))/365.0;
			cumul = cumul +	dt*expint(2.0*kap[i-1]*dt)/(gprev*gprev);
			gprev = gprev*exp(-dt*kap[i-1]);
			tprev = t[i];
		}
		
		if (i>numK)
			i=numK;
		dt = ((double)(tau[j]-tprev))/365.0;
		cumul = cumul +	dt*expint(2.0*kap[i-1]*dt)/(gprev*gprev);
		gprev = gprev*exp(-dt*kap[i-1]);
		tprev = tau[j];
		
		if (j>0)
			sig[j-1] = sqrt((zeta[j]-zeta[j-1])/cumul);
	}
	return;
}

/******************************************************************/
/* compute integral of exponential */
static double expint(double arg)
{	double ans;
	if (fabs(arg)< 1.e-6)
		ans = 1.0 + 0.5*arg + arg*arg/6.0;
	else
		ans = (exp(arg) - 1.0)/arg;
	return (ans);
}

/* compute the derivative of the above function */
static double expintder(double arg)
{	double ans;
	if (fabs(arg)< 1.e-6)
		ans = 0.5 + arg/3.0 + 0.125*arg*arg;
	else
		ans = ((arg-1.0)*exp(arg) + 1.0)/(arg*arg);
	return (ans);
}

/***************************************************************/
/* Gets kappa(t) = kappa[i] for t[i]<t<t[i+1] for i=0,..., n-1 */
/* Determines it from G[i] = G(t) at t=t[i] for i=0,..., n */
static void GetKapFromG(double b, long n, Date *t, double *G, double *kap)
{	double littleg, dt, ratio, theta;
	long i;

	littleg = b;
	for (i=n-1; i>=0; i--)
	{	dt = ((double)(t[i+1] - t[i]))/365.0;
		ratio = (G[i]-G[i+1])/(dt*littleg);
		theta =	FindTheta(ratio);
		kap[i] = theta/dt;
		littleg = littleg*exp(theta);
	}
	return;
}

/***************************************************************/
/* Use global Newton to solve (exp(theta)-1)/theta = ratio */
static double FindTheta(double ratio)
{	double theta, dtheta, fold, fnew;
	long iter, itermax, success;

/* initialize */
	itermax = 30;
	success = -2;		/* flag as no convergence */

    theta = 0.0;		/* intial guess */
	dtheta = 0.0;
    
/* global Newton */
    for (iter=1; iter<itermax && success<=0; iter++)
    {	fnew = expint(theta+dtheta) - ratio;
		if (iter>1 && (fabs(dtheta)<1.e-06 || fabs(fnew)<1.e-06))
			success++;							/* check for convergence */
		if (iter==1 || fabs(fnew)<fabs(fold))	/* see if step is successful */
		{	fold = fnew;
			theta = theta + dtheta;				/* if so, take another Newton step */
			dtheta = -fnew/expintder(theta);
		}
		else
			dtheta = dtheta/2.;					/* else cut step in half */

    }
	if (success<0)
		return (0.);
	else
		return (theta+dtheta);					/* take one more step & return */
}

/***********************************************************/
/* compute the time-average value of kappa(t) */
static double computeavgkap(long n, Date *t, double *kap)
{	double totalkap, totaltime, dt;
	long i;
	totalkap = 0.0;
	totaltime = 0.0;
	for (i=0; i<n; i++)
	{	dt = ((double)(t[i+1]-t[i]))/365.0;
		totaltime = totaltime + dt;
		totalkap = totalkap + kap[i]*dt;
	}
	return(totalkap/totaltime);
}

/***********************************************************/
/* find largest deviation from the average */
static long GetWorstDif(long n, double *Arr, double avg)
{	long i, iworst;
	double worst;

	worst = 0.;
	for (i=0; i<n; i++)
	{	if (fabs(Arr[i]-avg) >= worst)
		{	worst = fabs(Arr[i]-avg);
			iworst = i;
		}
	}
	return (iworst);
}

/***********************************************************/
/* find new (delta b)/b */
static double GetDeltaLogb(long n, double *Arr, Date *tArr, double avg, long iw)
{	double worst, dt, dtworst, cor, mincor;
	double isign, isignw;
	long i;	
			
	isignw = ((n-1-iw)%2) ? -1.0 : 1.0;
	worst = isignw*(Arr[iw] - avg);
	dtworst = ((double)(tArr[iw+1] - tArr[iw]))/365.0;
	mincor = 0.5*worst*dtworst;

	isign = -1.0;
	for (i=n-1; i>=0; i--)
	{	isign = -isign;
		dt = ((double)(tArr[i+1] - tArr[i]))/365.0;
		cor = worst + isign*(Arr[i] - avg);
		cor = 0.5 * cor / ((1.0/dt) + (1.0/dtworst));
		if (fabs(cor) < fabs(mincor))
			mincor = cor;
	}
	return (mincor);
}

/***********************************************************/
/******* Create & Free LGM Calibration Set *****************/
LGMCalSetPtr LGMCreateCalSet(long nPay, long nEx, long nLast)
{	LGMCalSet *ptr;
 
	if (nPay<1 || nPay>110 || nEx<1 || nEx>600)
		return (NULL);

	ptr = (LGMCalSet*) srt_calloc(1, sizeof(LGMCalSet));
	if (ptr==NULL)
		return (NULL);

/* common fixed leg pay dates */
	ptr->tPay = (Date*) srt_calloc(nPay+1, sizeof(Date));
	ptr->cvgpay = (double*) srt_calloc(nPay+1, sizeof(double));
	ptr->Dpay = (double*) srt_calloc(nPay+1, sizeof(double));

/* common exercise and start dates */
	ptr->tEx = (Date*) srt_calloc(nEx+1, sizeof(Date));
	ptr->tStart = (Date*) srt_calloc(nEx+1, sizeof(Date));
	ptr->DStart = (double*) srt_calloc(nEx+1, sizeof(double));

/* long swaptions */
	ptr->ifirst = (long*) srt_calloc(nEx+1, sizeof(long));
	ptr->nlong = (long*) srt_calloc(nEx+1, sizeof(long));
	ptr->cvgfirst = (double*) srt_calloc(nEx+1, sizeof(double));
	ptr->Rflong = (double*) srt_calloc(nEx+1, sizeof(double));
	ptr->Vlong = (double*) srt_calloc(nEx+1, sizeof(double));

/* short swaptions */
	ptr->nshort = (long*) srt_calloc(nEx+1, sizeof(long));
	ptr->Rfshort = (double*) srt_calloc(nEx+1, sizeof(double));
	ptr->Vshort = (double*) srt_calloc(nEx+1, sizeof(double));

/* caplets */
	ptr->tEnd = (Date*) srt_calloc(nEx+1, sizeof(Date));
	ptr->cvgcap = (double*) srt_calloc(nEx+1, sizeof(double));
	ptr->Dcap = (double*) srt_calloc(nEx+1, sizeof(double));
	ptr->Rfcap = (double*) srt_calloc(nEx+1, sizeof(double));
	ptr->Vcap = (double*) srt_calloc(nEx+1, sizeof(double));

/* 1 into k swaptions */
	ptr->Rfix = (double*) srt_calloc(nLast+1, sizeof(double));
	ptr->Vfix = (double*) srt_calloc(nLast+1, sizeof(double));

	if (ptr->tPay==NULL		|| ptr->cvgpay==NULL || ptr->Dpay==NULL     || 
		ptr->tEx==NULL		|| ptr->tStart==NULL || ptr->DStart==NULL   ||
		ptr->ifirst==NULL	|| ptr->nlong==NULL  || ptr->cvgfirst==NULL ||
		ptr->Rflong==NULL	|| ptr->Vlong==NULL  || ptr->nshort==NULL   ||
		ptr->Rfshort==NULL	|| ptr->Vshort==NULL || ptr->tEnd==NULL     ||
		ptr->cvgcap==NULL	|| ptr->Dcap==NULL   || ptr->Rfcap==NULL    ||
		ptr->Vcap==NULL		||  ptr->Rfix==NULL  || ptr->Vfix==NULL       )

	{	LGMFreeCalSet(&ptr);
		return(NULL);
	}

	ptr->nPay = nPay;
	ptr->nEx = nEx;
	return (ptr);
}

void LGMFreeCalSet(LGMCalSet **ptr)
{	if ((*ptr)==NULL) return;

	if ((*ptr)->tPay) srt_free((*ptr)->tPay);
	if ((*ptr)->cvgpay) srt_free((*ptr)->cvgpay);
	if ((*ptr)->Dpay) srt_free((*ptr)->Dpay); 
	if ((*ptr)->tEx) srt_free((*ptr)->tEx);
	if ((*ptr)->tStart) srt_free((*ptr)->tStart);
	if ((*ptr)->DStart) srt_free((*ptr)->DStart);
	if ((*ptr)->ifirst) srt_free((*ptr)->ifirst);
	if ((*ptr)->nlong) srt_free((*ptr)->nlong);
	if ((*ptr)->cvgfirst) srt_free((*ptr)->cvgfirst);
	if ((*ptr)->Rflong) srt_free((*ptr)->Rflong);
	if ((*ptr)->Vlong) srt_free((*ptr)->Vlong);
	if ((*ptr)->VegaLong) srt_free((*ptr)->VegaLong);
	if ((*ptr)->FrLong) srt_free((*ptr)->FrLong);
	if ((*ptr)->BetaLong) srt_free((*ptr)->BetaLong);
	if ((*ptr)->CEVLong) srt_free((*ptr)->CEVLong);
	if ((*ptr)->StDLong) srt_free((*ptr)->StDLong);
	if ((*ptr)->nshort) srt_free((*ptr)->nshort);
	if ((*ptr)->Rfshort) srt_free((*ptr)->Rfshort);
	if ((*ptr)->Vshort) srt_free((*ptr)->Vshort);
	if ((*ptr)->tEnd) srt_free((*ptr)->tEnd);
	if ((*ptr)->cvgcap) srt_free((*ptr)->cvgcap);
	if ((*ptr)->Dcap) srt_free((*ptr)->Dcap);	
	if ((*ptr)->Rfcap) srt_free((*ptr)->Rfcap);
	if ((*ptr)->Vcap) srt_free((*ptr)->Vcap);	
	if ((*ptr)->VegaCap) srt_free((*ptr)->VegaCap);	
	if ((*ptr)->Rfix) srt_free((*ptr)->Rfix);	
	if ((*ptr)->Vfix) srt_free((*ptr)->Vfix);

	(*ptr)->tPay = NULL;
	(*ptr)->cvgpay = NULL;
	(*ptr)->Dpay = NULL;
	(*ptr)->tEx = NULL;
	(*ptr)->tStart = NULL;
	(*ptr)->DStart = NULL;
	(*ptr)->ifirst = NULL;
	(*ptr)->nlong = NULL;
	(*ptr)->cvgfirst = NULL;
	(*ptr)->Rflong = NULL;
	(*ptr)->Vlong = NULL;
	(*ptr)->VegaLong = NULL;
	(*ptr)->nshort = NULL;
	(*ptr)->Rfshort = NULL;
	(*ptr)->Vshort = NULL;
	(*ptr)->tEnd = NULL;
	(*ptr)->cvgcap = NULL;
	(*ptr)->Dcap = NULL;
	(*ptr)->Rfcap = NULL;
	(*ptr)->Vcap = NULL;
	(*ptr)->VegaCap = NULL;
	(*ptr)->Rfix = NULL;
	(*ptr)->Vfix = NULL;

	srt_free((*ptr));
	*ptr = NULL;
	return;
}

/******************************************************/
/****** Create & Free LGM Calibration Parameters ******/
/* Create default calibration parameter structure */
LGMCalParmPtr LGMCreateCalParm(long numG, long numZ, long numR1, long numR2)
{	LGMCalParm *ptr = NULL;

	ptr = (LGMCalParm*) srt_calloc(1,sizeof(LGMCalParm));
	if (ptr == NULL)
		return NULL;
	
	ptr->calmeth = FixExp;
	ptr->Rmeth = dIRR;
	ptr->respectExer = 1;
	ptr->usecaps = 0;
	ptr->MinMonToEx = 6;
	ptr->keep1intok = 0;
	ptr->kap = 0.;
	ptr->usestarts = 1;

	if (numR1>0 || numR2>0)
		ptr->Rmeth = givenR;

/* G information */
	ptr->numG = 0;
	ptr->Gdate = NULL;
	ptr->G = NULL;
	if (numG>0)
	{	ptr->numG = numG;
		ptr->Gdate = (Date*) srt_calloc(numG, sizeof(Date));
		ptr->G = (double*) srt_calloc(numG, sizeof(double));
		if (ptr->Gdate == NULL || ptr->G == NULL)

		{	LGMFreeCalParm(&ptr);
			return(NULL);
		}

	}
/* zeta information */
	ptr->numZ = 0;
	ptr->zdate = NULL;
	ptr->zeta = NULL;
	if (numZ>0)
	{	ptr->numZ = numZ;
		ptr->zdate = (Date*) srt_calloc(numZ, sizeof(Date));
		ptr->zeta = (double*) srt_calloc(numZ, sizeof(double));
		if (ptr->zdate==NULL || ptr->zeta==NULL)
		{	LGMFreeCalParm(&ptr);
			return(NULL);
		}
	}
/* strike information */
	ptr->numR1 = 0;
	ptr->Rdate1 = NULL;
	ptr->Rdata1 = NULL;
	if (numR1>0)
	{	ptr->numR1 = numR1;
		ptr->Rdate1 = (Date*) srt_calloc(numR1, sizeof(Date));
		ptr->Rdata1 = (double*) srt_calloc(numR1, sizeof(double));
		if (ptr->Rdate1 == NULL || ptr->Rdata1 == NULL)
		{	LGMFreeCalParm(&ptr);
			return(NULL);
		}
	}
	ptr->numR2 = 0;
	ptr->Rdate2 = NULL;
	ptr->Rdata2 = NULL;
	if (numR2>0)
	{	ptr->numR2 = numR2;
		ptr->Rdate2 = (Date*) srt_calloc(numR2, sizeof(Date));
		ptr->Rdata2 = (double*) srt_calloc(numR2, sizeof(double));
		if (ptr->Rdate2==NULL || ptr->Rdata2==NULL)
		{	LGMFreeCalParm(&ptr);
			return(NULL);
		}
	}
	return(ptr);
}
 
/********************************************************/
/* Free calibration parameter structure */
void LGMFreeCalParm(LGMCalParm **ptr)
{	if ((*ptr)==NULL) return;

	if ((*ptr)->Gdate!=NULL) srt_free((*ptr)->Gdate);
	if ((*ptr)->G!=NULL) srt_free((*ptr)->G);
	if ((*ptr)->zdate!=NULL) srt_free((*ptr)->zdate);
	if ((*ptr)->zeta!=NULL) srt_free((*ptr)->zeta);
	
	if ((*ptr)->Rdate1!=NULL) srt_free((*ptr)->Rdate1);
	if ((*ptr)->Rdata1!=NULL) srt_free((*ptr)->Rdata1);
	if ((*ptr)->Rdate2!=NULL) srt_free((*ptr)->Rdate2);
	if ((*ptr)->Rdata2!=NULL) srt_free((*ptr)->Rdata2);

	(*ptr)->Gdate=NULL;
	(*ptr)->G=NULL;
	(*ptr)->zdate=NULL;
	(*ptr)->zeta=NULL;
	(*ptr)->Zeta1=NULL;

	(*ptr)->Rdate1=NULL;
	(*ptr)->Rdata1=NULL;
	(*ptr)->Rdate2=NULL;
	(*ptr)->Rdata2=NULL;

	srt_free((*ptr));
	*ptr = NULL;
	return;
}

/********************************************************/
/******* Swap/Cap Utilities******************************/
/********************************************************/
/* OVE: gets the market conventions:
		- either from the YC passed (if initialised)
		- or     from the default CCY parameters
*/
LGMErr LGMCcyDefaults(SrtCurvePtr yldcrv, LGMMarkConv *conventions)
{	String		ccy_str;
	SrtCcyParam	*ccy_param = NULL;
	Err			err;

	ccy_param = get_ccyparam_from_yldcrv(yldcrv);
	if (!ccy_param)
	{
		ccy_str = get_curve_ccy(yldcrv);
		err = swp_f_get_CcyParam_from_CcyStr(ccy_str, &ccy_param);
		if (err)
			return (err);
	}
	conventions->ccy = ccy_param->currency;
	conventions->lag = ccy_param->spot_lag;

/* swaption conventions */
	conventions->sbasis = ccy_param->swap_basis_code;
	conventions->sfreq = ccy_param->compd;
	conventions->sbdconv = ccy_param->swap_bus_day_conv;

/* caplet conventions */
	conventions->cbasis = ccy_param->cash_basis_code;
	conventions->cfreq = SRT_QUARTERLY;
	conventions->cbdconv = ccy_param->cash_bus_day_conv;

/*	minimum length of swap (in fixed leg periods) */
	conventions->minswap = 2;
	if (conventions->sfreq==SRT_ANNUAL)
		conventions->minswap = 1;

	return NULL;
}


LGMErr LGMMuniDefaults(char * szRefRate, LGMMarkConv *conventions)
{	
	SrtCcyParam	*ccy_param = NULL;
	Err			err;

	char csCurrency[] = "UTF";
	if ( err = swp_f_get_CcyParam_from_CcyStr(csCurrency, &ccy_param) )
		return err;

	conventions->ccy = ccy_param->currency;
	conventions->lag = ccy_param->spot_lag;

/* swaption conventions */
	conventions->sbasis = ccy_param->swap_basis_code;
	conventions->sfreq = ccy_param->compd;
	conventions->sbdconv = ccy_param->swap_bus_day_conv;

/* caplet conventions */
	conventions->cbasis = ccy_param->cash_basis_code;
	conventions->cfreq = SRT_QUARTERLY;
	conventions->cbdconv = ccy_param->cash_bus_day_conv;

/*	minimum length of swap (in fixed leg periods) */
	conventions->minswap = 2;
	if (conventions->sfreq==SRT_ANNUAL)
		conventions->minswap = 1;

	return NULL;
}

/********************************************************/
/* Construct the pay dates for a standard fixed leg with exercise
date tEx and last (theor) pay date tEnd, using market conventions in *conv.
If extraperiods>0, add extraperiod standard periods after tEnd */
/* Returns an array with tPay[0,1,...,nPay] where nPay = *nPayPtr, 
where tPay[0] is the start date, and where tPay[i] is the i-th pay date
On error, returns NULL */
Date *LGMFixLegSched(Date tEx, Date tEnd, long extraperiods,
					 LGMMarkConv *conv, long *nPayPtr)
{
	LGMErr error;
	long i, n, nPay;
	Date tStart, *tPay;
	long swmonth;
	CcyCode	ccy;

	error = NULL;
	ccy = conv->ccy;

	swmonth = 6;						/* swmonth is number of months per fixed leg period */
	if (conv->sfreq==SRT_ANNUAL)
		swmonth = 12;

/* Find spot of exercise */
	tStart = add_unit(tEx, conv->lag, SRT_BDAY, SUCCEEDING/*, ccy */);
	if (tEnd<tStart+10 || tEnd>tStart + LGMUTIL_MAXNB_PERIOD * 180)
		return(NULL);					/* serious date error */

/* Find out how many periods in the fixed leg */
	n=0;
	while ((n<=LGMUTIL_MAXNB_PERIOD) && (tStart < add_unit(tEnd, -n*swmonth, SRT_MONTH, conv->sbdconv)))
		n++;
	
	nPay = n;				
	if (extraperiods>0)
		nPay = nPay+extraperiods;
	if (nPay>=LGMUTIL_MAXNB_PERIOD || nPay<1)
		return (NULL);					/* serious date error */

/* Allocate space for fixed leg */
	tPay = (Date*) srt_calloc(nPay+1, sizeof(Date));
	if (tPay==NULL)
		return (NULL);					/* allocation failed */

/* Construct fixed leg */
	tPay[0] = tStart;
	for (i=1; i<=n; i++)
		tPay[i] = add_unit(tEnd, -(n-i)*swmonth, SRT_MONTH, conv->sbdconv/*, ccy */);

/* add extra periods if needed */
	for (i=1; i<=extraperiods; i++)
		tPay[n+i] = add_unit(tEnd, i*swmonth, SRT_MONTH, conv->sbdconv/*, ccy */);

	*nPayPtr = nPay;
	return (tPay);
}

/********************************************************/
/* Construct the coverages for a standard fixed leg with
start date tPay[0] and paydates tPay[i], i=1,...,nPay
   Returns an array with cvg[0,1,2,...,nPay] where
cvg[0] is UNDEFINED and cvg[i] is the cvg from tPay[i-1] to tPay[i]
On error, returns NULL */
double *LGMFixLegCvg(long nPay, Date *tPay, LGMMarkConv *conv)
{
	double *cvgPay;
	long i;

	if (nPay<1 || tPay==NULL)
		return (NULL);

	cvgPay = (double*) srt_calloc(nPay+1, sizeof(double));
	if (cvgPay==NULL)
		return(NULL);

	for (i=1; i<=nPay; i++)
		cvgPay[i] = coverage(tPay[i-1], tPay[i], conv->sbasis);

	return(cvgPay);
}

/********************************************************/
/* Construct the discount factors for the fixed leg with
start date tPay[0] and paydates tPay[i], i=1,...,nPay
/* Returns an array with DPay[0,1,2,...,nPay] where
DPay[0] is the discount factor to the start date, and
DPay[i] is the discount factor to tPay[i]
On error, returns NULL */
double *LGMFixLegDF(long nPay, Date *tPay, Date tNow, String ycName)
{
	double *DPay;
	long i;

	if (nPay<1 || tPay==NULL)
		return (NULL);

	DPay = (double*) srt_calloc(nPay+1, sizeof(double));
	if (DPay==NULL)
		return(NULL);

	for (i=0; i<=nPay; i++)
	{	DPay[i] = swp_f_df(tNow, tPay[i], ycName);
		if (DPay[i] == SRT_DF_ERROR)
		{	srt_free(DPay);
			return(NULL);
		}
	}
	return(DPay);
}

/********************************************************/
/* Construct the array of G(tPay) for the fixed leg with
start date tPay[0] and paydates tPay[i], i=1,...,nPay
/* Returns an array with GPay[0,1,2,...,nPay] where
GPay[0] is G(t) at the start date tPay[0], and
GPay[i] is G(t) at the pay date tPay[i]
On error, returns NULL */
double *LGMFixLegGs(long nPay, Date *tPay, LGM_TS *tsPtr)
{
	double *GPay;
	long i;

	if (nPay<1 || tPay==NULL || tsPtr==NULL)
		return (NULL);

	GPay = (double*) srt_calloc(nPay+1, sizeof(double));
	if (GPay==NULL)
		return(NULL);

	for (i=0; i<=nPay; i++)
		GPay[i] = LGMGFromTS(tPay[i], tsPtr);

	return(GPay);
}

/********************************************************/
/* Converts forward value of fixed leg into equivalent
fixed rate Rfix for a cash swaption */
LGMErr LGMGetFixRate(Date tEx, Date tEnd, double FwdVal,	/* definition of swaption */
					 String ycName,							/* yield curve name */
					 double *Rfix)							/* output */
{
	LGMErr error = NULL;
	LGMMarkConv conv;
	double DStart, DEnd, level;
	double *DPay = NULL, *cvgPay = NULL;
	Date tNow;
	Date *tPay = NULL;
	long i, nPay, extraperiods;
	SrtCurvePtr yldcrv = NULL;

	extraperiods = 0;
	yldcrv = lookup_curve(ycName);
	tNow = get_clcndate_from_yldcrv(yldcrv); /* real eval date irrelevent */
 
/* Get standard market conventions from market pointer */
	error = LGMCcyDefaults(yldcrv, &conv);
	if (error!=NULL)
		return(error);
	
/* Set up standard swap */
	tPay = LGMFixLegSched(tEx, tEnd, extraperiods, &conv, &nPay);

	if (tPay==NULL)
		return("Schedule failed in GetFixRate");

	cvgPay = LGMFixLegCvg(nPay, tPay, &conv);
	DPay = LGMFixLegDF(nPay, tPay, tNow, ycName);
	if (cvgPay==NULL || DPay==NULL)
	{	srt_free(tPay);
		srt_free(cvgPay);
		srt_free(DPay);
		return("Could not generate fixed leg in ExerInto");
	}

	DStart = DPay[0];		/* start date is tPay[0] */
	DEnd = DPay[nPay];

	level = 0.0;
	for (i=1; i<=nPay; i++)
		level = level + DPay[i]*cvgPay[i];
	
	*Rfix = (FwdVal*DStart-DEnd)/level;
	srt_free(tPay);
	srt_free(cvgPay);
	srt_free(DPay);
	return (error);
}

LGMErr LGMGetRfixFromX(Date tEx, Date tEnd, double x,			/* definition of swaption */ 
				  Date tNow, String ycName, LGM_TS *tsPtr, 	/* yield curve name */
				  double *Rfix)									/* output */
{
	LGMErr error;
	LGMMarkConv conv;
	double level, zeta, arg, expfac;
	double *DPay, *cvgPay, *GPay;
	Date *tPay;
	long i, nPay, extraperiods;
	SrtCurvePtr yldcrv = NULL;

	error = NULL;
	extraperiods = 0;
	yldcrv = lookup_curve(ycName);
 
/* Get standard market conventions from market pointer */
	error = LGMCcyDefaults(yldcrv, &conv);
	if (error!=NULL)
		return(error);
	
/* Set up standard swap */
	tPay = LGMFixLegSched(tEx, tEnd, extraperiods, &conv, &nPay);

	if (tPay==NULL)
		return("Schedule failed in GetFixRate");

/* get coverages, discount factors and G at pay dates */
	cvgPay = LGMFixLegCvg(nPay, tPay, &conv);
	DPay = LGMFixLegDF(nPay, tPay, tNow, ycName);
	GPay = LGMFixLegGs(nPay, tPay, tsPtr);
	if (cvgPay==NULL || DPay==NULL || GPay==NULL)
	{	srt_free(tPay);
		srt_free(cvgPay);
		srt_free(DPay);
		srt_free(GPay);
		return("Could not generate fixed leg for exer bdry");
	}
/* get zeta */
	zeta = LGMZetaFromTS(tEx, tNow, tsPtr);

/* calculate level at state variable x */
	level = 0.0;
	for (i=1; i<=nPay; i++)
	{	arg = (GPay[i]-GPay[0])*(x - 0.5*(GPay[i]+GPay[0])*zeta);
		expfac = exp(arg);
		level = level + DPay[i]*cvgPay[i]*expfac;
	}
	
	*Rfix = (DPay[0] - DPay[nPay]*expfac)/level;

/* free and return */
	srt_free(tPay);
	srt_free(cvgPay);
	srt_free(DPay);
	srt_free(GPay);

	return (error);
}



LGMErr LGMGetFloatRateFromX(Date tStart, Date tPay, double x,			
				  String ycName, LGM_TS *tsPtr, 	
				  double *FloatRate)
{

Err			 err=NULL;
SrtCurvePtr  yldcrv = NULL;
LGMMarkConv  conv;
Date		 tNow;
double		 dZetaStart;
double		 dGPay;
double		 dGStart;
double		 dDfStart,dDfPay,dDfStartPay;
double		 dCoverage;


	yldcrv = lookup_curve(ycName);

	err = LGMCcyDefaults(yldcrv, &conv);
	if (err!=NULL)
		return(err);
	
	yldcrv = lookup_curve(ycName);
	tNow = get_clcndate_from_yldcrv(yldcrv);
	
	dZetaStart = LGMZetaFromTS(tStart,tNow,tsPtr);
	dGPay =  LGMGFromTS(tPay,tsPtr);
	dGStart =  LGMGFromTS(tStart,tsPtr);

	dDfStart = swp_f_df(tNow,tStart,ycName);
	dDfPay= swp_f_df(tNow,tPay,ycName);

	dDfStartPay = (dDfPay/dDfStart)*exp((dGPay-dGStart)*x-0.5*dZetaStart*(dGPay*dGPay-dGStart*dGStart));

	dCoverage = coverage(tStart,tPay,conv.sbasis);
	*FloatRate = (1/dDfStartPay-1)/dCoverage;

	return err;

}

/********************************************************/

/****************************************************************/
/****** Math Utilities ******************************************/
/****************************************************************/
/* calculates the Gaussian probability density, protected against underflows */
double LGMsafeGauss(double z)
{	static double premult = .398942280401433;
	double arg, answer;

	arg = 0.5*z*z;
	if (arg > 50.0)
		answer = premult*exp(-50.);
	else
		answer = premult*exp(-arg);
	return (answer);
}

/****************************************************************/
/* calculates the cumulative normal distribution, protected against underflows */
double LGMsafeNorm(double z)
{	if (z<-10.)
		z = -10.0;
	else if (z>10.)
		z = 10.0;
	return(norm(z));
}

/***********************************************/
/* Gets the equivalent normal vol from CEV vol */
double LGMNormalVolFromCEV(
			double fwd,		/* today's forward price (forwarded to settlement date) */
			double strike,	/* strike */
			double mat,		/* time from today to expiry, in years */
			double CEVvol,	/* volatility parameter in CEV model */
			double beta)	/* exponenet in CEV model */
{	double fav, x, fac, normvol;

	if (beta > 0.98) beta = 1.0;			/* input is log normal vol */
	if (beta < 0.005) return (CEVvol);		/* input is normal vol */

	if (fwd<0.0001) fwd=0.0001;
	if (strike<0.0001) strike=0.0001;

	fav = (fwd+strike)/2.0;
	x = (fwd-strike)/(fwd+strike);

/* special case: near ATM */
	if (fabs(x)<0.001)								 
	{	fac = 1. + (2.0+beta)*(3.0+beta)*x*x/20.0;
		fac = 1. + beta*(1.0+beta)*x*x*fac/6.0;
		fac = fac/pow(fav, beta);
	}
/* special case: log normal	*/
	else if (beta>.98)						
		fac = (log(fwd/strike))/(fwd-strike);

/* general case */
	else
		fac = (pow(fwd, 1.0-beta) - pow(strike, 1.0-beta))/((1.0-beta)*(fwd-strike));

/* get the normal volatility */
	normvol = CEVvol/fac;
	normvol = normvol/(1.0 + beta*(2.0-beta)*normvol*normvol*mat/(24.0*fav*fav));
	/* SEE 4.21 A FORMULA  */
	return(normvol);	
}

/*************************************************************/
/* Returns the price of a European option for a normal model */
/* On error, returns -1 */
double LGMNormOptPrice(double fwd, double strike, double NormVol,
					   double mat, double dftoStart, 
					   SrtCallPutType callput)
{	
	double stddev, theta, price;

/* check */
	if (mat<0.)						/* exercise date has passed */
		return(0.);					
	if (callput>=SRT_LASTCALLPUTTYPE)	/* unknown option type */
		return (-1.);

/* get time value */
	stddev = sqrt(NormVol*NormVol*mat);
	if ((fabs(fwd-strike)>10.0*stddev) || stddev<=1.0e-12)
		theta = -10.;
	else
		theta = -fabs(fwd-strike)/stddev;
	
	price = stddev*(LGMsafeGauss(theta) + theta*LGMsafeNorm(theta));

/* add intrinsic value */
	if (callput==SRT_CALL && fwd>strike)
		price = price + fwd - strike;

	else if (callput==SRT_PUT && strike>fwd)
		price = price + strike - fwd;

	else if (callput==SRT_STRADDLE)
		price = price + price + fabs(fwd-strike);

	price = dftoStart*price;
	return(price);
}

/*************************************************************/
/* Computes the implied Normal Vol for a European option */
/* On error, returns -1 */
double LGMNormImpVol(double price, double fwd, double strike,
					 double mat, double dftoStart,
					 SrtCallPutType callput)
{
	double x, stddev, vol;
	double arg1, Nprice, Nder, dPrice, dstddev;
	long iter, itermax, success;

	if (mat<=0.0)
		return (-1.);
	if (callput>=SRT_LASTCALLPUTTYPE)	/* unknown option type */
		return (-1.);

/* strip out intrinsic value, and normalize df and strike to 1 */
	price = price/dftoStart;						/* strip out dis fac */
	
	if (callput==SRT_CALL && fwd>strike)						
		price = price - (fwd - strike);			/* strip out intrinsic value */
	if (callput==SRT_PUT && strike<fwd)						
		price = price - (strike - fwd);			/* strip out intrinsic value */
	if (callput==SRT_STRADDLE)
		price = 0.5*(price - fabs(fwd-strike));	/* strip out intrinsic value */

	if (price < 1.0e-08)
	{
		if (price > -1.0e-08)
		{
			return 1.0e-15;
		}
		else
		{
			return(-1.0);
		}
	}
	  
 /* Newton scheme */
	x = fabs(fwd-strike);
	stddev = (price + 0.5*x)/LGMsafeGauss(0);	/* guaranteed to converge from this initial guess */
	arg1 = -x/stddev;
	Nprice = stddev*LGMsafeGauss(arg1)-x*LGMsafeNorm(arg1);
	Nder = LGMsafeGauss(arg1); 

	success = -2;
	itermax = 25;
	for (iter=1; iter<itermax && success<0; iter++)
	{	dPrice = price-Nprice;
		dstddev = dPrice/Nder;
		stddev = stddev + dstddev;				/* new stddev */

/* check for convergence ... take two more steps after satisfactory convergence */		
		if (fabs(dPrice/price)<1.0e-7 || fabs (dPrice)<1.0e-11 || fabs(dstddev/stddev)<1.0e-7)
			success++;
		arg1 = -x/stddev;						/* calc new Black price & derivative */
		Nprice = stddev*LGMsafeGauss(arg1) - x*LGMsafeNorm(arg1);
		Nder = LGMsafeGauss(arg1);
	}

	if (success<-1)								/* failure! */
		return (-1.);
	stddev = stddev + (price-Nprice)/Nder;		/* last improvement */
	vol = stddev/sqrt(mat);
	return (vol);
}

/*****************************************************************/
/* Computes the price of a European option from Black's formula */
/* On error, returns -1 */
double LGMBlackOptPrice(double fwd, double strike, double BlackVol,
						double mat, double dftoStart, 
						SrtCallPutType callput)
 {
	double stddev, theta, arg, price;

	if (mat<0.)											/* exercise date has passed */
		return(0.);					
	if (callput>=SRT_LASTCALLPUTTYPE)					/* unknown option type */
		return (-1.);

	stddev = sqrt(BlackVol*BlackVol*mat);
	arg = log(fwd/strike);

	if ((fabs(arg)>10.0*stddev) || stddev<=1.0e-12)		/* nearly no time value */
	{
		theta = -10;
		price = 0.5*stddev*(fwd+strike)*(LGMsafeGauss(theta) + theta*LGMsafeNorm(theta));
		if (callput==SRT_CALL && fwd>strike)
			price = price + fwd - strike;

		else if (callput==SRT_PUT && strike>fwd)
			price = price + strike - fwd;

		else if (callput==SRT_STRADDLE)
			price = 2.0*price + fabs(fwd-strike);
		else
			;
		price = dftoStart*price;
		return(price);									/* return for "no time value" case */
	}

/* usual case */
 	theta = -fabs(arg)/stddev;

/* puts out of the money */
	if (fwd>=strike)
	{	price = strike*LGMsafeNorm(theta + 0.5*stddev)
							- fwd*LGMsafeNorm(theta - 0.5*stddev);
		if (callput==SRT_CALL)
			price = price + fwd - strike;
		if (callput==SRT_STRADDLE)
			price = price + price + fwd - strike;
		price = dftoStart*price;
		return(price);					/* return "puts out of money" case */
	}

/* calls out of the money */
	price = fwd*LGMsafeNorm(theta + 0.5*stddev)
						- strike*LGMsafeNorm(theta - 0.5*stddev);
	if (callput==SRT_PUT)
		price = price + strike - fwd;
	if (callput==SRT_STRADDLE)
		price = price + price + strike - fwd;
	price = dftoStart*price;
	return(price);
 }

/*****************************************************************/
/* Computes the implied Black Vol for a European option */
/* On error, returns -1 */
double LGMBlackImpVol(double price, double fwd, double strike,
					  double mat, double dftoStart,
					  SrtCallPutType callput)
{
	double x, theta, vol;
	double abslog, arg, Bprice, Bder, dPrice, dtheta;
	long iter, itermax, success;

	if (mat<0.)									/* exercise date has passed */
		return(0.);					
	if (callput>=SRT_LASTCALLPUTTYPE)			/* unknown option type */
		return (-1.);

/* strip out intrinsic value, and normalize df and strike to 1 */
	price = price/dftoStart;					/* strip out dis fac */
	
	if (fwd>=strike)
	{	if (callput==SRT_CALL)						
			price = price - (fwd - strike);		/* strip out intrin value */
		if (callput==SRT_STRADDLE)
			price = 0.5*(price - (fwd-strike));	/* strip out intrin value */
		price = price/fwd;						/* normalize */
		x = strike/fwd;
	}
	else
	{	if (callput==SRT_PUT)
			price = price - (strike - fwd);		/* strip out intrin value */
		if (callput==SRT_STRADDLE)
			price = 0.5*(price - (strike-fwd));	/* strip out intrin value */
		price = price/strike;					/* normalize */
		x = fwd/strike;
	}

	if (price>=x)
		return(-1.0);							/* negative  or infinite time value */

	if (price < 1.0e-08)
	{
		if (price > -1.0e-08)
		{
			return 1.0e-15;
		}
		else
		{
			return(-1.0);
		}
	}
	
 /* Newton scheme */
	abslog = fabs(log(x));
	theta = sqrt(2.0*abslog);					/* guaranteed to converge from */
	Bprice = 0.5*x-LGMsafeNorm(-theta);			/* this initial guess */
	Bder = x*LGMsafeGauss(0); 

	itermax = 25;
	success = -2;
	for (iter=1; iter<itermax && success<0; iter++)
	{	dPrice = price-Bprice;
		dtheta = dPrice/Bder;
		theta = theta + dtheta;					/* new theta */

/* check for convergence ... take two more steps after satisfactory convergence */		
			if (fabs(dPrice/price)<1.0e-7 || fabs (dPrice)<1.0e-11 || fabs(dtheta/theta)<1.0e-7)
			success++;

		arg = -abslog/theta + 0.5*theta;		/* calc new Black price & derivative */
		Bprice = x*LGMsafeNorm(arg) - LGMsafeNorm(arg - theta);
		Bder = x*LGMsafeGauss(arg);
	}
 
	if (success<-1)								/* failure! */
		return (-1.);
	theta = theta + (price-Bprice)/Bder;		/* last improvement */
	vol = theta/sqrt(mat);
	return (vol);
}

/* A BASIC DICHOTOMIE TO FIND THE CORRECT INDEX */
int BasicDichotomie(Date *t,long Indexmin,long Indexmax,Date tRef)
{
	long imin = Indexmin, imax =Indexmax , k;

	if (tRef >= t[Indexmax]) return Indexmax;
	while ((imax - imin) > 1)
	{
		k = (imax + imin) >> 1;
		if (tRef > t[k]) imin = k;
		else if (tRef < t[k]) imax = k;
		else return k;
	}
	return imin;
}


/* COMPUTE THE DIAG SWAPTIONS EQ. STRIKES FOR EACH EXERCISE DATES FOR EMKi METHODS 
	- STORE IN THE MIDAT STRUCTURE CvgFirst, CvgPay, HDates */

Err update_autocal_midat_struct(SrtSimMidAt		*MidAt,
								char			*YcName,
								Err				(*GetVol)(Date, Date, double, SRT_Boolean, double*),
								Err				(*GetBeta)(Date, Date, double*))
{
	Err					err = NULL;
	long				i, k, nCpn,  FirstPayIndex; 
	Date				*HDates, *tStart,*tEx,tNow, *tPay = (MidAt->tPay);
	double				*CvgFirst, *CvgPay, *Ks, *StDLong, DtPay, DtStart, DtEnd, 
						Payt, SwapLevel,SwapRate,YrtoExp, CEVVol, Beta,NORMVol;
	SrtCurvePtr			YldCrv;
	String				ccy_str;
	SrtCcyParam			*ccy_param = NULL;
	SrtBasisCode		FixedBasis;
	long				spot_lag , nPay = MidAt->nPay,
						nEx = MidAt->nEx;
	SrtDiffusionType	input_vol;

	YldCrv = lookup_curve(YcName);
	/* Get CCY conventions */
	ccy_param = get_ccyparam_from_yldcrv(YldCrv);
	
	if (!ccy_param)
	{
		ccy_str = get_curve_ccy(YldCrv);
		err = swp_f_get_CcyParam_from_CcyStr(ccy_str, &ccy_param);
		if (err) return (err);
	}

	FixedBasis = ccy_param->swap_basis_code;
	spot_lag  = get_spotlag_from_curve(YldCrv);
	tNow = get_clcndate_from_curve(YldCrv);

	CvgFirst = dvector(0, nEx-1);
	CvgPay	 = dvector(0, nPay-1);
	Ks		 = dvector(0, nEx-1);
	StDLong  = dvector(0, nEx-1);
	HDates   = lngvector( 1, nEx);
	tStart = lngvector(1,nEx);
	tEx = lngvector(1,nEx);

	for(i = 1; i <= nEx; i++) 
	{
		tStart[i] = MidAt->tStart[i-1];
		tEx[i] = add_unit(tStart[i],-spot_lag,SRT_BDAY, SUCCEEDING);
	}
	for( i = 1; i < nPay; i++) CvgPay[i-1]  = coverage(tPay[i-1], tPay[i], FixedBasis);

	for( i = 1; i <= nEx ; i++) 
	{
		FirstPayIndex = MidAt->FirstPay[i-1];
		CvgFirst[i-1] = coverage(tStart[i], tPay[FirstPayIndex], FixedBasis);
		HDates[i] = tStart[i];
	}
	
	for( i = 1; i <= nEx; i++)
	{	
		FirstPayIndex = MidAt->FirstPay[i-1];			
		nCpn = (nPay-FirstPayIndex-1);
		Payt = 0.0 ;
		SwapLevel = 0.;

		for(k = 0; k <= nCpn; k++) 
		{
			DtPay = swp_f_df(tNow, tPay[FirstPayIndex + k], YcName); 
			if (DtPay==SRT_DF_ERROR)   return(err);

			if(( k == 0) && ( nCpn > 0) ) 
			{
				SwapLevel += DtPay*CvgFirst[i-1];
				Payt += DtPay*MidAt->Payment[FirstPayIndex];
			}
			else if(( k > 0 ) && ( k < nCpn) )  
			{
				SwapLevel += DtPay*CvgPay[FirstPayIndex+k-1];
				Payt += DtPay*MidAt->Payment[FirstPayIndex+k];
			}
			
			if( k == nCpn) 
			{
				if(nCpn == 0)
				{
					SwapLevel += DtPay*CvgFirst[i-1];
					Payt += DtPay*(MidAt->Payment[FirstPayIndex+k]-1);
				}
				else 
				{
					SwapLevel += DtPay*CvgPay[FirstPayIndex+k-1];
					Payt += DtPay*(MidAt->Payment[FirstPayIndex+k]-1);
				}
			}
		}

		/* COMPUTE THE STDEV */
		DtStart = swp_f_df(tNow, tStart[i], YcName); 
		if (DtStart == SRT_DF_ERROR)   return(err);

		DtEnd = swp_f_df(tNow, tPay[FirstPayIndex + nCpn], YcName); 
		if (DtEnd == SRT_DF_ERROR)   return(err);

		Ks[i-1] = Payt/(MidAt->Strike[i-1]*SwapLevel) + DtEnd*(1-MidAt->Strike[i-1])/(MidAt->Strike[i-1]*SwapLevel); 

		SwapRate = (DtStart - DtEnd)/SwapLevel; /* CASH SWAP RATE */

		err = GetVol(tStart[i], tPay[FirstPayIndex + nCpn], SwapRate, SRT_FALSE, &CEVVol);	
		if(err) return err;

		err = GetBeta(tStart[i], tPay[FirstPayIndex + nCpn], &Beta);			
		if(err) return err;

		if(Beta == 0) input_vol = SRT_NORMAL;
		else input_vol = SRT_LOGNORMAL;

		YrtoExp = (double)(tEx[i]-tNow)*YEARS_IN_DAY;

		err = srt_f_optsarbvol(SwapRate,SwapRate,YrtoExp,CEVVol,0.0,0.0,0.0,input_vol,SRT_NORMAL,&NORMVol);
		if(err) return err;
		
		StDLong[i-1] = (Ks[i-1]-SwapRate)/(NORMVol*sqrt(YrtoExp)); 

	} 

	MidAt->MidAtCvgFirst = CvgFirst;
	MidAt->MidAtCvgPay = CvgPay;
	MidAt->MidAtStrike = Ks;
	MidAt->HDates = HDates;
	MidAt->StDLong = StDLong;

	if(tStart) free_lngvector(tStart,1,nEx); tStart = NULL;
	if(tEx) free_lngvector(tEx,1,nEx); tEx = NULL;

	return (err);
}
