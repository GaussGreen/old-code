/************************************************************************
 * Module:      driface
 * File:        cpis.c
 * Function:    CPI security with an option traded on AMEX
 * Author:     Vadim Borue

$ Header:$
 ************************************************************************/
#include "cpis.h"
#include "steps.h"

#include "dritkwrp.h"		/* TDrWrapperData routines */
#include "drlio.h"

#include "convert.h"
#include "yield.h"
#include "stub.h"

#include <math.h>
#include <errno.h>

#define xCPI_DEBUG
#define xCPI_TEST

#define CPI_INTERP_TYPE GTO_LINEAR_INTERP

#define CPI_BOND_FREQ 2
#define CPI_TREE_NODE_NUMBER 50
#define CPI_MAX_STDDEV	5.0

#define SHIFT_ZERO(x) ((x) = (IS_ALMOST_ZERO(x) ? 0.000001 : (x)))

#define	READ_DATA(type,ptr,str)	\
    { if (DrlFScanVType(fp, type, (void*) ptr) != SUCCESS) \
	  { GtoErrMsg("%s: can't read %s.\n", routine, str); \
								 goto done;}}
    
static FILE *fpTERM_PRN = NULL; 

typedef struct{ 	/* for a tree */
	double		fwdCPI;	/* cpi from curve */
	double		fwdAdd; /* additive term in CPIS MTM */
	double		dfMatur; /* disc at maturity */
	long		  nEx;
	double		stdDevs;				/* Number of stdDevs to spread the tree. */
	int			   maxNodes;		/*  Maximum number of nodes on a tree. */
	int			   midNode;				/* Index of the middle node. */
	double		*fwd;  /* forwards for an option*/
	double		*drift;
	double		*bpv;  /* payoffs for fwd */
	double		*factor;  /* factor to adjust for zero vesus swap rate spread */
	double		*vol;	/* accumulative vol  including sqrt(tau) */
	double		*std;  /* local transitional std */
	double		*rho;
	double		*vMean;		
	double		*vAssetEnd;  
	double		*vOldValue;  
	double		*vNewValue;	
	double		*bmStepSize;	/*  Distance between 2 nodes on the tree for BM.	*/
	double		*bmStepSizeAlloc;
} lclTreeCPI;

static int cpiSecurityLogInputs( TCPISecurity *p);
static int lclSetTreeCPI(lclTreeCPI *pTreeCPI, TCPISecurity *pcpi);
static void lclFreeTreeCPI(lclTreeCPI *pTreeCPI);
static int lclFwdYieldBond(TDate exd, TCurve *zc,double c, long ndate, TDate *dates, double *y);
static int lclPriceTreeCPI( lclTreeCPI *p, double *pv);

/*f---------------------------------------------------------------------
 * Pricing routine for CPI linked swaps.
 * Returns SUCCESS/FAILURE.
 */

DLL_EXPORT(int) DriCPISecurity(  TCPISecurity *pcpi){
    static	char	routine[] = "DriCPISecurity";
    int	status = FAILURE;
	lclTreeCPI treeCPI, *pTreeCPI = &treeCPI;

    TDate  valueDate = pcpi->discZC->fBaseDate;
    double spotCPI = pcpi->spotCPI, firstCPI = pcpi->firstCPI;

    SHIFT_ZERO(pcpi->cpiMR);

	if(lclSetTreeCPI(pTreeCPI, pcpi) != SUCCESS) goto done;

	if(lclPriceTreeCPI(pTreeCPI, &pcpi->priceOption) != SUCCESS) goto done;

    status = SUCCESS;
  done:
	pcpi->priceOption *= pcpi->origNotl;
	pcpi->priceSwap *= pcpi->origNotl;
	pcpi->price = pcpi->priceSwap + pcpi->priceOption;

    if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	lclFreeTreeCPI(pTreeCPI);
    return(status);
}

static int lclSetTreeCPI( lclTreeCPI *pTreeCPI, TCPISecurity *pcpi){
    static	char	routine[] = "lclSetTreeCPI";
	double freq = CPI_BOND_FREQ, acc = 1.0/freq;
    int	status = FAILURE, i, nEx;
	TDate today = pcpi->today;
	double fee = pcpi->fee, spread = pcpi->finSpread;
	double bpVol = pcpi->cpiBpVol, mr = pcpi->cpiMR;

	pTreeCPI->fwd = pTreeCPI->bpv = pTreeCPI->drift = 0;
	pTreeCPI->vol = pTreeCPI->std = pTreeCPI->rho = pTreeCPI->vMean = pTreeCPI->vAssetEnd = 0;  
	pTreeCPI->vOldValue = pTreeCPI->vNewValue = pTreeCPI->bmStepSizeAlloc = 0;

	pTreeCPI->maxNodes = CPI_TREE_NODE_NUMBER;
	pTreeCPI->stdDevs = CPI_MAX_STDDEV;
	pTreeCPI->midNode = pTreeCPI->maxNodes;

	if(pcpi->maturDate <= pcpi->discZC->fBaseDate){
		double x;
		pcpi->priceOption = 0.0;
		 x = pcpi->leverage*(pcpi->spotCPI - 1.0);
		 pcpi->priceSwap = x > 0.0 ? x : 0.0;
		status = SUCCESS;
		goto done;
	}

	/* MTM underlying CPIS */
	{
		double df, ndf, rdf, x;
		if(	GtoDiscountDate( pcpi->maturDate, pcpi->discZC, CPI_INTERP_TYPE, &df) != SUCCESS ||
			GtoDiscountDate( pcpi->maturDate, pcpi->nomZC, CPI_INTERP_TYPE, &ndf) != SUCCESS ||
			GtoDiscountDate( pcpi->maturDate, pcpi->realZC, CPI_INTERP_TYPE, &rdf) != SUCCESS)
			goto done;
		pTreeCPI->dfMatur = df;
		pTreeCPI->fwdCPI = df*pcpi->leverage*pcpi->spotCPI*rdf/ndf;
		pTreeCPI->fwdAdd = -df*pcpi->leverage;
		x = pTreeCPI->fwdAdd + pTreeCPI->fwdCPI;
		pcpi->priceSwap = x > 0.0 ? x : 0.0;
	}

	pTreeCPI->nEx = nEx= pcpi->numExDates;
	if(!nEx){
		status = SUCCESS;	goto done;
	}else
		if(!pcpi->trsNumDates || !pcpi->tipsNumDates){
			GtoErrMsg("%s: failed. No tips or treasury coupons left\n", routine);
			goto done;
		}

	if(!(pTreeCPI->fwd =  (double *) MALLOC(nEx*sizeof(double)))) goto done;
	if(!(pTreeCPI->drift =  (double *) MALLOC(nEx*sizeof(double)))) goto done;
	if(!(pTreeCPI->bpv =  (double *) MALLOC(nEx*sizeof(double)))) goto done;
	if(!(pTreeCPI->vol = (double *)MALLOC(nEx*sizeof(double)))) goto done;
	if(!(pTreeCPI->std = (double *)MALLOC(nEx*sizeof(double)))) goto done;
	if(!(pTreeCPI->rho = (double *)MALLOC(nEx*sizeof(double)))) goto done;
	if(!(pTreeCPI->bmStepSizeAlloc = (double *)MALLOC((nEx+1)*sizeof(double)))) goto done;
	{
		int n = 2*pTreeCPI->maxNodes + 1;
		if(!(pTreeCPI->vMean = (double *)MALLOC(n*sizeof(double)))) goto done;
		if(!(pTreeCPI->vAssetEnd = (double *)MALLOC(n*sizeof(double)))) goto done;
		if(!(pTreeCPI->vOldValue = (double *)MALLOC(n*sizeof(double)))) goto done;
		if(!(pTreeCPI->vNewValue = (double *)MALLOC(n*sizeof(double)))) goto done;
	}
	{ /* initialize */
		int m = pTreeCPI->maxNodes, mid = pTreeCPI->midNode;
		double *v = pTreeCPI->vNewValue;
		for(i=-m; i<=m; i++) v[mid + i] = 0.0;
	}
	pTreeCPI->bmStepSize = pTreeCPI->bmStepSizeAlloc + 1;
	pTreeCPI->bmStepSize[-1] = 0.0;

    for (i=0; i< nEx; i++) {
		TDate exd = pcpi->exDates[i];
		double ny, ry, dt, ndf, rdf,df, dff, z, fdt, totVol;
		double estAdd, exCPI, tau, rho, factor;
		tau = (double)(pcpi->ntDates[i] - today)/ 365.25;
		if(tau <= 0 ) tau= 1.0e-15;

		pTreeCPI->bmStepSize[i] = pTreeCPI->stdDevs*sqrt(tau)/pTreeCPI->maxNodes;
		factor = 0.5*(1.0 - exp(-2.0*mr*tau))/mr;
		totVol = bpVol*sqrt(factor);
		pTreeCPI->vol[i] = bpVol;
		if(i){
			double dtau = (double)(pcpi->ntDates[i] - pcpi->ntDates[i-1])/ 365.25;
			rho = exp(-mr*dtau);
			factor = rho;
		}else{
			rho = 0.0;
			factor = exp(-mr*tau);
		}
		pTreeCPI->rho[i] = rho;
		pTreeCPI->std[i] = sqrt(0.5*(1.0 - factor*factor)/mr);

		if(GtoDayCountFraction(exd,  pcpi->maturDate, GTO_B30_360, &dt) != SUCCESS) goto done;
		if( lclFwdYieldBond(exd, pcpi->nomZC, pcpi->trsCoupon, pcpi->trsNumDates, pcpi->trsDates, &ny) != SUCCESS ||
		    lclFwdYieldBond(exd, pcpi->realZC, pcpi->tipsCoupon, pcpi->tipsNumDates, pcpi->tipsDates, &ry) != SUCCESS)
			goto done;
		if(	GtoDiscountDate( exd, pcpi->discZC, CPI_INTERP_TYPE, &df) != SUCCESS ||
			GtoDiscountDate( exd, pcpi->nomZC, CPI_INTERP_TYPE, &ndf) != SUCCESS ||
			GtoDiscountDate( exd, pcpi->realZC, CPI_INTERP_TYPE, &rdf) != SUCCESS)
			goto done;

		/* this factor accounts for difference between swap rate and zero rate approximate acoounting of double shrinking */
		z = 2.0*mr*dt; 
		/* this an empirical formula and price of option crucially depends on this adjustment factor */
/*		factor = dt*(0.25+ (exp(-2.0*z)*(2.0*z*z*z + 3.0*z*z) + z*z - 4.0*z*z*exp(-z))/(z*z*(1-exp(-z))))/1.25; */
		factor = dt*(z + (z+2.0)*exp(-z) - 2.0)/z;

		pTreeCPI->vol[i] *= factor;
		totVol *= factor;
		pTreeCPI->drift[i] = exp(-0.5*totVol*totVol);

		dff = df/pTreeCPI->dfMatur;
		fdt = dt*freq;
		z = freq*(pow(dff, 1.0/fdt) - 1.0);
		z = 1.0 + acc*(z + spread);
		dff = 1.0/pow(z,fdt);
		estAdd = -dff*df*pcpi->leverage;
		exCPI = df*pcpi->leverage*pcpi->spotCPI*rdf/ndf;
		z = (1.0 + acc*ny)/(1.0 + acc*(ry + fee));
		exCPI *= pow(z, fdt)*dff;
		estAdd += exCPI - pTreeCPI->fwdAdd;
		pTreeCPI->bpv[i] = estAdd;
		pTreeCPI->fwd[i] = pTreeCPI->fwdCPI/estAdd;
	}

    GTO_IF_LOGGING({
		for(i=0; i<nEx; i++){
			DrlFPrintf(fpTERM_PRN, "Exdate=%s,   forward =%12.4f,  bpv=%12.4f, vol=%12.4f\n", 
					GtoFormatDate(pcpi->exDates[i]),pTreeCPI->fwd[i], pTreeCPI->bpv[i], pTreeCPI->vol[i]);
		}
    });

   status = SUCCESS;
  done:
    if (status != SUCCESS)
	GtoErrMsg("%s: failed.\n", routine);
    return(status);
}

static void lclFreeTreeCPI(lclTreeCPI *pTreeCPI){
	if(!pTreeCPI) return;
	if(pTreeCPI->fwd) FREE(pTreeCPI->fwd);
	if(pTreeCPI->drift) FREE(pTreeCPI->drift);
	if(pTreeCPI->bpv) FREE(pTreeCPI->bpv);
	if(pTreeCPI->vol) FREE(pTreeCPI->vol);
	if(pTreeCPI->std) FREE(pTreeCPI->std);
	if(pTreeCPI->rho) FREE(pTreeCPI->rho);
	if(pTreeCPI->bmStepSizeAlloc) FREE(pTreeCPI->bmStepSizeAlloc);
	if(pTreeCPI->vMean) FREE(pTreeCPI->vMean);
	if(pTreeCPI->vAssetEnd) FREE(pTreeCPI->vAssetEnd);
	if(pTreeCPI->vOldValue) FREE(pTreeCPI->vOldValue);
	if(pTreeCPI->vNewValue) FREE(pTreeCPI->vNewValue);
}

static int lclFwdYieldBond(TDate exd, TCurve *zc,double c, long ndate, TDate *dates, double *y){
    static	char	routine[] = "lclSetTreeCPI";
    int	status = FAILURE, i;
	long freq = CPI_BOND_FREQ;
	double df, price, cdf, matur;

	if(	GtoDiscountDate(exd, zc, CPI_INTERP_TYPE, &df) != SUCCESS) goto done;
	/* get the forward bond dirty price */
	for(i=0, price = 0.0; i<ndate; i++){
		TDate tc = dates[i];
		if(exd >= tc) continue;
		if(GtoDiscountDate(tc, zc, CPI_INTERP_TYPE, &cdf) != SUCCESS) goto done;
		price += cdf;
	}
	price *= c/freq;
	price += cdf;
	price /= df;

	if(GtoDayCountFraction(exd,  dates[ndate-1], GTO_B30_360, &matur) != SUCCESS) goto done;
	if(GtoBondYieldToMaturity(c, price, freq, matur, c, GTO_STUB_NONE, y) != SUCCESS) goto done;

   status = SUCCESS;
  done:
    if (status != SUCCESS)
	GtoErrMsg("%s: failed.\n", routine);
    return(status);
}

static int lclPriceTreeCPI( lclTreeCPI *p, double *pv){
    static	char	routine[] = "lclPriceTreeCPI";
    int	status = FAILURE,  i, nT, maxNode =  p->maxNodes;
	int totalNodes = 2*maxNode+1, midNode = p->midNode; 
	double *vOldValue = p->vOldValue, *vAssetEnd = p->vAssetEnd, *vNewValue = p->vNewValue;
	double *vMean = p->vMean, *step = p->bmStepSize,  *rho = p->rho, *bpv = p->bpv, *fwd = p->fwd;
	double *vol = p->vol, *std = p->std, *drift = p->drift;

	for(nT= p->nEx-1; nT>=0; nT--){
		double intrinsic, loLimit, hiLimit, exBdry,exPremium, oldExPremium, brownian, oldBrownian, mInv;
		double cF1,cF2;
		int nLowerSlice, exState, oldExState;

		for(i=0; i<totalNodes; i++){
			double B_T = step[nT]*(i- midNode);
			vOldValue[i] =  bpv[nT]*(1.0 - fwd[nT]*drift[nT]*exp(vol[nT]*B_T));
			vAssetEnd[i] = B_T;
		}
		loLimit = vAssetEnd[0];
		hiLimit = vAssetEnd[2*maxNode];

	/* Now compute the exercise bdry. */
		nLowerSlice = totalNodes;
		exBdry = hiLimit;
		exState= oldExState = 0;				/* sign of exPremium, or 0 for undecided. */
		exPremium = oldExPremium = 0.0;

		for(i=0; i<totalNodes; i++){
			intrinsic = vOldValue[i];
			exPremium = intrinsic - vNewValue[i];
			if(exPremium > 0.0){	/* Exercise the option. */
				vOldValue[i] =	intrinsic;
				exState = 1;
			}else{	/* Don't excercise */
				vOldValue[i] = vNewValue[i];
				exState = -1;
			}

		   /* Exercise boundary is where we change from exercised to non-exercised state (or vice versa).
		    We look for the most central exercise boundary - there is often more than one.
			*/
			if(exState*oldExState == -1	&&		/* Changed exercise state */
			   abs(i - midNode) < abs(nLowerSlice - midNode))	/* and we are more central than the last change. */
			{
				nLowerSlice = i; /* New exercise boundary. */
			/* Now estimate the value of the Brownian motion at the point where we exercise.
			This isn't necesarily on a node.
			*/
				if(exPremium * oldExPremium > 0.0){
						GtoErrMsg("%s: failed. Didn't change state at exercise boundary\n", routine);
						goto done;
				}
				brownian    = step[nT]*(i- midNode); 
				oldBrownian = step[nT]*(i - 1 - midNode);
				mInv = (oldBrownian - brownian)/(oldExPremium - exPremium);
				exBdry = brownian - exPremium*mInv;
			}

			oldExState = exState;
			oldExPremium = exPremium;
		}

		/* Integrate the value at the end over the distribution.
		 We have B(nT) = cF1*B(nT-1) + cF2*Z	is normally distributed.
		 */
		cF1 = rho[nT];
		cF2 = std[nT];

		for(i = 0; i < totalNodes; i++)
			vMean[i] = cF1 *step[nT-1]*(i- midNode); 
	
		if(DrTreesSplit1DSingleStep(
			vAssetEnd,				/* (I) underlying asset value at the end of the step */
			vOldValue,				/* (I) option payoff on underlying asset values at the end of the step */
			totalNodes,				
			nLowerSlice,		
			totalNodes,			
			vMean,					/* (I) underlying asset values at the start of the step */
			totalNodes,			/* (I) number of nodes at start of step */			
			cF2,					/* (I) standard deviation across step - NOT same as ChaseTools volatility! */
			hiLimit,				/* (I) if specified, value of high integration limit to use */	
			loLimit,				/* (I) if specified, value of low integration limit to use */
			exBdry,				
			hiLimit,			
			vNewValue) != SUCCESS){
					GtoErrMsg("%s: failed. ChaseTreesSplit1DSingleStep failed\n", routine);
					goto done;
			}
	}

	if(p->nEx)
		*pv = vNewValue[midNode];
	else
		*pv = 0.0;

   status = SUCCESS;
  done:
    if (status != SUCCESS)	GtoErrMsg("%s: failed.\n", routine);
    return(status);
}


/*f---------------------------------------------------------------------
 * Dr-Wrapper for {\tt DriCPISecurity}.
 * The argument {\tt dataFnam} specifies the name of the file
 * containing the data (if it is NULL, the default name
 * "cpiswap_w.dat" is used.
 * Returns SUCCESS/FAILURE.
 */

 DLL_EXPORT(int)  DriCpiSecurityW(char *dataFnam){
    static	char	routine[] = "DriCPISecurityW";
    int	status = FAILURE;
	char	priceType;
	TCPISecurity cpis, *pcpi = &cpis;

    FILE		*fp = NULL;
    static	char	defDataFnam[] = "cpis_w.dat";
    static	char	todayFnam[] = "today.dat";
    static	char	eqFnam[] = "equity.dyn";

    TDrWrapperData	*drWrap = NULL;
    double		pv;

	cpis.exDates = cpis.ntDates =  0;
	cpis.nomZC = cpis.realZC = cpis.discZC = 0;
	cpis.trsDates = cpis.tipsDates = 0;
	cpis.price = cpis.priceOption = cpis.priceSwap = -1.0e20;

	    /* Read today
     */
    if ((fp = fopen(todayFnam, "r")) == NULL) {
		GtoErrMsg("%s: can't open `%s' (%s).\n",
			  routine, todayFnam, strerror(errno));
		goto done;
    }
    READ_DATA(DRL_TDATE_T,  &(pcpi->today),"today");
    if (fp) fclose(fp);

    /* Read deal data   */
    if (dataFnam == NULL) dataFnam = defDataFnam;
    if ((fp = fopen(dataFnam, "r")) == NULL) {
		GtoErrMsg("%s: can't open `%s' (%s).\n",
			  routine, dataFnam, strerror(errno));
		goto done;
    }

    READ_DATA(DRL_TDATE_T, &(pcpi->maturDate),	"secuity maturity date");
	if(pcpi->maturDate < pcpi->today){
		GtoErrMsg("%s: the CPIS deal has expired.\n", routine);
		goto done;
    }
    READ_DATA(DRL_DOUBLE_T, &(pcpi->origNotl), "deal notional");
    READ_DATA(DRL_DOUBLE_T, &(pcpi->firstCPI), "deal initial CPI");
    READ_DATA(DRL_DOUBLE_T, &(pcpi->leverage), "deal CPI leverage");
    READ_DATA(DRL_PERCENT_T, &(pcpi->fee), "deal exercise fee");
    READ_DATA(DRL_PERCENT_T, &(pcpi->finSpread), "financing spread");

    READ_DATA(DRL_LONG_T, &(pcpi->numExDates), "num exercise dates");
    if(pcpi->numExDates < 1) {
		GtoErrMsg("%s: need at least 1 exercise date.\n", routine);
        goto done;
    }
	/* exercise notify dates */
    if(DrlLilVectArrayFpReadV(fp, 
			      pcpi->numExDates,
			      DRL_TDATE_T,  (void*) &pcpi->exDates,
			      DRL_TDATE_T,  (void*) &pcpi->ntDates,
			      DRL_NULL_T) == FAILURE){  
        GtoErrMsg("%s: Cannot read dates array.\n", routine);
        goto done;
    }
	{ /* cut dates */
		int i = 0, j =0;
		TDate tt= cpis.today;
		for(; i < pcpi->numExDates; i++){
			if(pcpi->exDates[i] < tt) continue;
			pcpi->exDates[j] = pcpi->exDates[i];
			pcpi->ntDates[j] = pcpi->ntDates[i];
			j++;
		}
		pcpi->numExDates = j;
	}

	/* set pegged bonds */
    READ_DATA(DRL_PERCENT_T, &(pcpi->trsCoupon), "treasury coupon");
    READ_DATA(DRL_LONG_T, &(pcpi->trsNumDates), "num exercise dates");
    if(pcpi->trsNumDates < 1) {
		GtoErrMsg("%s: need at least 1 treasury bond date.\n", routine);
        goto done;
    }
	if(DrlLilVectArrayFpReadV(fp, cpis.trsNumDates,
			      DRL_TDATE_T,  (void*) &pcpi->trsDates,
			      DRL_NULL_T) == FAILURE){  
			GtoErrMsg("%s: Cannot read dates array.\n", routine);
			goto done;
	}
	{ /* cut dates */
		int i = 0, j =0;
		TDate tt= cpis.today;
		for(; i < pcpi->trsNumDates; i++){
			if(pcpi->trsDates[i] <= tt) continue;
			pcpi->trsDates[j] = pcpi->trsDates[i];
			j++;
		}
		pcpi->trsNumDates = j;
	}
    if(pcpi->trsNumDates < 1) {
		GtoErrMsg("%s: Treasury bond has expiried.\n", routine);
        goto done;
    }
	
    READ_DATA(DRL_PERCENT_T, &(pcpi->tipsCoupon), "tips coupon");
    READ_DATA(DRL_LONG_T, &(pcpi->tipsNumDates), "num tips dates");
    if(pcpi->tipsNumDates < 1) {
		GtoErrMsg("%s: need at least 1 treasury bond date.\n", routine);
        goto done;
    }
	if(DrlLilVectArrayFpReadV(fp, cpis.tipsNumDates,
			      DRL_TDATE_T,  (void*) &pcpi->tipsDates,
			      DRL_NULL_T) == FAILURE){  
			GtoErrMsg("%s: Cannot read dates array.\n", routine);
			goto done;
	}
	{ /* cut dates */
		int i = 0, j =0;
		TDate tt= cpis.today;
		for(; i < pcpi->tipsNumDates; i++){
			if(pcpi->tipsDates[i] <= tt) continue;
			pcpi->tipsDates[j] = pcpi->tipsDates[i];
			j++;
		}
		pcpi->tipsNumDates = j;
	}
    if(pcpi->tipsNumDates < 1) {
		GtoErrMsg("%s: Tips bond have expired.\n", routine);
        goto done;
    }

    READ_DATA(DRL_CHAR_T, &priceType, "price type");

    /* Read model data 
     */   
    READ_DATA(DRL_PERCENT_T, &(pcpi->cpiBpVol),  "cpi vols");
    READ_DATA(DRL_PERCENT_T, &(pcpi->cpiMR),  "cpi MR");
    
    /* Close data file */ 
    if(fp) fclose(fp);

    /* Read market data
     */
    if (DriTDrWrapperDataGetFull(NULL, 
				 DRI_DRW_TYPE2_3CURVES, 
				 &drWrap) != SUCCESS)
	goto done;

	pcpi->nomZC =  drWrap->fZcCurve;
	pcpi->realZC =   drWrap->fRiskZcCurve;
	pcpi->discZC =	 drWrap->fDiscZcCurve;

    /* Read CPI spot value
     */
    if ((fp = fopen(eqFnam, "r")) == NULL) {
		GtoErrMsg("%s: can't open `%s' (%s).\n",
			  routine, eqFnam, strerror(errno));
		goto done;
    }
	{ 
		TDate dummyDate;
		READ_DATA(DRL_TDATE_T,  &dummyDate,	        "base date");
		READ_DATA(DRL_DOUBLE_T, &(pcpi->spotCPI),	"spot value");
		pcpi->spotCPI /= pcpi->firstCPI;
	}
    if (fp) fclose(fp);


    /* Open TERM.prn
     */
    fpTERM_PRN = fopen("TERM.PRN", "w");
    if (fpTERM_PRN IS NULL)
    {
        GtoErrMsg("%s: Cannot open TERM_PRN.\n",routine);
        goto done;
    }

#ifdef CPI_DEBUG
    GtoLoggingSet(1); 
#endif
    /* Log inputs in fpTERM_PRN     */
    GTO_IF_LOGGING({cpiSecurityLogInputs(pcpi);});

#ifdef CPI_TEST
	{
		int i;
		for(i=0; i< 20; i++){
/*
			double slope, lambda = 0.1;
			TCurve *zc = pcpi->discZC,  *nzc = pcpi->nomZC, *rzc = pcpi->realZC;
			int j;
			slope = 0.003*(i-10);
			for(j=0; j < zc->fNumItems; j++){
				double x = (double)(zc->fArray[j].fDate - zc->fBaseDate)/365.0;
				rzc->fArray[j].fRate = 0.04;
				zc->fArray[j].fRate =   0.05 + slope*(1.0-exp(-x*lambda));
				nzc->fArray[j].fRate =  0.04 + slope*(1.0-exp(-x*lambda));
			}
*/
/*			pcpi->cpiMR = 0.02 + 0.02*i; */
			pcpi->cpiBpVol = 0.002 + 0.001*i;
#endif
	    /* Call pricing routine   */
			if (DriCPISecurity( pcpi) == FAILURE)    goto done;
#ifdef CPI_TEST
			GtoLoggingSet(1); 
			   GTO_IF_LOGGING({
			 		 DrlFPrintf(fpTERM_PRN, "%g   %g  %g\n", pcpi->cpiBpVol, pcpi->priceOption, pcpi->priceSwap);
				});
			GtoLoggingSet(0); 
		}
	}
#endif

    GTO_IF_LOGGING({
		DrlFPrintf(fpTERM_PRN, "price=%12.4f,   underlying=%12.4f,  option=%12.4f\n", 
			pcpi->price, pcpi->priceSwap, pcpi->priceOption);
    });


	switch (priceType){
	case 'U':
		pv = cpis.priceSwap; 	break;
	case 'O':
		pv = cpis.priceOption; break;
	case 'S':
		pv = cpis.price; break;
	default:
        GtoErrMsg("%s: Unknown price type.\n",routine);
        goto done;
	}

    if (DriTDrWrapperDataPutPrice(pv) != SUCCESS)
		goto done;
    printf("Price:     %f \n", pv);
   
	status = SUCCESS;
  done:
    if (fp) fclose(fp);
    if(fpTERM_PRN != NULL) fclose(fpTERM_PRN);

    DriTDrWrapperDataFree(drWrap);
    
    if(cpis.exDates) FREE(cpis.exDates);
    if(cpis.ntDates) FREE(cpis.ntDates);
    if(cpis.trsDates) FREE(cpis.trsDates);
    if(cpis.tipsDates) FREE(cpis.tipsDates);

    if (status != SUCCESS)
	GtoErrMsg("%s: failed.\n", routine);
 
   return(status);
}

static int cpiSecurityLogInputs(TCPISecurity *pcpi){
    static	char	routine[] = "cpiSecurityLogInputs";
    int status = FAILURE;
    int timeStamp = GtoErrMsgTimeStamp (0); /* disable temporarily */
    FILE *prevFile= GtoErrMsgFilePointer(fpTERM_PRN);  /* to switch back 
							  to error.log after 
							  logging is done */
    long 		i;


    GtoErrMsg("\n%s INPUTS:\n", routine);
    GtoErrMsg("\n");
    GtoErrMsg("Today:                %s\n\n", GtoFormatDate(pcpi->today));
    GtoErrMsg("CPIS maturity:     %s\n\n", GtoFormatDate(pcpi->maturDate));
    GtoErrMsg("Original notional:   %10.4f\n\n", pcpi->origNotl);
    GtoErrMsg("Initial CPI:            %10.4f\n\n", pcpi->firstCPI);
    GtoErrMsg("Deal leverage:       %10.4f\n\n", pcpi->leverage);
    GtoErrMsg("Redemption fee:       %10.4f\n\n", pcpi->fee);
    GtoErrMsg("Credit spread:       %10.4f\n\n", pcpi->finSpread);

    GtoErrMsg("\n Notification and exercise dates\n\n");
    for (i=0; i<pcpi->numExDates; i++) {
	    GtoErrMsg("%3ld:   %s   %s \n", i+1, 
		      GtoFormatDate(pcpi->exDates[i]),
		      GtoFormatDate(pcpi->ntDates[i]));
    }		     

    GtoErrMsg("Treasury coupon:       %10.4f\n\n", pcpi->trsCoupon);
    GtoErrMsg("\n Treasury coupon dates\n\n");
    for (i=0; i<pcpi->trsNumDates; i++) {
	    GtoErrMsg("%3ld:   %s \n", i+1, 
		      GtoFormatDate(pcpi->trsDates[i]));
    }		     

    GtoErrMsg("Tips coupon:       %10.4f\n\n", pcpi->tipsCoupon);
    GtoErrMsg("\n Tips coupon dates\n\n");
    for (i=0; i<pcpi->tipsNumDates; i++) {
	    GtoErrMsg("%3ld:   %s \n", i+1, 
		      GtoFormatDate(pcpi->tipsDates[i]));
    }

    GtoErrMsg("\n\n");
    GtoErrMsg("CPI bp vol:       %10.4f\n\n", pcpi->cpiBpVol);
    GtoErrMsg("CPI MR:           %10.4f\n\n", pcpi->cpiMR);

    GtoErrMsg("Spot CPI fixing:      %10.4f\n\n", pcpi->spotCPI*pcpi->firstCPI);

    GtoPrintTCurve(pcpi->nomZC,  "Nominal Curve");
    GtoPrintTCurve(pcpi->realZC, "Real Curve");
    GtoPrintTCurve(pcpi->discZC, "Discount Curve");

    GtoErrMsg("\n\n\n");

    GtoErrMsgFilePointer(prevFile);   /* return to error.log */
    GtoErrMsgTimeStamp (timeStamp);

    status = SUCCESS;

    return(status);
}
