/****************************************************************
 * Module:	VNFM
 * Submodule:	Add-in Wrapper functions
 * Function:	
 * Author:	David Liu
 *****************************************************************/
#include "drlstd.h"			/* platform compatibility */
#include <math.h>

#include "date_sup.h"
#include "convert.h"
#include "yearfrac.h"
#include "gtomat.h"

#include "drlmem.h"		/* DrlDoubleMatrAlloc() */
#include "drlvtype.h"		/* DrlLilStructGet() */
#include "drlts.h"		/* DrlTCurveWrapRead() */
#include "drlmatrixo.h"

#include "macros.h"	
#include "bastypes.h"

#define	_vnfm_SOURCE
#include "vnfmanly.h"
#include "vnfmwrap.h"
#include "vnfmwrpobj.h"
#include "vnfmdata.h"
#include "volcurvo.h"

#if defined(_WINDLL) || !defined(TESTLIB)
# undef __DEBUG__
#endif


/*f-----------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_GENFWDVOL_O.
 *                                                             
 * <br><br>
 * Add-in wrapper for function <i> VnfmGenerateFwdVolL</i>.  \\
 *
 * The inputs are the same as <i> VnfmGenerateFwdVolL</i>. \\
 * The output <i> fwdVolO</i> returns an object handle of <b> TVolCurve</b>
 * class which contains following fields:                                      
 * <br>
 * 	<br> <i> baseDate</i>: <i> refDateL</i>.
 * 	<br> <i> dates</i>: <i> obsEndDatesL</i>.
 * 	<br> <i> volatilities</i>: <i> volsL</i>.
 * <br>
 * Here we use obsEndDatesL instead of obsStartDatesL to 
 * avoid the case where all the start dates are the 
 * same (today), which is not allowed in GtoVolCurveNew 
 * function.

 */
DLL_EXPORT(int)
VnfmGenerateFwdVolO(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDateL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRateL,	/*  3 'F' (I) zero coupon rates */

	FloatL *backboneqL,	/*  4 'L' (I) dist type (0=LN, 0.5=N) */
	FloatL *betaL,		/*  5 'F' (I) array of mr coeff */
	FloatL *alphaL,		/*  6 'F' (I) array of weight coeff */
	TDateL *dateL,		/*  7 'D' (I) array of dates */
	FloatL *sigmaL,		/*  8 'F' (I) array of spot volatilities */
	FloatL *rhoL,		/*  9 'F' (I) array of correlation */

	TDateL *obsStartDatesL, /* 10 'D' (I) observation start dates */
        TDateL *obsEndDatesL,   /* 11 'D' (I) observation end dates */
        TDateL *resetDatesL,    /* 12 'D' (I) rate reset dates */
        double *maturitiesL,    /* 13 'F' (I) rate maturities */
        long   *frequenciesL,   /* 14 'L' (I) rate frequencies */
        TVolCurve **fwdVolO)    /* 15 'F' (O) Vols output */
{

    	static	char routine[] = "VnfmGenerateFwdVolO";
	int status = FAILURE;
	FloatL *volsL = NULL;

	volsL = NEW_ARRAY(double, 
			  (long)obsEndDatesL[0]+1); 
	if (volsL IS NULL) 
		goto done;
	
        volsL[0] = (long)obsEndDatesL[0];

	if(VnfmGenerateFwdVolL( refDateL,
				zcDateL,
				zcRateL,	
				backboneqL,
				betaL,	
				alphaL,
				dateL,
				sigmaL,
				rhoL,
				obsStartDatesL, 
				obsEndDatesL,
				resetDatesL,
				maturitiesL,
				frequenciesL,
				volsL) == FAILURE)
	goto done;

	/* Construct VolCurve 
	 * Use obsEndDatesL instead of obsStartDatesL to 
	 * avoid the case where all the start dates are the 
	 * same (today), which is not allowed in GtoVolCurveNew 
	 * function.
	 */
        *fwdVolO = GtoVolCurveNew( (long)obsEndDatesL[0],
                                   refDateL[1],
                                   obsEndDatesL+1,
                                   volsL+1,
                                   "V");

        if (*fwdVolO IS NULL)
        {
            GtoErrMsg ("%s:Failed to construct fwdVolO object.\n",routine);
            goto done;
        }

        status = SUCCESS;

   done:
        if (status IS FAILURE)
            GtoErrMsg ("%s: Failed.\n",routine);

        FREE(volsL);

        return status;
}



/*f-----------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_GENBVOL_O.
 *                                                             
 * <br><br>
 * Add-in wrapper for function <i> VnfmGenerateVolCurveL</i>.  \\
 *
 * The inputs are the same as <i> VnfmGenerateVolCurveL</i>. \\
 * The output <i> bVolO</i> returns an object handle of <b> TVolCurve</b>
 * class which contains following fields:                                      
 * <br>
 *   <br> <i> baseDate</i>: <i> bvRefDateL</i>.
 *   <br> <i> dates</i>: <i> bvDatesL</i>.
 *   <br> <i> volatilities</i>: <i> bvRatesL</i>. 
 * <br>
 */
DLL_EXPORT(int)
VnfmGenerateVolCurveO(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDateL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRateL,	/*  3 'F' (I) zero coupon rates */

	FloatL *backboneqL,	/*  4 'L' (I) dist type (0=LN, 0.5=N) */
	FloatL *betaL,		/*  5 'F' (I) array of mr coeff */
	FloatL *alphaL,		/*  6 'F' (I) array of weight coeff */
	TDateL *dateL,		/*  7 'D' (I) array of dates */
	FloatL *sigmaL,		/*  8 'F' (I) array of spot volatilities */
	FloatL *rhoL,		/*  9 'F' (I) array of correlation */

	TDateIntervalL *bvMatL,	/* 10 'F' (I) vol maturity */
	IntL   *volFreqL,	/* 11 'L' (I) vol frequency */
	IntL   *volTypeL,	/* 12 'L' (I) 0=base, 1=fwd */
	TDateL *bvRefDateL,	/* 13 'D' (I) vol ref date */
	TDateL *bvDatesL,	/* 14 'D' (I) array of vol dates */
	TVolCurve **bVolO)	/* 15 'F' (O) base vol */
{

    	static	char routine[] = "VnfmGenerateBaseVolO";
	int status = FAILURE;
	FloatL *bvRatesL = NULL;

	bvRatesL = NEW_ARRAY(double,
			     (long)bvDatesL[0]+1); 
	if (bvRatesL IS NULL) 
		goto done;

        bvRatesL[0] = (long)bvDatesL[0];

	if ( VnfmGenerateVolCurveL( refDateL,
				    zcDateL,
				    zcRateL,
				    backboneqL,
				    betaL,	
				    alphaL,
				    dateL,
				    sigmaL,	
		    		    rhoL,
				    bvMatL,
		 		    volFreqL,	
				    volTypeL,
				    bvRefDateL,
				    bvDatesL,
				    bvRatesL) == FAILURE)
	goto done;

	/* Construct VolCurve */
        *bVolO = GtoVolCurveNew((int)bvDatesL[0],
                                bvRefDateL[1],
                                bvDatesL+1,
                                bvRatesL+1,
                                "Linear");

        if (*bVolO IS NULL)
        {
            GtoErrMsg ("%s:Failed to construct bVolO object.\n",routine);
            goto done;
        }

        status = SUCCESS;

   done:
        if (status IS FAILURE)
            GtoErrMsg ("%s:Failed to create a VolCurve.\n",routine);

        FREE(bvRatesL);

        return status;
}


/*f-----------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_GENSWMAT_O.
 *                                                             
 * <br><br>
 * Add-in wrapper for function <i> VnfmGenerateSwaptionMatrixL</i>.  \\
 *
 * The inputs are the same as <i> VnfmGenerateSwaptionMatrixL</i>. \\
 * The output <i> swVolO</i> converts <i> swVolL</i> to a [<i> swExpL</i> $x$ 
 * <i> swMatL</i>] matrix object of <b> MAT</b> class, which can be accessed using 
 * the following fields:
 * <br>
 *    <br> Use <i> row1, row2,</i> \ldots, to access each 
 *  			 row of the swaption volatility matrix <i> swVolO</i>.
 *    <br> Use <i> col1, col2,</i> \ldots, to access each 
 *			 column of the swaption volatility matrix <i> swVolO</i>.
 * <br>
 */
DLL_EXPORT(int)
VnfmGenerateSwaptionMatrixO(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDateL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRateL,	/*  3 'F' (I) zero coupon rates */

	FloatL *backboneqL,	/*  4 'L' (I) dist type (0=LN, 0.5=N) */
	FloatL *betaL,		/*  5 'F' (I) array of mr coeff */
	FloatL *alphaL,		/*  6 'F' (I) array of weight coeff */
	TDateL *dateL,		/*  7 'D' (I) array of dates */
	FloatL *sigmaL,		/*  8 'F' (I) array of spot volatilities */
	FloatL *rhoL,		/*  9 'F' (I) array of correlation */

	TDateL *swRefDateL,	/* 10 'D' (I) vol ref date */
	IntL   *swTypeL,	/* 11 'L' (I) array of matrix param [0..2]: */
				/*        [0] matr type (0=vertical, 1=diag) */
				/*        [1] vol frequency (1,2,4,12) */
				/*        [2] NOT USED (pass 0) */
	TDateIntervalL *swMatL,	/* 12 'F' (I) array of mat intervals */
	TDateIntervalL *swExpL,	/* 13 'F' (I) array of exp intervals */
	TMatrix2D      **swVolO)/* 14 'F' (O) matrix of swaption vols */
{
	static	char		routine[] = "VnfmGenerateSwaptionMatrixO";
	int status = FAILURE;
	FloatL *swVolL = NULL;

	swVolL = NEW_ARRAY(double,
			   (long)(swMatL[0]*swExpL[0])+1);
	if (swVolL IS NULL) {
		GtoErrMsg ("malloc of swVolL failed\n");
		goto done;
	}

	swVolL[0] = (long)(swMatL[0]*swExpL[0]);

	if ( VnfmGenerateSwaptionMatrixL( refDateL,
					  zcDateL,
					  zcRateL,
					  backboneqL,
					  betaL,	
					  alphaL,
					  dateL,
					  sigmaL,		
					  rhoL,	
					  swRefDateL,
					  swTypeL,
					  swMatL,		
					  swExpL,	
					  swVolL) == FAILURE)
        goto done;;

	*swVolO = GtoMatrix1DTo2DNew ((long)swExpL[0], 
			                      (long)swMatL[0], 
			                      swVolL +1 );

	status = SUCCESS;

   done:
        if (status IS FAILURE)
            GtoErrMsg ("%s:Failed to create a TMatrix2D.\n",routine);

        FREE(swVolL);

        return status;
}



/*f-----------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_GENACORR_O.
 *                                                             
 * <br><br>
 * Add-in wrapper for function <i> VnfmGenerateCorrelationMatrixL</i>.  \\
 *
 * The inputs are the same as <i> VnfmGenerateCorrelationMatrixL</i>. \\
 * The output <i> corrO</i> converts <i> corrL</i> to a [<i> tMatL</i> $x$ 
 * <i> tMatL</i>] matrix object of <b> MAT</b> class, which can be accessed using 
 * the following fields:                                      
 * <br>
 * 	<br> Use <i> row1, row2,</i> \ldots, to access each 
 *			row of the correlation matrix <i> corrO</i>. 
 * 	<br> Use <i> col1, col2,</i> \ldots, to access each 
 *			column of the correlation matrix <i> corrO</i>.
 * <br>
 */
DLL_EXPORT(int)
VnfmGenerateCorrelationMatrixO(
        TDateL *refDateL,       /*  1 'D' (I) reference date */
        TDateL *zcDateL,        /*  2 'D' (I) zero coupondates */
        FloatL *zcRateL,        /*  3 'F' (I) zero coupon rates */

        FloatL *backboneqL,     /*  4 'F' (I) [1] backbone q (0=LN, 0.5=N) */
        FloatL *betaL,          /*  5 'F' (I) array of mr coeff */
        FloatL *alphaL,         /*  6 'F' (I) array of weight coeff */
        TDateL *dateL,          /*  7 'D' (I) array of dates */
        FloatL *sigmaL,         /*  8 'F' (I) array of spot volatilities */
        FloatL *rhoL,           /*  9 'F' (I) array of correlation */

        TDateIntervalL *tExpL,  /* 10 'F' (I) expiration time */
        TDateIntervalL *tMatL,  /* 11 'F' (I) array of mat times */
        TMatrix2D      **corrO) /* 12 'F' (O) matrix of correlations */
{
	static  char            routine[] = "VnfmGenerateCorrelationMatrixO";
        int             status = FAILURE;
	FloatL 		*corrL = NULL;

	corrL = NEW_ARRAY(double,
			  (long)(tMatL[0]*tMatL[0])+1);
	if (corrL IS NULL) {
		GtoErrMsg ("malloc of corrL failed\n");
		goto done;
	}

	corrL[0] = (long)(tMatL[0]*tMatL[0]);

	status = VnfmGenerateCorrelationMatrixL( refDateL,    
        					 zcDateL,    
        					 zcRateL,   
        					 backboneqL,  
        					 betaL,      
        					 alphaL,    
        					 dateL,    
        					 sigmaL,  
        					 rhoL,   
						 tExpL,	     
						 tMatL,	
						 corrL);

	if (status IS FAILURE)
		goto done;

	*corrO = GtoMatrix1DTo2DNew ((long)tMatL[0], 
			                     (long)tMatL[0], 
			                     corrL + 1 );

   done:
        if (status IS FAILURE)
            GtoErrMsg ("%s:Failed to create a TMatrix2D.\n",routine);

        FREE(corrL);

        return status;
}



/*f-----------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_CAL2F1SPVOL_GEN_O.
 *                                                             
 * <br><br>
 * Add-in wrapper for function <i> VnfmCalib1V2FGeneralL</i>.
 *
 * The inputs are the same as <i> VnfmCalib1V2FGeneralL</i>. \\
 * The outputs <i> spotVolO</i>, <i> bVolO</i>, and <i> swVolO</i> 
 * return three objects of <b> TVolCurve</b>, <b> TVolCurve</b>, 
 * and <b> MAT</b> classes respectively, which contain the following 
 * fields:
 * \begin{list}%
 * {}{\setlength{\leftmargin}{\leftmargin}}
 *   <br> <i> spotVolO:</i>
 *      <br>
 *  	   <br> <i> baseDate</i>: <i> refDateL</i>.
 * 	   <br> <i> dates</i>: <i> nfDatesL</i>.
 * 	   <br> <i> volatilities</i>: <i> spotVolsL</i>.
 *      <br>
 *   <br> <i> bVolO:</i>
 *      <br>
 *         <br> <i> baseDate</i>: <i> refDateL</i>. 
 *         <br> <i> dates</i>: <i> bvDatesL</i>.
 *         <br> <i> volatilities</i>: <i> bvRatesL</i>.
 *      <br>
 *   <br> <i> swVolO:</i> convert <i> swVolL</i> into a [<i> swExpL</i> $x$ 
 *	   <i> swMatL</i>] matrix.
 *	<br>
 *	   <br> Use <i> row1, row2,</i> \ldots, to access each row 
 *		 of the swaption volatility matrix <i> swVolO</i>.
 *    	   <br> Use <i> col1, col2,</i> \ldots, to access each column
 *		 of the swaption volatility matrix <i> swVolO</i>.
 *      <br>
 * \end{list}
 */
DLL_EXPORT(int)
VnfmCalib1V2FGeneralO(
	TDateL         *refDateL,   /*  1 'D' (I) reference date */
	TDateL         *zcDateL,    /*  2 'D' (I) zero coupondates */
	FloatL         *zcRateL,    /*  3 'F' (I) zero coupon rates */

	FloatL         *backboneqL, /*  4 'L' (I) num scalars:
				     *        [1] distType
				     *        [2] generate base vol 
				     *        [3] generate swaption matrix 
				     *        [4] minimum volatility rate mat */
	FloatL         *nfParamsL,  /*  5 'F' (I) N-fact params */
	TDateL         *nfDatesL,   /*  6 'D' (I) array of dates */
	FloatL         *rateMatL,   /*  7 'F' (I) rate mat [0..nDates-1] */
	IntL           *rateFreqL,  /*  8 'L' (I) rate freq [0..nDates-1] */
	FloatL         *rateVolL,   /*  9 'F' (I) rate vol [0..nDates-1] */
	TDateIntervalL *bvMatL,	    /* 10 'F' (I) vol maturity */
	TDateL         *bvDatesL,   /* 11 'D' (I) vol dates */
	IntL           *swTypeL,    /* 12 'L' (I) [0]=type, [1]=freq */
	TDateIntervalL *swMatL,	    /* 13 'F' (I) array of mat intervals */
	TDateIntervalL *swExpL,	    /* 14 'F' (I) array of exp intervals */
	TVolCurve      **spotVolO,  /* 15 'F' (O) spot vol */
	TVolCurve      **bVolO,	    /* 16 'F' (O) base vol */
	TMatrix2D      **swVolO)    /* 17 'F' (O) swaption vol matrix */
{
    	static	char routine[] = "VnfmCalib1V2FGeneralO";
	int		status = FAILURE;
	FloatL *spotVolsL= NULL;
	FloatL *bvRatesL= NULL;
	FloatL *swVolL= NULL;

	spotVolsL = NEW_ARRAY(double,
			     (long)nfDatesL[0]+1);
	if (spotVolsL IS NULL) 
                goto done;

	bvRatesL  = NEW_ARRAY(double,
                             (long)bvDatesL[0]+1);
        if (bvRatesL IS NULL) 
                goto done;

	swVolL    = NEW_ARRAY(double,
                             (long)(swMatL[0]*swExpL[0])+1);
        if (swVolL IS NULL) 
                goto done;

        spotVolsL[0] = (long)nfDatesL[0];
        bvRatesL[0]  = (long)bvDatesL[0];
        swVolL[0]    = (long)(swMatL[0]*swExpL[0]);

	if ( VnfmCalib1V2FGeneralL( refDateL,
				    zcDateL,
				    zcRateL,
				    backboneqL,
				    nfParamsL,
				    nfDatesL,
				    rateMatL,
				    rateFreqL,
				    rateVolL,
				    bvMatL,	
				    bvDatesL,
				    swTypeL,
				    swMatL,
				    swExpL,
				    spotVolsL,
				    bvRatesL,
				    swVolL)  == FAILURE)
        goto done;

        /* Construct spot Vol */
        *spotVolO = GtoVolCurveNew((int)nfDatesL[0],
                                   refDateL[1],
                                   nfDatesL+1,
                                   spotVolsL+1,
                                   "V");

        if (*spotVolO IS NULL)
        {
            GtoErrMsg ("%s:Failed to construct spotVolO object.\n",routine);
            goto done;
        }

        /* Construct base Vol */
        *bVolO = GtoVolCurveNew((int)bvDatesL[0],
                                refDateL[1],
                                bvDatesL+1,
                                bvRatesL+1,
                                "Linear");

        if (*bVolO IS NULL)
        {
            GtoErrMsg ("%s:Failed to construct bVolO object.\n",routine);
            goto done;
        }

	*swVolO =  GtoMatrix1DTo2DNew ((long)swExpL[0], 
			                       (long)swMatL[0], 
			                       swVolL + 1);

	if (*swVolO IS NULL)
        {
            GtoErrMsg ("%s:Failed to construct swVol object.\n",routine);
            goto done;
        }

	status = SUCCESS;

    done:
	if (status IS FAILURE)
            GtoErrMsg ("%s: Failed.\n",routine);

        FREE(spotVolsL);
        FREE(bvRatesL);
        FREE(swVolL);

	return (status);
}



/*f-----------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_CAL2F1SPVOL_GEN_NEW_O.
 *                                                             
 * <br><br>
 * Add-in wrapper for function <i> VnfmCalib1V2FGeneralNewL</i>.
 *
 * The inputs are the same as <i> VnfmCalib1V2FGeneralNewL</i>. \\
 * The outputs <i> spotVolO</i>, <i> bVolO</i>, and <i> swVolO</i> 
 * return three objects of <b> TVolCurve</b>, <b> TVolCurve</b>, 
 * and <b> MAT</b> classes respectively, which contain the following 
 * fields:
 * \begin{list}%
 * {}{\setlength{\leftmargin}{\leftmargin}}
 *   <br> <i> spotVolO:</i>
 *      <br>
 *  	   <br> <i> baseDate</i>: <i> refDateL</i>.
 * 	   <br> <i> dates</i>: <i> nfDatesL</i>.
 * 	   <br> <i> volatilities</i>: <i> spotVolsL</i>.
 *      <br>
 *   <br> <i> bVolO:</i>
 *      <br>
 *         <br> <i> baseDate</i>: <i> refDateL</i>. 
 *         <br> <i> dates</i>: <i> bvDatesL</i>.
 *         <br> <i> volatilities</i>: <i> bvRatesL</i>.
 *      <br>
 *   <br> <i> swVolO:</i> convert <i> swVolL</i> into a [<i> swExpL</i> $x$ 
 *	   <i> swMatL</i>] matrix.
 *	<br>
 *	   <br> Use <i> row1, row2,</i> \ldots, to access each row 
 *		 of the swaption volatility matrix <i> swVolO</i>.
 *    	   <br> Use <i> col1, col2,</i> \ldots, to access each column
 *		 of the swaption volatility matrix <i> swVolO</i>.
 *      <br>
 * \end{list}
 */
DLL_EXPORT(int)
VnfmCalib1V2FGeneralNewO(
	TDateL         *refDateL,   /*  1 'D' (I) reference date */
	TDateL         *zcDateL,    /*  2 'D' (I) zero coupondates */
	FloatL         *zcRateL,    /*  3 'F' (I) zero coupon rates */

	FloatL         *floatScalarsL, /*  4 'L' (I) num scalars:
				     *        [1] back bone q
				     *        [2] qL (1=normal,0=lognormal)
				     *        [3] qR (1=normal,0=lognormal)
				     *        [4] Fsh
				     *        [5] generate base vol 
				     *        [6] generate swaption matrix 
				     *        [7] minimum volatility rate mat */
	FloatL         *nfParamsL,  /*  5 'F' (I) N-fact params */
	TDateL         *nfDatesL,   /*  6 'D' (I) array of dates */
	CharBlockL     *inVTypeL,   /*  7 'C' (I) input vol type:LOG, NORM */
	FloatL         *rateMatL,   /*  8 'F' (I) rate mat [0..nDates-1] */
	IntL           *rateFreqL,  /*  9 'L' (I) rate freq [0..nDates-1] */
	FloatL         *rateVolL,   /* 10 'F' (I) rate vol [0..nDates-1] */
	CharBlockL     *outVTypeL,  /* 11 'L' (I) output vol type:LOG, NORM */
	TDateIntervalL *bvMatL,	    /* 12 'F' (I) vol maturity */
	TDateL         *bvDatesL,   /* 13 'D' (I) vol dates */
	IntL           *swTypeL,    /* 14 'L' (I) [0]=type, [1]=freq */
	TDateIntervalL *swMatL,	    /* 15 'F' (I) array of mat intervals */
	TDateIntervalL *swExpL,	    /* 16 'F' (I) array of exp intervals */
	TVolCurve      **spotVolO,  /* 17 'F' (O) spot vol */
	TVolCurve      **bVolO,	    /* 18 'F' (O) base vol */
	TMatrix2D      **swVolO)    /* 19 'F' (O) swaption vol matrix */
{
    	static	char routine[] = "VnfmCalib1V2FGeneralNewO";
	int		status = FAILURE;
	FloatL *spotVolsL= NULL;
	FloatL *bvRatesL= NULL;
	FloatL *swVolL= NULL;

	spotVolsL = NEW_ARRAY(double,
			     (long)nfDatesL[0]+1);
	if (spotVolsL IS NULL) 
                goto done;

	bvRatesL  = NEW_ARRAY(double,
                             (long)bvDatesL[0]+1);
        if (bvRatesL IS NULL) 
                goto done;

	swVolL    = NEW_ARRAY(double,
                             (long)(swMatL[0]*swExpL[0])+1);
        if (swVolL IS NULL) 
                goto done;

        spotVolsL[0] = (long)nfDatesL[0];
        bvRatesL[0]  = (long)bvDatesL[0];
        swVolL[0]    = (long)(swMatL[0]*swExpL[0]);

	if ( VnfmCalib1V2FGeneralNewL( 
				    refDateL,
				    zcDateL,
				    zcRateL,
				    floatScalarsL,
				    nfParamsL,
				    nfDatesL,
				    inVTypeL,
				    rateMatL,
				    rateFreqL,
				    rateVolL,
				    outVTypeL,
				    bvMatL,	
				    bvDatesL,
				    swTypeL,
				    swMatL,
				    swExpL,
				    spotVolsL,
				    bvRatesL,
				    swVolL)  == FAILURE)
        goto done;

        /* Construct spot Vol */
        *spotVolO = GtoVolCurveNew((int)nfDatesL[0],
                                   refDateL[1],
                                   nfDatesL+1,
                                   spotVolsL+1,
                                   "V");

        if (*spotVolO IS NULL)
        {
            GtoErrMsg ("%s:Failed to construct spotVolO object.\n",routine);
            goto done;
        }

        /* Construct base Vol */
        *bVolO = GtoVolCurveNew((int)bvDatesL[0],
                                refDateL[1],
                                bvDatesL+1,
                                bvRatesL+1,
                                "Linear");

        if (*bVolO IS NULL)
        {
            GtoErrMsg ("%s:Failed to construct bVolO object.\n",routine);
            goto done;
        }

	*swVolO =  GtoMatrix1DTo2DNew ((long)swExpL[0], 
			                       (long)swMatL[0], 
			                       swVolL + 1);

	if (*swVolO IS NULL)
        {
            GtoErrMsg ("%s:Failed to construct swVol object.\n",routine);
            goto done;
        }

	status = SUCCESS;

    done:
	if (status IS FAILURE)
            GtoErrMsg ("%s: Failed.\n",routine);

        FREE(spotVolsL);
        FREE(bvRatesL);
        FREE(swVolL);

	return (status);
}



/*f-----------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_CAL1SPVOL_ARB_O.
 *                                                             
 * <br><br>
 * Add-in wrapper for function <i> VnfmCalib1VVolCurveArbL</i>.  \\
 *
 * The inputs are the same as <i> VnfmCalib1VVolCurveArbL</i>. \\
 * The output <i> spotVolO</i> returns an object handle of <b> TVolCurve</b>
 * class which contains following fields:                                      
 * <br>
 * 	<br> <i> baseDate</i>: <i> refDateL</i>.
 * 	<br> <i> dates</i>: <i> datesL</i>.
 * 	<br> <i> volatilities</i>: <i> spotVolL</i>. 
 * <br>
 */
DLL_EXPORT(int)
VnfmCalib1VVolCurveArbO(
        TDateL *refDateL,       /*  1 'D' (I) zero coupon value date */
        TDateL *zcDatesL,       /*  2 'D' (I) array of zero coupon dates */
        FloatL *zcRatesL,       /*  3 'F' (I) array of zero coupon rates */
                                /*        Input Model: */
        FloatL *backboneqL,     /*  4 'L' (I) dist type (0=LN, 0.5=N) */
        FloatL *betaL,          /*  5 'F' (I) array of mr coeff */
        FloatL *alphaL,         /*  6 'F' (I) array of weight coeff */
        TDateL *dateL,          /*  7 'D' (I) array of dates */
        FloatL *sigmaL,         /*  8 'F' (I) array of volatilities */
        FloatL *rhoL,           /*  9 'F' (I) array of correlations */
                                /*        Input benchmark vol curve: */
	TDateL *rateResetL,	/* 10,'F' <I> array of rate reset[0..nDates-1]*/
	double *rateMatL,	/* 11,'F' <I> array of vol mat [0..nDates-1] */
	IntL   *rateFreqL,	/* 12,'L' <I> array of vol freq [0..nDates-1] */
	FloatL *rateVolL,	/* 13,'F' <I> array of vol [0..nDates-1] */
        FloatL *tStartL,        /* 14,'F' (I) start/end time for calibration */
                                /*        Output Calibrated Model: */
	TVolCurve **spotVolO)	/* 15 'F' (O) Spot Vol */
{
static  char            routine[] = "VnfmCalib1VVolCurveArbO";
        int             status = FAILURE;
	FloatL 		*spotVolRateL = NULL;

	spotVolRateL = NEW_ARRAY(double,
				 (long)dateL[0]+1);

	if (spotVolRateL IS NULL) {
		GtoErrMsg ("malloc of spotVolRateL failed\n");
		goto done;
	}

	spotVolRateL[0] = (long)dateL[0];

	if ( VnfmCalib1VVolCurveArbL( refDateL,       
        			      zcDatesL,      
        			      zcRatesL,     
        			      backboneqL,   
        			      betaL,       
        			      alphaL,     
        		              dateL,     
        			      sigmaL,   
        			      rhoL,    
				      rateResetL,
				      rateMatL,
				      rateFreqL,
				      rateVolL,
        			      tStartL,  
				      spotVolRateL) == FAILURE ) 
	goto done;

	/* Construct VolCurve */
	*spotVolO = GtoVolCurveNew( (int)dateL[0],
                              	    refDateL[1],
			     	    dateL+1,
                             	    spotVolRateL+1,
                             	    "V");

	if (*spotVolO IS NULL)
	{
	    GtoErrMsg ("%s:Failed to construct spotVolO object.\n",routine);
	    goto done;
	}

	status = SUCCESS;

   done:
        if (status IS FAILURE)
            GtoErrMsg ("%s:Failed to create a volCurve.\n",routine);

        FREE(spotVolRateL);

        return status;
}




/*f-----------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_CAL1SPVOL_O.
 *                                                             
 * <br><br>
 * Add-in wrapper for function <i> VnfmCalib1VVolCurveNewL</i>.  \\
 *
 * The inputs are the same as <i> VnfmCalib1VVolCurveNewL</i>. \\
 * The output <i> spotVolO</i> returns an object handle of <b> TVolCurve</b>
 * class which contains following fields:                                      
 * <br>
 * 	<br> <i> baseDate</i>: <i> refDateL</i>.
 * 	<br> <i> dates</i>: <i> datesL</i>.
 * 	<br> <i> volatilities</i>: <i> spotVolL</i>. 
 * <br>
 */
DLL_EXPORT(int)
VnfmCalib1VVolCurveNewO(
        TDateL *refDateL,       /*  1 'D' (I) zero coupon value date */
        TDateL *zcDatesL,       /*  2 'D' (I) array of zero coupon dates */
        FloatL *zcRatesL,       /*  3 'F' (I) array of zero coupon rates */
                                /*        Input Model: */
        FloatL *backboneqL,     /*  4 'L' (I) dist type (0=LN, 0.5=N) */
        FloatL *betaL,          /*  5 'F' (I) array of mr coeff */
        FloatL *alphaL,         /*  6 'F' (I) array of weight coeff */
        TDateL *dateL,          /*  7 'D' (I) array of dates */
        FloatL *sigmaL,         /*  8 'F' (I) array of volatilities */
        FloatL *rhoL,           /*  9 'F' (I) array of correlations */
                                /*        Input base vol curve: */
	double *rateMatL,	/* 10,'F' <I> array of vol mat [0..nDates-1] */
	IntL   *rateFreqL,	/* 11,'L' <I> array of vol freq [0..nDates-1] */
	FloatL *rateVolL,	/* 12,'F' <I> array of vol [0..nDates-1] */
        FloatL *tStartL,        /* 13,'F' (I) start/end time for calibration */
                                /*        Output Calibrated Model: */
	TVolCurve **spotVolO)	/* 14 'F' (O) Spot Vol */
{
static  char            routine[] = "VnfmCalib1VVolCurveNewO";
        int             status = FAILURE;
	FloatL 		*spotVolRateL = NULL;

	spotVolRateL = NEW_ARRAY(double,
				 (long)dateL[0]+1);

	if (spotVolRateL IS NULL) {
		GtoErrMsg ("malloc of spotVolRateL failed\n");
		goto done;
	}

	spotVolRateL[0] = (long)dateL[0];

	if ( VnfmCalib1VVolCurveNewL( refDateL,       
        			      zcDatesL,      
        			      zcRatesL,     
        			      backboneqL,   
        			      betaL,       
        			      alphaL,     
        		              dateL,     
        			      sigmaL,   
        			      rhoL,    
				      rateMatL,
				      rateFreqL,
				      rateVolL,
        			      tStartL,  
				      spotVolRateL) == FAILURE ) 
	goto done;

	/* Construct VolCurve */
	*spotVolO = GtoVolCurveNew( (int)dateL[0],
                              	    refDateL[1],
			     	    dateL+1,
                             	    spotVolRateL+1,
                             	    "V");

	if (*spotVolO IS NULL)
	{
	    GtoErrMsg ("%s:Failed to construct spotVolO object.\n",routine);
	    goto done;
	}

	status = SUCCESS;

   done:
        if (status IS FAILURE)
            GtoErrMsg ("%s:Failed to create a volCurve.\n",routine);

        FREE(spotVolRateL);

        return status;
}



/*f-----------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_CAL1SPVOL_OLD_O.
 *                                                             
 * <br><br>
 * Add-in wrapper for function <i> VnfmCalib1VVolCurveOldL</i>.  \\
 *
 * The inputs are the same as <i> VnfmCalib1VVolCurveOldL</i>. \\
 * The output <i> spotVolO</i> returns an object handle of <b> TVolCurve</b>
 * class which contains following fields:                                   
 * <br>
 * 	<br> <i> baseDate</i>: <i> refDateL</i>.
 * 	<br> <i> dates</i>: <i> datesL</i>.
 * 	<br> <i> volatilities</i>: <i> spotVolL</i>. 
 * <br>
 */
DLL_EXPORT(int)
VnfmCalib1VVolCurveOldO(
        TDateL *refDateL,       /*  1 'D' (I) zero coupon value date */
        TDateL *zcDatesL,       /*  2 'D' (I) array of zero coupon dates */
        FloatL *zcRatesL,       /*  3 'F' (I) array of zero coupon rates */
                                /*        Input Model: */
        FloatL *backboneqL,     /*  4 'L' (I) dist type (0=LN, 0.5=N) */
        FloatL *betaL,          /*  5 'F' (I) array of mr coeff */
        FloatL *alphaL,         /*  6 'F' (I) array of weight coeff */
        TDateL *dateL,          /*  7 'D' (I) array of dates */
        FloatL *sigmaL,         /*  8 'F' (I) array of volatilities */
        FloatL *rhoL,           /*  9 'F' (I) array of correlations */
                                /*        Input base vol curve: */
        IntL   *calTypeL,       /* 10 'L' (I) calibration type [0] */
        double *volMatL,        /* 11 'F' (I) vol maturity [0] */
        IntL   *volFreqL,       /* 12 'L' (I) volatility frequency [0] */
        TDateL *volDatesL,      /* 13 'D' (I) array of vol dates */
        FloatL *volRates1L,     /* 14 'F' (I) array of vol values index # 1*/
        FloatL *tStartL,        /* 15 'F' (I) start/end time for calibration */
                                /*        Output Calibrated Model: */
	TVolCurve **spotVolO)	/* 16 'F' (O) spot volatility object */
{
static  char            routine[] = "VnfmCalib1VVolCurveOldO";
        int             status = FAILURE;
	FloatL 		*spotVolL = NULL;

	spotVolL = NEW_ARRAY(double,
			     (long)dateL[0]+1);

	if (spotVolL IS NULL) {
		GtoErrMsg ("malloc of spotVolL failed\n");
		goto done;
	}

	spotVolL[0] = (long)dateL[0];

	if (  VnfmCalib1VVolCurveOldL(  refDateL,    
        			     	zcDatesL,       
        				zcRatesL,      
        				backboneqL,  
        				betaL,      
        				alphaL,    
        				dateL,    
        				sigmaL,  
        				rhoL,   
        				calTypeL,    
        				volMatL,    
        				volFreqL,  
        				volDatesL,
        				volRates1L, 
        				tStartL,   
					spotVolL) == FAILURE)
	goto done;
	
	/* Construct VolCurve */
        *spotVolO = GtoVolCurveNew( (int)dateL[0],
                                    refDateL[1],
                                    dateL+1,
                                    spotVolL+1,
                                    "V");

        if (*spotVolO IS NULL)
        {
            GtoErrMsg ("%s:Failed to construct spotVolO object.\n",routine);
            goto done;
        }

        status = SUCCESS;

   done:
        if (status IS FAILURE)
            GtoErrMsg ("%s:Failed to create a volCurve.\n",routine);

        FREE(spotVolL);

        return status;
}



/*f-----------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_CALPAR_SHORT_O.
 *                                                             
 * <br><br>
 * Add-in wrapper for function <i> VnfmCalibParamShortTermL</i>.  \\
 *
 * The inputs are the same as <i> VnfmCalibParamShortTermL</i>. \\
 * The outputs <i> paOptO</i>, <i> optModO</i>, and <i> constModO</i> 
 * convert <i> paOptL</i>, <i> optModL</i>, and <i> constModL</i> into
 * three object handles of <b> MAT</b> class.  Use field <i> col1</i> 
 * to access the return values.  
 */
DLL_EXPORT(int)
VnfmCalibParamShortTermO(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDatesL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRatesL,	/*  3 'F' (I) zero coupon rates */

	TDateL *datesL,		/*  4 'D' (I) array of dates */
				/*        Input Model: */
	IntL   *paOptFlagL,	/*  5 'L' (I) optimize flag */
	FloatL *paMinL,		/*  6 'F' (I) min value */
	FloatL *paMaxL,		/*  7 'F' (I) max value */
	FloatL *paMidL,		/*  8 'F' (I) initial value */
	
	IntL   *nOptBenchL,	/*  9 'L' (I) # of optimized benchmarks */
	char   *optBenchL,	/* 10 'C' (I) optim. benchmarks */
	FloatL *optWeightL,	/* 11 'F' (I) optim. market weight */
	FloatL *optMidL,	/* 12 'F' (I) optim. market value */

	IntL   *nConstBenchL,	/* 13 'L' (I) # of constraints benchmarks */
	char   *constBenchL,	/* 14 'C' (I) constr. benchmark */
	FloatL *constMinL,	/* 15 'F' (I) constr. min value */
	FloatL *constMaxL,	/* 16 'F' (I) constr. max value */

	FloatL *floatScalarsL,	/* 17 'F' (I) dist type (0=LN, 0.5=N) */

	TMatrix2D **paOptO,	/* 18 'F' (O) optimal parameters */
	TMatrix2D **optModO,	/* 19 'F' (O) optim. model value */
	TMatrix2D **constModO)	/* 20 'F' (O) constr. model value */

{
static	char		routine[] = "VnfmCalibParamShortTermO";
	int		status = FAILURE;
	FloatL *paOptL = NULL;
	FloatL *optModL = NULL;
	FloatL *constModL = NULL;

	paOptL    = NEW_ARRAY(double,
			   (long)paOptFlagL[0]+1);
	if (paOptL IS NULL) {
		GtoErrMsg ("malloc of paOptL failed\n");
		goto done;
	}

	optModL   = NEW_ARRAY(double,
                           (long)optBenchL[0]+1);
	if (optModL IS NULL) {
		GtoErrMsg ("malloc of optModL failed\n");
		goto done;
	}

	constModL = NEW_ARRAY(double,
                           (long)constBenchL[0]+1);
	if (constModL IS NULL) {
		GtoErrMsg ("malloc of constModL failed\n");
		goto done;
	}

	paOptL[0] = (long)paOptFlagL[0];
	optModL[0] = (long)optBenchL[0];
	constModL[0] = (long)constBenchL[0];

	status = VnfmCalibParamShortTermL( refDateL,	
					   zcDatesL, 
					   zcRatesL,
		                           datesL,	
					   paOptFlagL,	
				  	   paMinL,
					   paMaxL,	
				  	   paMidL,		
     	 				   nOptBenchL,	
					   optBenchL,
					   optWeightL,
					   optMidL,	
					   nConstBenchL,	
					   constBenchL,
					   constMinL,	
					   constMaxL,
					   floatScalarsL,	
					   paOptL,
					   optModL,
					   constModL);

	if (status IS FAILURE)
		goto done;


	*paOptO    = GtoMatrix1DTo2DNew ((long)paOptFlagL[0], 
				                     1, 
				                     paOptL +1 );
	*optModO   = GtoMatrix1DTo2DNew ((long)optBenchL[0], 
				                     1, 
				                     optModL +1 );
	*constModO = GtoMatrix1DTo2DNew ((long)constModL[0], 
				                     1, 
				                     constModL +1 );

   done:
        if (status IS FAILURE)
            GtoErrMsg ("%s:Failed to create a TMatrix2D.\n",routine);

        FREE(paOptL);
        FREE(optModL);
        FREE(constModL);

        return status;
}



/*f-----------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_CALPAR_SQ1V_O.
 *                                                             
 * <br><br>
 * Add-in wrapper for function <tt>VnfmCalibParamSquareL</tt>.  \\
 * The inputs are the same as <tt>VnfmCalibParamSquareL</tt>. \\
 * The outputs <tt> paOptO</tt>, <tt> optModO</tt>, and <tt> constModO</tt> 
 * convert <tt> paOptL</tt>, <tt> optModL</tt>, and <tt> constModL</tt> into
 * three object handles of <b>MAT</b> class.  Use field {\it col1</tt> 
 * to access the return values.  
 */


DLL_EXPORT(int)
VnfmCalibParamSquareO(
	TDateL *refDateL,	/*  1 'D' (I) reference date */
	TDateL *zcDatesL,	/*  2 'D' (I) zero coupondates */
	FloatL *zcRatesL,	/*  3 'F' (I) zero coupon rates */

	TDateL *datesL,		/*  4 'D' (I) array of dates */
				/*        Input Model: */
	IntL   *paOptFlagL,	/*  5 'L' (I) optimize flag */
	FloatL *paMinL,		/*  6 'F' (I) min value */
	FloatL *paMaxL,		/*  7 'F' (I) max value */
	FloatL *paMidL,		/*  8 'F' (I) initial value */
	
	IntL   *nOptBenchL,	/*  9 'L' (I) # of optimized benchmarks */
	char   *optBenchL,	/* 10 'C' (I) optim. benchmarks */
	FloatL *optWeightL,	/* 11 'F' (I) optim. market weight */
	FloatL *optMidL,	/* 12 'F' (I) optim. market value */

	IntL   *nConstBenchL,	/* 13 'L' (I) # of constraints benchmarks */
	char   *constBenchL,	/* 14 'C' (I) constr. benchmark */
	FloatL *constMinL,	/* 15 'F' (I) constr. min value */
	FloatL *constMaxL,	/* 16 'F' (I) constr. max value */

	FloatL *floatScalarsL,	/* 17 'F' (I) dist type (0=LN, 0.5=N) */

	FloatL *rateMatL,	/* 18 'F' (I) rate mat [0..nDates-1] */
	IntL *rateFreqL,	/* 19 'L' (I) rate freq [0..nDates-1] */
	FloatL *rateVolL,	/* 20 'F' (I) rate vol [0..nDates-1] */

	TMatrix2D **paOptO,	/* 21 'F' (O) optimal parameters */
	TMatrix2D **optModO,	/* 22 'F' (O) optim. model value */
	TMatrix2D **constModO)	/* 23 'F' (O) constr. model value */

{
static	char		routine[] = "VnfmCalibParamSquareO";
	int		status = FAILURE;
	FloatL *paOptL = NULL;
	FloatL *optModL = NULL;
	FloatL *constModL = NULL;

	paOptL    = NEW_ARRAY(double,
			   (long)paOptFlagL[0]+1);
	if (paOptL IS NULL) {
		GtoErrMsg ("malloc of paOptL failed\n");
		goto done;
	}

	optModL   = NEW_ARRAY(double,
                           (long)optBenchL[0]+1);
	if (optModL IS NULL) {
		GtoErrMsg ("malloc of optModL failed\n");
		goto done;
	}

	constModL = NEW_ARRAY(double,
                           (long)constBenchL[0]+1);
	if (constModL IS NULL) {
		GtoErrMsg ("malloc of constModL failed\n");
		goto done;
	}

	paOptL[0]    = (long)paOptFlagL[0];
	optModL[0]   = (long)optBenchL[0];
	constModL[0] = (long)constBenchL[0];

	status = VnfmCalibParamSquareL(
			refDateL,	
			zcDatesL,
			zcRatesL,
			datesL,	
			paOptFlagL,
			paMinL,	
			paMaxL,
			paMidL,		
			nOptBenchL,	
			optBenchL,
			optWeightL,
			optMidL,
			nConstBenchL,
			constBenchL,
			constMinL,
			constMaxL,
			floatScalarsL,	
			rateMatL,
			rateFreqL,
			rateVolL,
			paOptL,	
			optModL,
			constModL);

	if (status IS FAILURE)
		goto done;

	*paOptO    = GtoMatrix1DTo2DNew ((long)paOptFlagL[0], 
				                     1, 
				                     paOptL +1 );
	*optModO   = GtoMatrix1DTo2DNew ((long)optBenchL[0], 
				                     1, 
				                     optModL +1 );
	*constModO = GtoMatrix1DTo2DNew ((long)constModL[0], 
				                     1, 
				                     constModL +1 );

   done:
        if (status IS FAILURE)
            GtoErrMsg ("%s:Failed to create a TMatrix2D.\n",routine);

        FREE(paOptL);
        FREE(optModL);
        FREE(constModL);

        return status;
}



/*f-----------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_SMOOTH_SWMAT_O.
 *                                                             
 * <br><br>
 * Add-in wrapper for function <i> VnfmSmoothSwaptionMatrixL</i>. 
 *
 * The inputs are the same as <i> VnfmSmoothSwaptionMatrixL</i>. \\
 * The outputs <i> nfParamsOutO</i>, <i> outputMktO</i>, and <i> spvolMatO</i> 
 * return three objects of <b> MAT</b> class, which have the following 
 * fields:
 * \begin{list}%
 * {}{\setlength{\leftmargin}{\leftmargin}}
 *   <br> <i> nfParamsOutO</i>
 *      <br>
 *  	   <br> <i> col1</i>: <i> nfParamsOutL</i>.
 *      <br>
 *   <br> <i> outputMktO:</i> convert <i> outputMktL</i> into a
 *                     [<i> swExpL</i> $x$ <i> swMatL</i>] matrix.
 *      <br>
 *         <br> Use <i> row1, row2,</i> \ldots, to access 
 *		 each row of the swaption volatility matrix <i> outputMktO</i>.
 *         <br> Use <i> col1, col2,</i> \ldots, to access 
 *	 	 each column of the swaption volatility matrix <i> outputMktO</i>.
 *      <br>
 *   <br> <i> spvolMatO:</i> convert <i> spvolMatL</i> into a 
 *	               [<i> swExpL</i> $x$ <i> swMatL</i>] matrix.
 *	<br>
 *	   <br> Use <i> row1, row2,</i> \ldots, to access each 
 *		 row of the spot volatility matrix <i> spvolMatO</i>.
 *    	   <br> Use <i> col1, col2,</i> \ldots, to access each 
 *		 column of he spot volatility matrix <i> spvolMatO</i>.
 *      <br>
 * \end{list}
 */
DLL_EXPORT(int)
VnfmSmoothSwaptionMatrixO(
	TDateL *refDateL,	 /* 01 'D' (I) reference date */
	TDateL *zcDateL,	 /* 02 'D' (I) zero coupondates */
	FloatL *zcRateL,	 /* 03 'F' (I) zero coupon rates */

	IntL   *swTypeL,	 /* 04 'L' (I) array of matrix param [2]: */
				 /*        [0] matr type (0=vertical, 1=diag) */
				 /*        [1] vol frequency (1,2,4,12) */
	double *swMatL,		 /* 05 'F' (I) array of mat intervals */
	double *swExpL,		 /* 06 'F' (I) array of exp intervals */
	double *midMktL,	 /* 07 'F' (I) mid market matrix */
	double *bidToMidL,	 /* 08 'F' (I) bid to mid matrix (>=0) */
	double *matWeightL,	 /* 09 'F' (I) maturity weights */
	double *expWeightL,	 /* 10 'F' (I) expiration weights */
	long   *integerScalarsL, /* 11 'L' (I) numeric scalars */
				 /*        [1] optimType */
				 /*        [2] normType */
	double *doubleScalarsL,	 /* 12 'F' (I) float scalars */
				 /*        [1] smoothParam */
	double    *nfParamsInL,	 /* 13 'F' (I) array of parameters */
	TMatrix2D **nfParamsOutO,/* 14 'F' (O) output parameters */
	TMatrix2D **outputMktO,	 /* 15 'F' (O) output swaption matrix */
	TMatrix2D **spvolMatO)	 /* 16 'F' (O) output spot vol matrix */
{
        static	char  routine[] = "VnfmSmoothSwaptionMatrixO";
	int	status = FAILURE;
	FloatL *nfParamsOutL = NULL;	
	FloatL *outputMktL = NULL;
	FloatL *spvolMatL = NULL;

	nfParamsOutL = NEW_ARRAY(double,
                                 (long)nfParamsInL[0]+1);
        if (nfParamsOutL IS NULL) 
                goto done;

        outputMktL   = NEW_ARRAY(double,
                                 (long)(swMatL[0]*swExpL[0])+1);
        if (outputMktL IS NULL) 
                goto done;

        spvolMatL    = NEW_ARRAY(double,
                                 (long)(swMatL[0]*swExpL[0])+1);
        if (spvolMatL IS NULL) 
                goto done;

        nfParamsOutL[0] = (long)nfParamsInL[0];
        outputMktL[0]   = (long)(swMatL[0]*swExpL[0]);
        spvolMatL[0]    = (long)(swMatL[0]*swExpL[0]);

        status = VnfmSmoothSwaptionMatrixL ( 
                         refDateL,
					     zcDateL,
		                 zcRateL,
					     swTypeL,
					     swMatL,	
					     swExpL,
					     midMktL,	
					     bidToMidL,
					     matWeightL,	
					     expWeightL,
					     integerScalarsL,	 
					     doubleScalarsL,	
					     nfParamsInL,
					     nfParamsOutL,
					     outputMktL,
					     spvolMatL);	


	if (status IS FAILURE)
		goto done;

	*nfParamsOutO = GtoMatrix1DTo2DNew ((long)nfParamsOutL[0], 
				                        1, 
				                        nfParamsOutL+1);

	*outputMktO   = GtoMatrix1DTo2DNew ((long)swExpL[0], 
				                        (long)swMatL[0], 
				                        outputMktL+1);

	*spvolMatO    = GtoMatrix1DTo2DNew ((long)swExpL[0], 
				                        (long)swMatL[0], 
		                     		    spvolMatL+1);

    done:
	if (status IS FAILURE)
            GtoErrMsg ("%s: Failed.\n",routine);

        FREE(nfParamsOutL);
        FREE(outputMktL);
        FREE(spvolMatL);

	return (status);
}



/*f-----------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_SWMAT_CHECK_O.
 *                                                             
 * <br><br>
 * Add-in wrapper for function <i> VnfmCheckSwaptionCalibrationL</i>.
 *
 * The inputs are the same as <i> VnfmCheckSwaptionCalibrationL</i>. \\
 * The outputs <i> tExpFailedO</i> and <i> tMatFailedO</i> 
 * return two objects of <b> MAT</b> class, which contains the 
 * following fields:
 * \begin{list}%
 * {}{\setlength{\leftmargin}{\leftmargin}}
 *   <br> <i> tExpFailedO:</i> convert <i> tExpFailedL</i> into a 
 *  	   [<i> swExpL</i>$\cdot$ <i> tExpToCheckL</i> $x$ 1] matrix.  \\
 *    	   Use <i> col1</i> to access the failed expirations. 
 *   <br> <i> tMatFailedO:</i> convert <i> tMatFailedL</i> into a 
 *  	   [<i> swExpL</i>$\cdot$ <i> tExpToCheckL</i> $x$ 1] matrix.  \\
 *    	   Use <i> col1</i> to access the failed maturities. 
 *   <br> <i> outputSWVolO:</i> convert <i> outputSWVolL</i> into a 
 *  	   [<i> swExpL</i> $x$ <i> swMatL</i>] matrix.  \\
 *    	   Use <i> col1</i> to access the failed maturities. 
 * \end{list}
 */
DLL_EXPORT(int)
VnfmCheckSwaptionCalibrationO(
        TDateL *refDateL,       /*  1 'D' (I) reference date */
        TDateL *zcDateL,        /*  2 'D' (I) zero coupondates  */
        FloatL *zcRateL,        /*  3 'F' (I) zero coupon rates */

        FloatL *nfParamsL,      /*  4 'F' (I) N-F params (b1,b2,r,a) */
        FloatL *floatScalarsL,  /*  5 'F' (I) numerical scalars */
                                /*        [1] vol twk size */
                                /*        [2] distribution type (0=LN, 0.5=N) */
                                /*        [3] min spot vol  */
        IntL *swTypeL,          /*  6 'L' (I) mattype, freq  */
        TDateIntervalL *swMatL, /*  7 'F' (I) array of mat intervals */
        TDateIntervalL *swExpL, /*  8 'F' (I) array of exp intervals */
        FloatL *swVolL,         /*  9 'F' (I) swaption volatilities */
 
	IntL   *volCheckFlagL,  /* 10 'L' (I) Vol check flags  */
				/*        [1] check spot vol ratio flag */
				/*            1=Yes, 0=No               */
				/*        [2] vol adjutment flag      */
				/*            1=Yes, 0=No        */
				/*        [3] 1 = output adjusted swap vol*/
				/*            0 = output vol correction */
	FloatL *numScalarsL,    /* 11 'L' (I) numerical scalars  */
				/*        [1] min vol adjust amount  */
				/*        [2] max spot vol raio >1.0     */ 

        FloatL *tExpToCheckL,   /* 12 'F' (I) array of exp to check */
        IntL   *finalMatL,      /* 13 'L' (I) array of final flags  */
                                /*            TRUE = final; FALSE = CMS  */

        TMatrix2D **tExpFailedO, /* 14 'F' (O) output failed expiration list */
        TMatrix2D **tMatFailedO, /* 15 'F' (O) output failed maturity list  */
	TMatrix2D **outputSWVolO)/* 16 'F' (O) output mod swaption vols  */
{
    	static	char routine[] = "VnfmCheckSwaptionCalibrationO";
	int		status = FAILURE;

	long 	numMaxFailedPts = (long)(swExpL[0]*tExpToCheckL[0]);
	long	numSwpMat = (long)(swExpL[0]*swMatL[0]);

	FloatL *tExpFailedL = NULL;
	FloatL *tMatFailedL = NULL;
	FloatL *outputSWVolL = NULL;

	if((tExpFailedL = NEW_ARRAY(double, numMaxFailedPts+1)) == NULL ||
	   (tMatFailedL = NEW_ARRAY(double, numMaxFailedPts+1)) == NULL ||
	   (outputSWVolL = NEW_ARRAY(double, numSwpMat+1)) == NULL)
                goto done;

	tExpFailedL[0] = numMaxFailedPts;
	tMatFailedL[0] = numMaxFailedPts;
	outputSWVolL[0] = numSwpMat;

	if (VnfmCheckSwaptionCalibrationL(refDateL,    
        				  zcDateL,    
        				  zcRateL,   
        				  nfParamsL, 
        				  floatScalarsL,  
        				  swTypeL,         
        				  swMatL,
        				  swExpL,
        				  swVolL,       
					      volCheckFlagL,
					      numScalarsL,
        				  tExpToCheckL,  
        				  finalMatL,
        				  tExpFailedL,  
        				  tMatFailedL,
					      outputSWVolL) == FAILURE)
        	goto done;

        /* Construct tExpFailedO */
        *tExpFailedO = GtoMatrix1DTo2DNew (numMaxFailedPts,
				                           1,
				                           tExpFailedL + 1);
	if (*tExpFailedO IS NULL)
        {
            GtoErrMsg ("%s: failed to construct tExpFailedO object.\n",routine);
            goto done;
        }

        /* Construct tMatFailedO */
        *tMatFailedO = GtoMatrix1DTo2DNew (numMaxFailedPts,
				                           1,
				                           tMatFailedL + 1);
	if (*tMatFailedO IS NULL)
        {
            GtoErrMsg ("%s: failed to construct tMatFailedO object.\n",routine);
            goto done;
        }

        /* Construct outputSWVolO */
        *outputSWVolO = GtoMatrix1DTo2DNew ((long)swExpL[0],
				                            (long)swMatL[0],
				                            outputSWVolL + 1);
	if (*outputSWVolO IS NULL)
        {
            GtoErrMsg ("%s: failed to construct outputSWVolO object\n",routine);
            goto done;
        }

	status = SUCCESS;

 done:
	if (status IS FAILURE)
            GtoErrMsg ("%s: failed.\n",routine);

        FREE(tExpFailedL);
        FREE(tMatFailedL);
        FREE(outputSWVolL);

	return (status);
}





/*f-------------------------------------------------------------
 * <b> Add-in Function Name:</b> TF_CREATE_VNFM_O.
 *                                                             
 * <br><br>
 * Create the VnfmData structure given the zero curve and model
 * parameters:
 * <br>
 * <br>[refDateL] zero curve value date.
 * <br>[zcDateL] zero curve coupon dates.
 * <br>[zcRateL] zero curve coupon rates.
 * <br>[backboneqL] array of 1 double:\\
 *	(1) backboneq (0 for lognormal, 0.5 for normal).
 * <br>[betaL] array of mean-reversion coefficients (as may as factors).
 * <br>[alphaL] array of factor weights (as may as factors).
 * <br>[dateL] array of dates used for the spot volatility
 * (the spot volatility is assumed to be constant between
 * two consecutive dates).\\
 * {\bf WARNING: Unlike the ALIB routines, the convention is that
 * the spot volatilities and correlations of index $i$
 * apply between date $i$ and date $i+1$.
 * The first date in the array MUST be the volatility reference date (today)}.
 * <br>[sigmaL] range of spot volatilities, indexed by dates in row
 * and factors in column. {\bf These volatilities are
 * currently all overwritten, the argument being here for future use.
 * Pass a dummy input value of 0 for all elements in the array}.
 * <br>[rhoL] range of correlations, indexed by dates in in row
 * and by factor in columns in the following way:
 * the correlations $(\rho_{ij}(t))_{1&lt;= i<j &lt;= N}$
 * between factors $i$ and $j$ is stored in the
 * column of index $k$ (where $k=0,\dots\,N-1$) where
 * $k = {N(N-1)/ 2} - {(N-i)(N-i-1)/ 2} + (j-i-1)$
 * with  $0&lt;= k &lt;= {N(N-1)/ 2}-1$.
 * <br>[that] Output the VnfmData structure.
 * <br>
 */
 
DLL_EXPORT(int) VnfmDataCreateO(
	TDateL *refDateL,       /* 1 'D' (I) reference date */
	TDateL *zcDatesL,       /* 2 'D' (I) array of zero coupon dates */     
	FloatL *zcRatesL,       /* 3 'F' (I) array of zero coupon rates */    
 
	FloatL *backboneqL,     /* 4 'F' (I) back bone Q */
	FloatL *betaL,          /* 5 'F' (I) array of mr coeff */
	FloatL *alphaL,         /* 6 'F' (I) array of weight coeff */
	TDateL *dateL,          /* 7 'D' (I) array of dates */
	FloatL *sigmaL,         /* 8 'F' (I) volatility arrays */
	FloatL *rhoL,           /* 9 'F' (I) correlation arrays */
	
	VnfmData **thatO)	/* 10 '' (O) VnfmData structure */
{ 
        static char routine[]="VnfmDataCreateO";
        int   status = FAILURE;


	*thatO = VnfmDataCreate(refDateL,     
				zcDatesL,    
				zcRatesL,   
 
				backboneqL,
				betaL,    
				alphaL,  
				dateL,  
				sigmaL,
				rhoL); 

	if (*thatO == NULL) goto done;

	status = SUCCESS;

done:
	if (status != SUCCESS) {
	    GtoErrMsg("%s: failed.\n", routine);
	}

	return(status);
}


