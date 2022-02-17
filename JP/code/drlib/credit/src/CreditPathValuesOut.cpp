//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CreditPathValuesOut.cpp
//
//   Description : Value grids to be returned by the credit engine
//
//   Author      : Jay Blumenstein
//
//   Date        : 10 June 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Nrfns.hpp"
#include "edginc/CreditPathValuesOut.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/Format.hpp"
#include "edginc/Maths.hpp"

#define TOLERANCE 1E-10

DRLIB_BEGIN_NAMESPACE

CClassConstSP const CreditPathValuesOut::TYPE = CClass::registerClassLoadMethod(
                                                                                "CreditPathValuesOut", typeid(CreditPathValuesOut), load);

DEFINE_TEMPLATE_TYPE(CreditPathValuesOutArray);


void CreditPathValuesOut::load(CClassSP& clazz)
{
    REGISTER(CreditPathValuesOut, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultCreditPathValuesOut);
    clazz->setPublic(); // make visible to EAS/spreadsheet
    FIELD(pathDates, "path dates"); 
    FIELD_MAKE_TRANSIENT(pathDates);
    FIELD(pathDateTypes, "0=exposure date, 1=mtm date, 2=both");
    FIELD_MAKE_TRANSIENT(pathDateTypes);
    FIELD(pathValues, "values on each path date");
    FIELD_MAKE_TRANSIENT(pathValues);
    FIELD(peakExposure, "peak exposure"); // output
    FIELD_MAKE_OPTIONAL(peakExposure);
    FIELD(aveExposure, "average exposure"); // output
    FIELD_MAKE_OPTIONAL(aveExposure);
    FIELD(DRE, "Derivative Risk Equivalent"); // output
    FIELD_MAKE_OPTIONAL(DRE);
    FIELD(cvr, "Capital Valuation Reserve"); // output
    FIELD_MAKE_OPTIONAL(cvr);
    FIELD(pathExposureDates, "exposure dates"); // output
    FIELD_MAKE_OPTIONAL(pathExposureDates);
    FIELD(peakProfile, "array of peak exposure for each default date"); // output
    FIELD_MAKE_OPTIONAL(peakProfile);
    FIELD(aveProfile, "array of average exposure for each default date"); // output
    FIELD_MAKE_OPTIONAL(aveProfile);
    FIELD(stdProfile, "array of exposure stdev for each default date"); // output
    FIELD_MAKE_OPTIONAL(stdProfile);
    FIELD(dreProfile, "array of DRE for each default date"); // output
    FIELD_MAKE_OPTIONAL(dreProfile);
    FIELD(gemBandDates, "dates for GEMS21 banding"); // note: not optional! required in order to compute peak and avg
    FIELD(gemPeak, "GEMS21 peak values");
    FIELD_MAKE_OPTIONAL(gemPeak);
    FIELD(gemAve, "GEMS21 avg values");
    FIELD_MAKE_OPTIONAL(gemAve);
    FIELD(gemStd, "GEMS21 std values");
    FIELD_MAKE_OPTIONAL(gemStd);
    FIELD(gemDRE, "GEMS21 DRE values");
    FIELD_MAKE_OPTIONAL(gemDRE);
    
    // transient fields
    FIELD(exposureCcy, "exposure ccy");
    FIELD_MAKE_TRANSIENT(exposureCcy);
    
    FIELD(lastExposureDate, "last exposure date");
    FIELD_MAKE_TRANSIENT(lastExposureDate);

    FIELD(creditDebug,"credit debug object");
    FIELD_MAKE_OPTIONAL(creditDebug);
}

/** access methods */
DateTimeArraySP& CreditPathValuesOut::getPathDates() {
    return pathDates;
}
CIntArraySP& CreditPathValuesOut::getPathDateTypes() {
    return pathDateTypes;
}
DoubleArrayArraySP& CreditPathValuesOut::getPathValues() {
    return pathValues;
}

CreditPathValuesOut::CreditPathValuesOut() : CObject(TYPE){
    pathDates = DateTimeArraySP(   );
    pathDateTypes = CIntArraySP(   ); 
    pathValues = DoubleArrayArraySP(   ); 
    peakExposure = aveExposure = cvr = DRE = -1.0; // init to invalid results
    lastExposureDate = DateTime(0,0);
    creditDebug = CreditDebugSP(new CreditDebug()); 
    pathExposureDates = DateTimeArraySP(   );
    peakProfile = DoubleArraySP(   );
    aveProfile = DoubleArraySP(   );
    stdProfile = DoubleArraySP(   );
    dreProfile = DoubleArraySP(   );
    gemPeak = DoubleArraySP(   );
    gemAve = DoubleArraySP(   );
    gemStd = DoubleArraySP(   );
    gemDRE = DoubleArraySP(   );
};

/** add another CreditPathValuesOut object's values to this one*/
void CreditPathValuesOut::add(const CreditPathValuesOut* rhs)
{
    const string method = "CreditPathValuesOut::add";
    
    if (exposureCcy != rhs->exposureCcy)
        throw ModelException(method, "cannot add path values of different ccy: " +
                             exposureCcy + " vs " + rhs->exposureCcy);
    
    // need to take the max of the 2 lastExposureDates.
    resetLastExposureDate(rhs);

    const DoubleArrayArraySP source = rhs->pathValues;
    
    int numPaths = pathValues->size();
    /** check outer dimension */
    if (numPaths != source->size())
    {
        throw ModelException(method, "Arrays must be of equal size.");
    }
    if (numPaths == 0)
        return; // nothing to do
    
    /** check inner dimensions */
    int i;  
    for (i = 0; i < numPaths; i++)
    {
        if ((*pathValues)[i].size() != (*source)[i].size())
        {
            throw ModelException(method, "Arrays must be of equal size.  Failed at array index" +
                                 Format::toString(i));
        }
    }
    
    /** add the two arrays and the result is in pathValues */
    int numDates = (*pathValues)[0].size(); // assuming all sizes the same
    for (i = 0; i < numPaths; i++)
    {
        for (int j = 0; j < numDates; j++) 
        {
            (*pathValues)[i][j] += (*source)[i][j];
        }
    }
}

/** calculate peak, average exposures and dre/stdev from path values */
/** also compute GEMS21 banded ave/std/peak/dre */
void CreditPathValuesOut::calcPeakAverageStdAndDRE(
    const CMarketDataSP&        market,
    const MarginAcctSP&        marginAcct,	
    double		                percentile,     
    const string&              creditCcyName,
    const CreditSpreadCurveSP&	creditSpread,
    const IModelSP              model)
{
    const string method = "CreditPathValuesOut::calcPeakAndAverage";
    int j;

    // first, extract the yield curve from the market.
    YieldCurveSP discount;
    try {
                CClassConstSP clazz = CClass::forName("YieldCurve");
        if (!clazz) {
            throw ModelException("Invalid object type", method);
        }
              
        if (!!model)
        {
            MarketObjectSP marketObj = market->GetData(creditCcyName, clazz);
            if ( YieldCurve::TYPE->isInstance(marketObj.get()))
            {
                discount = YieldCurveSP(dynamic_cast<YieldCurve*>(marketObj.get()));
                discount->getMarket(model.get(), market.get());
            }
            else
                throw ModelException("Cannot find yield curve for currency", method);
        }
        else
        {
            MarketObjectSP marketObj = market->GetData(creditCcyName, clazz);
            if ( YieldCurve::TYPE->isInstance(marketObj.get()))
            {
                discount = YieldCurveSP(dynamic_cast<YieldCurve*>(marketObj.get()));
                // can't get market if there's no model specified
                //discount->getMarket(model.get(), market.get());
            }
        }

        
    } catch (exception& e){
        throw ModelException(e, method, "Failed when trying to get yield curve from market.");
    }
    
    
    int numPaths = pathValues->size();
    if (numPaths == 0)
        throw ModelException(method, "no path value data.");
    
    int numDates = (*pathValues)[0].size(); // assuming all sizes the same
    
    ///////////////////////////////////////////////////////////////////////////////
    //
    //  Initialize exposure dates info
    //
    
    CIntArray exposureInds;                                    // keep track of which dates are exposure dates
    pathExposureDates = DateTimeArraySP(new DateTimeArray(0)); // member variable (used for GEMS as well)
    for (j = 0; j < numDates; j++)
    {
        if ((*pathDateTypes)[j] != 1)
        {
            pathExposureDates->push_back((*pathDates)[j]);
            exposureInds.push_back(j);
        }
    }
    int numExposureDates = pathExposureDates->size(); 
    
    // place to store the new exposures, after we've analyzed margin.  Will later be put in
    // creditDebug object.  These numbers are what is used for the profiles computations.
    DoubleArrayArraySP& exposureWithMargin = creditDebug->exposureWithMargin;
    exposureWithMargin = DoubleArrayArraySP(new DoubleArrayArray(numPaths));
    for (j = 0; j < numPaths; j++) { (*exposureWithMargin)[j].resize(numExposureDates); }

    /////////////////////////////////////////////////////////////////////////////// 
    //
    //  Initialize MTM dates info
    //

    // init member pathMTMDates 
    pathMTMDates = DateTimeArraySP(new DateTimeArray(0));
    CIntArray mtmInds;
    for (j = 0; j < numDates; j++)
    {
        if ((*pathDateTypes)[j] != 0)
        {
            pathMTMDates->push_back((*pathDates)[j]);
            mtmInds.push_back(j);
        }
    }

    int numMTMDates = pathMTMDates->size();

    peakProfile = DoubleArraySP(new DoubleArray(numExposureDates));
    aveProfile	= DoubleArraySP(new DoubleArray(numExposureDates));
    dreProfile	= DoubleArraySP(new DoubleArray(numExposureDates));
    stdProfile	= DoubleArraySP(new DoubleArray(numExposureDates));
    
    /** compute average and peak profiles */
    vector<double> v(numPaths);
    
    for (int iPath = 0; iPath < numPaths; iPath++)
    {
        // need to initialize margin account.
        double      vmAcct = marginAcct->upfrontMargin; /* keeps track of variation account  */
        DateTime    vmDate;                             /* date variation account last grown */
        double      vmDiscount;                         /* discount factor for var margin    */
        double      imAcct = marginAcct->initialMargin; /* shorthand for initial margin acct */
        DateTime    imDate = (*pathDates)[0];			/* date for initial margin. use first exposure date, ie value date */
        double      imDiscount;                         /* discount factor for init margin   */
        DateTime    sampleDate;                         /* date assoc'd with current sample*/
        int			marginCallDelay = marginAcct->marginCallDelay;
        
        /*
         * Find date to base margin account on.
         * The variation margin start at upfrontMargin.
         */
        vmDate = market->GetReferenceDate();
        
        DateTime currMarginDate;

        int iMTM = 0; // keep track of what mtm date we are at
        // Iterate over all exposure dates.
        for (int iExp = 0; iExp < numExposureDates; iExp++)
        {   
            DateTime exposureDate = (*pathExposureDates)[iExp];
			DateTime lastMTMDate = marginAcct->marginCallHol->addBusinessDays(exposureDate, -marginCallDelay-1);
            double valueOnExposureDate =(*pathValues)[iPath][ exposureInds[iExp] ];
            // need to update running MTM.
            
            /*
             * Determine collateral available on default date.
             * This requires simulating the variation margin account.
             * Check if we need to consider any more MTM dates before the iExp 
             * exposure date for the possibility of a margin call (which would
             * update the value of vmAccount, vmDate
             */
            
            DateTimeArray& MTMDates = *pathMTMDates;
            while ( iMTM < numMTMDates &&
                    MTMDates[iMTM] <= lastMTMDate )
            {
                double valueOnMTMDate =(*pathValues)[iPath][ mtmInds[iMTM] ];
                DateTime vmDateTmp = MTMDates[iMTM];
                double  vmAcctTmp = valueOnMTMDate;
                double  marginCall;
                
                /*
                 * Grow variation margin, and calc margin call.
                 */
                
                vmDiscount = discount->pv(vmDate,vmDateTmp);                
                
                marginCall = vmAcctTmp - vmAcct/vmDiscount;
                
                /*
                 * Apply appropriate unsecured threshold.
                 */
                
                if (vmAcctTmp >= 0.0)
                {
                    vmAcctTmp = Maths::max(0.0, vmAcctTmp - marginAcct->thresholdMarginCpty);
                }
                else
                {
                    vmAcctTmp = Maths::min(0.0, vmAcctTmp + marginAcct->thresholdMarginJPM);
                }
                
                /*
                 * Check if margin call big enough to bother with.
                 */
                
                // to do: add Maths::abs
                if ( Maths::max(marginCall,0.0) + Maths::max(-marginCall,0.0) > marginAcct->minMarginCall )
                {
                    /*
                     * Margin call was made.
                     */
                    vmAcct = vmAcctTmp;
                    vmDate = vmDateTmp;
                    
                }
                
                iMTM++;
            }
            
            /*
             * Now vmAcct holds the variation margin available at
             * the sample date associated with sd.
             */
            
            /*
             * Grow initial and variation margin to sample date
             */
            
            vmDiscount = discount->pv(vmDate,exposureDate);
            imDiscount = discount->pv(imDate,exposureDate);
            
            /*
             * Finally calculate exposure.
             * = replacement cost - initial margin + variation margin
             * Note: discount factors used to FV margin accounts.
             */

            double exposure = Maths::max(0.0, valueOnExposureDate - 
                                         (imAcct/imDiscount + vmAcct/vmDiscount)); 

            (*exposureWithMargin)[iPath][iExp] = exposure;
            
        } // exposure dates
    } // paths
        
    // now set average/peak/stdev profiles.
    computeExposureProfiles(numExposureDates, percentile, exposureWithMargin);
    
	// now compute DRE profile
	computeDreProfile(numExposureDates, discount, creditSpread);

	// now compute average and peak.
    EdgCreditBandedStatistics(discount, creditSpread);

    // now compute GEMS21 bands.
    calcExposureProfileBandForGEMS(discount);    
}

/** Computes peak/ave/std exposure, given an exposureWithMargin profile 
    for each exposure date and path */
void CreditPathValuesOut::computeExposureProfiles(
                        int numExposureDates,
                        double percentile,
                        const DoubleArrayArraySP& exposureWithMargin)
                        
{
    int numPaths = pathValues->size();
    int percentileCount = (int)(percentile*numPaths + 0.5);

    vector<double> v(numPaths); // will hold all exposures.  Re-used for each date.

    
    /** compute average exposure = average over all dates in profile */
    for (int iExp = 0; iExp < numExposureDates; iExp++) 
    {   
        double expDateSum = 0;
        double varDateSum = 0;
        for (int iPath = 0; iPath < numPaths; iPath++) 
        {
            v[iPath] = (*exposureWithMargin)[iPath][iExp];
            expDateSum += v[iPath];
            varDateSum += v[iPath] * v[iPath];
        }

		expDateSum /= (double)(numPaths);
		varDateSum = varDateSum/(double)(numPaths) - expDateSum * expDateSum;

		if ( varDateSum < - TOLERANCE * max(1.0, expDateSum * expDateSum) )
			throw ModelException("CreditPathValuesOut::computeExposureProfiles",
							"negative variance " + Format::toString(varDateSum) +
							" at " + (*pathExposureDates)[iExp].toString());
		else if ( varDateSum < 0 )
			varDateSum = 0;

        (*aveProfile)[iExp] = expDateSum;
		(*stdProfile)[iExp] = sqrt( varDateSum );
        (*peakProfile)[iExp] = nrselect(percentileCount, numPaths, &*v.begin()-1);
        if ((*peakProfile)[iExp] <0.0) { (*peakProfile)[iExp] = 0.0; }
    }
}
             
         
/*f
* Buckets an exposure profile into a new set of dates.
*
* Note: A bucket exposure profile is uses the same type of data elements
* in the time series, but the data has different meanings.  In an exposure
* profile, the expected exposure is centered on the sample date, and "extends"
* half way to each neighbor.  In standard GEMS-21 definitions, the exposure
* starts at the sample date and runs until next one.  
* 
* For example, a 3M point on the exposure profile is exposure centered around
* 3Months, while the GEMS banded 3M point contains all the exposure between
* the 1-Month and 3-Month points.
*
* Also, code could be optimized into a single loop, but I thought looping
* over bands X exposures would be much simpler.
*
*/

// to do
// put bandDates, peaks, and avgs into a class.


void CreditPathValuesOut::calcExposureProfileBandForGEMS(const YieldCurveSP &discount)
{
    string method="EdgExposureProfileBandForGEMS";
    int i;
    /* 
     * Check that dates are monotonic
     */
    
    int nBands = gemBandDates.size();
    gemPeak = DoubleArraySP(new DoubleArray(nBands));
    gemAve = DoubleArraySP(new DoubleArray(nBands));
    gemStd = DoubleArraySP(new DoubleArray(nBands));
    gemDRE = DoubleArraySP(new DoubleArray(nBands));
    
    if (nBands == 0) return; // nothing to do
    for (i=1; i<nBands; i++)
    {
        if (gemBandDates[i].getDate() <= gemBandDates[i-1].getDate()) 
        {
            throw ModelException(method, "Band date " + gemBandDates[i].toString() + 
                                 " is out of sequence");
        }
    }
    
    /*
     * Loop over each band, calculating it.
     */
    int nExpos = pathExposureDates->size();
    for (i=0; i<nBands; i++)
    {
        double  avgNumerator =   0.0;
        double  dreNumerator =   0.0;
        double  avgDenominator = 0.0;
        DateTime   startOfBand;        /* band starts accumulating exposure */
        DateTime   endOfBand;          /* last date band accumulates exposure */
        
        if (i==0)
        {
            startOfBand = gemBandDates[i].rollDate(-36500);  /* back 100 years */
        }
        else
        {
            startOfBand = gemBandDates[i-1];
        }
        if (i==nBands-1) 
        {
            endOfBand = gemBandDates[i].rollDate(36500);   /* forward 100 years */
        }
        else
        {
            endOfBand = gemBandDates[i];
        }
        
        /*
         * Check for any overlap from each exposure.
         * Note: Some of the date variables are made doubles.  This is so
         * they can accurately take into account half-days.  Otherwise
         * there is an occasional round off error, and double counting
         * of some of the days.
         */
        for (int j=0; j<nExpos; j++)
        {
            DateTime      expThis;                /* this exposure date        */
            DateTime      expBefore;              /* date of exposure before   */
            DateTime      expAfter;               /* date of exposure after    */
            double     expStarts;                 /* date this exposure starts */     
            double     expEnds;                   /* date this exposure ends   */     
            double     startOfOverlap;            /* first date of overlap     */
            double     endOfOverlap;              /* last date of overlap      */
            
            expThis =   (*pathExposureDates)[j];
            expBefore = (*pathExposureDates)[Maths::max(j-1,0)];
            expAfter =  (*pathExposureDates)[Maths::min(j+1,nExpos-1)];
            
            expStarts = ((double)expBefore.getDate() + (double)expThis.getDate()) / 2.0;
            expEnds =   ((double)expAfter.getDate() + (double)expThis.getDate()) / 2.0;
            
            /*
             *  If any of this exposure falls in the band:
             *      whichever starts last or ends first
             */
            startOfOverlap = Maths::max((double)startOfBand.getDate(), expStarts);
            endOfOverlap =   Maths::min((double)endOfBand.getDate(),   expEnds);
            
            double startofBandDouble = (double)startOfBand.getDate();
            double endOfBandDouble = (double)endOfBand.getDate();
            if (startofBandDouble <  expEnds &&
                expStarts <= endOfBandDouble)
            {
                double  df;
                double  t = (endOfOverlap - startOfOverlap);
                
                if (t < 1.0)                /* exposure is at least one day */
                {
                    t = 1.0;
                }
                
                t /= 365.0;                 /* convert from days to years */
                
                                            
                /* Calculate discount factor to middle of overlap */
                df = discount->pv( DateTime( (int)((startOfOverlap + endOfOverlap)/2.0) , 0 ) ); 

                /*
                 * Peak Exposure is the peak of the peaks
                 */
                
                // to do: check that when Tom says maximum here he means "peak" maximum.
                if ((*gemPeak)[i] < (*peakProfile)[j])
                {
                    (*gemPeak)[i] = (*peakProfile)[j];
                }
                
                if ((*gemStd)[i] < (*stdProfile)[j])
                {
                    (*gemStd)[i] = (*stdProfile)[j];
                }
                
                /*
                 *                        sum of PV(exposure)*exposureTime
                 * Expected exposure =    --------------------------------
                 *                        sum of PV(1.0)*exposureTime
                 *
                 * Which is the average exposure weighted by time,
                 * discounted to present value
                 */
                
                avgNumerator +=   df * t * (*aveProfile)[j];
                dreNumerator +=   df * t * (*dreProfile)[j];
                /* 
                 * If this is the last exposure, should pad exposure with
                 * zero.  Do this by extending the time to the end.
                 * Note: Don't do this for the last band, as that's
                 * extends to infinity, and is simply the average after it.
                 */
                if (j==nExpos-1 && i<nBands-1)
                {
                    t = (endOfBandDouble - startOfOverlap) / 365.0;
                }
                
                avgDenominator += df * t;
                
            }
        }
        
        if (Maths::isZero(avgDenominator))
        {
            (*gemAve)[i] = 0.0;
            (*gemDRE)[i] = 0.0;
        }
        else
        {
            (*gemAve)[i] = avgNumerator / avgDenominator;
            (*gemDRE)[i] = dreNumerator / avgDenominator;
        }
    } // i

/*
 * This is really quite lazy.  Instead of being intelligent about
 * creating the banded profile, we'll just remove trailing zeros
 * in a brute force way.  
 *
 * Use the last exposure date in the profile as a proxy for the 
 * instrument maturity date.  This is reasonable, as the credit
 * model will always include that date in the output profile.
 */

//    {
//        DateTime maDateTime = 0;
//
//        if ((*pathExposureDates)->numItems > 0)
//        {
//            maDateTime = (*pathExposureDates)->sampleDate[ (*pathExposureDates)->numItems - 1];
//        }
//
//        EdgExposure(*pathExposureDates)ileTruncateZeros((*pathExposureDates)Out,maDateTime);
//    }

// to do: don't do anything in above if past maturity date

}

void CreditPathValuesOut::resetLastExposureDate(const CreditPathValuesOut* rhs)
{
    if (this->lastExposureDate < rhs->lastExposureDate)
    {
        this->lastExposureDate = rhs->lastExposureDate;
    }
};


void CreditPathValuesOut::flushGrids()
{
    pathValues = DoubleArrayArraySP(   );
    creditDebug->exposureWithMargin = DoubleArrayArraySP(   );
}


/*f
 *
 * Calculates various statistics desired for a banded exposure profile.
 *
 * Expected = average of the exposures, weighted by their lengths, 
 *            (each exposure being present valued).
 *
 * Peak = maximum of all the PEAK_EXPOSURE_CONFIDENCE percentile exposures.  Not present valued.
 *
 * CVR (Capital Valuation Reserve):
 *
 * This does NOT use risk cap.  Rather the following methodology:
 *
 * All the expected banded exposures are summed.
 * Each is discounted to today, and multiplied by the probability
 * that it will occur.  This is determined by multiplying the credit
 * spread times the year fraction covered by the exposure.
 *
 * The justification for this is that the credit spread is like
 * an annual probability number (e.g. if 0, the counter party is
 * like the US Feds, and hence risk free.  If 2%, there is a 2%
 * discount for the likelihood of default, which implies the default
 * is expected 2% of the time.
 *
 * Note: All three statistics are coded in a single routine, as CVR and 
 * Average share a lot of code.  It is also likely that all will be
 * needed at the same time.  
 * 
 * If CVR is not needed, simply pass in a creditSpread curve as NULL.
 *
 */

void CreditPathValuesOut::EdgCreditBandedStatistics(
   YieldCurveSP        discount,     /* (I) Currency to use for discounting  */
   CreditSpreadCurveSP creditSpread) /* (I) Credit spreads  (NULL if no CVR) */

//   double            *expected,     /* (O) Expected exposure                */
//   double            *peak,         /* (O) Peak loan equivalent             */
//   double            *cvr)          /* (O) Capital Valuation Reserve        */
{
    const string method = "EdgCreditBandedStatistics";


    long           i;                   /* iterates over exposures           */
    double         avgNumerator;        /* Numerator for expected exposure   */
    double         avgDenominator;      /* Denominator for expected exposure */
    IYieldCurveSP   riskyCurve;          /* Discount + Credit Spread Yields   */
    bool computingCVR = false;			/* if true we are computing cvr as well */
    int n = numNonZeroExposureDates();	/* number of exposures in profile    */

    /*
     * Zero out the sums
     */
//	*expected =		 0.0;
    avgNumerator =   0.0;
    avgDenominator = 0.0;
    cvr =           0.0;    
//    *peak =          0.0;

    

    if (creditSpread.get() != 0)
    {
        computingCVR = true;
        // create a risky curve
        riskyCurve = creditSpread.get()->makeRiskyCurve(*discount.get());
    }

    /*
     * Iterate over each exposure, calculating each statistic
     */
    for (i=0; i<n; i++)
    {
        double     df;                            /* PV at exposure date     */
        double     riskyDF;                       /* Risky PV at exposure    */
        double     s;                             /* interpd credit spread   */
        double     t;                             /* time weighting factor   */
        DateTime      expThis;                       /* this exposure date      */
        DateTime      expBefore;                     /* date of exposure before */
        DateTime      expAfter;                      /* date of exposure after  */

        /*
         * Date of this exposure 
         */
        expThis = (*pathExposureDates)[i];

        /*
         * Date of exposure before this one
         * or this date if it's the first.
         */
        if (i==0)
        {
            expBefore = expThis;
        }
        else
        {
            expBefore = (*pathExposureDates)[i-1];
        }

        /*
         * Date of exposure after this one
         * or this date if it's the last.
         */
        if (i==n-1)
        {
            expAfter = expThis;
        }
        else
        {
            expAfter = (*pathExposureDates)[i+1];
        }

        /*
         * The exposure is weighted by the time between
         * the midpoints of it and its neighbors.
         *
         * t = midpoint(a,b) - midpoint(b,c)
         *   = (a+b)/2 - (b+c)/2
         *   = (a/2 + b/2) - (b/2 + c/2)
         *   = (a/2 - c/2)
         *   = (a-c)/2
         */

        t = (expAfter.getDate() - expBefore.getDate()) / 2.0;

        if (t < 1.0)                          /* exposure is at least one day*/
        {
            t = 1.0;
        }

        t /= 365.0;                           /* convert from days to years */


        /*
         * Calculate the default probability by calculating the forward
         * spread of (USD curve + credit spread) - USD curve.
         */
        if (!computingCVR)
        {
            s = 0.0;                          /* no curve, use zero */
            riskyDF = 1.0;                    /* no curve, no PV    */
        }
        else
        {
            DateTime   startOfExposure = DateTime( (expBefore.getDate() + expThis.getDate()) / 2 , 0 );
            DateTime   endOfExposure =  DateTime( (expAfter.getDate()  + expThis.getDate()) / 2 , 0 );
            double  usdForward;
            double  riskyForward;

            /*
             * Special case: Exposure over a single date.
             * Simply fudge it to include an extra day.
             */
            if (startOfExposure == endOfExposure)
            {
                endOfExposure.rollDate(1);
            }

            DayCountConventionSP actual365f(new Actual365F());			

            usdForward = discount->fwd(startOfExposure, endOfExposure, actual365f.get(), 0);

            riskyForward = riskyCurve->fwd(startOfExposure, endOfExposure, actual365f.get(), 0);
                              
            s = riskyForward - usdForward;
        }


        /*
         * Calculate discount factors to PV exposure and CVR 
         */
        df = discount->pv(expThis);

        if (computingCVR)
        {
            riskyDF = riskyCurve->pv(expThis);
        }

        /*
         *                        sum of PV(exposure)*exposureTime
         * Expected exposure =    --------------------------------
         *                        sum of PV(1.0)*exposureTime
         *
         * Which is the average exposure weighted by time,
         * discounted to present value
         */

        avgNumerator   +=   df * t * (*aveProfile)[i];
        avgDenominator +=   df * t;

        /*
         * The idea is to multiply the expected exposure (which is in $)
         * by the probability of default.  As the credit spread
         * is an annualized probability, we weight it by the time
         * "covered" by the exposure.  Basically this stretches
         * half-way to each of the exposure's neighbors.
         *
         * We PV at the risky rate.
         */
        

        cvr += riskyDF * t * s * (*aveProfile)[i];
    }

    if (Maths::isZero(avgDenominator))
    {
        aveExposure = 0.0;  // member variable
    }
    else
    {
        aveExposure = avgNumerator / avgDenominator;    // member variable
    }

	int j;

	/*
	 * Peak Exposure is the peak of the peaks
	 */
	peakExposure = 0.0; // member variable
	for (j = 0; j < peakProfile->size(); j++)
	{
		if ((*peakProfile)[j] > peakExposure)
		peakExposure = (*peakProfile)[j];
	}

	/* 
	 * DRE is the peak of DREs
	 */
	DRE = 0.0; // member variable
	for (j = 0; j < dreProfile->size(); j++)
	{
		if ((*dreProfile)[j] > DRE)
		DRE = (*dreProfile)[j];
	}
}

int CreditPathValuesOut::numNonZeroExposureDates() const
{
    return CreditPathValuesOut::numHistoricalDates(
        pathExposureDates, lastExposureDate);
}


/** calculate DRE profile */
void CreditPathValuesOut::computeDreProfile(
		int							numExposureDates,
		const YieldCurveSP			&discount,
		const CreditSpreadCurveSP	&creditSpread)
{
	const string method = "computeDreProfile";
	double ndProb, avg, std;
    IYieldCurveSP   riskyCurve;          /* Discount + Credit Spread Yields   */

	// create a risky curve
    riskyCurve = creditSpread.get()->makeRiskyCurve(*discount.get());

	for(int i=0; i<numExposureDates; i++)
	{
		ndProb = riskyCurve->pv((*pathExposureDates)[i])/discount->pv((*pathExposureDates)[i]);
		avg = (*aveProfile)[i];
		std = (*stdProfile)[i];
		(*dreProfile)[i] = min( sqrt( avg*avg + std*std/ndProb ), 
							max( avg, (*peakProfile)[i] ) );
	}

}
DRLIB_END_NAMESPACE
