/***************************************************************************
 *      SCCS Keyword Information
 *      ------------------------
 *      Module name     :  frmpp.c
 *      Company Name    :  JP Morgan Securities Inc.
 *      Authors         :  Davis (Shuenn-Tyan) Lee
 *                         (Derivatives Research)
 *      Code version    :  1.18
 *      Extracted       :  3/12/97 at 09:54:19
 *      Last Updated    :  3/10/97 at 15:28:11
 ***************************************************************************
 *      FRMS prepay model (author: Forrest Quinn)
 *      Model version:  3.1.18
 *
 *      Copyright 1996 J.P. Morgan & Co. Incorporated. All rights reserved.
 ***************************************************************************/


#include <math.h>


/* GTO hdrs */
extern "C" {
#include "cdate.h"
#include "ldate.h"
#include "cerror.h"
#include "convert.h"
#include "macros.h"
}

#include "mbsbase.h"
#include "ppconst.h"
#include "ppbase.h"
#include "mbsconst.h"
#include "frmpp.h"

#define VERSION_3


#define NUM_DISCRETE_POINTS 401
#define LOWER_RATE_BOUND 0.01
#define UPPER_RATE_BOUND 0.20

ostream& operator<<(ostream& s, const DRSmmVecMatrix& a)
{
    s << *(a.m_item);
    return s;
}

ostream& operator<<(ostream& s, const TMbsSmmVecMatrix& a)
{
    s << a.RatesGrid ;
    s << a.smmVecMatrix ;
    return s;
}

ostream& operator<<(ostream& s, const DRSmmVecMatrixItem& a)
{
    s << a.m_ratePts;
    s << a.m_matrix;
    return s;
}
//TMbsDouble2	smmVecMatrix;
//TMbsDouble	RatesGrid;


DRSmmVecMatrixInput::DRSmmVecMatrixInput(long startMthIx,
                                         long endMthIx,
                                         DArray& trueWacs,
                                         DArray& effWacs,
                                         long effOrigTerm,
                                         DArray& wams,
                                         long ix0Moy,
                                         long inclSchedAmort,
                                         double  twkAgeRamp,
                                         DArray& seasonality,
                                         DArray& seasLogist,
                                         DArray& scurveLogist,
                                         DArray& ltAges,
                                         double wac,
                                         double critRat,
                                         double multBurnedSteep)
                                             : m_startMthIx(startMthIx),
                                               m_endMthIx(endMthIx),
                                               m_trueWacs(trueWacs),
                                               m_effWacs(effWacs),
                                               m_effOrigTerm(effOrigTerm),
                                               m_wams(wams),
                                               m_ix0Moy(ix0Moy),
                                               m_inclSchedAmort(inclSchedAmort),
                                               m_twkAgeRamp(twkAgeRamp),
                                               m_seasonality(seasonality),
                                               m_seasLogist(seasLogist),
                                               m_scurveLogist(scurveLogist),
                                               m_ltAges(ltAges),
                                               m_wac(wac),
                                               m_critRat(critRat),
                                               m_multBurnedSteep(multBurnedSteep),
                                               DRCacheInput("SmmMatrix") {}

DRCacheItem* DRSmmVecMatrixInput::create ()
{
    return new DRSmmVecMatrixItem(*this);
}

DRCacheInput* DRSmmVecMatrixInput::copy ()
{
    return new DRSmmVecMatrixInput(*this);
}

void DRSmmVecMatrixInput::Print(ostream& s) const
{
    s << "Start " << m_startMthIx << endl;
    s << "End " << m_endMthIx << endl;
}

bool DRSmmVecMatrixInput::ok_to_use (DRCacheInput& a)
{
    DRSmmVecMatrixInput& input = (DRSmmVecMatrixInput&) a;

    bool ans = (m_startMthIx == input.m_startMthIx) &&
               (m_endMthIx == input.m_endMthIx) &&
               (m_trueWacs == input.m_trueWacs) &&
               (m_effWacs == input.m_effWacs) &&
               (m_effOrigTerm == input.m_effOrigTerm) &&
               (m_wams == input.m_wams) &&
               (m_ix0Moy == input.m_ix0Moy) &&
               (m_inclSchedAmort == input.m_inclSchedAmort) &&
               (m_seasonality == input.m_seasonality) &&
               (m_seasLogist == input.m_seasLogist) &&
               (m_scurveLogist == input.m_scurveLogist) &&
               (m_ltAges == input.m_ltAges) &&
               (m_wac == input.m_wac) &&
               (m_critRat == input.m_critRat) &&
               (m_multBurnedSteep == input.m_multBurnedSteep);

    return ans;
};

DRSmmVecMatrixItem::DRSmmVecMatrixItem(DRSmmVecMatrixInput& input)
{
    int numDates = input.m_endMthIx - input.m_startMthIx+1;
    m_ratePts.resize(NUM_DISCRETE_POINTS);
    m_matrix.resize (numDates, NUM_DISCRETE_POINTS);

    TMbsDouble grFrmLagRates0(input.m_endMthIx + 1),
                              grFrmLagRates1(input.m_endMthIx+1);

    DArray smmVect0 (numDates);
    DArray smmVect1 (numDates);

    int j;
    for (j=0;j<NUM_DISCRETE_POINTS;j++)
    {
        m_ratePts[j]=LOWER_RATE_BOUND +
            j*(UPPER_RATE_BOUND-LOWER_RATE_BOUND)/(NUM_DISCRETE_POINTS-1);
    }

    for (j=0;j<NUM_DISCRETE_POINTS;j++)
    {
	grFrmLagRates0 = m_ratePts[j];

	/* Get the SMM's for the model before the impact of triggered
	 * prepayments */
	if (SmmVector(input.m_startMthIx,
		input.m_endMthIx,
		input.m_trueWacs,
		input.m_effWacs,
		input.m_effOrigTerm,
		input.m_wams,
		grFrmLagRates0,
		input.m_ix0Moy,
		input.m_inclSchedAmort,
                input.m_twkAgeRamp,
		input.m_seasonality,
		input.m_seasLogist,
		input.m_scurveLogist,
		input.m_ltAges,
		smmVect0) IS FAILURE)
	{
		GtoErrMsg("Error in SmmVector\n");
		return;
	}

	grFrmLagRates1 = MAX(m_ratePts[j], input.m_wac/input.m_critRat);

	if (SmmVector(input.m_startMthIx,
		input.m_endMthIx,
		input.m_trueWacs,
		input.m_effWacs,
		input.m_effOrigTerm,
		input.m_wams,
		grFrmLagRates1,
		input.m_ix0Moy,
		input.m_inclSchedAmort,
		input.m_twkAgeRamp,
		input.m_seasonality,
		input.m_seasLogist,
		input.m_scurveLogist,
		input.m_ltAges,
		smmVect1) IS FAILURE)
	{
		GtoErrMsg("Error in SmmVector\n");
		return;
	}
		
	for (int k = 0; k < numDates; k++) {
            if (smmVect0[k] > smmVect1[k]) {
		double smmPostCall = 1 - (1 - smmVect0[k]) / (1 -  smmVect1[k]);
		smmVect0[k] = 1 - (1 - smmVect1[k]) * pow((1 - smmPostCall), input.m_multBurnedSteep);
	    }
	}			
	m_matrix.store_array(smmVect0, j, 2);
    }
}


int DRSmmVecMatrixItem::bytes()
{
    return sizeof(*this) + m_matrix.size(1) * m_matrix.size(2) * 8 + m_ratePts.size() * 8;
}

DRSmmVecMatrix::DRSmmVecMatrix (long    startMthIx,
			  	long	endMthIx,	
				DArray& trueWacs,
				DArray& effWacs,
				long	effOrigTerm,
				DArray& wams,
				long	ix0Moy,
				long	inclSchedAmort,
                                double  twkAgeRamp,
				DArray& seasonality,
				DArray& seasLogist, 
				DArray& scurveLogist,
				DArray& ltAges, 	
				double	wac,
				double	critRat,
				double  multBurnedSteep)
				       :m_input(startMthIx,
						endMthIx,	
						trueWacs,
						effWacs,
						effOrigTerm,
						wams,
						ix0Moy,
						inclSchedAmort,
                                                twkAgeRamp,
						seasonality,
						seasLogist, 
						scurveLogist,
						ltAges, 	
						wac,
						critRat,
						multBurnedSteep)
{
    m_cacheItemPtr = theCache.get (m_input);
    m_item = (DRSmmVecMatrixItem*) m_cacheItemPtr.c_ptr();
}


double DRSmmVecMatrix::GetSmmVector (double rate, int month)
{
    if (rate < LOWER_RATE_BOUND) rate = LOWER_RATE_BOUND;
    if (rate > UPPER_RATE_BOUND) rate = UPPER_RATE_BOUND;

    double LowerBoundRate = m_item->m_ratePts[0];
    double UpperBoundRate = m_item->m_ratePts[NUM_DISCRETE_POINTS-1];
    double step	= (UpperBoundRate-LowerBoundRate)/(NUM_DISCRETE_POINTS-1) ;
    int    lower_cell = (int) floor((rate-LowerBoundRate)/step);
    int    upper_cell = MIN(lower_cell+1, NUM_DISCRETE_POINTS - 1);
    double lowerSMM = m_item->m_matrix[month][lower_cell];
    double upperSMM = m_item->m_matrix[month][upper_cell];
    double lower_rate = m_item->m_ratePts[lower_cell];
    double upper_rate = m_item->m_ratePts[upper_cell];	
    double SMM;

    if (lower_cell == upper_cell)
	SMM = lowerSMM;
    else
        SMM = lowerSMM+((rate-lower_rate)/(upper_rate-lower_rate))*(upperSMM-lowerSMM);

    return(SMM);
}


// pass upper bd, lower bd
bool TMbsSmmVecMatrix::NeedToReBuild(
    long    startMthIx,	 /* (I) Starting month for which smm's are
                          * to be generated */
    long    endMthIx,	 /* (I) Last month for which smm's are to
                          * be generated */
    double *trueWacs,	 /* (I) Array of gross coupons, starting
                          * with that which applies to prepays in
                          * month startMthIx */
    double *effWacs,	 /* (I) effective WACs (possibly shifted),
                          * starting with that which applies to
                          * prepays in month startMthIx */
    long    effOrigTerm, /* (I) WARM+WALA */
    double *wams,        /* (I) Array of monthly WAM's (not used
                          * for group effects);  start with value
                          * for prepays in month startMthIx; wam
                          * and wala are valid as of the the month
                          * prior to startMthIx (should be in line
                          * with standard quotes) */
    long    ix0Moy,      /* (I) onth of year for month indexed by zero
                          * (integer between 1 and 12) */
    long    inclSchedAmort,/* (I) zero means ignore scheduled
                          * amortization; 1 means include it;
                          * inclSchedAmort indicates whether to add
                          * scheduled amortization to the prepayment
                          * rates */
    double *seasonality, /* (I) Seasonal multiplier for turnover
                          * component */
    double *seasLogist,  /* (I) Array of 4 parameters to get age of
                          * seasoning as a logistic */
    double *scurveLogist,/* (I) array of 4 s-curve logistic parameters
                          * to use (may have been altered by prepay
                          * tweaks) */
    double *ltAges, 	 /* (I) array of six longer term aging
                          * parameters to use adjust refinance ratio */
    double  wac,
    double  critRat,
    double  multBurnedSteep,
    double  twkAgeRamp)
{
    int i;
    if (oldwac != wac) return true;
    if (oldcritRat != critRat) return true;
    if (oldmultBurnedSteep != multBurnedSteep) return true;
    if (oldtwkAgeRamp != twkAgeRamp) return true;

    if (oldSeasonality.size() == 0) return true;
    if (oldStartMthIx != startMthIx) return true;
    if (oldEndMthIx != endMthIx) return true;

    if (oldEffOrigTerm != effOrigTerm) return true;
    if (oldix0Moy != ix0Moy) return true;
    if (oldinclSchedAmort != inclSchedAmort) return true;

    int numMonths = endMthIx - startMthIx +1;
    for (i = 0; i < numMonths; i++) {
        if (oldWams[i] != wams[i]) return true;
        if (oldTrueWacs[i] != trueWacs[i]) return true;
        if (oldEffWacs[i] != effWacs[i]) return true;
    }
    for (i=0; i<12; i++) if (oldSeasonality[i]  != seasonality[i]) return true;
    for (i=0; i<4;  i++) if (oldSeasLogist[i]   != seasLogist[i]) return true;
    for (i=0; i<4;  i++) if (oldSCurveLogist[i] != scurveLogist[i]) return true;
    for (i=0; i<6;  i++) if (oldltAges[i]       != ltAges[i]) return true;

    return false;
}

void TMbsSmmVecMatrix::Build(
    double  Lower_Rate,
    double  Upper_Rate,
    long    startMthIx,     /* (I) Starting month for which smm's
                             * are to be generated */
    long    endMthIx,       /* (I) Last month for which smm's are to be
                             * generated */
    double *trueWacs,       /* (I) Array of gross coupons, starting with
                             * that which applies to prepays in month
                             * startMthIx */
    double *effWacs,        /* (I) effective WACs (possibly shifted),
                             * starting with that which applies to
                             * prepays in month startMthIx */
    long    effOrigTerm,    /* (I) WARM+WALA */
    double *wams,           /* (I) Array of monthly WAM's (not used for
                             * group effects);  start with value for
                             * prepays in month startMthIx; wam and wala
                             * are valid as of the the month prior to
                             * startMthIx (should be in line with
                             * standard quotes) */
    long    ix0Moy,         /* (I) onth of year for month indexed by
                             * zero (integer between 1 and 12) */
    long    inclSchedAmort, /* (I) zero means ignore scheduled
                             * amortization; 1 means include it;
                             * inclSchedAmort indicates whether to add
                             * scheduled amortization to the prepayment
                             * rates */
    double *seasonality,    /* (I) Seasonal multiplier for turnover
                             * component */
    double *seasLogist,     /* (I) Array of 4 parameters to get age of
                             * seasoning as a logistic */
    double *scurveLogist,   /* (I) array of 4 s-curve logistic parameters
                             * to use (may have been altered by prepay
                             * tweaks) */
    double *ltAges,         /* (I) array of six longer term aging
                             * parameters to use adjust refinance ratio */
    double  wac,
    double  critRat,
    double  multBurnedSteep,
    double  twkAgeRamp)
{
    oldSeasonality.resize(12);
    oldSeasLogist.resize(4);
    oldSCurveLogist.resize(4);
    oldltAges.resize(6);

    oldwac = wac;
    oldcritRat = critRat;
    oldmultBurnedSteep = multBurnedSteep;
    oldtwkAgeRamp = twkAgeRamp;

    int numMonths = endMthIx - startMthIx + 1;
    oldTrueWacs.resize(numMonths);
    oldEffWacs.resize(numMonths);
    oldWams.resize(numMonths);

    oldEffOrigTerm = effOrigTerm;
    oldix0Moy = ix0Moy;
    oldinclSchedAmort = inclSchedAmort;
    oldStartMthIx = startMthIx;
    oldEndMthIx	= endMthIx;

    int i;
    for (i=0; i<12; i++) oldSeasonality[i]  = seasonality[i];
    for (i=0; i<4;  i++) oldSeasLogist[i]   = seasLogist[i];
    for (i=0; i<4;  i++) oldSCurveLogist[i] = scurveLogist[i];
    for (i=0; i<6;  i++) oldltAges[i]       = ltAges[i];
    for (i=0; i<numMonths; i++) {
        oldTrueWacs[i] = trueWacs[i];
        oldEffWacs[i] = effWacs[i];
        oldWams[i] = wams[i];
    }

    int j,k;
    double* smmVect0;
    double* smmVect1 ;
    TMbsDouble GridOfRates;
    TMbsDouble grFrmLagRates0(endMthIx + 1), grFrmLagRates1(endMthIx+1);

    GridOfRates.resize(NUM_DISCRETE_POINTS);
    smmVecMatrix.resize (endMthIx-startMthIx+1, NUM_DISCRETE_POINTS);

    if ((smmVect0 = NEW_ARRAY(double,endMthIx-startMthIx+1)) IS NULL)
    {
        GtoErrMsg("Allocation failure for smmVect0\n");
        goto done;
    }

    if ((smmVect1 = NEW_ARRAY(double,endMthIx-startMthIx+1)) IS NULL)
    {
        GtoErrMsg("Allocation failure for smmVect1\n");
        goto done;
    }

    for (j=0;j<NUM_DISCRETE_POINTS;j++)
    {
        GridOfRates[j]=Lower_Rate+j*(Upper_Rate-Lower_Rate)/(NUM_DISCRETE_POINTS-1);
    }

    RatesGrid=GridOfRates;

    for (j=0;j<NUM_DISCRETE_POINTS;j++)
    {
        grFrmLagRates0 = GridOfRates[j];

        /* Get the SMM's for the model before the impact of triggered
         * prepayments */
        if (SmmVector(startMthIx,
                      endMthIx,
                      trueWacs,
                      effWacs,
                      effOrigTerm,
                      wams,
                      grFrmLagRates0,
                      ix0Moy,
                      inclSchedAmort,
                      twkAgeRamp,
                      seasonality,
                      seasLogist,
                      scurveLogist,
                      ltAges,
                      smmVect0) IS FAILURE)
        {
            GtoErrMsg("Error in SmmVector\n");
            goto done;
        }

        grFrmLagRates1 = MAX(GridOfRates[j], wac/critRat);

        if (SmmVector(startMthIx,
                      endMthIx,
                      trueWacs,
                      effWacs,
                      effOrigTerm,
                      wams,
                      grFrmLagRates1,
                      ix0Moy,
                      inclSchedAmort,
                      twkAgeRamp,
                      seasonality,
                      seasLogist,
                      scurveLogist,
                      ltAges,
                      smmVect1) IS FAILURE)
        {
            GtoErrMsg("Error in SmmVector\n");
            goto done;
        }

        for (k = 0; k < endMthIx-startMthIx+1; k++) {
            if (smmVect0[k] > smmVect1[k]) {
                double smmPostCall = 1 - (1 - smmVect0[k]) / (1 -  smmVect1[k]);
                smmVect0[k] = 1 - (1 - smmVect1[k]) * pow((1 - smmPostCall),multBurnedSteep);
            }
        }

        TMbsDouble mySmmVec (smmVect0, endMthIx-startMthIx+1);
        smmVecMatrix.store_array(mySmmVec,j, 2);
    }
done:
    FREE(smmVect0);
    FREE(smmVect1);
}


bool TMbsAvgLifeVec::NeedToReBuild(
    double  baseCPR,   /* (I) a continuous CPR */
    double  refWAC)    /* (I) continuous coupon */
{
    if (oldbaseCPR != baseCPR) return true;
    if (oldrefWAC != refWAC) return true;

    return false;
}

void TMbsAvgLifeVec::Build(double  baseCPR, /* (I) a continuous CPR */
                           double  refWAC)  /* (I) continuous coupon */
{
    oldbaseCPR = baseCPR;
    oldrefWAC = refWAC;

    avgLife.resize (361);

    for (int i = 0; i <= 360; i++) {
        AvgLife(double(i)/12., baseCPR, refWAC, &avgLife[i]);
    }
}

double TMbsAvgLifeVec::GetAvgLife (int month)
{
    if (month < 0) return(0);
    if (month > 360) return(avgLife[360]);

    return(avgLife[month]);
}


// interpolates if necessary
double TMbsSmmVecMatrix::GetSmmVector (double rate, int month)
{
//    if (rate < LOWER_RATE_BOUND) {cout << rate << endl; throw "Too Low"; }
//    if (rate > UPPER_RATE_BOUND) { cout << rate << endl; throw "Too High" ;}
    if (rate < LOWER_RATE_BOUND) rate = LOWER_RATE_BOUND;
    if (rate > UPPER_RATE_BOUND) rate = UPPER_RATE_BOUND;

    double LowerBoundRate = RatesGrid[0];
    double UpperBoundRate = RatesGrid[NUM_DISCRETE_POINTS-1];
    double step = (UpperBoundRate-LowerBoundRate)/(NUM_DISCRETE_POINTS-1) ;
    int    lower_cell = (int) floor((rate-LowerBoundRate)/step);
    int    upper_cell = MIN(lower_cell+1, NUM_DISCRETE_POINTS - 1);
    double lowerSMM = smmVecMatrix[month][lower_cell];
    double upperSMM = smmVecMatrix[month][upper_cell];
    double lower_rate = RatesGrid[lower_cell];
    double upper_rate = RatesGrid[upper_cell];	

    double SMM;
    if (lower_cell == upper_cell)
        SMM = lowerSMM;
    else SMM = lowerSMM+((rate-lower_rate)/(upper_rate-lower_rate))*(upperSMM-lowerSMM);

    return(SMM);
}



/*************************************************************************
 *  FWD DECLARATION OF STATIC (private) FUNCS
 *************************************************************************/


static int
GpDensity(
    double  incRatio,    /* (I) incentive ratio, i.e., gross WAC divided by
                            relevant rate for refinance opp'y */
    double *groupPs,     /* (I) is an array with the following six parameters
                            for the group distribution, in order:
                                gpMinIncent        : Incentive ratio at which
                                                     there is first a trigger
                                                     effect.
                                pPeakLoc           : Incentive ratio at which
                                                     there is a peak in trigger
                                                     density per remaining (non-
                                                     triggered) balance.
                                gpAsymDenPerPeakDen: Ratio of density per
                                                     remaining balance for
                                                     limiting extreme premium to
                                                     peak density.
                                gpMultiplier       : Speed multiplier to the
                                                     multi-group effect (always
                                                     as fraction of remaining,
                                                     non-triggered balance).
                                gpSkewPwr          : Skewness (internals use a
                                                     power function) for the
                                                     trigger distribution (1 for
                                                     no skew adj).
                                gpSigmaDistrib     : The sigma of the
                                                     distribution. */
    double *densityMid); /* (O) Output the desity of the group distribution at
                            the given incentive level for the group parameters
                            specified. */

static int
GroupFactors(
    long    origMthIx,     /* (I) Index number for the month of origination */
    long    endMthIx,      /* (I) Last month for which smm's are to be
                              generated */
    double *wacs,          /* (I) Gross weighted average coupons (this model
                              only allows a single value for all months),
                              where wac[i] is that in month w/lagged refi
                              rate = grFrmLagRates[i] */
    double *grFrmLagRates, /* (I) Array of mortgage commitment rates (lagged);
                              grFrmLagRates is the array of market refi
                              interest rates (appropriately lagged) */
    double  gpAddlSmm,     /* (I) Additional smm for a group that is triggered
                              (works as a hazard, e.g., an additional 10% smm
                              means that after the non-triggered prepayments,
                              10% of the survivors would additionally prepay) */
    double *groupPs,       /* (I) Array of NUM_GROUPPS group params */
    double *gpFactors);    /* (O) Fills the array gpFactors with factors that
                              represent only triggered prepayments (net impact
                              across all groups). */

static int
OrigGroups(
    double *groupPs,       /* (I) is an array with the following six
                              parameters for the group distribution, in order:
                              gpMinIncent        : Incentive ratio at which
                                                   there is first a trigger
                                                   effect.
                              pPeakLoc           : Incentive ratio at which
                                                   there is a peak in trigger
                                                   density per remaining (non-
                                                   triggered) balance.
                              gpAsymDenPerPeakDen: Ratio of density per
                                                   remaining balance for
                                                   limiting extreme premium to
                                                   peak density.
                              gpMultiplier       : Speed multiplier to the
                                                   multi-group effect (always
                                                   as fraction of remaining,
                                                   non-triggered balance).
                              gpSkewPwr          : Skewness (internals use a
                                                   power function) for the
                                                   trigger distribution (1
                                                   for no skew adj).
                              gpSigmaDistrib     : The sigma of the
                                                   distribution. */
    double *origMinRatios, /* (O) Original minimum incentive ratios for each
                              group */
    double *origGpWghts,   /* (O) Multi-group weights for new production
                              (before any burnout) */
    double *maxRatios);    /* (O) Maximum incentive ratios for each group */

static int
HazardNormal(
    double  xVal,       /* (I) Incentive ratio at which there is first a
                           trigger effect */
    double  xPeak,      /* (I) Incentive ratio at which there is a peak in
                           trigger density per remaining (non-triggered)
                           balance */
    double  asymDensity,/* (I) Ratio of density per remaining balance for
                           limiting extreme premium to peak density */
    double  multiplier, /* (I) Speed multiplier to the multi-group effect
                           (always as fraction of remaining, non-triggered
                           balance) */
    double  skewPwr,    /* (I) Skewness (internals use a power function) for
                           the trigger distribution (1 for no skew adj) */
    double  sigma,      /* (I) controls the sigma of the distribution */
    double *gp_Density);/* (O) A density function to be used as distribution
                           of groups.  It is called a "hazard" because it is
                           actually the fraction of remaining (non-triggered)
                           distribution.  Thus, the values do not sum to 1
                           and in fact are usually of infinite sum, with the
                           function not converging to zero for large xVal
                           (in fact it goes to asymDensity). */

int
HistSmmFix(
    double  grCpnUse,      /* (I) The gross coupon used for prepayment
                              incentive calculation:  sometimes one wants to
                              adjust the actual coupon used by (for example)
                              tilting the coupon towards that expected after
                              a reset */
    double  grCpnForSchAm, /* (I) The gross coupon used for calculation of
                              scheduled amortization */
    double  age,           /* (I) Age at the end of the monthly period */
    double  longIxLagged,  /* (I) Long index (e.g., 10-year CMT), appropriately
                              lagged for impact on desired month */
    long    moy,           /* (I) Month of year of prepayment (1 for those
                              reported in January, ..., 12 for those reported
                              in Dec).  Can input 0 to omit any seasonal
                              adjustment */
    double  grFrmSprd,     /* (I) Spread of par "gross" coupon in secondary ARM
                              market (3rd month, TBA) over the long index */
    long    origTermInMths,/* (I) WARM+WALA */
    long    inclSchedAmort,/* (I) zero means ignore scheduled amortization;
                              1 means include it;  inclSchedAmort indicates
                              whether to add scheduled amortization to the
                              prepayment rates */
    double  twkAgeRamp,    /* (I) tweak age ramp */
    double *seasonality,   /* (I) Array of 12 values: seasonal multiples for
                              speeds up to current coupon speed */
    double *seasLogist,    /* (I) Array of 4 parameters to get age of
                              seasoning as a logistic */
    double *logisticPs,    /* (I) First 5 values are used as inputs to
                              LogisticFix for the "S"-curve; 6th value is used
                              as parameter to further slow down discounts (so
                              deeper discounts continue to prepay slower even
                              after "S"-curve is flat */
    double *histSmm_Fix);

double
ParamsToCpr2(
    double refiInc,
    double *logisticPs);

double
ParamsToCpr(
    double refiInc,
    double *logisticPs);

double
LogisticFix(
    double x,              /* (I) Logistic takes x as the variable input */
    double rightAsymptote, /* (I) Logistic asymptotes to this value for very
                              large positive x values */
    double leftAsymptote,  /* (I) Logistic asymptotes to this value for very
                              large negative x values */
    double xWidth,         /* (I) The width (in x-space) of the logistic
                              transition */
    double inflection,     /* (I) The x-value at which the logistic has its
                              inflection point (i.e., its steepest
                              derivative) */
    double skewPwr);       /* (I) Used as exponent of power transformation used
                              to skew the logistic about inflection point */

static int
GetSmmSchFix(
    double  age,            /* (I) Age at the end of the monthly period */
    long    origTermInMths, /* (I) WARM+WALA */
    double  grCpn,          /* (I) The gross coupon used for calculation of
                               scheduled amortization */
    double *smmSch);

static int
SeasRampFix(
    double  age,         /* (I) age in months */
    double  twkAgeRamp,  /* (I) tweak age ramp (months) */
    double  incRatio,    /* (I) incentive ratio, i.e., gross WAC divided by
                            relevant rate for refinance opp'y */
    double *seasLogist,  /* (I) Array of 4 parameters to get age of seasoning
                            as a logistic */
    double *seasRamp);   /* (O) Age of seasoning calculated from a logistic
                            function */

static int
IndexAtm(
    long    ixOld,     /* (I) the guess of the "ATM" index (usually just use
                          last month's index here) */
    double  incRatio,  /* (I) a "global" input: stores the max incentive ratio
                          for each of the (multi-) groups */
    double *maxRatios, /* (I) Maximum incentive ratios for each group */
    long   *ixAtm);    /* (O) Returns the index for which the current
                          incentive ratio is within the group of that index
                          (or nearest possible).  For example, if the current
                          ratio is outside the range of all groups, then the
                          appropriate extreme is returned */


static int
CombineLagRates
   (TDate   startDate,
    TDate   amortIndexStartDate, /* (I) start date of rate array */
    long    numAmortIndexRates,  /* (I) # rates in array */
    double *adjAmortIndexRates,  /* (I) array of 2-pt commitment rates */
    long    amortIndexType,      /* (I) type of FHLMC commitment rate used;
                                  * should be either MBS_PP_REFI_TYPE_FH30 
                                  * or MBS_PP_REFI_TYPE_FH15 */
    long    histAmortIndexStartDate, /* (I) starting month (day-of-mon ignored)
                                      * for adjHistAmortIndexRates[] rates */
    long    numHistAmortIndexRates,  /* (I) # of rates in adj_commit[] */
    double *adjHistAmortIndexRates,  /* (I) array of average monthly
                                      * 2-point commitment rates (assumed same
                                      * type as adjAmortIndexRates[]) */
    double *wtlags,
    TDate  *userEarliestCommitDate,
    double *grFrmLagRates);      /* (O) array of lagged rates */

void FillInHistRates (int numRates, double* adjHistRates, int numRatesToFill, double* fillRates);

void PerformLagging (int numRates, double* grFrmLagRates, double *wtLags);


static int
TweakInc
   (long    startMthIx,       /* (I) Starting month number in array (for
                                 grFrmLagRates) */
    long    endMthIx,         /* (I) Ending   month number in array (for
                                 grFrmLagRates) */
    double  twkCallInc,       /* (I) Amount by which to tweak incentive (as
                                 fraction of WAC; e.g., .01 means calls delayed
                                 by 1% of WAC) */
    double  critRat,          /* (I) Critical ratio at which callability is
                                 approximated to begin;  Shift refinance
                                 incentive according to critRat */
    double *wacs,             /* WACs for given security */
    double *grFrmLagRates,    /* (I) Array of mortgage commitment rates
                                 (lagged);  grFrmLagRates is the array of
                                 market refi interest rates (appropriately
                                 lagged) */
    double *grFrmLagRatTran); /* (O) Transformed grFrmLagRates that have the
                                 effect of shifting callability (see
                                 twkCallInc input) */

static
int SetDefaultMgrpTweaks
   (long    mbsAgency,         /* (I) GNMA, etc. */
    long    mbsTerm,           /* (I) MBS_PP_MBSTERM_30, etc */
    double  turnoverTwk,       /* (I) user defined */
    double  incentPtTwk,       /* (I) user defined */
    double  totSteepTwk,       /* (I) user defined */
    double  baseSteepTwk,      /* (I) user defined */
    double  ageRampTwk,        /* (I) user defined */
    double *twkTurn,           /* (O) */
    double *twkCallInc,        /* (O) */
    double *multBurnedSteep,   /* (O) Adjusts steepness of callability beyond
                                * "elbow" (excluding group effects handled by
                                * multTrigSteep */
    double *multTrigSteep,     /* (O) */
    double *twkAgeRamp);       /* (O) */

/***************************************************************************
 *  PRIVATE FUNCTION
 *  Compute the ratio after adjusting for a roll down the mortgage curve and
 *  a long term aging effect. The "roll" captures the fact that 15-year
 *  product should have a lower market coupon than 30-year product (typically
 *  about 50 b.p. in recent history) and extends this relationship for any
 *  given term by an average life interpolation.
 *  The long term aging allows for older product to be less sensitive to the
 *  ratio than newer product; this is modeled by always applying a minimum
 *  incentive multiplier with the remainder subject to exponential decay with
 *  age.
 *
 *  Parameters used in adjusting the refinance ratio. Values assumed or
 *  derived in model estimation were:
 *  baseCPR = .06 (assumed)
 *  ratio15Year = 1.0625 (assumed)
 *  refWAC = .085 (assumed)
 *  premMinMult = .01 (fit) - for premium
 *  discMinMult = .99 (fit) - for discount
 *  decayRateAnn = .14 (fit)
 ***************************************************************************/
static
int AdjRatio
   (double  ratio,        /* (I) Standard incentive measure (should be
                           * consistent with the original term for the given
                           * product). */
    double  ageInMths,    /* (I) Age in months */
    double  origTermMths, /* (I) Original term in months (usually should use
                           * exactly 360 or 180). */
    double  baseCPR,      /* (I) */
    double  ratio15Year,  /* (I) */
    double  refWAC,       /* (I) */
    double  premMinMult,  /* (I) */
    double  discMinMult,  /* (I) */
    double  decayRateAnn, /* (I) */
    double *adjRatio);    /* (O) */

/***************************************************************************
 *  PRIVATE FUNCTION
 *  Routine to calculate average life given a continuous CPR at prepay and
 *  continuous coupon at coupon, with final maturity given by final
 *  (scheduled amortization is included).  All units must be consistent
 *  (e.g., final in years and annualized contin prepays and coupons).
 ***************************************************************************/
static
int AvgLife
   (double  finalMat,  /* (I) final maturity in months */
    double  speed,     /* (I) a continuous CPR */
    double  coupon,    /* (I) continuous coupon */
    double *avgLife);  /* (O) average life */

/***************************************************************************
 *  PRIVATE FUNCTION
 *  Routine to calculate the integral of t*exp(ct) from zero to "final".
 ***************************************************************************/
static
int IntegralTexpCT
   (double  final,
    double  c,
    double *integralTexpCT);

/***************************************************************************
 *  PUBLIC FUNCTIONS
 *  INPUT:
 *  <mbsDeal>
 *  mbsAgency - issuing mbsAgency (e.g., MBS_AGENCY_FNMA, ...).
 *  mbsTerm - Either MBS_PP_REFI_TYPE_FH30 or MBS_PP_REFI_TYPE_FH15;
 *            determines whether model will treat collateral as 30-
 *            or 15-year product (regardless of warm+wala).
 *  <mbsPrepayAssump>
 *  startDate - starting month for prepays.
 *  numAmortMons - # months of prepays needed.
 *  amortForm - form for output prepays--one of: MBS_PP_SPD_CPR,
 *              MBS_PP_SPD_SMM, MBS_PP_SPD_PSA.
 *  wala - WALA of MBS, in months.
 *  warm - WARM of MBS, in months.
 *  grossCpn - Single WAC assumed for group effects; wac is gross coupon
 *             (7 for 7%).
 *  inclSchedAmort - zero means ignore scheduled amortization; 1 means
 *                   include it; inclSchedAmort indicates whether to add
 *                   scheduled amortization to the prepayment rates.
 *  <mbsRateEnv>
 *  amortIndexType - type of FHLMC commitment rate used; should be either
 *                   MBS_PP_REFI_TYPE_FH30 or MBS_PP_REFI_TYPE_FH15.
 *  === user-supplied hist/future 2-pt fh rates ===
 *  amortIndexStartDate - month of first adjAmortIndexRates rate
 *                        (day-of-month ignored).
 *  numAmortIndexRates - number of adjAmortIndexRates[] rates supplied.
 *  adjAmortIndexRates - past & future monthly-average 2-point commitment
 *                       rates (in decimal)
 *  === our hard-coded hist 2-pt fh rates ===
 *  histAmortIndexStartDate - starting month (day-of-mon ignored) for
 *                            adjHistAmortIndexRates[] rates.
 *  numHistAmortIndexRates - # of rates in adjHistAmortIndexRates[].
 *  adjHistAmortIndexRates - array of average monthly 2-point commitment
 *                           rates (assumed same type as adjAmortIndexRates[]).
 ***************************************************************************/
int
frm_mgrp_prepays
   (TMbsDeal           *mbsDeal,           /* (I) */
    TMbsRateEnv        *mbsRateEnv,        /* (I) */
    TMbsPrepayDefaults *mbsPrepayDefaults, /* (I) */
    double *smmArr)     /* (O) UTPUT of smmFixed: Fills smmFixed with array
                         * of smm values; these combine the standard S and
                         * aging-curves with multi-group effects. The first
                         * entry to smmVect corresponds to startMthIx. */
{
    static char routine[] = "frm_mgrp_prepays";
    int status = FAILURE;
    long i;
    long iSeas;
    long imon;
    long endMthIx;
    long startMthIx;

    long ix0moy;         /* (I) onth of year for month indexed by zero
                           (integer between 1 and 12); it is the month of
                           year for the month indexed by 0 */
    long numTotalRates;
    long startWalaForGp; /* Single WALA as of startMthIx used for the
                          * group impact */
    /* 4 parameters for tweaking (two are multiplicative, so base value is 1):
     * twkTurn: Adjusts turnover by this amount (e.g., .01 means 1% increase
     * as hazard in turnover).  twkCallInc: Shifts point at which callability
     * starts without impact to discounts (linear btn c.c. and the "elbow").
     * multBurnedSteep: Adjusts steepness of callability beyond "elbow"
     * (excluding group effects handled by multTrigSteep).  multTrigSteep:
     * Multiplier to steepness of multi-group prepayments (1 means no change;
     * 2 doubles group pace as hazard). */
    double twkTurn;
    double twkCallInc;
    double multBurnedSteep;
    double multTrigSteep;
    double twkAgeRamp;
    double newPrepay;
	
    //**
    static TMbsDouble smmVect0;
    static TMbsDouble wams;
    static TMbsDouble gpFactors; /* For factors tracking only multi-group
                                  * effects (ignoring sch + non-trig prepays) */
    static TMbsDouble trueWacs;    /* actual WACs of bond */
    static TMbsDouble effWacs;     /* effective WACs (may be shifted) */
    static TMbsDouble effAllWacs;  /* effective WACs (hist. + fut.) */
    static TMbsDouble grFrmLagRates;    /* (I) Array of mortgage commitment
                                         * rates (lagged);  grFrmLagRates is
                                         * the array of market refi interest
                                         * rates (appropriately lagged) */
    static long oldNumTotalRates = 0;
    static TMbsDouble grFrmLagRatesTran;
    //**
    double seasonality[DEF_MBS_PP_NUM_SEASONALITY];
    double groupPs[DEF_MBS_PP_MGRP_NUM_GROUPPS];
    TMonthDayYear mdy;
    TDate  userEarliestCommitDate;

    long iCnt = 0;
    long iFromStart;

        static int ii =0;
    /* if not already done, init mbs base lib to run "silently" */
    MbsBaseInit("JPM_PREPAY",FALSE,FALSE,NULL);

    /* safe size for rate arrays */
    numTotalRates = mbsRateEnv->numHistAmortIndexRates +
                    mbsRateEnv->numHistAmortIndexRates + 360;

    if (numTotalRates != oldNumTotalRates) {
	oldNumTotalRates = numTotalRates;
	grFrmLagRates.resize(numTotalRates);
	wams.resize(numTotalRates);
	trueWacs.resize(numTotalRates);
	effWacs.resize(numTotalRates);
	effAllWacs.resize(numTotalRates);
	grFrmLagRatesTran.resize(numTotalRates);
    }

    if (CombineLagRates((mbsDeal->mbsPrepayAssump)->startDate,
                        mbsRateEnv->amortIndexStartDate,
                        mbsRateEnv->numAmortIndexRates,
                        mbsRateEnv->adjAmortIndexRates,
                        mbsRateEnv->amortIndexType,
                        mbsRateEnv->histAmortIndexStartDate,
                        mbsRateEnv->numHistAmortIndexRates,
                        mbsRateEnv->adjHistAmortIndexRates,
                        mbsPrepayDefaults->mgrpLagWgts,
                       &userEarliestCommitDate,
                        grFrmLagRates) IS FAILURE)
    {
        GtoErrMsg("%s : error in CombineLagRates\n",routine);
        return status;
    }

//    if (mbsRateEnv->numAmortIndexRates == 1) {
//
//        for (int h = 0; h<numTotalRates; h++) {
//            std::cerr << h << " " << grFrmLagRates[h] << std::endl;
//        }
//     } 

    /* CombineLagRates makes sure that userEarliestCommitDate
     * is the earlier start date of the user-supplied rates
     * and our stored historical rates */
    if (GtoDateToMDY(userEarliestCommitDate,
                     &mdy) IS FAILURE)
    {
        GtoErrMsg("%s : error in GtoMDYToDate\n",routine);
        return status;
    }
    ix0moy = mdy.month;

    if( MonsDiff(userEarliestCommitDate,
                 (mbsDeal->mbsPrepayAssump)->startDate,
                &startMthIx) ISNT SUCCESS )
    {
        GtoErrMsg("%s : MonsDiff failed\n",routine);
        return status;
    }

	endMthIx = startMthIx + (mbsDeal->mbsPrepayAssump)->numAmortMons - 1;
    startWalaForGp = (mbsDeal->mbsPrepayAssump)->wala + 1;

    if (startMthIx < startWalaForGp+DEF_MBS_PP_MGRP_NUM_WTLAGS)
    {
        GtoErrMsg("%s : Not enough hist. refi. rates (eariest start date %s) WALA (%d)\n",
        routine,GtoFormatDate(userEarliestCommitDate),
        (mbsDeal->mbsPrepayAssump)->wala);
        return status;
    }

    if ((mbsDeal->mbsPrepayAssump)->numAmortMons >= DEF_MBS_PP_MGRP_MAX_PREPAYS)
    {
        GtoErrMsg("%s : Cannot generate more than 400 months of prepays\n",
            routine);
        return status;
    }

    /* Compute true and "effective" (possibly shifted) WACS & WARMs
     * for each month from today into future;
     * Note that wac[i] is WAC of month in which lagged refi rate
     * is grFrmLagRates[i+startMthIx]
     */
    for (i = 0; i < endMthIx-startMthIx+1; i++)
    {
        trueWacs[i] = (mbsDeal->mbsPrepayAssump)->grossCpn;
        /* we no longer shift wacs, so effective WAC is same as true WAC */
        effWacs[i] = trueWacs[i];

        wams[i] = (mbsDeal->mbsPrepayAssump)->warm - i -1;
    }

    /* Also compute array of effective WACs (possibly shifted)
     * from start of historical info, into future;
     * Note that mapping of WAC/rate arrays is different:
     * wac[i] is that of month in which lagged refi rate
     * is grFrmLagRates[i]
     */
    for (i = 0; i <= endMthIx; i++)
    {
        /* we no longer shift WACs, so these effective WACs
         * are just the true WACs */
        effAllWacs[i] = (mbsDeal->mbsPrepayAssump)->grossCpn;
    }

    /* Before using user-supplied tweaks, combine with default 
     * tweaks for this type of collateral */
    if (SetDefaultMgrpTweaks(mbsDeal->mbsAgency, 
                             mbsDeal->mbsTerm,
                             (mbsDeal->mbsPrepayAssump)->turnoverTwk,
                             (mbsDeal->mbsPrepayAssump)->incentPtTwk,
                             (mbsDeal->mbsPrepayAssump)->totSteepTwk,
                             (mbsDeal->mbsPrepayAssump)->baseSteepTwk,
                             (mbsDeal->mbsPrepayAssump)->ageRampTwk,
                             &twkTurn,
                             &twkCallInc,
                             &multBurnedSteep,
                             &multTrigSteep,
                             &twkAgeRamp) ISNT SUCCESS)
    {
        GtoErrMsg("%s: Failed to set net tweaks\n",routine);
        return status;
    }

    if (multTrigSteep IS 0)
    {
        GtoErrMsg("%s : multTrigSteep cannot equal to ZERO\n",routine);
        return status;
    }

    for (iSeas = 0; iSeas < DEF_MBS_PP_NUM_SEASONALITY; iSeas++)
    {
        seasonality[iSeas] = mbsPrepayDefaults->seasonality[iSeas]*(1.+twkTurn);
    }

    /* Shift S-Curve to reflect a change in housing turnover */
/*  Removed following statements and instead of above section 07-08-96
    scurveLogist[0] = 1.-(1.-scurveLogist[0])*(1.-twkTurn);
    scurveLogist[1] = 1.-(1.-scurveLogist[1])*(1.-twkTurn);
*/

    /* Apply multiplier for steepness due to multi-group prapayments */
    for (i = 0; i < DEF_MBS_PP_MGRP_NUM_GROUPPS; i++)
    {
        groupPs[i] = mbsPrepayDefaults->groupPs[i];
    }
    groupPs[3] = groupPs[3] * multTrigSteep;


    if (TweakInc(startMthIx,
                 endMthIx,
                 twkCallInc,
                 mbsPrepayDefaults->critRat,
                 effAllWacs,      
                 grFrmLagRates,
                 grFrmLagRatesTran) IS FAILURE)
    {
        GtoErrMsg("%s : error in TweakInc\n",routine);
        return status;
    }


    smmVect0.resize(endMthIx-startMthIx+1);
	
    static DArray s_seasonality;
    static DArray s_seasLogist;
    static DArray s_scurveLogist;
    static DArray s_ltAges;

    s_seasonality.resize(12);
    s_seasLogist.resize(4);
    s_scurveLogist.resize(9);
    s_ltAges.resize(6);

    int q;
    for (q = 0; q < 12; q++) s_seasonality[q] = seasonality[q];
    for (q = 0; q < 4; q++) s_seasLogist[q] = (mbsPrepayDefaults->seasLogist)[q];
    for (q = 0; q < 9; q++) s_scurveLogist[q] = (mbsPrepayDefaults->scurveLogist)[q];
    for (q = 0; q < 6; q++) s_ltAges[q] = (mbsPrepayDefaults->ltAges)[q];

    static TMbsSmmVecMatrix vecMatrix;
    if (vecMatrix.NeedToReBuild (			
	startMthIx,
	endMthIx,
	trueWacs,//
	effWacs,//
	(mbsDeal->mbsPrepayAssump)->warm+
	(mbsDeal->mbsPrepayAssump)->wala,
	wams,//
	ix0moy,
	(mbsDeal->mbsPrepayAssump)->inclSchedAmort,
	seasonality,//
	mbsPrepayDefaults->seasLogist,//
	mbsPrepayDefaults->scurveLogist,//
	mbsPrepayDefaults->ltAges,//
	(mbsDeal->mbsPrepayAssump)->grossCpn,
	mbsPrepayDefaults->critRat,
	multBurnedSteep,
        twkAgeRamp))
    {
	vecMatrix.Build (
		LOWER_RATE_BOUND,
		UPPER_RATE_BOUND,
		startMthIx,
		endMthIx,
		trueWacs,
		effWacs,
		(mbsDeal->mbsPrepayAssump)->warm+
		(mbsDeal->mbsPrepayAssump)->wala,
		wams,
		ix0moy,
		(mbsDeal->mbsPrepayAssump)->inclSchedAmort,
		seasonality,
		mbsPrepayDefaults->seasLogist,
		mbsPrepayDefaults->scurveLogist,
		mbsPrepayDefaults->ltAges,
		(mbsDeal->mbsPrepayAssump)->grossCpn,
		mbsPrepayDefaults->critRat,
		multBurnedSteep,
                twkAgeRamp);
    }

    for (iCnt = 0; iCnt < endMthIx - startMthIx + 1; iCnt++)
    {
	smmVect0[iCnt] = vecMatrix.GetSmmVector(grFrmLagRatesTran[iCnt+startMthIx], iCnt);
    }
    

    /* Check that a positive group multiplier was requested
     * (othw ignore groups) */
    if (groupPs[3] > 0)
    {
	gpFactors.resize(endMthIx);

        /* Get the factors resulting from triggered (multi-group)
         * prepayments only, ignoring all others */
	 if (GroupFactors(startMthIx-startWalaForGp,
                         endMthIx,
                         effAllWacs,
                         grFrmLagRatesTran,
                         mbsPrepayDefaults->gpAddlSmm,
                         groupPs,
                         gpFactors) IS FAILURE)
        {
            GtoErrMsg("%s : error in GroupFactors\n",routine);
            goto done;
        }
    }

 
    /* Combine triggered and non-triggered prepayments to get the
     * resulting prepayment vector */
    for (iCnt = startMthIx-1; iCnt < endMthIx; iCnt++)
    {
        iFromStart = iCnt-startMthIx+1;
        if (groupPs[3] > 0)
        {
            /* Combine with group effects into vector of speeds */
            smmArr[iFromStart] = 1 - (1 - smmVect0[iFromStart]) *
                                 gpFactors[iCnt] / gpFactors[iCnt-1];
        }
        else
        {  
            smmArr[iFromStart] = smmVect0[iFromStart];
        }
    }
 
//bmc


//    if (ii==0) {
//        ii=1;
//    
//    std::cerr << "start date " << mbsDeal->mbsPrepayAssump->startDate << std::endl
//        << "numAmortMons "<< mbsDeal->mbsPrepayAssump->numAmortMons<< std::endl
//        << "grossCpn "<< mbsDeal->mbsPrepayAssump->grossCpn<< std::endl
//        << "warm "<< mbsDeal->mbsPrepayAssump->warm<< std::endl
//        << "wala "<< mbsDeal->mbsPrepayAssump->wala<< std::endl
//        << "amortForm "<< mbsDeal->mbsPrepayAssump->amortForm<< std::endl
//        << "turnoverTwk "<< mbsDeal->mbsPrepayAssump->turnoverTwk<< std::endl
//        << "incentPtTwk "<< mbsDeal->mbsPrepayAssump->incentPtTwk<< std::endl
//        << "totSteepTwk "<< mbsDeal->mbsPrepayAssump->totSteepTwk
//        << "baseSteepTwk "<< mbsDeal->mbsPrepayAssump->baseSteepTwk<< std::endl
//        << "ageRampTwk "<< mbsDeal->mbsPrepayAssump->ageRampTwk<< std::endl
//        << "desiredRefiPoints "<< mbsPrepayDefaults->desiredRefiPoints<< std::endl 
//        << "ioMultiplier "<< mbsPrepayDefaults->ioMultiplier<< std::endl
//        << "numGroups "<< mbsPrepayDefaults->numGroups<< std::endl
//        << "gpAddlSmm "<< mbsPrepayDefaults->gpAddlSmm << std::endl
//        << "critRat "<< mbsPrepayDefaults->critRat<< std::endl;
//    }

    /* adjust prepay form: currently in smm */
    if ((mbsDeal->mbsPrepayAssump)->amortForm ISNT MBS_PP_SPD_SMM)
    {
        for (imon = 0; imon < endMthIx-startMthIx+1; imon++)
        {
            if (ConvertPrepays(MBS_PP_SPD_SMM,
                              (mbsDeal->mbsPrepayAssump)->amortForm,
                              smmArr[imon],
                              (mbsDeal->mbsPrepayAssump)->wala+imon,
                             &newPrepay) IS FAILURE)
            {
                GtoErrMsg("%s: ERROR returned in ConvertPrepays\n",routine);
                goto done;
            }
            smmArr[imon] = newPrepay;
        }
    }
    //std::cerr << " cpr " << smmArr[0] << std::endl;

    status = SUCCESS;

done:

    if (status IS FAILURE)
    {
        GtoErrMsg("%s : Failed\n",routine);
    }
    return(status);
}


/*************************************************************************
 *  INPUTS to GroupFactors:
 *  origMthIx    : Index number for the month of origination.
 *  endMthIx     : Index of the last month number for which results will
 *                 be generated.
 *  wacs         : Gross weighted average coupons (this model
 *                 only allows a single value for all months),
 *                 where wac[i] is that in month w/lagged refi
 *                 rate = grFrmLagRates[i]
 *   grFrmLagRates: Array of mortgage commitment rates (lagged);
 *                 grFrmLagRates is the array of market refi
 *                 interest rates (appropriately lagged) 
 *  wac          : Gross weighted average coupon (this model only allows
 *                 a single value for all months).
 *  grFrmLagRates: Gross coupon rates available as refinance rates in
 *                 each period (adjusted for lag).
 *  groupPs      : is an array with the following six parameters for the
 *                 group distribution, in order:
 *    gpMinIncent        : Incentive ratio at which there is first a
 *                         trigger effect.
 *    pPeakLoc           : Incentive ratio at which there is a peak in
 *                         trigger density per remaining (non-triggered)
 *                         balance.
 *    gpAsymDenPerPeakDen: Ratio of density per remaining balance for
 *                         limiting extreme premium to peak density.
 *    gpMultiplier       : Speed multiplier to the multi-group effect
 *                         (always as fraction of remaining,
 *                         non-triggered balance).
 *    gpSkewPwr          : Skewness (internals use a power function)
 *                         for the trigger distribution (1 for no skew adj).
 *    gpSigmaDistrib     : The sigma of the distribution.
 *  gpAddlSmm    : Additional smm for a group that is triggered (works
 *                 as a hazard, e.g., an additional 10% smm means that
 *                 after the non-triggered prepayments, 10% of the
 *                 survivors would additionally prepay).
 *
 *  OUTPUT of GroupFactors:
 *  Fills the array gpFactors with factors that represent only triggered
 *  prepayments (net impact across all groups).
 *************************************************************************/
static int
GroupFactors(
    long    origMthIx,
    long    endMthIx,
    double *wacs,
    double *grFrmLagRates,
    double  gpAddlSmm,
    double *groupPs,    /* (I) Array of NUM_GROUPPS group params */
    double *gpFactors)
{
    static char routine[] = "GroupFactors";
    int status = FAILURE;
    long ixAtm;
    long newixAtm;
    long iCount;
    long jCount;
    double sumWghts;
    double incRatio;
    double wghtLoss;
    double fractOfTrig;

    static DArray oldGroupPs;//(DEF_MBS_PP_MGRP_NUM_GROUPPS);
    static double* oldGpFactors = NULL;
    static TMbsDouble gpWghts;       /* Multi-group weights (updated for burnout
                                      * as appropriate) */
    static TMbsDouble origMinRatios; /* Original minimum incentive ratios for
                                      * each group */
    static TMbsDouble origGpWghts;   /* Multi-group weights for new production
                                      * (before any burnout) */
    static TMbsDouble minRatios;     /* Updated minimum incentive ratios for
                                      * each group */
    static TMbsDouble maxRatios;

    bool ans = true;

    if (oldGroupPs.size() == 0) ans = false;
    else {
	for (int i = 0; i < DEF_MBS_PP_MGRP_NUM_GROUPPS; i++)
		ans &= (oldGroupPs[i] == groupPs[i]);
    }

    if (gpWghts.size() == 0 || !ans || oldGpFactors != gpFactors) {
	oldGroupPs.resize(DEF_MBS_PP_MGRP_NUM_GROUPPS);
	gpWghts.resize(DEF_MBS_PP_MGRP_NUM_GROUPS);
	origMinRatios.resize(DEF_MBS_PP_MGRP_NUM_GROUPS);
	origGpWghts.resize(DEF_MBS_PP_MGRP_NUM_GROUPS);
	minRatios.resize(DEF_MBS_PP_MGRP_NUM_GROUPS);
	maxRatios.resize(DEF_MBS_PP_MGRP_NUM_GROUPS);
	oldGroupPs = DArray(groupPs, DEF_MBS_PP_MGRP_NUM_GROUPPS);
	oldGpFactors = gpFactors;
	if (OrigGroups(groupPs,
		origMinRatios,
		origGpWghts,
		maxRatios) IS FAILURE)
	{
		GtoErrMsg("%s : error in OrigGroups\n",routine);
		goto done;
	}	
    }
	
    for (iCount = 0; iCount < DEF_MBS_PP_MGRP_NUM_GROUPS; iCount++)
    {
        minRatios[iCount] = origMinRatios[iCount];
        gpWghts[iCount] = origGpWghts[iCount];
    }
    sumWghts = 1;
 
    /* Initial guess for the search routine */
    ixAtm = 0;
 
    if (origMthIx > 1)
    {
        for (iCount = 0; iCount < origMthIx-1; iCount++)
        {
            /* Starting factor (generally 1) before there is a chance
             * for prepays */
            gpFactors[iCount] = sumWghts;
        }
    }

    for (iCount = origMthIx-1; iCount < endMthIx; iCount++)
    {
        incRatio = wacs[iCount+1] / grFrmLagRates[iCount+1];
 
        /* Find the index for which the trigger is "at-the-money"
         * (or the closest one) */
        if (IndexAtm(ixAtm,
                     incRatio,
                     maxRatios,
                    &newixAtm) IS FAILURE)
        {
            GtoErrMsg("%s : error in indexAtm\n",routine);
            goto done;
        }
        ixAtm = newixAtm;
 
        if (ixAtm > 0)
        {
            for (jCount = 0; jCount < ixAtm; jCount++)
            {
                wghtLoss = gpWghts[jCount] * gpAddlSmm;
                gpWghts[jCount] = gpWghts[jCount] - wghtLoss;
                sumWghts = sumWghts - wghtLoss;
            }
        }
 
        /* For the "at-the-money" group, calculate the portion of the
         * group that is triggered, and adjust the minimum incentive.
         * Note that this approximation was designed to perform well
         * (smoothly and consistently) for up/down shift scenarios */
        fractOfTrig = (incRatio - minRatios[ixAtm]) /
                (maxRatios[ixAtm] - minRatios[ixAtm]);

        if (fractOfTrig > 1)
        {
            fractOfTrig = 1;
        }
        else if (fractOfTrig < 0)
        {
            fractOfTrig = 0;
        }

        wghtLoss = gpWghts[ixAtm] * gpAddlSmm * fractOfTrig;

        gpWghts[ixAtm] = gpWghts[ixAtm] - wghtLoss;
        sumWghts = sumWghts - wghtLoss;

        /* Update the minimum incentive for this group by assuming
         * that the earliest exercisers went first. */

        minRatios[ixAtm] = gpAddlSmm * fractOfTrig * maxRatios[ixAtm] +
                (1 - gpAddlSmm * fractOfTrig) * minRatios[ixAtm];
        gpFactors[iCount] = sumWghts;
    }

    status = SUCCESS;
 
done:
    if (status IS FAILURE)
    {
        GtoErrMsg("%s : Failed\n",routine);
    }
    return(status);
}


/*************************************************************************
 *  INPUTS to indexAtm
 *  ixOld     : the guess of the "ATM" index (usually just use last
 *              month's index here).
 *  maxRatios : an input from OrigGroups function: stores the max incentive
		ratio for each of the (multi-)groups.
 *
 *  OUTPUT of indexAtm:
 *  Returns the index for which the current incentive ratio is within the
 *  group of that index (or nearest possible).
 *  For example, if the current ratio is outside the range of all groups,
 *  then the appropriate extreme is returned.
 *************************************************************************/
static int
IndexAtm(
    long    ixOld,
    double  incRatio,
    double *maxRatios,
    long    *indexAtm)
{
    static char routine[] = "indexAtm";
    int status = FAILURE;
    long iStep;
    long iCount;
 
    if (incRatio > maxRatios[ixOld])
    {
        iStep = 1;
        for (iCount = ixOld; iCount < DEF_MBS_PP_MGRP_NUM_GROUPS; iCount+=iStep)
        { 
            if (incRatio*(double)(iStep) <= maxRatios[iCount]*(double)(iStep))
            {
                (*indexAtm) = iCount;
                return(SUCCESS);
            }
        }
    }     
    else
    {
        iStep = -1;
        for (iCount = ixOld; iCount >= 0; iCount+=iStep)
        { 
            if (incRatio*(double)(iStep) <= maxRatios[iCount]*(double)(iStep))
            {
                (*indexAtm) = iCount;
                if (iStep IS -1)
                {
                    /* When decrementing, identify the change after the
                     * desired value is passed. */
                    (*indexAtm) = (*indexAtm) + 1;
                }
                return(SUCCESS);
            }
        }
    }

    /* If we fall through to here, then result is on one of the extremes. */
    if (iStep IS -1)
    {
        (*indexAtm) = 1 - 1;
    }
    else
    {
        (*indexAtm) = DEF_MBS_PP_MGRP_NUM_GROUPS - 1;
    }

    status = SUCCESS;
 
//done:
    if (status IS FAILURE)
    {
        GtoErrMsg("%s : Failed\n",routine);
    }
    return(status);
}


/*************************************************************************
 *  INPUTS to OrigGroups:
 *  groupPs is an array with the following six parameters, in order:
 *  (1): Incentive ratio at which there is first a trigger effect.
 *  (2): Incentive ratio at which there is a peak in trigger density
 *       per remaining (non-triggered) balance.
 *  (3): Ratio of density per remaining balance for limiting extreme
 *       premium to peak density.
 *  (4): Speed multiplier to the multi-group effect (always as fraction
 *       of remaining, non-triggered balance).
 *  (5): Skewness (internals use a power function) for the trigger
 *       distribution (1 for no skew adj).
 *  (6): The sigma of the distribution.
 *
 *  OUTPUTS of OrigGroups:
 *  Fills several multigroup arrays with their starting (pre-burnout)
 *  values:
 *    origGpWghts   : Multi-group weights for new production (before
 *                    any burnout).
 *    maxRatios     : Maximum incentive ratios for each group.
 *    origMinRatios : Original minimum incentive ratios for each group.
 *    numGroups     : Number of groups used in multi-group.
 *************************************************************************/
static int
OrigGroups(
    double *groupPs,
    double *origMinRatios,
    double *origGpWghts,
    double *maxRatios)
{
    static char routine[] = "OrigGroups";
    int status = FAILURE;
    long   iCount;
    double oldMax;
    double remainingWght;
    double densityMid;
    double ratioMid;
    double decrInv;   /* Used as decrement (step size) of inverse of ratio. */
 
    remainingWght = 1;

    /* We will use almost 100 groups (at least for now). */

    /* Smaller decrements used to avoid a blow-up at zero. */
    decrInv = 1. / (double)(DEF_MBS_PP_MGRP_NUM_GROUPS+1);

    oldMax = 1;
    for (iCount = 0; iCount < DEF_MBS_PP_MGRP_NUM_GROUPS; iCount++)
    {
        origMinRatios[iCount] = oldMax;
        maxRatios[iCount] = 1 / ((1 / oldMax) - decrInv);
        ratioMid = 0.5 * (origMinRatios[iCount] + maxRatios[iCount]);

        /* Density at midpoint. */
        if (ratioMid <= groupPs[0])
        {
            densityMid = 0;
        }
        else
        {   
            GpDensity(ratioMid, groupPs, &densityMid);
        }

        /* Impact of this (hazard-like) density across the range of ratios
         * spanned by this group; applied only to fraction remaining. */
        origGpWghts[iCount] = remainingWght * (1 - exp(densityMid *
                (origMinRatios[iCount] - maxRatios[iCount])));

        /* Update the remaining weight (not in any group so far). */
        remainingWght = remainingWght - origGpWghts[iCount];

        oldMax = maxRatios[iCount];
    }

    /* Remaining fraction of homeowners not triggered within any of
     * our groups is bucketed in the least rational group. */
    origGpWghts[DEF_MBS_PP_MGRP_NUM_GROUPS-1] =
        origGpWghts[DEF_MBS_PP_MGRP_NUM_GROUPS-1] + remainingWght;

    status = SUCCESS;
 
//done:
    if (status IS FAILURE)
    {
        GtoErrMsg("%s : Failed\n",routine);
    }
    return(status);
}


/*************************************************************************
 *  INPUTS to GpDensity:
 *  incRatio           : the ratio of WAC to market coupon used in the
 *                       prepayment funtion.
 *  groupPs            : an array with the following six parameters, in
 *                       order:
 *     gpMinIncent        : Incentive ratio at which there is first a
 *                          trigger effect.
 *     gpPeakLoc          : Incentive ratio at which there is a peak
 *                          in trigger density per remaining
 *                          (non-triggered) balance.
 *     gpAsymDenPerPeakDen: Ratio of density per remaining balance for
 *                          limiting extreme premium to peak density.
 *     gpMultiplier       : Speed multiplier to the multi-group effect
 *                          (always as fraction of remaining,
 *                          non-triggered balance).
 *     gpSkewPwr          : Skewness (internals use a power function)
 *                          for the trigger distribution (1 for no
 *                          skew adj).
 *     gpSigmaDistrib     : The sigma of the distribution.
 *  Sixth entry controls the sigma of the distribution.
 *  WARNING: incRatio must be greater than groupPs(0) (the gpMinIncent
 *           value).
 *
 *  OUTPUT of GpDensity:
 *  Returns the desity of the group distribution at the given incentive
 *  level for the group parameters specified.
 *************************************************************************/
static int
GpDensity(
    double  incRatio,
    double *groupPs,
    double *densityMid)
{
    static char routine[] = "GpDensity";
    int status = FAILURE;
    double logInc;
    double logPeak;
    double incRatioDiff;
    double gp_Density;
 
    incRatioDiff = incRatio - groupPs[0];
    logInc = log(incRatioDiff);
    logPeak = log(groupPs[1] - groupPs[0]);
 
    /* Density of the log-distribution (which is normal). */
    if (HazardNormal(logInc,
                     logPeak,
                     groupPs[2],
                     groupPs[3],
                     groupPs[4],
                     groupPs[5],
                    &gp_Density) IS FAILURE)
    {
        GtoErrMsg("%s : error in hazardNormal\n",routine);
        goto done;
    }
 
    /* Transform log density into a density for original (non-log)
     * measure. */
    (*densityMid) = gp_Density / incRatioDiff;
 
    status = SUCCESS;
 
done:
    if (status IS FAILURE)
    {
        GtoErrMsg("%s : Failed\n",routine);
    }
    return(status);
}


/*************************************************************************
 *  INPUTS to hazardNormal:
 *  See GpDensity above: inputs are an incentive level, and groupPs
 *  parameters. However, this routine is usually
 *  called in log space, where a more normal distribution is used.
 *
 *  OUTPUT of hazardNormal:
 *  A density function to be used as distribution of groups. It is
 *  called a "hazard" because it is actually the
 *  fraction of remaining (non-triggered) distribution. Thus, the values
 *  do not sum to 1 and in fact are
 *  usually of infinite sum, with the function not converging to zero
 *  for large xVal (in fact it goes to asymDensity).
 *************************************************************************/
static int
HazardNormal(
    double  xVal,
    double  xPeak,
    double  asymDensity,
    double  multiplier,
    double  skewPwr,
    double  sigma,
    double *hazardNormal)
{
    static char routine[] = "hazardNormal";
    int status = FAILURE;

    (*hazardNormal) = exp(-(pow(fabs((xVal-xPeak)/sigma),(2.*skewPwr))));
    if (xVal > xPeak)
    {
        /* Adjustment prevents the hazard rate to ever drop lower
         * than asymDensity (from the peak of 1). Note that we have
         * effectively rescaled the volatility used in order to match
         * not only level and slope but also curvature where the left
         * and right sides meet (at xPeak). */
        (*hazardNormal) = asymDensity + (1-asymDensity) * pow((*hazardNormal),
                        (pow((1./(1.-asymDensity)), (2.*skewPwr))));
    }

    /* Rescale according to the multiplier. */
    (*hazardNormal) = (*hazardNormal) * multiplier / (2.507 * sigma);

    status = SUCCESS;

//done:
    if (status IS FAILURE)
    {
        GtoErrMsg("%s : Failed\n",routine);
    }
    return(status);
}


/*************************************************************************
 *  logisticPs: First 5 values are used as inputs to LogisticFix for the
 *  "S"-curve; 6th value is used as parameter to further slow down
 *  discounts (so deeper discounts continue to prepay slower even after
 *  "S"-curve is flat).
 *************************************************************************/
int
SmmVector(
    long    startMthIx,
    long    endMthIx,
    double *trueWacs,       /* (I) Array of gross coupons, starting with that
                               which applies to prepays in month startMthIx */
    double *effWacs,        /* (I) effective WACs (possibly shifted), starting
                             * with that which applies to prepays in 
                             * month startMthIx */
    long    effOrigTerm,
    double *wams,
    double *grFrmLagRates,
    long    ix0Moy,
    long    inclSchedAmort,
    double  twkAgeRamp,
    double *seasonality,
    double *seasLogist,
    double *scurveLogist,
    double *ltAges,
    double *smmVect)
{
    static char routine[] = "SmmVector";
    int status = FAILURE;
    long monthNum;
    long monthIx;
    long moy;
    long numYears;
    double ratio1;
    double ratio2;
    double age;
    double histSmm_Fix;
    double baseCPR;
    double ratio15Year;
    double refWAC;
    double premMinMult;
    double discMinMult;
    double decayRateAnn;
    double lagRateUse;
 
    for (monthNum = startMthIx; monthNum <= endMthIx; monthNum++)
    {
        monthIx = monthNum - startMthIx;
        numYears = (long)((ix0Moy + monthNum - 1) / 12);

        /* Calculates mod(ix0Moy+monthNum, 12), using 12 instead of 0. */
        moy = ix0Moy + monthNum - 12 * numYears;
        age = (double) effOrigTerm - wams[monthIx];

#ifdef VERSION_3
        /* @# Need to modify S-logistic, aging, and multi-group parameters.
         * Need to change interface which added parameters used in adjusting
         * refinance incentive ratio.
         * We'll do this later because the result should be fully tested.
         * Davis Lee, 11/07/1996
         */
        /* The refinance incentive ratio before adjustment for roll down
         * the curve and long term aging.
         */
        ratio1 = effWacs[monthIx] / grFrmLagRates[monthNum];

        /* The ratio after adjusting for roll down the curve and long term
         * aging.
         */
        baseCPR = ltAges[0];
        ratio15Year = ltAges[1];
        refWAC = ltAges[2];
        premMinMult = ltAges[3];
        discMinMult = ltAges[4];
        decayRateAnn = ltAges[5];

        if (AdjRatio(ratio1,
                     age,
                     effOrigTerm,
                     baseCPR,
                     ratio15Year,
                     refWAC,
                     premMinMult,
                     discMinMult,
                     decayRateAnn,
                    &ratio2) IS FAILURE)
        {
            GtoErrMsg("%s: ERROR returned in AdjRatio\n",routine);
            goto done;
        }

        /* Adjust the refinance rate to reflect the desired ratio. */
        lagRateUse = grFrmLagRates[monthNum] * ratio1 / ratio2;
#else
        lagRateUse = grFrmLagRates[monthNum];
#endif

        if (HistSmmFix(effWacs[monthIx],
                       trueWacs[monthIx],
                       age,
                       lagRateUse,
                       moy,
                       0,  
                       effOrigTerm,
                       inclSchedAmort,
                       twkAgeRamp,
                       seasonality,
                       seasLogist,
                       scurveLogist,
                      &histSmm_Fix) IS FAILURE)
        {
            GtoErrMsg("%s: ERROR returned in HistSmmFix\n",routine);
            goto done;
        }
        smmVect[monthIx] = histSmm_Fix;
    }

    status = SUCCESS;

done:
    if (status IS FAILURE)
    {
        GtoErrMsg("%s : Failed\n",routine);
    }
    return(status);
}


/*************************************************************************
 *  INPUTS to HistSmmFix:
 *  WARNING: Change name (to something like fixSmm) to avoid confusion
 *           with the ARM model.
 *
 *  grCpnUse          : The gross coupon used for prepayment incentive
 *                      calculation:  sometimes one wants to adjust the
 *                      actual coupon used by (for example) tilting the
 *                      coupon towards that expected after a reset.
 *  grCpnForSchAm     : The gross coupon used for calculation of scheduled
 *                      amortization.
 *  age               : Age at the end of the monthly period.
 *  longIxLagged      : Long index (e.g., 10-year CMT), appropriately
 *                      lagged for impact on desired month.
 *  moy               : Month of year of prepayment (1 for those reported
 *                      in January, ..., 12 for those reported in Dec).
 *                      Can input 0 to omit any seasonal adjustment.
 *  grArmSprd         : Spread of par "gross" coupon in secondary ARM
 *                      market (3rd month, TBA) over the long index.
 *  origTermInMths    : WARM+WALA
 *  inclSchedAmort    : nonzero means to use scheduled amort; zero means
 *                      to omit scheduled amort.
 *  absoluteRateEffect: Magnitude of "absolute" rate effect vs. relative
 *                      rate impact of coupon vs. refi rate.   As of
 *                      first writing, this value was estimated as .33.
 *  longIxHist        : Historic level of long index, used to benchmark
 *                      the absolute rate effect (faster for lower rates,
 *                      slower for higher rates than this).
 *  grArmSprdHist     : Like grArmSprd, except this is the historic level
 *                      that prevailed when model was estimated.
 *  logisticPs        : Logistic parameters for the "S"-curve (first 5
 *                      values), then 6th value is used to further slow
 *                      discounts.  The 6th value should be less than 1,
 *                      and the negative of its natural log is also used
 *                      as a minimal incentive sensitivity for premiums.
 *
 *  OUTPUT of HistSmmFix:
 *  Returns the model's speed in SMM for the desired single month.
 *************************************************************************/
int
HistSmmFix(
    double  grCpnUse,
    double  grCpnForSchAm,
    double  age,
    double  longIxLagged,
    long    moy,
    double  grFrmSprd,
    long    origTermInMths,
    long    inclSchedAmort,
    double  twkAgeRamp,
    double *seasonality,
    double *seasLogist,
    double *logisticPs,
    double *histSmm_Fix)
{
    static char routine[] = "HistSmmFix";
    int status = FAILURE;
    double refiInc;
    double refiInc0;
    double mthlyCpr0;
    double mthlyCpr;
    double cprApply;
    double smmNotAppl;
    double smmTotal;
    double smmApply;
    double smmSch;
    double seasMult;
	
    /* No incentive when the WAC is the commit rate. */
    refiInc0 = 1;
 
    /* CPR from logistic function @ no incentive. */
#ifdef VERSION_3
    mthlyCpr0 = ParamsToCpr2(refiInc0,
                             logisticPs);
#else
    mthlyCpr0 = ParamsToCpr(refiInc0,
                            logisticPs);
#endif
         
    refiInc = grCpnUse / (longIxLagged + grFrmSprd);
 
    /* Apply logistic parameters to get CPR before age,
     * seasonality adjustments. */
#ifdef VERSION_3
    mthlyCpr = ParamsToCpr2(refiInc,
                            logisticPs);
#else
    mthlyCpr = ParamsToCpr(refiInc,
                           logisticPs);
#endif
         
    if (mthlyCpr > mthlyCpr0)
    {
        /* Prepays are faster than in the case of no incentive: apply
         * seasonality only to the no incentive portion. */
        cprApply = mthlyCpr0;
    }
    else
    {   
        cprApply = mthlyCpr;
    }
 
    /* Change units from CPR to SMM. */
    smmTotal = 1 - pow((1.-mthlyCpr), (1./12.));

    if (moy > 0)
    {
        /* moy > 0 means that the user wants to apply seasonality for
         * given month of year.  Change units from % to pure and from
         * CPR to SMM, here for portion of prepays that use seasonality. */
        smmApply = 1 - pow((1.-cprApply), (1./12.));
 
        /* The prepays that are not dependent on seasonality. */
        smmNotAppl = 1-(1-smmTotal)/(1-smmApply);
 
        /* Total prepayments after applying seasonality only to those
         * up to "0-incentive". */
        smmTotal = 1-(1-smmNotAppl)*(1-seasonality[moy-1]*smmApply);
    }
 
    /* Convert from SMM to CPR. */
    mthlyCpr = 1. - pow((1.-smmTotal), 12.);
 
    /* Seasoning Ramp component */
    if (SeasRampFix(age,
                    twkAgeRamp,
                    refiInc,
                    seasLogist,
                   &seasMult) IS FAILURE)
    {
        GtoErrMsg("%s: ERROR returned in SeasRampFix\n",routine);
        goto done;
    }
 
    /* Adjust CPR for seasoning. */
    mthlyCpr = mthlyCpr * seasMult;
 
    /* Return smm for this month. */
    (*histSmm_Fix) = 1. - pow((1.-mthlyCpr), (1./12.));
 
    /*  Include scheduled amortization if user has selected this. */
    if (inclSchedAmort ISNT 0)
    {
        if (GetSmmSchFix(age,
                         origTermInMths,
                         grCpnForSchAm,
                        &smmSch) IS FAILURE)
        {
            GtoErrMsg("%s: ERROR returned in GetSmmSchFix\n",routine);
            goto done;
        }
        (*histSmm_Fix) = 1 - (1 - (*histSmm_Fix)) * (1 - smmSch);
    }
 
    status = SUCCESS;

done:
    if (status IS FAILURE)
    {
        GtoErrMsg("%s : Failed\n",routine);
    }
    return(status);
}

/*************************************************************************
 *  Function to transform refinance incentive and logistic parameters into
 *  the S-curve CPR (before seas, aging, etc.).
 *  Same as ParamsToCpr, except this is expanded to use logisticPs(6-8),
 *  which allow for an S-curve to be used.
 *  For premiums while discounts (below some incentive) use an exponentially
 *  decaying speed (with some limiting speed).
 *  The additional parameters are:
 *  logisticPs[6] : Incentive ratio at which the exponential behavior begin
 *                  (applies below the ratio).
 *  logisticPs[7] : Minimum speed (applies at a zero incentive).
 *  logisticPs[8] : Power of ratio that is used in the exponential decay.
 *************************************************************************/
double
ParamsToCpr2(
    double  refiInc,
    double *logisticPs)
{
    double expRatio;
    double minSpeed;
    double pwr;
    double modSpeed;
    double cpr;

    expRatio = logisticPs[6];
    if (refiInc < expRatio)
    {
        minSpeed = logisticPs[7];
        pwr = logisticPs[8];

        /* Model speed at transition point between logistic and exponential -
           not at current incentive */
        modSpeed = ParamsToCpr(expRatio,logisticPs);
        cpr = minSpeed + (modSpeed-minSpeed) *pow(refiInc/expRatio,pwr);
    }
    else
    {
        cpr = ParamsToCpr(refiInc,logisticPs);
    }

    return(cpr);
}    

/*************************************************************************
 *  Function to transform refinance incentive and logistic parameters into
 *  the S-curve CPR (before seas, aging, etc.).
 *************************************************************************/
double
ParamsToCpr(
    double  refiInc,
    double *logisticPs)
{
    double speedMult;
    double incOverInfl;
    double cpr;

    /* CPR from logistic function. */
    cpr = LogisticFix(refiInc,
                      logisticPs[0],
                      logisticPs[1],
                      logisticPs[2],
                      logisticPs[3],
                      logisticPs[4]);

/*
    incOverInfl = refiInc - logisticPs[4];
*/
    /* bug found by Forrest    07-08-96 */
    incOverInfl = refiInc - logisticPs[3];

    /* These adjustments guarantee that there is always at least some
     * minimal prepayment response to incentive changes (Unlike a pure
     * "S"-curve, which might totally flatten in extreme premium and/or
     * discount regions). */
    if (incOverInfl < 0)
    {
        /* Adjustment for very lower coupons (so that deeper discounts
         * continue to slow even after "S"-curve has flattened). */
        speedMult = pow(logisticPs[5], (-1.*incOverInfl));
    }
    else
    {   
        /* Use a similar adjustment for high coupons, but here we
         * guarantee that larger premiums continue to pay faster. */
        speedMult = 1 - log(logisticPs[5]) * incOverInfl;
    }

    /* Do not keep all of speed adjustment for high coupons (a flat
     * burned out response is generally better there);  this
     * adjustment logistically phases in the magnitude of speed
     * adjustment, with the full effect for extreme discounts. */
    speedMult = 1 + (speedMult - 1) * (logisticPs[0] - cpr) /
                (logisticPs[0] - logisticPs[1]);

    cpr = cpr * speedMult;

    return(cpr);
}


/*************************************************************************
 *************************************************************************/
static int
GetSmmSchFix(
    double  age,
    long    origTermInMths,
    double  grCpn,
    double *smmSch)
{
    static char routine[] = "GetSmmSchFix";
    int
        status = FAILURE;
    double
        cpnMult,
        cpnMultPwr;

    if (age < (double) origTermInMths)
    {
        cpnMult = 1 + grCpn / 12;
        cpnMultPwr = pow(cpnMult, ((double) origTermInMths-age+1));
        (*smmSch) = (cpnMult - 1) / (cpnMultPwr - 1);
    }
    else
    {   
        /* Schedule calls for complete paydown by this age. */
        (*smmSch) = 1;
    }

    status = SUCCESS;

//done:
    if (status IS FAILURE)
    {
        GtoErrMsg("%s : Failed\n",routine);
    }
    return(status);
}



/*************************************************************************
 *  Calculates the seasoning Ramp. Age of seasoning calculated from a
 *  logistic function.
 *  INPUTS to seasRampFix:
 *  age         : age in months;
 *  incRatio    : incentive ratio, i.e., gross WAC divided by relevant rate
 *                for refinance opp'y.
 *  seasLogist : Array of 4 parameters to get age of seasoning as a logistic
 *            function of incentive:
 *           Maximum age of seasoning
 *           Minimum age of seasoning
 *           Width of transition (in ratios)
 *           Inflection point of transition.
 *************************************************************************/
static int
SeasRampFix(
    double  age,
    double  twkAgeRamp,
    double  incRatio,
    double *seasLogist,
    double *seasRampFix)
{
    static char routine[] = "SeasRampFix";
    int
        status = FAILURE;
    double
        ageOfSeas;
 
    ageOfSeas = LogisticFix(-1.*incRatio,
                            seasLogist[0],
                            seasLogist[1],
                            seasLogist[2],
                            -1.*seasLogist[3],
                            1);
   
    if (fabs(ageOfSeas) < 0.0001)
    {
        ageOfSeas = 0.0001;
    }

    /* Are we on the seasoning Ramp? */
    if ((age+twkAgeRamp) < ageOfSeas)
    {
        (*seasRampFix) = (age+twkAgeRamp) / ageOfSeas;
    }
    else
    {   
        (*seasRampFix) = 1;
    }

    status = SUCCESS;

//done:
    if (status IS FAILURE)
    {
        GtoErrMsg("%s : Failed\n",routine);
    }
    return(status);
}


/*************************************************************************
 *  INPUTS to LogisticFix:
 *  x: Logistic takes x as the variable input.
 *  rightAsymptote: Logistic asymptotes to this value for very large
 *  positive x values.
 *  leftAsymptote: Logistic asymptotes to this value for very large
 *  negative x values.
 *  xWidth: The width (in x-space) of the logistic transition.
 *  inflection: The x-value at which the logistic has its inflection
 *  point (i.e., its steepest derivative).
 *  skewPwr: Used as exponent of power transformation used to skew the
 *  logistic about inflection point.
 *
 *  OUTPUT of LogisticFix:
 *  Returns the logistic function of x with the given input parameters
 *  (see code for exact formula).
 *************************************************************************/
double
LogisticFix(
    double x,
    double rightAsymptote,
    double leftAsymptote,
    double xWidth,
    double inflection,
    double skewPwr)
{
    double
        a,
        bExp,
        c,   
        d,
        xSkew,
        xWidthSkew,
        inflectionSkew,
        logistic_Fix,
	expAdj;

    /* Transformation for skewness is to replace x with xSkew and
     * xWidth with xWidthSkew. */
    if (x/inflection > 0)
    {
        xSkew = pow((x/inflection), skewPwr);
    }
    else
    {   
        xSkew = pow(-1.*(-x/inflection), skewPwr);
    }

    xWidthSkew = skewPwr * xWidth / inflection;
    inflectionSkew = 1;
    expAdj = MAX(xSkew, inflectionSkew) / xWidthSkew;

    a = rightAsymptote;
    bExp = exp(xSkew / xWidthSkew - expAdj);
    d = exp(inflectionSkew / xWidthSkew - expAdj);
    c = leftAsymptote * d;

    logistic_Fix = (a * bExp + c) / (bExp + d);

    return(logistic_Fix);
}

/***************************************************************************
 *  PRIVATE FUNCS
 ***************************************************************************/


static int
CombineLagRates
   (TDate   startDate,
    TDate   amortStartDate, /* (I) start date of rate array */
    long    numAmortRates,  /* (I) # rates in array */
    double *adjAmortRates,  /* (I) array of 2-pt commitment rates */
    long    amortType,      /* (I) type of FHLMC commitment rate used;
                             * should be either MBS_PP_REFI_TYPE_FH30 
                             * or MBS_PP_REFI_TYPE_FH15 */
    long    histStartDateInYYYYMMDD, /* (I) starting month (day-of-mon ignored)
                                      * for adjHistAmortIndexRates[] rates */
    long    numHistRates,   /* (I) # of rates in adj_commit[] */
    double *adjHistRates,   /* (I) array of average monthly
                             * 2-point commitment rates (assumed
                             * same type as adjAmortIndexRates[]) */
    double *wtLags,
    TDate  *userEarliestCommitDate,
    double *grFrmLagRates)  /* (O) array of lagged rates */
{
    int amortStartIndex;
    TDate histStartDate;

    TDateOf(histStartDateInYYYYMMDD, &histStartDate);

    long numMonthsHistToStart;
    MonsDiff (histStartDate, amortStartDate, &numMonthsHistToStart);

    if (numMonthsHistToStart > 0) {
	amortStartIndex = numMonthsHistToStart;
	*userEarliestCommitDate = histStartDate;
	FillInHistRates (numHistRates, adjHistRates, (int) numMonthsHistToStart, grFrmLagRates);
    }
    else {
	amortStartIndex = 0;
	*userEarliestCommitDate = amortStartDate;
    }

    int numRates = amortStartIndex + numAmortRates;

    for (int i = 0; i < numAmortRates; i++)
	grFrmLagRates [i+amortStartIndex] = adjAmortRates[i];

   // int j;
   // for ( j=0; j<numRates;j++) 
   //     std::cerr << "j " << j << " rate " << grFrmLagRates[j]<<std::endl;

    PerformLagging (numRates, grFrmLagRates, wtLags);

//    std::cerr <<"post lag\n";
 //   for ( j=0; j<numRates;j++) 
 //       std::cerr << "j " << j << " rate " << grFrmLagRates[j]<<std::endl;

    return SUCCESS;
}

void FillInHistRates (int numRates, double* adjHistRates, int numRatesToFill, double* fillRates)
{
    int i;
    for (i = 0; i < numRates && i <numRatesToFill; i++)
	fillRates[i] = adjHistRates[i];

    for (i = numRates; i < numRatesToFill; i++)
	fillRates[i] = fillRates[numRates - 1];
}

void PerformLagging (int numRates, double* grFrmLagRates, double *wtLags)
{
    /* calculate lag rates from commintment rates */
    for (int i = numRates-1; i >= 0; i--)
    {
        grFrmLagRates[i] = 0.;
        if (i >= DEF_MBS_PP_MGRP_NUM_WTLAGS)
        {
            for (int j = 0; j < DEF_MBS_PP_MGRP_NUM_WTLAGS; j++)
            {
                grFrmLagRates[i] += grFrmLagRates[i-1-j] * wtLags[j];
            }
        }    
        else
        {   
            grFrmLagRates[i] = grFrmLagRates[i+1];
        }
    }

}


/**********************************************************************
 *  Routine to transform lagged rates (for given WAC) to reflect a
 *  shift in the incentive required for callability.
 *  INPUTS to TweakInc:
 *      startMthIx   : Starting month number in array (for grFrmLagRates).
 *      endMthIx     : Ending   month number in array (for grFrmLagRates).
 *      twkCallInc   : Amount by which to tweak incentive (as fraction of
 *                     WAC;  e.g., .01 means calls delayed by 1% of WAC).
 *      critRat      : Critical ratio at which callability is approximated
 *                     to begin;  Shift refinance incentive according to
 *                     critRat.
 *      wacs         : WACs for given security.
 *      grFrmLagRates: Rates for use in prepayment analysis (they should
 *                     already be lagged).
 *
 *  OUTPUT of TweakInc : Transformed grFrmLagRates that have the effect of
 *                       shifting callability (see twkCallInc input).
 **********************************************************************/

static int
TweakInc
   (long    startMthIx,
    long    endMthIx,
    double  twkCallInc,
    double  critRat,
    double *wacs,
    double *grFrmLagRates,
    double *grFrmLagRatTran)
{
    static char routine[] = "TweakInc";
    int status = FAILURE;
    long iCount;
    double thRatio;
    double thFract;
    double xJunk;

    for (iCount = 0; iCount <= endMthIx; iCount++)
    {
        xJunk = grFrmLagRates[iCount];

        /* Original incentive ratio. */
        thRatio = wacs[iCount] / grFrmLagRates[iCount];

        if (thRatio >= critRat + twkCallInc)
        {
            /* Beyond critRat incentive, use the full shift in incentive. */
            thFract = 1;
        }
        else
        {
            if (thRatio > 1)
            {
                /* Intermed (btn c.c. and critRat) shifted linearly. */
                thFract = (thRatio - 1) / (critRat + twkCallInc - 1);
            }
            else
            {
                /* Discounts do not experience a prepayment shift. */
                thFract = 0;
            }
        }
        thRatio = thRatio - twkCallInc * thFract;

        /* Translate new incentive ratio into market interest
         * rate (for refi).
         */
        grFrmLagRatTran[iCount] = wacs[iCount] / thRatio;
    }
    status = SUCCESS;

    return(status);
}



/*******************************************************
 *  SetDefaultMgrpTweaks()
 *  Determines the default (hard-coded) MGrp prepay tweaks
 *  (4 values) as a function of the type of collateral
 *******************************************************/
static
int SetDefaultMgrpTweaks
   (long    mbsAgency,       /* (I) GNMA, etc. */
    long    mbsTerm,         /* (I) MBS_PP_MBSTERM_30, etc */
    double  turnoverTwk,     /* (I) user defined */
    double  incentPtTwk,     /* (I) user defined */
    double  totSteepTwk,     /* (I) user defined */
    double  baseSteepTwk,    /* (I) user defined */
    double  ageRampTwk,      /* (I) user defined */
    double *twkTurn,         /* (O) */
    double *twkCallInc,      /* (O) */
    double *multBurnedSteep, /* (O) */
    double *multTrigSteep,   /* (O) */
    double *twkAgeRamp)      /* (O) */
{
    static char routine[] = "SetDefaultMgrpTweaks";
    int status = FAILURE;

    /* For terms <= 180, treat as 15yr */
    if( mbsTerm IS MBS_PP_MBSTERM_15 )
    {
        switch (mbsAgency)
        {
        case MBS_AGENCY_FNMA:
            (*twkTurn) = 0.;
            (*twkCallInc) = 0.;
            (*multBurnedSteep) = 1.;
            (*multTrigSteep) = 1.;
            (*twkAgeRamp) = 0.;
            break;
        case MBS_AGENCY_FHLMC:
            (*twkTurn) = 0.;
            (*twkCallInc) = 0.;
            (*multBurnedSteep) = 1.;
            (*multTrigSteep) = 1.;
            (*twkAgeRamp) = 0.;
            break;
        case MBS_AGENCY_GOLD:
            (*twkTurn) = 0.;
            (*twkCallInc) = 0.;
            (*multBurnedSteep) = 1.;
            (*multTrigSteep) = 1.;
            (*twkAgeRamp) = 0.;
            break;
        case MBS_AGENCY_GNMAI:
            (*twkTurn) = 0.;
            (*twkCallInc) = 0.;
            (*multBurnedSteep) = 1.;
            (*multTrigSteep) = 1.;
            (*twkAgeRamp) = 0.;
            break;
        case MBS_AGENCY_GNMAII:
            (*twkTurn) = 0.;
            (*twkCallInc) = 0.;
            (*multBurnedSteep) = 1.;
            (*multTrigSteep) = 1.;
            (*twkAgeRamp) = 0.;
            break;
        case MBS_AGENCY_WHOLE:
            (*twkTurn) = 0.;
            (*twkCallInc) = 0.;
            (*multBurnedSteep) = 1.7;
            (*multTrigSteep) = 1.7;
            (*twkAgeRamp) = 0.;
            break;
        default:
            GtoErrMsg("%s: Invalid mbsAgency: %ld\n",routine,mbsAgency);
            goto done;
        }
    }
    /* Assume all other terms to be 30yr */
    else
    {
        switch (mbsAgency)
        {
        case MBS_AGENCY_FNMA:
            (*twkTurn) = 0.;
            (*twkCallInc) = 0.;
            (*multBurnedSteep) = 1.;
            (*multTrigSteep) = 1.;
            (*twkAgeRamp) = 0.;
            break;
        case MBS_AGENCY_FHLMC:
            (*twkTurn) = 0.;
            (*twkCallInc) = 0.;
            (*multBurnedSteep) = 1.;
            (*multTrigSteep) = 1.;
            (*twkAgeRamp) = 0.;
            break;
        case MBS_AGENCY_GOLD:
            (*twkTurn) = 0.;
            (*twkCallInc) = 0.;
            (*multBurnedSteep) = 1.;
            (*multTrigSteep) = 1.;
            (*twkAgeRamp) = 0.;
            break;
        case MBS_AGENCY_GNMAI:
            (*twkTurn) = 0.;
            (*twkCallInc) = 0.;
            (*multBurnedSteep) = 1.;
            (*multTrigSteep) = 1.;
            (*twkAgeRamp) = 0.;
            break;
        case MBS_AGENCY_GNMAII:
            (*twkTurn) = 0.;
            (*twkCallInc) = 0.;
            (*multBurnedSteep) = 1.;
            (*multTrigSteep) = 1.;
            (*twkAgeRamp) = 0.;
            break;
        case MBS_AGENCY_WHOLE:
            (*twkTurn) = 0.;
            (*twkCallInc) = 0.;
            (*multBurnedSteep) = 1.7;
            (*multTrigSteep) = 1.7;
            (*twkAgeRamp) = 0.;
            break;
        default:
            GtoErrMsg("%s: Invalid mbsAgency: %ld\n",routine,mbsAgency);
            goto done;
        }
    }

    /* combine user-supplied tweaks with default tweaks
     * for this type of collateral */
    (*twkTurn) = (1+turnoverTwk)*(1+(*twkTurn))-1.;
    (*twkCallInc) += incentPtTwk;
    (*multBurnedSteep) *= totSteepTwk;
    (*multTrigSteep) *= baseSteepTwk;
    (*twkAgeRamp) += ageRampTwk;

    status = SUCCESS;

done:
    return status;
}

/***************************************************************************
 *  PRIVATE FUNCTION
 *  Compute the ratio after adjusting for a roll down the mortgage curve and
 *  a long term aging effect. The "roll" captures the fact that 15-year
 *  product should have a lower market coupon than 30-year product (typically
 *  about 50 b.p. in recent history) and extends this relationship for any
 *  given term by an average life interpolation.
 *  The long term aging allows for older product to be less sensitive to the
 *  ratio than newer product; this is modeled by always applying a minimum
 *  incentive multiplier with the remainder subject to exponential decay with
 *  age.
 *
 *  Parameters used in adjusting the refinance ratio. Values assumed or
 *  derived in model estimation were:
 *  baseCPR = .06 (assumed)
 *  ratio15Year = 1.0625 (assumed)
 *  refWAC = .085 (assumed)
 *  premMinMult = .01 (fit) - for premium
 *  discMinMult = .99 (fit) - for discount
 *  decayRateAnn = .14 (fit)
 ***************************************************************************/
static
int AdjRatio
   (double  ratio,        /* (I) Standard incentive measure (should be
                           * consistent with the original term for the given
                           * product). */
    double  ageInMths,    /* (I) Age in months */
    double  origTermMths, /* (I) Original term in months (usually should use
                           * exactly 360 or 180). */
    double  baseCPR,      /* (I) */
    double  ratio15Year,  /* (I) */
    double  refWAC,       /* (I) */
    double  premMinMult,  /* (I) */
    double  discMinMult,  /* (I) */
    double  decayRateAnn, /* (I) */
    double *adjRatio)     /* (O) */
{
    static char routine[] = "AdjRatio";
    int status = FAILURE;
    double thAL;
    double intermedRatio;
    double origAL;
    double multUse;
    static double avgLife[361];

    /* First time the array of average life need to be pre-computed
     * or the parameters has been changed.
     * It will save a lots of computation time when call AvgLife()
     * function many times.   Davis Lee   03/02/2000
     */
    static TMbsAvgLifeVec avgLifeVec;
    if (avgLifeVec.NeedToReBuild(baseCPR, refWAC)) {
        avgLifeVec.Build(baseCPR, refWAC);
    }

    /* First adjust the incentive for the roll down the MBS curve
     * (lower coupon for shorter term).
     *  Note that at present this adjustment is quite approximate: might
     * refine this in the future.
     */
    thAL = avgLifeVec.GetAvgLife(origTermMths-ageInMths);
//    AvgLife((origTermMths-ageInMths)/12., baseCPR, refWAC, &thAL);

    static 	double AL30=0;
    if (AL30==0)
    {
        AL30 = avgLifeVec.GetAvgLife(360);
//	AvgLife(30, baseCPR, refWAC, &AL30);
    }

    static 	double AL15=0;
    if (AL15==0)
    {
        AL15 = avgLifeVec.GetAvgLife(180);
//	AvgLife(15, baseCPR, refWAC, &AL15);
    }


    if (origTermMths IS 360)
    {
        /* Usually this is the case; save time by not recalculating
         * it needlessly.
         */
        origAL = AL30;
    }
    else
    {
        if (origTermMths IS 180)
        {
            origAL = AL15;
        }
        else
        {
            origAL = avgLifeVec.GetAvgLife(origTermMths);
//            AvgLife(origTermMths/12., baseCPR, refWAC, &origAL);
        }
    }
 
    /* Approximation: assumes a linear relationship for new loans between
     * AL at baseCPR and market gross cpn.
     * Also, rather than using true WAC, for now we use a single reference
     * WAC for all cases. Thus, this is in reality an aging adjustment that
     * behaves the same as a roll down MBS curve in some circumstances.
     */
    /* For 15 years model we simply approximate by the the same slope as
     * 30 years model (extrapolation) */
    intermedRatio = ratio * (1 + (ratio15Year-1) * (origAL-thAL) / (AL30-AL15));
 
    /* Also adjust for the long term aging effect (portion given by minMult
     * is always present; remainder of incentive suffers from exponential
     * decay with age).
     */

    if (intermedRatio >= 1)
    {
        multUse = pow(premMinMult+(1-premMinMult)/(1+decayRateAnn),
            (ageInMths/12.));
        *adjRatio = 1 + (intermedRatio - 1) * multUse;
    }
    else
    {
        multUse = pow(discMinMult+(1-discMinMult)/(1+decayRateAnn),
            (ageInMths/12.));
        *adjRatio = pow(intermedRatio,multUse);
    }

    status = SUCCESS;

//done:
    return status;
}


/***************************************************************************
 *  PRIVATE FUNCTION
 *  Routine to calculate average life given a continuous CPR at prepay and
 *  continuous coupon at coupon, with final maturity given by final
 *  (scheduled amortization is included).  All units must be consistent
 *  (e.g., final in years and annualized contin prepays and coupons).
 ***************************************************************************/
static
int AvgLife
   (double  finalMat,  /* (I) final maturity in months */
    double  speed,     /* (I) a continuous CPR */
    double  coupon,    /* (I) continuous coupon */
    double *avgLife)   /* (O) average life */
{
    static char routine[] = "AvgLife";
    int status = FAILURE;
    double thExp;
    double const0;
    double const1;
    double integral1;
    double integral2;
   
    if (finalMat <= 0)
    {
        /* No life left in this case. */
        *avgLife = 0;
        goto done;
    }

    thExp = exp(coupon * finalMat);

    const0 = 1. / (1 - 1 / thExp);
    const1 = 1. / (thExp - 1);

    IntegralTexpCT(finalMat, coupon-speed, &integral1);
    IntegralTexpCT(finalMat, -speed, &integral2);
    *avgLife = const1 * (coupon - speed) * integral1 +
               const0 * speed * integral2;

    status = SUCCESS;

done:
    return status;
}

/***************************************************************************
 *  PRIVATE FUNCTION
 *  Routine to calculate the integral of t*exp(ct) from zero to "final".
 ***************************************************************************/
static
int IntegralTexpCT
   (double  final,
    double  c,
    double *integralTexpCT)
{
    static char routine[] = "IntegralTexpCT";
    int status = FAILURE;
    double thExp;

    thExp = exp(c * final);

    *integralTexpCT = (final * thExp - (thExp - 1) / c) / c;

    status = SUCCESS;

//done:
    return status;
}

