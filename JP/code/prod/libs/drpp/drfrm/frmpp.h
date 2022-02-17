/***************************************************************************
 *	SCCS Keyword Information
 *	------------------------
 *	Module name	:  frmpp.h
 *	Company Name	:  JP Morgan Securities Inc.
 *	Authors  	:  Forrest Quinn
 *			   Davis (Shuenn-Tyan) Lee
 *			   (Derivatives Research)
 *	Code version    :  1.4
 *	Extracted	:  5/19/97 at 13:54:15
 *	Last Updated	:  5/19/97 at 13:54:00
 ***************************************************************************
 *      Public functions of GNMA II ARM prepay model (author: Forrest Quinn)
 *      Model version:  3.1.4
 *
 *      Copyright 1996 J.P. Morgan & Co. Incorporated. All rights reserved.
 ***************************************************************************/
#ifndef __frmpp_h
#define __frmpp_h

#include "drstruct.h"
typedef DArray TMbsDouble;
typedef DMatrix TMbsDouble2;

extern "C" {
#include "cdate.h"
#include "bastypes.h"
}
#include "ppconst.h"
#include "mbsconst.h"
#include "drcache.h"

class DRSmmVecMatrixInput : public DRCacheInput {
protected:
	DRSmmVecMatrixInput () : DRCacheInput("SmmVecMatrix") {}
	DRSmmVecMatrixInput (
		long startMthIx,
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
		double multBurnedSteep);
	
	friend class DRSmmVecMatrixItem;
	friend class DRSmmVecMatrix;
	virtual DRCacheItem* create ();
	virtual DRCacheInput* copy ();
	virtual void Print(ostream& s) const;
	virtual bool ok_to_use (DRCacheInput& a);

	long m_startMthIx;
	long m_endMthIx;
	DArray m_trueWacs;
	DArray m_effWacs;
	long m_effOrigTerm;
	DArray m_wams;
	long m_ix0Moy;
	long m_inclSchedAmort;
        double m_twkAgeRamp;
	DArray m_seasonality;
	DArray m_seasLogist;
	DArray m_scurveLogist;
	DArray m_ltAges;
	double m_wac;
	double m_critRat;
	double m_multBurnedSteep;
};

class DRSmmVecMatrixItem : public DRCacheItem {
protected:
	DRSmmVecMatrixItem(DRSmmVecMatrixInput&);
	friend class DRSmmVecMatrix;
	friend class DRSmmVecMatrixInput;
	virtual int bytes();
	DMatrix m_matrix;
	DArray m_ratePts;	
	friend ostream& operator<<(ostream&, const DRSmmVecMatrixItem&);
};

class DRSmmVecMatrix {
public:
    DRSmmVecMatrix(long    startMthIx,    /* (I) Starting month for which
                                           * smm's are to be generated */
                   long    endMthIx,      /* (I) Last month for which
                                           * smm's are to be generated */
                   DArray& trueWacs,      /* (I) Array of gross coupons,
                                           * starting with that which applies
                                           * to prepays in month startMthIx */
                   DArray& effWacs,       /* (I) effective WACs (possibly
                                           * shifted), starting with that
                                           * which applies to prepays in
                                           * month startMthIx */
                   long    effOrigTerm,   /* (I) WARM+WALA */
                   DArray& wams,          /* (I) Array of monthly WAM's (not
                                           * used for group effects);  start
                                           * with value for prepays in month
                                           * startMthIx; wam and wala are valid
                                           * as of the the month prior to
                                           * startMthIx (should be in line
                                           * with standard quotes) */
                   long    ix0Moy,        /* (I) onth of year for month indexed
                                           * by zero (integer between 1 and 12)
                                           */
                   long    inclSchedAmort,/* (I) zero means ignore scheduled
                                           * amortization;  1 means include it;
                                           * inclSchedAmort indicates whether
                                           * to add scheduled amortization to
                                           * the prepayment rates */
                   double  twkAgeRamp,    /* (I) tweak age ramp (months) */
                   DArray& seasonality,	  /* (I) Seasonal multiplier for
                                           * turnover component */
                   DArray& seasLogist,    /* (I) Array of 4 parameters to get
                                           * age of seasoning as a logistic */
                   DArray& scurveLogist,  /* (I) array of 4 s-curve logistic
                                           * parameters to use (may have been
                                           * altered by prepay tweaks) */
                   DArray& ltAges,        /* (I) array of six longer term
                                           * aging parameters to use adjust
                                           * refinance ratio */
                   double  wac,
                   double  critRat,
                   double  multBurnedSteep);

	double GetSmmVector (double rate, int month);
	friend ostream& operator<<(ostream&, const DRSmmVecMatrix&);

protected:
	DRCacheItemPtr m_cacheItemPtr;
	DRSmmVecMatrixItem* m_item;
	DRSmmVecMatrixInput m_input;
};

class TMbsSmmVecMatrix {
	TMbsDouble2		smmVecMatrix;
	TMbsDouble		RatesGrid;
	double oldwac, oldcritRat, oldmultBurnedSteep, oldtwkAgeRamp;

	long oldStartMthIx;
	long oldEndMthIx;
	long oldEffOrigTerm;
	long oldix0Moy;
	long oldinclSchedAmort;
	TMbsDouble oldTrueWacs;
	TMbsDouble oldEffWacs;
	TMbsDouble oldWams;
	TMbsDouble oldSeasonality;
	TMbsDouble oldSeasLogist;
	TMbsDouble oldSCurveLogist;
	TMbsDouble oldltAges;
public:
	bool NeedToReBuild(
                 long    startMthIx,	    /* (I) Starting month for which
                                             * smm's are to be generated */
                 long    endMthIx,          /* (I) Last month for which smm's
                                             * are to be generated */
                 double* trueWacs,          /* (I) Array of gross coupons,
                                             * starting with that which
                                             * applies to prepays in month
                                             * startMthIx */
                 double* effWacs,           /* (I) effective WACs (possibly
                                             * shifted), starting with that
                                             * which applies to prepays in
                                             * month startMthIx */
                 long    effOrigTerm,       /* (I) WARM+WALA */
                 double* wams,              /* (I) Array of monthly WAM's
                                             * (not used for group effects);
                                             * start with value for prepays
                                             * in month startMthIx; wam and
                                             * wala are valid as of the the
                                             * month prior to startMthIx
                                             * (should be in line with
                                             * standard quotes) */
                 long    ix0Moy,            /* (I) onth of year for month
                                             * indexed by zero (integer
                                             * between 1 and 12) */
                 long    inclSchedAmort,    /* (I) zero means ignore scheduled
                                             * amortization; 1 means include
                                             * it;  inclSchedAmort indicates
                                             * whether to add scheduled
                                             * amortization to the prepayment
                                             * rates */
                 double* seasonality,       /* (I) Seasonal multiplier for
                                             * turnover component */
                 double* seasLogist,        /* (I) Array of 4 parameters to
                                             * get age of seasoning as a
                                             * logistic */
                 double* scurveLogist,      /* (I) array of 4 s-curve logistic
                                             * parameters to use (may have been
                                             * altered by prepay tweaks) */
                 double* ltAges,            /* (I) array of six longer term
                                             * aging parameters to use adjust
                                             * refinance ratio */
                 double  wac,
                 double  critRat,
                 double  multBurnedSteep,
                 double  twkAgeRamp);

void Build(double  Upper_Rate,
           double  Lower_Rate,
           long    startMthIx,	    /* (I) Starting month for which smm's
                                     * are to be generated */
           long    endMthIx,        /* (I) Last month for which smm's are
                                     * to be generated */
           double* trueWacs,        /* (I) Array of gross coupons, starting
                                     * with that which applies to prepays in
                                     * month startMthIx */
           double* effWacs,         /* (I) effective WACs (possibly shifted),
                                     * starting with that which applies to
                                     * prepays in month startMthIx */
           long    effOrigTerm,     /* (I) WARM+WALA */
           double* wams,            /* (I) Array of monthly WAM's (not used
                                     * for group effects);  start with value
                                     * for prepays in month startMthIx; wam
                                     * and wala are valid as of the the month
                                     * prior to startMthIx (should be in line
                                     * with standard quotes) */
           long    ix0Moy,          /* (I) onth of year for month indexed by
                                     * zero (integer between 1 and 12) */
           long    inclSchedAmort,  /* (I) zero means ignore scheduled
                                     * amortization; 1 means include it;
                                     * inclSchedAmort indicates whether to
                                     * add scheduled amortization to the
                                     * prepayment rates */
           double* seasonality,     /* (I) Seasonal multiplier for turnover
                                     * component */
           double* seasLogist,      /* (I) Array of 4 parameters to get age
                                     * of seasoning as a logistic */
           double* scurveLogist,    /* (I) array of 4 s-curve logistic
                                     * parameters to use (may have been
                                     * altered by prepay tweaks) */
           double* ltAges,          /* (I) array of six longer term aging
                                     * parameters to use adjust refinance
                                     * ratio */
           double  wac,
           double  critRat,
           double  multBurnedSteep,
           double  twkAgeRamp);
													
	double GetSmmVector (double rate, int month);
	friend ostream& operator<<(ostream&, const TMbsSmmVecMatrix&);
};


class TMbsAvgLifeVec {
	TMbsDouble  avgLife;
	double      oldbaseCPR, oldrefWAC;

public:
	bool NeedToReBuild(
                   double  baseCPR,    /* (I) a continuous CPR */
                   double  refWAC);    /* (I) continuous coupon */

        void Build(double  baseCPR,    /* (I) a continuous CPR */
                   double  refWAC);    /* (I) continuous coupon */

	double GetAvgLife (int month);
};


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
    double *smmArr);    /* (O) UTPUT of smmFixed: Fills smmFixed with array
                         * of smm values; these combine the standard S and
                         * aging-curves with multi-group effects. The first
                         * entry to smmVect corresponds to startMthIx. */


int
HistSmmFix(
    double  grCpnUse,
    double  grCpnForSchAm,
    double  age,
    double  longIxLagged,
    long     moy,
    double  grFrmSprd,
    long     origTermInMths,
    long     inclSchedAmort,
    double *seasonality,
    double *seasLogist,
    double *logisticPs,
    double *histSmm_Fix);


int 
SmmVector(
    long    startMthIx,     /* (I) Starting month for which smm's are to be
                               generated */
    long    endMthIx,	    /* (I) Last month for which smm's are to be
                               generated */
    double *trueWacs,       /* (I) Array of gross coupons, starting with that
                               which applies to prepays in month startMthIx */
    double *effWacs,        /* (I) effective WACs (possibly shifted), starting
                             * with that which applies to prepays in 
                             * month startMthIx */
    long    effOrigTerm,    /* (I) WARM+WALA */
    double *wams,           /* (I) Array of monthly WAM's (not used for group
                               effects);  start with value for prepays in month
                               startMthIx; wam and wala are valid as of the the
                               month prior to startMthIx (should be in line
                               with standard quotes) */
    double *grFrmLagRates,  /* (I) Array of mortgage commitment rates (lagged);
                               grFrmLagRates is the array of market refi
                               interest rates (appropriately lagged) */
    long    ix0Moy,         /* (I) onth of year for month indexed by zero
                               (integer between 1 and 12) */
    long    inclSchedAmort, /* (I) zero means ignore scheduled amortization;
                               1 means include it;  inclSchedAmort indicates
                               whether to add scheduled amortization to the
                               prepayment rates */
    double  twkAgeRamp,     /* (I) tweak age ramp (months) */
    double *seasonality,    /* (I) Seasonal multiplier for turnover component */
    double *seasLogist,     /* (I) Array of 4 parameters to get age of
                               seasoning as a logistic */
    double *scurveLogist,   /* (I) array of 4 s-curve logistic parameters to
                               use (may have been altered by prepay tweaks) */
    double *ltAges,         /* (I) array of six longer term aging parameters to
                               use adjust refinance ratio */
    double *smmVect);       /* (O) Returns the model's speed in SMM for the
                               desired single month */	  
static int 
AvgLife
   (double  finalMat,  /* (I) final maturity in months */
    double  speed,     /* (I) a continuous CPR */
    double  coupon,    /* (I) continuous coupon */
    double *avgLife);  /* (O) average life */

#endif
