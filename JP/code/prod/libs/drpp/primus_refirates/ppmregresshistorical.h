/* ////////////////////////////////////////////////////////////////////////
// J.P. Morgan Securities Inc.		     Fixed Income Market Strategies
//					     Extracted: %H% @ %T%
// Module:     %M%
// Version:    %I%
// Updated:    %G% @ %U%
// Author:     Paul Honig
// Modified:   Leah G Brown
//
// Description: structure definitions for historical rates for use by
//              prepayment models
//
//////////////////////////////////////////////////////////////////////// */

#ifndef __ppm_regress_historical_h
#define __ppm_regress_historical_h

#include "enumeration.h"
#include "extern.h"

/* For blended mortgage spreads, these macros define the fraction of spread
   off of the one year relative to the ten year treasury to use.
   These macros are currently used by PpmRegress.c and PpmRegressHistWrap.C */
#define MSPRD_15YR_TSY1YR_FRAC 0.464
#define MSPRD_30YR_TSY1YR_FRAC 0.215

/* A description of historical rates */
typedef struct
{
    int    nRates;             /* number of monthly 30yr rates and 30-15yr spreads */
    int    firstDateYear;      /* year of date of first rate in the arrays */
    int    firstDateMonth;     /* month of date of first rate in the arrays */
    int    firstDateDay;       /* day of date of first rate in the arrays */
    double *mtgRates30yr;      /* array of 30yr mortgage rates */
    double *mtgRates15yr;      /* array of 15yr mortgage rates */
    double *mtgRates30yrPts;   /* array of points on 30yr mortgage rates */
    double *mtgRates15yrPts;   /* array of points on 15yr mortgage rates */
    double *tsy10yrRates;      /* array of tsy 10yr rates */
    double *tsy1yrRates;       /* array of tsy 1yr rates */
    double *ssprd10yr;         /* array of 10yr swap spread */
    double *lib1yrRates;       /* array of cms 1yr rates */
    
    double fn30yrServicing;    /* 30 yr FNMA CC / FHLMC Wkly Commitment Spread */
    double fn15yrServicing;    /* 15 yr FNMA CC / FHLMC Wkly Commitment Spread */
    double gn30yrServicing;    /* 30 yr GNMA CC / FHLMC Wkly Commitment Spread */
    double gn15yrServicing;    /* 15 yr GNMA CC / FHLMC Wkly Commitment Spread */
    double fh30yrServicing;    /* 30 yr FHLMC CC / FHLMC Wkly Commitment Spread */
    double fh15yrServicing;    /* 15 yr FHLMC CC / FHLMC Wkly Commitment Spread */
    double fn30yrNetCurrCpn;   /* 30 yr FNMA CC */
    double fn15yrNetCurrCpn;   /* 15 yr FNMA CC */
    double gn30yrNetCurrCpn;   /* 30 yr GNMA CC */
    double gn15yrNetCurrCpn;   /* 15 yr GNMA CC */
    double fh30yrNetCurrCpn;   /* 30 yr FHLMC CC */
    double fh15yrNetCurrCpn;   /* 15 yr FHLMC CC */
    double tsy1yrForSpread;    /* one year Tsy rate from current environment */
    double tsy10yrForSpread;   /* ten year Tsy rate from current environment */
    double lib1yrForSpread;    /* one year Cms rate from current environment */
    double ssprd10yrForSpread; /* ten year swap spread from current environment */
    double cmtHump1yrForSpread;   /* 10yr CMT1yr hump from current environment */
    double cmtHump10yrForSpread;  /* 10yr CMT10yr hump from current environment */
    double ARMFixedCCSpread;   /* ARM / Fixed Rate Current Coupon Spread */
    double latestMtgRate30yrPts;/* 30 Yr Points going forward */
    double latestMtgRate15yrPts;/* 15 Yr Points going forward */

    double mtgSpreadAdditiveFrac;         /* additive mtg spread factor. */
    MtgSpreadType  mtgSpreadType;   /* TSY_MTG_SPREAD or SWAP_MTG_SPREAD .*/


    unsigned int isShallowCopy:1;      /* Is this a shallow copy? */
    unsigned int isLastHistRateIncomplete:1; /* Is the last rate in the    */
					     /* historical array for month */
					     /* without a complete month's */
					     /* worth of fhlmc commitment data */

    int partialMonthDateYear;  /* year of last date in a partial data month for which */
			       /* historical information is valid */
    int partialMonthDateMonth; /* month of last date in a partial data month for which */
			       /* historical information is valid */
    int partialMonthDateDay;   /* day of last date in a partial data month for which */
			       /* historical information is valid */

    double  offRunMtgSpread15YrAdj;
    double  offRunMtgSpread30YrAdj;
} PpmRegressHistoricalRates;

/***************************************************************************
 * Constructor function for PpmRegressHistoricalRates structures.
 ***************************************************************************/
EXTERNC PpmRegressHistoricalRates* PpmRegressHistRatesNew();

/***************************************************************************
 * Function to create a shallow copy of a PpmRegressHistoricalRates structure.
 ***************************************************************************/
EXTERNC PpmRegressHistoricalRates* PpmRegressHistRatesShallowCopy(
    const PpmRegressHistoricalRates*);

/***************************************************************************
 * Destructor function for PpmRegressHistoricalRates structures.
 * NOTE: This handles destruction of shallow copies as well.
 ***************************************************************************/
EXTERNC void PpmRegressHistRatesFree(PpmRegressHistoricalRates*);

EXTERNC int
ppmRegressReadHistFile(
    char *userPathDaily,
    char *userPathWeekly,
    PpmRegressHistoricalRates *histRates,
    int asOfDateYear, int asOfDateMonth, int asOfDateDay);

#define XONFAILURE(x, y) { 		\
	if ((x) == PRIMUS_FAILURE) {	  	\
		goto ERROR_EXIT;  	\
	}			  	\
}

#define XONFALSE2(x, y) { 		\
        if (!(x)) {	         	\
		goto ERROR_EXIT;  	\
	}			  	\
}

#endif
