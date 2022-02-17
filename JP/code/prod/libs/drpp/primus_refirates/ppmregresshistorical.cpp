/* ////////////////////////////////////////////////////////////////////////
// J.P. Morgan Securities Inc.		     U.S. Fixed Income Analytics
//					     Extracted: %H% @ %T%
// Module:     %M%
// Version:    %I%
// Updated:    %G% @ %U%
// Author:     Jerry Cohen
//
// Description:
//   Functions for manipulating PpmRegressHistoricalRates structures.
//
////////////////////////////////////////////////////////////////////////*/

/*Copyright 1996, Morgan Guaranty Trust Company of New York.  All rights reserved*/

#include "ppmregresshistorical.h"
extern "C" {
#include "ldate.h"
#include "convert.h"
}
#include <string.h>

/*
ppmRegressBuildFilePath() : build a file path to find flat data files
in the appropriate directory
*/
static int
ppmRegressBuildFilePath(
    int  pathLength,
    char *path,
    char *fileName)
{
    static char *overridePath;
    static char *dataPath;
    static int  gotEnv = FALSE;
    char        *thePath, *latest = "";

    if (!gotEnv)
    {
	overridePath = getenv("PPMREGRESSPATH");
	dataPath = getenv("DATA");
	gotEnv = TRUE;
    }
    if (overridePath)
	thePath = overridePath;
    else if (dataPath)
    {
	thePath = dataPath;
	latest = "LATEST/";
    }
    else
	XONFALSE2(FALSE,
		  "ppmRegressBuildFilePath():No path for regression prepayment model data files");
    
    XONFALSE2((int)strlen(thePath) + (int)strlen(fileName) + (int)strlen(latest) + 1 < pathLength,
	      "ppmRegressBuildFilePath() : Path too long");
    strcpy(path, thePath);
    strcat(path, "/");
    strcat(path, latest);
    strcat(path, fileName);
    
    return SUCCESS;
ERROR_EXIT:
    return PRIMUS_FAILURE;
}

/***************************************************************************
 *
 * PpmRegressHistoricalRates constructor function.
 *
 ***************************************************************************/
PpmRegressHistoricalRates* PpmRegressHistRatesNew()
{
    PpmRegressHistoricalRates* objp =
	(PpmRegressHistoricalRates*)malloc(sizeof(PpmRegressHistoricalRates));

    if (objp != NULL)
    {
	objp->isShallowCopy = 0;
	objp->mtgRates30yr = NULL;
	objp->mtgRates15yr = NULL;
	objp->mtgRates30yrPts = NULL;
	objp->mtgRates15yrPts = NULL;
	objp->tsy10yrRates = NULL;
	objp->tsy1yrRates = NULL;
	objp->ssprd10yr = NULL;
	objp->lib1yrRates = NULL;

	objp->mtgSpreadAdditiveFrac = 1.0;
	objp->mtgSpreadType = TSY_MTG_SPREAD;
	objp->offRunMtgSpread15YrAdj = 0.;
	objp->offRunMtgSpread30YrAdj = 0.;
    }
    
    return objp;
}

/***************************************************************************
 *
 * Function to make a shallow copy of a PpmRegressHistoricalRates
 * object.  The pointer (i.e. array) elements of the source
 * object are not deep-copied so we mark the copy as shallow to prevent
 * the destructor function from freeing memory not owned by the copy.
 *
 ***************************************************************************/
PpmRegressHistoricalRates* PpmRegressHistRatesShallowCopy(
    const PpmRegressHistoricalRates* srcp)
{
    PpmRegressHistoricalRates* copyp = NULL;

    if (srcp != NULL &&	(copyp = PpmRegressHistRatesNew()) != NULL)
    {
	memcpy((void*)copyp, (void*)srcp,
	       sizeof(PpmRegressHistoricalRates));
	copyp->isShallowCopy = 1;
    }

    return copyp;
}

/***************************************************************************
 *
 * Destructor function for PpmRegressHistoricalRates objects which is
 * smart enough to distinguish between normal objects and shallow
 * copies.
 *
 ***************************************************************************/
void PpmRegressHistRatesFree(PpmRegressHistoricalRates* objp)
{
    if (objp != NULL)
    {
	if (!objp->isShallowCopy)
	{
	    free(objp->mtgRates30yr);
	    free(objp->mtgRates15yr);
	    free(objp->mtgRates30yrPts);
	    free(objp->mtgRates15yrPts);
	    free(objp->tsy10yrRates);
	    free(objp->tsy1yrRates);
	    free(objp->ssprd10yr);
	    free(objp->lib1yrRates);
	}
	free(objp);
    }
}

/*
completePeriod() : determines whether there is a complete monthly period in
the last month of the historical rates file.  If not, the previous
month becomes the last entry and rate projections are used going
forward.
*/
static int
TDateCompare(TDate newDate, TDate endDate)
{
    long days = 0;
    GtoDaysDiff(newDate, endDate, GTO_ACT_365, &days);
    return -days;
}

static void
addWeek(TDate startDate, TDate *endDate)
{
     TDateInterval  weekInterval;
 
     GtoMakeDateInterval(7, 'D', &weekInterval);
     GtoDtFwdAny(startDate, &weekInterval, endDate);
}

static int
completePeriod(TDate endDate, TDate rateDate)
{
    int retval = TRUE;
    TDate newDate;

    
	addWeek(rateDate, &newDate);
    if (TDateCompare(newDate, endDate) < 0)
	retval = FALSE;
 
    return retval;
}

/*
ppmRegressReadHistCCFile() : read in the historical current coupon
file.  This file, updated daily, contains the spread between the
treasury 10 year and the JPM current net coupon for GN, FN, FH 15 and
30 year mortgages.
*/
#define PPM_REGRESS_FILE_PATH_SIZE 255
static int
ppmRegressReadHistCCFile(
    char *userPath,
    PpmRegressHistoricalRates *histRates,
    TDate asOfDate)
{
    FILE *fp = NULL;
    char buf[200];
    int i;
    TMonthDayYear mdy;
    TDate aDate;
    int mon, day, year;
    double tsy10yr,
	   fn30yrCCRate, fn15yrCCRate,
	   gn30yrCCRate, gn15yrCCRate,
	   fh30yrCCRate, fh15yrCCRate,
	   ARMFixedCCSpread, tsy1yr,
	   ssprd10yr, lib1yr,
	   cmtHump1yr, cmtHump10yr;

	fp = fopen(userPath, "r");
	if (!fp) goto ERROR_EXIT;
    
    for (i = 0; i < 2; i++)
    {
	XONFALSE2(fgets(buf, 200, fp), 
		  "ppmRegressReadCCHistFile() : Error reading history header lines");
    }

    i = 0;
    while (fgets(buf, 200, fp) != NULL)
    {
	XONFALSE2(sscanf(buf, "%d/%d/%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
			 &mon, &day, &year,
			 &tsy10yr,
			 &gn30yrCCRate, &gn15yrCCRate,
			 &fn30yrCCRate, &fn15yrCCRate,
			 &fh30yrCCRate, &fh15yrCCRate,
			 &ARMFixedCCSpread, &tsy1yr,
			 &ssprd10yr, &lib1yr,
			 &cmtHump1yr, &cmtHump10yr)
		 == 16,
		 "ppmRegressReadHistCCFile() : Error parsing rate line");

	mdy.month = mon;
	mdy.day = day;
	mdy.year = year + ((year > 70) ? 1900 : 2000);
	GtoMDYToDate(&mdy, &aDate);
	
	if (i > 0 && TDateCompare(aDate, asOfDate) > 0)
	    break;
	histRates->fn30yrNetCurrCpn = fn30yrCCRate;
	histRates->fn15yrNetCurrCpn = fn15yrCCRate;
	histRates->gn30yrNetCurrCpn = gn30yrCCRate;
	histRates->gn15yrNetCurrCpn = gn15yrCCRate;
	histRates->fh30yrNetCurrCpn = fh30yrCCRate;
	histRates->fh15yrNetCurrCpn = fh15yrCCRate;
	histRates->tsy1yrForSpread = tsy1yr;
	histRates->tsy10yrForSpread = tsy10yr;
	histRates->ARMFixedCCSpread = ARMFixedCCSpread;
	histRates->lib1yrForSpread = lib1yr;
	histRates->ssprd10yrForSpread = ssprd10yr;
	histRates->cmtHump1yrForSpread = cmtHump1yr;
	histRates->cmtHump10yrForSpread = cmtHump10yr;
	i++;
    }
    XONFALSE2(i > 0, "ppmRegressReadHistCCFile(): No current coupon rates in file");

 /* the Off-Run mtg spread adjustment is the blended difference */
 /* between hot-run rates and off-run rates */
    histRates->offRunMtgSpread15YrAdj =
	MSPRD_15YR_TSY1YR_FRAC *
	        ( histRates->tsy1yrForSpread - (histRates->lib1yrForSpread +
						histRates->cmtHump1yrForSpread/100.)) +
	(1 - MSPRD_15YR_TSY1YR_FRAC) *
	        (histRates->tsy10yrForSpread - (histRates->tsy10yrForSpread +
						histRates->cmtHump10yrForSpread/100.));
    
    histRates->offRunMtgSpread30YrAdj =
	MSPRD_30YR_TSY1YR_FRAC *
	        ( histRates->tsy1yrForSpread - (histRates->lib1yrForSpread +
						histRates->cmtHump1yrForSpread/100.)) +
	(1 - MSPRD_30YR_TSY1YR_FRAC) *
	        (histRates->tsy10yrForSpread - (histRates->tsy10yrForSpread +
						histRates->cmtHump10yrForSpread/100.));
    fclose(fp);
    return SUCCESS;
ERROR_EXIT:
    if (fp)
	fclose(fp);
    return PRIMUS_FAILURE;
}
/*
ppmRegressReadHistFile() : reads in the historical rates file and
converts it from weekly data into monthly data.  
*/
#define NUM_WEEKS_IN_AVG 8

int
ppmRegressReadHistFile(
    char *userPathDaily,
    char *userPathWeekly,
    PpmRegressHistoricalRates *histRates,
    int asOfYear, int asOfMonth, int asOfDay)
{
    FILE *fp = NULL;
    char buf[200];
    int nRatesAllocated, nRates, i, j, useLagFull, servicingIdx;
    double *mtgRates30yr = NULL, *mtgRates15yr = NULL,
	   *tsy10yrRates = NULL, *tsy1yrRates = NULL,
	   *fn30yrCCRates = NULL, *fn15yrCCRates = NULL,
	   *gn30yrCCRates = NULL, *gn15yrCCRates = NULL,
	   *fh30yrCCRates = NULL, *fh15yrCCRates = NULL,
	   *mtgRates30yrPts = NULL, *mtgRates15yrPts = NULL,
	   *ssprd10yrRates = NULL, *lib1yrRates = NULL;
    double spreadWeight[NUM_WEEKS_IN_AVG],sumSpreadWeight;

    TDate *rateDates=NULL, startDate, endDate, aDate;
    int mon, day, year;
    double mtgRate30yr, mtgRate15yr, cmt1yr, cmt10yr,
	   fn30yrCCRate, fn15yrCCRate,
	   gn30yrCCRate, gn15yrCCRate,
	   fh30yrCCRate, fh15yrCCRate,
	   mtgRate30yrPts, mtgRate15yrPts,
	   ssprd10yr, lib1yr;
    TDate asOfDate;
    TMonthDayYear mdy;
    TDateInterval  monthInterval, minusMonthInterval;

    GtoMakeDateInterval(-1, 'M', &minusMonthInterval);
    GtoMakeDateInterval(1, 'M', &monthInterval);

    mdy.month = asOfMonth;
    mdy.year = asOfYear;
    mdy.day = asOfDay;
    GtoMDYToDate(&mdy, &asOfDate);

    
    fp = fopen(userPathWeekly, "r");
	if (!fp) goto ERROR_EXIT;

    for (i = 0; i < 2; i++)
    {
	XONFALSE2(fgets(buf, 200, fp), 
		  "ppmRegressReadHistFile() : Error reading history header lines");
    }
    if (strstr(buf, "Use Lag Full"))
	useLagFull = TRUE;
    else
	useLagFull = FALSE;
    nRates = 0;
    nRatesAllocated = 3 * 52;
    XONFALSE2(mtgRates30yr = (double *) malloc(sizeof(double) * nRatesAllocated),
	      "ppmRegressReadHistFile() : Error allocating rates array");

    XONFALSE2(mtgRates15yr = (double *) malloc(sizeof(double) * nRatesAllocated),
	      "ppmRegressReadHistFile() : Error allocating rates array");

    XONFALSE2(tsy10yrRates = (double *) malloc(sizeof(double) * nRatesAllocated),
	      "ppmRegressReadHistFile() : Error allocating rates array");

    XONFALSE2(tsy1yrRates = (double *) malloc(sizeof(double) * nRatesAllocated),
	      "ppmRegressReadHistFile() : Error allocating rates array");

    XONFALSE2(rateDates = (TDate *) malloc(sizeof(TDate) * nRatesAllocated),
	      "ppmRegressReadHistFile() : Error allocating dates array");

    XONFALSE2(fn30yrCCRates = (double *) malloc(sizeof(double) * nRatesAllocated),
	      "ppmRegressReadHistFile() : Error allocating rates array");

    XONFALSE2(fn15yrCCRates = (double *) malloc(sizeof(double) * nRatesAllocated),
	      "ppmRegressReadHistFile() : Error allocating rates array");

    XONFALSE2(gn30yrCCRates = (double *) malloc(sizeof(double) * nRatesAllocated),
	      "ppmRegressReadHistFile() : Error allocating rates array");

    XONFALSE2(gn15yrCCRates = (double *) malloc(sizeof(double) * nRatesAllocated),
	      "ppmRegressReadHistFile() : Error allocating rates array");

    XONFALSE2(fh30yrCCRates = (double *) malloc(sizeof(double) * nRatesAllocated),
	      "ppmRegressReadHistFile() : Error allocating rates array");

    XONFALSE2(fh15yrCCRates = (double *) malloc(sizeof(double) * nRatesAllocated),
	      "ppmRegressReadHistFile() : Error allocating rates array");

    XONFALSE2(mtgRates30yrPts = (double *) malloc(sizeof(double) * nRatesAllocated),
	      "ppmRegressReadHistFile() : Error allocating rates array");

    XONFALSE2(mtgRates15yrPts = (double *) malloc(sizeof(double) * nRatesAllocated),
	      "ppmRegressReadHistFile() : Error allocating rates array");

    XONFALSE2(ssprd10yrRates = (double *) malloc(sizeof(double) * nRatesAllocated),
	      "ppmRegressReadHistFile() : Error allocating rates array");

    XONFALSE2(lib1yrRates = (double *) malloc(sizeof(double) * nRatesAllocated),
	      "ppmRegressReadHistFile() : Error allocating rates array");

    while (fgets(buf, 200, fp) != NULL)
    {
	if (nRatesAllocated <= nRates)
	{
	    nRatesAllocated *= 2;
	    XONFALSE2(mtgRates30yr = (double *)
		      realloc(mtgRates30yr, sizeof(double) * nRatesAllocated),
		      "ppmRegressReadHistFile() : Error allocating rates array");
		      
	    XONFALSE2(mtgRates15yr = (double *)
		      realloc(mtgRates15yr, sizeof(double) * nRatesAllocated),
		      "ppmRegressReadHistFile() : Error allocating rates array");
		      
	    XONFALSE2(tsy10yrRates = (double *)
		      realloc(tsy10yrRates, sizeof(double) * nRatesAllocated),
		      "ppmRegressReadHistFile() : Error re-allocating rates array");

	    XONFALSE2(tsy1yrRates = (double *)
		      realloc(tsy1yrRates, sizeof(double) * nRatesAllocated),
		      "ppmRegressReadHistFile() : Error re-allocating rates array");

	    XONFALSE2(rateDates = (TDate *)
		      realloc(rateDates, sizeof(TDate) * nRatesAllocated),
		      "ppmRegressReadHistFile() : Error re-allocating rates array");

	    XONFALSE2(fn30yrCCRates = (double *)
		      realloc(fn30yrCCRates, sizeof(double) * nRatesAllocated),
		      "ppmRegressReadHistFile() : Error allocating rates array");
		      
	    XONFALSE2(fn15yrCCRates = (double *)
		      realloc(fn15yrCCRates, sizeof(double) * nRatesAllocated),
		      "ppmRegressReadHistFile() : Error allocating rates array");
		      
	    XONFALSE2(gn30yrCCRates = (double *)
		      realloc(gn30yrCCRates, sizeof(double) * nRatesAllocated),
		      "ppmRegressReadHistFile() : Error allocating rates array");
		      
	    XONFALSE2(gn15yrCCRates = (double *)
		      realloc(gn15yrCCRates, sizeof(double) * nRatesAllocated),
		      "ppmRegressReadHistFile() : Error allocating rates array");
		      
	    XONFALSE2(fh30yrCCRates = (double *)
		      realloc(fh30yrCCRates, sizeof(double) * nRatesAllocated),
		      "ppmRegressReadHistFile() : Error allocating rates array");
		      
	    XONFALSE2(fh15yrCCRates = (double *)
		      realloc(fh15yrCCRates, sizeof(double) * nRatesAllocated),
		      "ppmRegressReadHistFile() : Error allocating rates array");
		      
	    XONFALSE2(mtgRates30yrPts = (double *)
		      realloc(mtgRates30yrPts, sizeof(double) * nRatesAllocated),
		      "ppmRegressReadHistFile() : Error allocating rates array");
		      
	    XONFALSE2(mtgRates15yrPts = (double *)
		      realloc(mtgRates15yrPts, sizeof(double) * nRatesAllocated),
		      "ppmRegressReadHistFile() : Error allocating rates array");

	    XONFALSE2(tsy10yrRates = (double *)
		      realloc(tsy10yrRates, sizeof(double) * nRatesAllocated),
		      "ppmRegressReadHistFile() : Error re-allocating rates array");

	    XONFALSE2(tsy1yrRates = (double *)
		      realloc(tsy1yrRates, sizeof(double) * nRatesAllocated),
		      "ppmRegressReadHistFile() : Error re-allocating rates array");

	    XONFALSE2(ssprd10yrRates = (double *)
		      realloc(ssprd10yrRates, sizeof(double) * nRatesAllocated),
		      "ppmRegressReadHistFile() : Error allocating rates array");

	    XONFALSE2(lib1yrRates = (double *)
		      realloc(lib1yrRates, sizeof(double) * nRatesAllocated),
		      "ppmRegressReadHistFile() : Error allocating rates array");
	}
	XONFALSE2(sscanf(buf, "%d/%d/%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
			 &mon, &day, &year,
			 &mtgRate30yr,
			 &mtgRate15yr,
			 &cmt1yr, &cmt10yr,
			 &gn30yrCCRate, &gn15yrCCRate,
			 &fn30yrCCRate, &fn15yrCCRate,
			 &fh30yrCCRate, &fh15yrCCRate,
			 &mtgRate30yrPts,
			 &mtgRate15yrPts,
			 &ssprd10yr, &lib1yr)
		 == 17,
		 "ppmRegressReadHistFile() : Error parsing rate line");
		 
	mtgRates30yr[nRates] = mtgRate30yr;
	mtgRates15yr[nRates] = mtgRate15yr;
	tsy10yrRates[nRates] = cmt10yr;
	tsy1yrRates[nRates] = cmt1yr;
	mdy.month = mon;
	mdy.day = day;
	mdy.year = year + ((year > 70) ? 1900 : 2000);
	GtoMDYToDate(&mdy, rateDates + nRates);
	mtgRates30yrPts[nRates] = mtgRate30yrPts;
	mtgRates15yrPts[nRates] = mtgRate15yrPts;

	fn30yrCCRates[nRates] = fn30yrCCRate;
	fn15yrCCRates[nRates] = fn15yrCCRate;
	gn30yrCCRates[nRates] = gn30yrCCRate;
	gn15yrCCRates[nRates] = gn15yrCCRate;
	fh30yrCCRates[nRates] = fh30yrCCRate;
	fh15yrCCRates[nRates] = fh15yrCCRate;

	ssprd10yrRates[nRates] = ssprd10yr;
	lib1yrRates[nRates] = lib1yr;
	
	nRates++;
    }

    if (useLagFull)
    {
	GtoDateToMDY(rateDates[0], &mdy);
	mdy.day = 1;
	GtoMDYToDate(&mdy, &startDate);

	histRates->firstDateYear = mdy.year;
	histRates->firstDateMonth = mdy.month;
	histRates->firstDateDay = mdy.day;
    }
    else
    {
	GtoDateToMDY(rateDates[0], &mdy);
	mdy.day = 15;
	GtoMDYToDate(&mdy, &startDate);

	histRates->firstDateYear = mdy.year;
	histRates->firstDateMonth = mdy.month;
	histRates->firstDateDay = 1;

	
	GtoDateToMDY(rateDates[0], &mdy);
	if (mdy.day <= 15)
	{
	    GtoDtFwdAny(startDate, &minusMonthInterval, &startDate);
	}
	else
	{
	    mdy.year = histRates->firstDateYear;
	    mdy.month = histRates->firstDateMonth;
	    mdy.day = histRates->firstDateDay;
	    GtoMDYToDate(&mdy, &aDate);
	    GtoDtFwdAny(aDate, &monthInterval, &aDate);
	    GtoDateToMDY(aDate, &mdy);
	    
	    histRates->firstDateYear = mdy.year;
	    histRates->firstDateMonth = mdy.month;
	    histRates->firstDateDay = mdy.day;
	}
    }
    GtoDtFwdAny(startDate, &monthInterval, &endDate);

    histRates->nRates = 0;
    nRatesAllocated = 3 * 12;
    XONFALSE2(histRates->mtgRates30yr = (double *) malloc(sizeof(double) * nRatesAllocated),
	      "ppmRegressReadHistFile() : Error allocating rates array");

    XONFALSE2(histRates->mtgRates15yr = (double *) malloc(sizeof(double) * nRatesAllocated),
	      "ppmRegressReadHistFile() : Error allocating rates array");

    XONFALSE2(histRates->tsy10yrRates = (double *) malloc(sizeof(double) * nRatesAllocated),
	      "ppmRegressReadHistFile() : Error allocating rates array");

    XONFALSE2(histRates->tsy1yrRates = (double *) malloc(sizeof(double) * nRatesAllocated),
	      "ppmRegressReadHistFile() : Error allocating rates array");

    XONFALSE2(histRates->mtgRates30yrPts = (double *) malloc(sizeof(double) * nRatesAllocated),
	      "ppmRegressReadHistFile() : Error allocating rates array");

    XONFALSE2(histRates->mtgRates15yrPts = (double *) malloc(sizeof(double) * nRatesAllocated),
	      "ppmRegressReadHistFile() : Error allocating rates array");

    XONFALSE2(histRates->ssprd10yr = (double *) malloc(sizeof(double) * nRatesAllocated),
	      "ppmRegressReadHistFile() : Error allocating rates array");

    XONFALSE2(histRates->lib1yrRates = (double *) malloc(sizeof(double) * nRatesAllocated),
	      "ppmRegressReadHistFile() : Error allocating rates array");

    i = servicingIdx = 0;
    while (i < nRates && TDateCompare(rateDates[i], asOfDate) <= 0)
    {
	double sumMtgRate30yr = 0., sumMtgRate15yr = 0.,
	       sumTsy10yr = 0., sumTsy1yr = 0.,
	       sumMtgRate30yrPts = 0., sumMtgRate15yrPts = 0.,
	       sumSsprd10yr = 0., sumLib1yr = 0.;
	int ratesInPeriod = 0;
	long busDays=0;
	
	while(((i < nRates && TDateCompare(rateDates[i], asOfDate) <= 0)) &&
	      (TDateCompare(rateDates[i], startDate) < 0))
	    i++;
	while(((i < nRates && TDateCompare(rateDates[i], asOfDate) <= 0)) &&
	      (TDateCompare(rateDates[i], endDate) < 0))
	{
		servicingIdx = i;
		if ((i == 0) || (TDateCompare(rateDates[i-1], startDate) < 0))
		{
			long isHoliday, isWeekend;
			GtoBusinessDaysDiff(startDate, rateDates[i], "NONE", &busDays);
			GtoIsHoliday(startDate, "NONE", &isHoliday);
			GtoIsWeekend(startDate, &isWeekend);
			if (!isHoliday && !isWeekend)
				busDays += 1;
		}
		else
			GtoBusinessDaysDiff(rateDates[i-1], rateDates[i], "NONE", &busDays);
	    
		// busDays *= -1;
		
		sumMtgRate30yr += busDays * mtgRates30yr[i];
 		sumMtgRate15yr += busDays * mtgRates15yr[i];
 		sumTsy10yr += busDays * tsy10yrRates[i];
 		sumTsy1yr += busDays * tsy1yrRates[i];
 		sumMtgRate30yrPts += busDays * mtgRates30yrPts[i];
 		sumMtgRate15yrPts += busDays * mtgRates15yrPts[i];
 		sumSsprd10yr += busDays * ssprd10yrRates[i];
 		sumLib1yr += busDays * lib1yrRates[i];
 
		i++;
 		ratesInPeriod += busDays;
	}
	if ((i < nRates) && (TDateCompare(rateDates[i], asOfDate) <= 0))
	{
		// KLUDGE to handle single refi rate for a month
		if (ratesInPeriod <= 1) 
		{
			sumMtgRate30yr = mtgRates30yr[i-1];
 			sumMtgRate15yr = mtgRates15yr[i-1];
 			sumTsy10yr = tsy10yrRates[i-1];
 			sumTsy1yr = tsy1yrRates[i-1];
 			sumMtgRate30yrPts = mtgRates30yrPts[i-1];
 			sumMtgRate15yrPts = mtgRates15yrPts[i-1];
 			sumSsprd10yr = ssprd10yrRates[i-1];
 			sumLib1yr = lib1yrRates[i-1];
			ratesInPeriod = 1;
		}
		else
		{
			GtoBusinessDaysDiff(rateDates[i-1], endDate, "NONE", &busDays);
			long isHoliday, isWeekend;
			GtoIsHoliday(endDate, "NONE", &isHoliday);
			GtoIsWeekend(endDate, &isWeekend);
			if (!isHoliday && !isWeekend)
				busDays -= 1; 
		
			// busDays *= -1;
		
			sumMtgRate30yr += busDays * mtgRates30yr[i];
 			sumMtgRate15yr += busDays * mtgRates15yr[i];
 			sumTsy10yr += busDays * tsy10yrRates[i];
 			sumTsy1yr += busDays * tsy1yrRates[i];
 			sumMtgRate30yrPts += busDays * mtgRates30yrPts[i];
 			sumMtgRate15yrPts += busDays * mtgRates15yrPts[i];
 			sumSsprd10yr += busDays * ssprd10yrRates[i];
 			sumLib1yr += busDays * lib1yrRates[i];
 
			ratesInPeriod += busDays;
		}
		histRates->isLastHistRateIncomplete = 0;
	}
	else
	{
		GtoBusinessDaysDiff(rateDates[i-1], endDate, "NONE", &busDays);
		// busDays *= -1;
		if (busDays <= 1)
			histRates->isLastHistRateIncomplete = 0;
		else
		{
			histRates->isLastHistRateIncomplete = 1;
			GtoDateToMDY(rateDates[i-1], &mdy);
			histRates->partialMonthDateYear = mdy.year;
			histRates->partialMonthDateMonth = mdy.month;
			histRates->partialMonthDateDay = mdy.day;
		}
	}
	XONFALSE2(ratesInPeriod > 0,
		  "ppmRegressReadHistFile() : Missing some rates in ppm_regress_historical_rates.dat");
	if (nRatesAllocated <= histRates->nRates)
	{
	    nRatesAllocated *= 2;
	    XONFALSE2(histRates->mtgRates30yr = (double *)
		      realloc(histRates->mtgRates30yr,
			      sizeof(double) * nRatesAllocated),
		      "ppmRegressReadHistFile() : Error allocating rates array");
		      
	    XONFALSE2(histRates->mtgRates15yr = (double *)
		      realloc(histRates->mtgRates15yr,
			      sizeof(double) * nRatesAllocated),
		      "ppmRegressReadHistFile() : Error allocating rates array");
		      
	    XONFALSE2(histRates->tsy10yrRates = (double *)
		      realloc(histRates->tsy10yrRates,
			      sizeof(double) * nRatesAllocated),
		      "ppmRegressReadHistFile() : Error re-allocating rates array");

	    XONFALSE2(histRates->tsy1yrRates = (double *)
		      realloc(histRates->tsy1yrRates,
			      sizeof(double) * nRatesAllocated),
		      "ppmRegressReadHistFile() : Error re-allocating rates array");

	    XONFALSE2(histRates->mtgRates30yrPts = (double *)
		      realloc(histRates->mtgRates30yrPts,
			      sizeof(double) * nRatesAllocated),
		      "ppmRegressReadHistFile() : Error allocating rates array");
		      
	    XONFALSE2(histRates->mtgRates15yrPts = (double *)
		      realloc(histRates->mtgRates15yrPts,
			      sizeof(double) * nRatesAllocated),
		      "ppmRegressReadHistFile() : Error allocating rates array");

	    XONFALSE2(histRates->ssprd10yr = (double *)
		      realloc(histRates->ssprd10yr,
			      sizeof(double) * nRatesAllocated),
		      "ppmRegressReadHistFile() : Error re-allocating rates array");

	    XONFALSE2(histRates->lib1yrRates = (double *)
		      realloc(histRates->lib1yrRates,
			      sizeof(double) * nRatesAllocated),
		      "ppmRegressReadHistFile() : Error re-allocating rates array");

	}

	histRates->mtgRates30yr[histRates->nRates] = sumMtgRate30yr / (double) ratesInPeriod;
	histRates->mtgRates15yr[histRates->nRates] = sumMtgRate15yr / (double) ratesInPeriod;
	histRates->tsy10yrRates[histRates->nRates] = sumTsy10yr / (double) ratesInPeriod;
	histRates->tsy1yrRates[histRates->nRates] = sumTsy1yr / (double) ratesInPeriod;
	histRates->mtgRates30yrPts[histRates->nRates] = sumMtgRate30yrPts / (double) ratesInPeriod;
	histRates->mtgRates15yrPts[histRates->nRates] = sumMtgRate15yrPts / (double) ratesInPeriod;
	histRates->ssprd10yr[histRates->nRates] = sumSsprd10yr / (double) ratesInPeriod;
	histRates->lib1yrRates[histRates->nRates] = sumLib1yr / (double) ratesInPeriod;

	(histRates->nRates)++;

	// if (i < nRates && TDateCompare(rateDates[i], asOfDate) < 0)
	//    histRates->isLastHistRateIncomplete = 0;
	// else if (completePeriod(endDate, rateDates[i-1]))
	//    histRates->isLastHistRateIncomplete = 0;
	// else
	//    histRates->isLastHistRateIncomplete = 1;

	GtoDtFwdAny(startDate, &monthInterval, &startDate);
	GtoDtFwdAny(endDate, &monthInterval, &endDate);
    }

    mtgRate30yr = mtgRate15yr = mtgRate30yrPts = mtgRate15yrPts = 0.;
    fn30yrCCRate = fn15yrCCRate = gn30yrCCRate = gn15yrCCRate = fh30yrCCRate = fh15yrCCRate = 0.;
    spreadWeight[0] = 1.;
    sumSpreadWeight = 1.;
    for (i = 1; (i < NUM_WEEKS_IN_AVG) && (i < servicingIdx); i++)
    {
	spreadWeight[i] = spreadWeight[i-1]*.48;
	sumSpreadWeight +=spreadWeight[i];
    }
    for (j = 0; j < i; j++)
	spreadWeight[j] = spreadWeight[j]/sumSpreadWeight;
    for (i = servicingIdx, j = 0; (i >= 0) && (i > servicingIdx - NUM_WEEKS_IN_AVG); i--, j++)
    {
	    fn30yrCCRate += fn30yrCCRates[i]*spreadWeight[j];
	    fn15yrCCRate += fn15yrCCRates[i]*spreadWeight[j];
	    gn30yrCCRate += gn30yrCCRates[i]*spreadWeight[j];
	    gn15yrCCRate += gn15yrCCRates[i]*spreadWeight[j];
	    fh30yrCCRate += fh30yrCCRates[i]*spreadWeight[j];
	    fh15yrCCRate += fh15yrCCRates[i]*spreadWeight[j];

	    mtgRate30yr  += mtgRates30yr[i]*spreadWeight[j];
	    mtgRate15yr  += mtgRates15yr[i]*spreadWeight[j];

	    mtgRate30yrPts  += mtgRates30yrPts[i]*spreadWeight[j];
	    mtgRate15yrPts  += mtgRates15yrPts[i]*spreadWeight[j];
    }
    histRates->fn30yrServicing = (mtgRate30yr - fn30yrCCRate);
    histRates->fn15yrServicing = (mtgRate15yr - fn15yrCCRate);
    histRates->gn30yrServicing = (mtgRate30yr - gn30yrCCRate);
    histRates->gn15yrServicing = (mtgRate15yr - gn15yrCCRate);
    histRates->fh30yrServicing = (mtgRate30yr - fh30yrCCRate);
    histRates->fh15yrServicing = (mtgRate15yr - fh15yrCCRate);

    histRates->latestMtgRate30yrPts = mtgRate30yrPts;
    histRates->latestMtgRate15yrPts = mtgRate15yrPts;
    // histRates->partialMonthDateYear = asOfYear;
    // histRates->partialMonthDateMonth = asOfMonth;
    // histRates->partialMonthDateDay = asOfDay;
    XONFAILURE(ppmRegressReadHistCCFile(userPathDaily, histRates, asOfDate),
	      "ppmRegressReadHistFile() : Error reading current net coupon file");


    if (fp)
	fclose(fp);
    if (mtgRates30yr)
	free(mtgRates30yr);
    if (mtgRates15yr)
	free(mtgRates15yr);
    if (tsy10yrRates)
	free(tsy10yrRates);
    if (tsy1yrRates)
	free(tsy1yrRates);
    if (rateDates)
	free(rateDates);
    if (fn30yrCCRates)
	free(fn30yrCCRates);
    if (fn15yrCCRates)
	free(fn15yrCCRates);
    if (gn30yrCCRates)
	free(gn30yrCCRates);
    if (gn15yrCCRates)
	free(gn15yrCCRates);
    if (fh30yrCCRates)
	free(fh30yrCCRates);
    if (fh15yrCCRates)
	free(fh15yrCCRates);
    if (mtgRates30yrPts)
	free(mtgRates30yrPts);
    if (mtgRates15yrPts)
	free(mtgRates15yrPts);
    if (ssprd10yrRates)
	free(ssprd10yrRates);
    if (lib1yrRates)
	free(lib1yrRates);
    return SUCCESS;
ERROR_EXIT:
    if (fp)
	fclose(fp);
    if (mtgRates30yr)
	free(mtgRates30yr);
    if (mtgRates15yr)
	free(mtgRates15yr);
    if (tsy10yrRates)
	free(tsy10yrRates);
    if (tsy1yrRates)
	free(tsy1yrRates);
    if (rateDates)
	free(rateDates);
    if (fn30yrCCRates)
	free(fn30yrCCRates);
    if (fn15yrCCRates)
	free(fn15yrCCRates);
    if (gn30yrCCRates)
	free(gn30yrCCRates);
    if (gn15yrCCRates)
	free(gn15yrCCRates);
    if (fh30yrCCRates)
	free(fh30yrCCRates);
    if (fh15yrCCRates)
	free(fh15yrCCRates);
    if (mtgRates30yrPts)
	free(mtgRates30yrPts);
    if (mtgRates15yrPts)
	free(mtgRates15yrPts);
    if (ssprd10yrRates)
	free(ssprd10yrRates);
    if (lib1yrRates)
	free(lib1yrRates);
    return PRIMUS_FAILURE;
}

