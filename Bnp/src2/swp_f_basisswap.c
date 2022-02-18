/**********************************************************************************************
Basis Swap using CP rates.
Author : Yann Samuelides
*********************************************************************************************/

#include "opfnctns.h"
#include "srt_h_all.h"
#include "swp_h_all.h"
#include "swp_h_vol.h"
#include "utallhdr.h"

Err CPBasisSwap(
    long    StartDate,
    double  Maturity,
    char*   cshortRefRateCode, /* Rate corresponding to the maturity of the CP */
    char*   clongRefRateCode,  /* Rate corresponding to the Libor leg */
    char*   szYieldCurveName,
    double* Vol, /* Vol used to take into account the delayed payment effect */
    double* ConstSpread,
    double* SpikeLibor,
    double* SpikeCP,
    double* Spread, /* Spread CP-Libor for each open day */
    double* Dates,  /* Open days during the life of the swap */
    double* TotalMargins,
    double* MarginsSpreads /* Part of the margin coming from the spread CP-Libor */
)
{
    Err   err = NULL;
    int   i, period, periodCP;
    int   nbsperiod;           /* number of periods in the basis swap */
    long* tPeriodDates = NULL; /*start dates of the different periods in the basis swap with the
                                  last date at the end*/
    long           TheoEndDate;
    long           daycounter;
    double*        SpreadPeriodCP = NULL;
    long           EndFixingPeriodCP, EndFixingPeriodLibor;
    long           today;
    int            spotlag;
    double*        dshortFra             = NULL;
    long*          tEveryDayDates        = NULL;
    long*          tDatesCP              = NULL;
    double*        dlongFra              = NULL;
    double*        PeriodDFs             = NULL;
    double*        DFs                   = NULL;
    long*          ShortEndDates         = NULL;
    int*           numdaysopenperiodCP   = NULL;
    double*        MeanspreadperiodCP    = NULL;
    double*        MeanLiborperiodCP     = NULL;
    double*        MeanconvexityperiodCP = NULL;
    long*          TheoDates             = NULL;
    double         totalmargin, marginspread, marginslope, marginconvexity, marginsecondorder;
    int            numdaysopenperiod = 0;
    int            numperiodCP;
    int            maxnumofdaysperiodCP;
    SwapDP         sdpFra;
    SrtCompounding CPFrequency, BasisSwapFrequency;
    SrtBasisCode   CPbasis, BasisSwapbasis;
    SrtCrvPtr      yldcrv;
    yldcrv  = lookup_curve(szYieldCurveName);
    today   = get_today_from_curve(yldcrv);
    spotlag = get_spotlag_from_curve(yldcrv);

    /* get info from the RefRate Code of the CP */
    if (err = swp_f_get_ref_rate_details(cshortRefRateCode, &CPbasis, &CPFrequency))
    {
        return err;
    }

    /* get info from the RefRate Code of the Swap */
    if (err = swp_f_get_ref_rate_details(clongRefRateCode, &BasisSwapbasis, &BasisSwapFrequency))
    {
        return err;
    }

    nbsperiod   = (int)(Maturity * BasisSwapFrequency);    /* Number of periods in the basis swap */
    numperiodCP = (int)(CPFrequency / BasisSwapFrequency); /* Number of CP periods in a period of
                                                              the basis swap */
    maxnumofdaysperiodCP =
        (int)(360 / CPFrequency); /* Maximum number of open days in a CP period */

    /*Memory allocation*/
    tPeriodDates          = lngvector(0, nbsperiod);
    tEveryDayDates        = lngvector(0, maxnumofdaysperiodCP);
    dshortFra             = dvector(0, maxnumofdaysperiodCP);
    SpreadPeriodCP        = dvector(0, maxnumofdaysperiodCP);
    ShortEndDates         = lngvector(0, maxnumofdaysperiodCP);
    DFs                   = dvector(0, maxnumofdaysperiodCP);
    PeriodDFs             = dvector(0, nbsperiod - 1);
    dlongFra              = dvector(0, nbsperiod - 1);
    MeanspreadperiodCP    = dvector(0, numperiodCP);
    MeanLiborperiodCP     = dvector(0, numperiodCP);
    MeanconvexityperiodCP = dvector(0, numperiodCP);
    numdaysopenperiodCP   = ivector(0, numperiodCP);
    tDatesCP              = lngvector(0, numperiodCP);
    TheoDates             = lngvector(0, nbsperiod);

    /* Initialise the Schedule */
    TheoEndDate = add_unit(
        StartDate, (int)(nbsperiod * 12 / BasisSwapFrequency), SRT_MONTH, NO_BUSDAY_CONVENTION);

    /* Calculate the start dates of each period in the basis swap */
    for (i = 0; i <= nbsperiod; i++)
    {
        TheoDates[i] = add_unit(
            TheoEndDate,
            -(int)((nbsperiod - i) * 12 / BasisSwapFrequency),
            SRT_MONTH,
            NO_BUSDAY_CONVENTION);
        tPeriodDates[i] = add_unit(TheoDates[i], 0, SRT_BDAY, MODIFIED_SUCCEEDING);
    }

    daycounter = 0;

    /* Loop on the periods of the basis swap */
    for (period = 0; period < nbsperiod; period++)
    {
        /* first initialises the sdpFra */
        err = swp_f_setSwapDP(
            add_unit(tPeriodDates[period], spotlag, SRT_BDAY, MODIFIED_SUCCEEDING),
            tPeriodDates[period + 1],
            BasisSwapFrequency,
            BasisSwapbasis,
            &sdpFra);
        if (err)
            return err;

        /* input the spot lag  */
        sdpFra.spot_lag = spotlag;

        /* computation of the Fra of the Period*/
        err = swp_f_ForwardRate_SwapDP(
            &sdpFra, szYieldCurveName, clongRefRateCode, &(dlongFra[period]));
        if (err)
            return err;

        /* computes the DF between today and THE LAST DAY OF THE PERIOD */
        PeriodDFs[period] = swp_f_df(
            today,
            add_unit(tPeriodDates[period + 1], -1, SRT_BDAY, MODIFIED_SUCCEEDING),
            szYieldCurveName);

        /* Calculates the dates of CP periods within the swap period */
        for (i = 0; i <= numperiodCP; i++)
        {
            tDatesCP[i] = add_unit(
                TheoDates[period + 1],
                -(int)((numperiodCP - i) * 12 / CPFrequency),
                SRT_MONTH,
                NO_BUSDAY_CONVENTION);
            tDatesCP[i] = add_unit(tDatesCP[i], 0, SRT_BDAY, MODIFIED_SUCCEEDING);
        }

        /* Loop on the number of days of the CP period */
        for (periodCP = 0; periodCP < numperiodCP; periodCP++)
        {
            /* Determines the NY open days in the CP period */

            tEveryDayDates[0]             = tDatesCP[periodCP];
            i                             = 1;
            numdaysopenperiodCP[periodCP] = 0;
            while (i < maxnumofdaysperiodCP)
            {
                tEveryDayDates[i] =
                    add_unit(tEveryDayDates[i - 1], 1, SRT_BDAY, MODIFIED_SUCCEEDING);
                if (tEveryDayDates[i] >= tDatesCP[periodCP + 1])
                {
                    i = maxnumofdaysperiodCP;
                }
                i++;
                numdaysopenperiodCP[periodCP]++;
            }

            MeanspreadperiodCP[periodCP] = 0;

            /* Calculates the spread between the CP and the Libor for all open days of the period*/
            for (i = 0; i < numdaysopenperiodCP[periodCP]; i++)
            {
                EndFixingPeriodCP = add_unit(tEveryDayDates[i], 30, SRT_DAY, NO_BUSDAY_CONVENTION);
                EndFixingPeriodCP = add_unit(EndFixingPeriodCP, 0, SRT_BDAY, SUCCEEDING);
                EndFixingPeriodLibor = add_unit(
                    add_unit(tEveryDayDates[i], spotlag, SRT_BDAY, MODIFIED_SUCCEEDING),
                    1,
                    SRT_MONTH,
                    NO_BUSDAY_CONVENTION);
                EndFixingPeriodLibor =
                    add_unit(EndFixingPeriodLibor, 0, SRT_BDAY, MODIFIED_SUCCEEDING);
                SpreadPeriodCP[i] = ConstSpread[period];
                if (year(EndFixingPeriodCP) != year(tEveryDayDates[i]))
                {
                    SpreadPeriodCP[i] += SpikeCP[period];
                }
                if (year(EndFixingPeriodLibor) != year(tEveryDayDates[i]))
                {
                    SpreadPeriodCP[i] -= SpikeLibor[period];
                }
                MeanspreadperiodCP[periodCP] += SpreadPeriodCP[i];
            }
            MeanspreadperiodCP[periodCP] =
                MeanspreadperiodCP[periodCP] /
                (numdaysopenperiodCP[periodCP] + 0.0); /* Mean of the spread during the CP period */

            /* Copies the values of the dates of the period and the spreads of the period */
            for (i = 0; i < numdaysopenperiodCP[periodCP]; i++)
            {
                Spread[daycounter] = SpreadPeriodCP[i];
                Dates[daycounter]  = tEveryDayDates[i];
                daycounter++;
            }

            /* Compute the short Term Libors and the DF needed.*/
            MeanLiborperiodCP[periodCP]     = 0;
            MeanconvexityperiodCP[periodCP] = 0;
            for (i = 0; i < numdaysopenperiodCP[periodCP]; i++)
            {
                /* first initialise the sdpFra */
                ShortEndDates[i] = add_unit(
                    add_unit(tEveryDayDates[i], spotlag, SRT_BDAY, MODIFIED_SUCCEEDING),
                    (int)(CPFrequency / 12.0),
                    SRT_MONTH,
                    NO_BUSDAY_CONVENTION);
                err = swp_f_setSwapDP(
                    add_unit(tEveryDayDates[i], spotlag, SRT_BDAY, MODIFIED_SUCCEEDING),
                    ShortEndDates[i],
                    CPFrequency,
                    CPbasis,
                    &sdpFra);
                if (err)
                    return err;

                /* input the spot lag  */
                sdpFra.spot_lag = spotlag;

                /* computation of the Fra */
                err = swp_f_ForwardRate_SwapDP(
                    &sdpFra, szYieldCurveName, cshortRefRateCode, &(dshortFra[i]));
                if (err)
                    return err;

                /* Computes the mean of the short Libor within the CP period */
                MeanLiborperiodCP[periodCP] += dshortFra[i];

                /*	Get the Discount Factors */
                DFs[i] = swp_f_df(
                    today,
                    add_unit(ShortEndDates[i], 0, SRT_BDAY, MODIFIED_SUCCEEDING),
                    szYieldCurveName);

                /* Computes the mean of the convexity effects within the CP period */
                MeanconvexityperiodCP[periodCP] +=
                    dshortFra[i] * (DFs[i] / PeriodDFs[period] - 1) *
                    (add_unit(ShortEndDates[i], 0, SRT_BDAY, MODIFIED_SUCCEEDING) - today - 0.0) /
                    360.0;
            }

            /* Normalizes by the number of days to get the mean */
            MeanLiborperiodCP[periodCP] =
                MeanLiborperiodCP[periodCP] / (numdaysopenperiodCP[periodCP] + 0.0);
            MeanconvexityperiodCP[periodCP] = Vol[period] * Vol[period] *
                                              MeanconvexityperiodCP[periodCP] /
                                              (numdaysopenperiodCP[periodCP] + 0.0);
        }

        /********************* Calculates the margin of the basis swap for the period.***********************************/

        marginspread      = 0;
        marginslope       = 0;
        marginconvexity   = 0;
        marginsecondorder = 0;

        for (periodCP = 0; periodCP < numperiodCP; periodCP++)
        {
            /* Contribution of the spread */
            marginspread += MeanspreadperiodCP[periodCP];
            /* Contribution of the slope */
            marginslope += MeanLiborperiodCP[periodCP];
            /* Contribution of the convexity adjustment */
            marginconvexity += MeanconvexityperiodCP[periodCP];
            /* Contribution of the second order term */
            marginsecondorder +=
                ((numdaysopenperiodCP[periodCP] + 0.0) / 360.0 - 1.0 / (2.0 * CPFrequency)) *
                (MeanspreadperiodCP[periodCP] + MeanLiborperiodCP[periodCP]) *
                (MeanspreadperiodCP[periodCP] + MeanLiborperiodCP[periodCP]);
        }

        /* Compounding adjustment */
        marginsecondorder =
            (BasisSwapFrequency + 0.0) / (CPFrequency + 0.0) *
            (marginsecondorder + 1.0 / (2.0 * CPFrequency) * (marginspread + marginslope) *
                                     (marginspread + marginslope));
        marginspread = marginspread * (BasisSwapFrequency + 0.0) / (CPFrequency + 0.0);
        marginslope =
            marginslope * (BasisSwapFrequency + 0.0) / (CPFrequency + 0.0) - dlongFra[period];
        marginconvexity = marginconvexity * (BasisSwapFrequency + 0.0) / (CPFrequency + 0.0);

        /* Total margin of the period */
        totalmargin = marginspread + marginslope + marginconvexity + marginsecondorder;

        /* Copies the values of the margins of the period */
        MarginsSpreads[period] = marginspread;
        TotalMargins[period]   = totalmargin;
    }

    /*Memory freeage*/
    free_lngvector(tPeriodDates, 0, nbsperiod);
    free_lngvector(tEveryDayDates, 0, maxnumofdaysperiodCP);
    free_dvector(dshortFra, 0, maxnumofdaysperiodCP);
    free_lngvector(ShortEndDates, 0, maxnumofdaysperiodCP);
    free_dvector(DFs, 0, maxnumofdaysperiodCP);
    free_dvector(SpreadPeriodCP, 0, maxnumofdaysperiodCP);
    free_dvector(PeriodDFs, 0, nbsperiod - 1);
    free_dvector(dlongFra, 0, nbsperiod - 1);
    free_dvector(MeanspreadperiodCP, 0, numperiodCP);
    free_dvector(MeanLiborperiodCP, 0, numperiodCP);
    free_dvector(MeanconvexityperiodCP, 0, numperiodCP);
    free_ivector(numdaysopenperiodCP, 0, numperiodCP);
    free_lngvector(tDatesCP, 0, numperiodCP);
    free_lngvector(TheoDates, 0, nbsperiod);

    return err;
}

Err get_fra(long Fradate, int spotlag, char* cRefRateCode, char* szYieldCurveName, double* dFra);

Err srt_f_CPCap(
    long           StartDate,
    long           EndDate,
    double         Strike,
    double         chi,
    char*          cCPRefRateCode, /* Rate corresponding to the CP */
    char*          cRefRateCode,   /* Rate corresponding to the Libor 3M */
    char*          szYieldCurveName,
    char*          szVolCurveName,
    double*        price,
    SrtCallPutType CallPut)
{
    Err            err = NULL;
    long*          fixingdays;
    long           nfixingdays = 0;
    double*        vols;
    double*        Libors;
    double*        spreads;
    long           days;
    double         lognorm = 1.0;
    double         sumCP, sumnormalvol, vol;
    double         maturity;
    long           today, spotlag;
    long           enddate;
    int            i;
    SrtCompounding float_compounding;
    SrtBasisCode   float_basis;
    SrtCrvPtr      yldcrv;
    yldcrv  = lookup_curve(szYieldCurveName);
    today   = get_today_from_curve(yldcrv);
    spotlag = get_spotlag_from_curve(yldcrv);

    err = swp_f_get_ref_rate_details(cRefRateCode, &float_basis, &float_compounding);

    if (err)
    {
        return err;
    }

    /*first get the number of fixingdays*/

    days = StartDate;

    while (days < EndDate)
    {
        nfixingdays++;
        days = add_unit(days, 1, SRT_BDAY, MODIFIED_SUCCEEDING);
    }

    /*Memory Allocation*/
    fixingdays = lngvector(0, nfixingdays - 1);
    vols       = dvector(0, nfixingdays - 1);
    Libors     = dvector(0, nfixingdays - 1);
    spreads    = dvector(0, nfixingdays - 1);

    fixingdays[0] = StartDate;
    spreads[0]    = swp_f_spread(fixingdays[0], fixingdays[0], cCPRefRateCode);

    for (i = 1; i < nfixingdays; i++)
    {
        fixingdays[i] = add_unit(fixingdays[i - 1], 1, SRT_BDAY, MODIFIED_SUCCEEDING);
        spreads[i]    = swp_f_spread(fixingdays[i], fixingdays[i], cCPRefRateCode);
    }

    sumCP        = 0;
    sumnormalvol = 0;

    for (i = 0; i < nfixingdays; i++)
    {
        get_fra(fixingdays[i], spotlag, cRefRateCode, szYieldCurveName, &(Libors[i]));
        enddate = add_unit(
            add_unit(fixingdays[i], spotlag, SRT_BDAY, MODIFIED_SUCCEEDING),
            12 / float_compounding,
            SRT_MONTH,
            MODIFIED_SUCCEEDING);
        err = swp_f_vol(
            szVolCurveName,
            add_unit(fixingdays[i], spotlag, SRT_BDAY, MODIFIED_SUCCEEDING),
            enddate,
            Strike - spreads[i],
            &(vols[i]),
            &lognorm);
        vols[i] = chi * Libors[i] * vols[i];
        sumCP += Libors[i] + spreads[i];
        sumnormalvol += vols[i];
    }

    vol = sumnormalvol / sumCP;

    maturity = (EndDate - today) / 365.25;

    *price = coverage(StartDate, EndDate, float_basis) *
             srt_f_optblksch(sumCP / nfixingdays, Strike, vol, maturity, 1.0, CallPut, PREMIUM);
    *price *= swp_f_df(today, EndDate, szYieldCurveName);

    /*Free Memory*/
    free_lngvector(fixingdays, 0, nfixingdays - 1);
    free_dvector(vols, 0, nfixingdays - 1);
    free_dvector(Libors, 0, nfixingdays - 1);
    free_dvector(spreads, 0, nfixingdays - 1);

    return err;
}
