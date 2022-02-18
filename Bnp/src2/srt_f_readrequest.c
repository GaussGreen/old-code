/**********************************************************************
 *      Name: srt_f_readrequest.c                                     *
 *  Function:                                                         *
 * Copyright: (C) Paribas Capital Markets Ltd.                        *
 *--------------------------------------------------------------------*
 *    Author: Arnaud Sahaguet                                         *
 *      Date: 19/11/94                                                *
 *--------------------------------------------------------------------*
 *    Inputs:                                                         *
 *   Returns:                                                         *
 *   Globals:                                                         *
 *--------------------------------------------------------------------*
 * Modification Record                                                *
 * Date     Inits   Comments                                          *
 * dd/mm/yy                                                           *
 **********************************************************************/

#include "srt_h_all.h"
#include "srt_h_readrequest.h"

static SrtLst* get_next_sigma_date(SrtLst*);
static SrtLst* get_next_tau_date(SrtLst*);

/* When you want to shift the bucket defined by dates d1 and d2, you really
want to move the the sigma or tau given at the date d2.
i.e. sigma(d2) is correct between d1 and d2. If you want to shift sigma(d2),
you don't want to shift what is before d1. -> you don't want to shift d1.
To do so, we define the bucket by date d1+epsilon and d2 */

Err srt_f_readrequest(String* requests, int* buckets, int nr, TermStruct ts, SrtIOStruct* iolist)
{
    int  j;
    Date date1, date2;
    Err  err = NULL;

    for (j = 0; j < nr; j++)
    {
        strupper(requests[j]);
        strip_white_space(requests[j]);

        /* -------------------- Vega request -------------------- */

        if (strcmp(requests[j], "VEGA") == 0)
        {
            /* we convert a bucket number into 2 dates that define the bucket */

            if (err = srt_f_sig_ts_dates(buckets[j], &date1, &date2, ts))
            {
                fprintf(stderr, "INCORRECT SIGMA REQUEST\n");
                continue;
            }

            /* we create a first default request - since no shift value has been specified by the
             * user */
            /* default shift = +0.01 */

            if (err = srt_f_IOstructsetsigmashift(
                    iolist, buckets[j], date1, date2, 0.01, SH_ABSOLUTE, SRT_NO, 0, ""))
            {
                return err;
            }

            continue;
        }

        /* ----------------------- Tau request ------------------- */

        if (strcmp(requests[j], "TAU") == 0)
        {
            /* we convert a bucket number into 2 dates that define the bucket */

            if (err = srt_f_tau_ts_dates(buckets[j], &date1, &date2, ts))
            {
                fprintf(stderr, "INCORRECT TAU REQUEST\n");
                continue;
            }

            /* we create a default request - since no shift value has been specified by the user */
            /* default shift = +1 */

            if (err = srt_f_IOstructsettaushift(
                    iolist, buckets[j], date1, date2, 0.01, SH_ABSOLUTE, SRT_NO, 0, ""))
            {
                return err;
            }

            continue;
        }

        /* ------------------------ RATE request ------------------- */

        if (strcmp(requests[j], "PARALLEL") == 0)
        {
            if (err = srt_f_IOstructsetrshift(
                    iolist, "PARALLEL", 0.01, SH_ABSOLUTE, SRT_NO, 0, "", NULL))
            {
                return err;
            }
            /* we create a rate request for a parallel shift whose shift is defined by the user */
        }

    } /* end of for loop */

    return err;
}
/* ------------------------------------------------------------------------- */

Err srt_f_readresult(String greek, int bucket, double* ret_val)
{
    Err          err;
    SrtIOStruct* iolist;
    double       premium, vega1, vega2, tau1, tau2, value;

    iolist = get_request_list();

    if (!strcmp(greek, "PREMIUM"))
    {
        if (err = srt_f_IOstructgetpremiumval(*iolist, &premium))
        {
            return err;
        }
        value = premium;
    }
    else if (!strcmp(greek, "VEGA"))
    {
        if (err = srt_f_IOstructgetsigmashiftval(*iolist, bucket, 0.01, SH_ABSOLUTE, &vega1))
        {
            return err;
        }

        if (err = srt_f_IOstructgetpremiumval(*iolist, &vega2))
        {
            return err;
        }

        value = 100 * (vega1 - vega2);
    }
    else if (!strcmp(greek, "TAU"))
    {
        if (err = srt_f_IOstructgettaushiftval(*iolist, bucket, 0.01, SH_ABSOLUTE, &tau1))
        {
            return err;
        }

        if (err = srt_f_IOstructgetpremiumval(*iolist, &tau2))
        {
            return err;
        }

        value = 100 * (tau1 - tau2);
    }
    else
        return ("no such greek");

    *ret_val = value;

    return NULL;
}

static SrtLst* get_next_sigma_date(SrtLst* ts)
{
    SrtLst* ls;
    ls = ts;
    while ((ls != NULL) && (((IrmTermStructVal*)ls->element->val.pval)->val_origin != SIGMA_DATE) &&
           (((IrmTermStructVal*)ls->element->val.pval)->val_origin != BOTH_DATE))
        ls = ls->next;

    return ls;
}

/* ------------------------------------------------------------------------- */

static SrtLst* get_next_tau_date(SrtLst* ts)
{
    SrtLst* ls;

    ls = ts;
    while ((ls != NULL) && (((IrmTermStructVal*)ls->element->val.pval)->val_origin != TAU_DATE) &&
           (((IrmTermStructVal*)ls->element->val.pval)->val_origin != BOTH_DATE))
        ls = ls->next;

    return ls;
}

/* ------------------------------------------------------------------------- */

Err srt_f_sig_ts_dates(int bucket_index, Date* date1, Date* date2, TermStruct ts)
{
    Err     err       = NULL;
    int     index     = bucket_index;
    int     index_max = 0;
    SrtLst* ls;
    int     i;

    /*****
       Determine the number of sigma dates in the Term Structure:
            i dates	-> i+1	buckets
    *****/

    for (ls = get_next_sigma_date(ts.head); ls != NULL; ls = get_next_sigma_date(ls))
    {
        if (ls != NULL)
        {
            index_max++;
            ls = ls->next;
        }
    }

    /* we check the bucket defined by bucket_index is correct */

    if (index_max == 0)
    {
        return "No sigma dates";
    }
    if ((bucket_index <= 0) || (bucket_index > index_max + 1))
    {
        return "Bad bucket number...";
    }

    ls = ts.head;
    switch (index)
    {
        /* for the first bucket */
    case 1:
        ls = get_next_sigma_date(ts.head);
        /** the starting date is unknown, so -1 will be ok **/
        *date1 = -1;
        *date2 = (Date)((IrmTermStructVal*)ls->element->val.pval)->date;
        break;

        /* for the other buckets */
    default:

        /* for the last bucket */
        if (index == index_max + 1)
        {
            for (ls = get_next_sigma_date(ts.head), i = index_max - 1; ((ls != NULL) && (i > 0));
                 ls = get_next_sigma_date(ls), i--)
            {
                if (ls != NULL)
                    ls = ls->next;
            }
            /* the last sigma date */
            *date1 = (Date)(((IrmTermStructVal*)(ls->element->val.pval))->date + 0.01);
            /* a date very far in the future, so 100,000 will be ok */
            *date2 = 100000;
        }

        /* for the other buckets */
        else
        {
            for (ls = get_next_sigma_date(ts.head), i = index - 1; ((ls != NULL) && (i > 0));
                 ls = get_next_sigma_date(ls), i--)
            {
                if (ls != NULL)
                    ls = ls->next;
            }
            *date1 = (Date)(((IrmTermStructVal*)(ls->previous->element->val.pval))->date + 0.01);
            *date2 = (Date)((IrmTermStructVal*)ls->element->val.pval)->date;
        }
        break;
    }
    return err;
}
/* ------------------------------------------------------------------------- */

Err srt_f_tau_ts_dates(int bucket_index, Date* date1, Date* date2, TermStruct ts)
{
    Err     err       = NULL;
    int     index     = bucket_index;
    int     index_max = 0;
    SrtLst* ls;
    int     i;

    /*****
       Determine the number of tau dates in the Term Structure:
            i dates	-> i+1	buckets
    *****/

    for (ls = get_next_tau_date(ts.head); ls != NULL; ls = get_next_tau_date(ls))
    {
        if (ls != NULL)
        {
            index_max++;
            ls = ls->next;
        }
    }

    /* we check the bucket defined by bucket_index is correct */

    if (index_max == 0)
    {
        return "No tau dates";
    }
    if ((bucket_index <= 0) || (bucket_index > index_max + 1))
    {
        return "Bad bucket number...";
    }

    ls = ts.head;
    switch (index)
    {
        /* for the first bucket */
    case 1:
        ls = get_next_tau_date(ts.head);
        /** the starting date is unknown, so -1 will be ok **/
        *date1 = -1;
        *date2 = (Date)((IrmTermStructVal*)ls->element->val.pval)->date;
        break;

        /* for the other buckets */
    default:

        /* for the last bucket */
        if (index == index_max + 1)
        {
            for (ls = get_next_tau_date(ts.head), i = index_max - 1; ((ls != NULL) && (i > 0));
                 ls = get_next_tau_date(ls), i--)
            {
                if (ls != NULL)
                    ls = ls->next;
            }
            /* the last tau date */
            *date1 = (Date)(((IrmTermStructVal*)(ls->element->val.pval))->date + 0.01);
            /* a date very far in the future, so 100,000 will be ok */
            *date2 = 100000;
        }

        /* for the other buckets */
        else
        {
            for (ls = get_next_tau_date(ts.head), i = index - 1; ((ls != NULL) && (i > 0));
                 ls = get_next_tau_date(ls), i--)
            {
                if (ls != NULL)
                    ls = ls->next;
            }
            *date1 = (Date)(((IrmTermStructVal*)(ls->previous->element->val.pval))->date + 0.01);
            *date2 = (Date)((IrmTermStructVal*)ls->element->val.pval)->date;
        }
        break;
    }
    return err;
}
