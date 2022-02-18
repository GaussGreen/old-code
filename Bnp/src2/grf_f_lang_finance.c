/******************************************************************************\
*                Copyright (c) 1995 PARIBAS Capital Markets Group              *
********************************************************************************
*
*   SYSTEM          :   GRF
*
*   MODULE NAME     :   GRF_F_LANG_FINANCE
*
*   PURPOSE         :   Interpretation of GRFN financial functions.
*
*   AUTHOR          :
*
*   DATE            :
*
*   VERSION         :
*
*   DESCRIPTION     :
*
*   FUNCTIONS USED  :
*
*   EXT FUNCTIONS   :   coverage (srt_f_utldates.c)
*
*
********************************************************************************


/******************************************************************************\
*                               Import Include Files                           *
\******************************************************************************/

#include "grf_h_all.h"
#include "grf_h_public.h"
#include "math.h"
#include "srt_h_product_list.h"
#include "srtundutils.h"
#include "swp_h_swap_generic.h"

/******************************************************************************\
*                           Private Function Definitions                       *
\******************************************************************************/

/*****************************************************************************

    FUNCTION    :   grfn_make_swapdates

    DESCRIPTION :   Allocates dates according to sdp stored in dfind,
                    cvg, dfindlen which will find their way into a comll.
                    These are the dates of a swap, and may be used to compute
                    a swap rate, level payment, or cap.

                    dfind points to array of swap dates (type date);

                    cvg points to array of fractions of years in the
                    appropriate basis.

                    dfindlen will be length of these arrays.

******************************************************************************/

static Err grfn_make_swapdates(
    SwapDP* swapdp, Date today, StructType type, Date** df_date, double** cvg, int* num_dates)
{
    Err   err = NULL;
    int   num_start_dates;
    Date* start_dates;
    Date* end_dates;

    /* Create the dates and coiverages of the fixed leg of swap and extracts them */
    err = swp_f_make_FixedLegDatesAndCoverages(
        swapdp, today, df_date, num_dates, &start_dates, &end_dates, cvg, &num_start_dates);

    srt_free(start_dates);
    srt_free(end_dates);

    return err;

} /* END static Err grfn_make_swapdates(...) */

/* ---------------------------------------------------------------------------- */

/* Function to compute right away all the pay dates , the spreads (with their coverage)
   on the floating leg  */
static Err grfn_make_floating_leg(
    SwapDP*    float_sdp,
    Date       today,
    StructType type,
    Date**     float_df_date,
    double**   float_cvg,
    double**   spread,
    int*       float_dfindlen,
    String     ref_rate_code)
{
    Err   err = NULL;
    int   num_start_dates;
    Date* fixing_dates;
    Date* start_dates;
    Date* end_dates;

    /* If there is no reference rate mentionned, return an error */
    if (!ref_rate_code)
        return serror("No reference rate specified in grfn_make_floating_leg");

    /* If dealing with CASH, there is no need to do this ( float leg = df(f) - df(i) )*/
    if (!strcmp(ref_rate_code, "CASH"))
    {
        *float_dfindlen = 0;
        return NULL;
    }

    /* Make the floating leg: dates, times, cvg, spreads ...*/
    err = swp_f_make_FloatLegDatesCoveragesAndSpreads(
        float_sdp,
        today,
        ref_rate_code,
        float_df_date,
        float_dfindlen,
        &fixing_dates,
        &start_dates,
        &end_dates,
        float_cvg,
        spread,
        &num_start_dates);

    if (start_dates)
        srt_free(start_dates);
    if (end_dates)
        srt_free(end_dates);
    if (fixing_dates)
        srt_free(fixing_dates);

    return err;

} /* END grfn_make_floating_leg(...) */

/* ----------------------------------------------------------------------------- */

/* -----------------------------------------------------------------------------

    FUNCTION        :   check_financial_args

    CALLED BY       :   grf_interp_fin_func

    DESCRIPTION     :   This function checks the arguments of a GRFN financial
                        function against the contents of a coded string
                        comprising of special characters defined in the
                        GRFN symbol table.

                        Example: The string "secb" will be decoded as types
                        of 4 args: (s)tart, (e)nd, (c)ompounding, (b)asis.

                        If there is a "start" and an "end", direction is set.

   ----------------------------------------------------------------------------- */

Err check_financial_args(
    String          name,
    COMLLType       type,
    COMLL_PTR       lastarg,
    GrfnDeal*       gd,
    Date*           s,
    long*           e_nfp,
    BasisCode*      b,
    SrtCompounding* c,
    long*           aux,
    String          short_ref_rate,
    String          long_ref_rate,
    double**        dstore,
    long*           dstorelen,
    comllchar**     sstore,
    long*           sstorelen)
{
    /*..Local variables..*/

    int    i, code_len;
    int    n = 0;
    String code; /* Code for the argument types */
    Err    err;

    /*..Main Body..*/

    code = NULL;

    /* Set the expected sequence of arguments for the function */
    /* s = start, e = end, c = compounding, b = basis, d = double, t = string */
    switch (type)
    {
    case COMLL_F_PAYOFF:
        code = "dt";
        break;

    case COMLL_F_CMS:
        code = "secbddtr";
        break;

    case COMLL_F_SWAP:
        code = "secblr";
        break;

    case COMLL_F_LEVEL:
    case COMLL_F_CLEAN:
    case COMLL_F_YIELD:

        code = "secb";
        break;

    case COMLL_F_CVG:
        code = "seb";
        break;

    case COMLL_F_FRA:

        code = "sebr";
        break;

    case COMLL_F_DF:

        code = "se";
        break;

    case COMLL_F_ACCINT:

        code = "aasb";
        break;

    } /* end switch(type) */

    /* Return if code was not set up above */

    if (code != NULL)
        code_len = strlen(code);
    else
        return NULL;

    /* Makes the link between the code letter and the COMLL_ */
    for (i = 0; i < code_len; i++)
    {
        switch (code[i])
        {
        case 'g': /* Go */

            break;

        case 'd': /* Double */

            (*dstorelen)++;

            (*dstore) = (double*)realloc((*dstore), (*dstorelen) * sizeof(double));
            if ((*dstore) == NULL)
                return serror("%s %s", GRERR_MEMORY_ALLOC, name);

            (*dstore)[(*dstorelen) - 1] = lastarg->dval;

            break;

        case 't': /* sTring */

            (*sstorelen)++;

            (*sstore) = (comllchar*)realloc((*sstore), (*sstorelen) * sizeof(comllchar));
            if ((*sstore) == NULL)
                return serror("%s %s", GRERR_MEMORY_ALLOC, name);

            strcpy((*sstore)[(*sstorelen) - 1], lastarg->sval);

            break;

        case 'a': /* Aux range */

            aux[n++] = DTOL(lastarg->dval);

            if (n == 1)
            {
                CheckAuxRange1(aux[n - 1], gd);
            }
            else
            {
                CheckAuxRange2(aux[n - 1], aux[n - 2], gd);
            }

            break;

        case 'r': /* short reference Rate name */

            if (!lastarg)
            {
                /* No reference rate is mentionned: set the default one : CASH */
                strncpy(short_ref_rate, "CASH", strlen("CASH"));
                short_ref_rate[strlen("CASH")] = '\0';
            }
            else
            {
                if (lastarg->type != COMLL_STRING)
                    return serror("%s: %s", GRERR_BAD_REF_RATE, lastarg->sval);
                else
                {
                    strncpy(short_ref_rate, lastarg->sval, strlen(lastarg->sval));
                    short_ref_rate[strlen(lastarg->sval)] = '\0';
                }

                lastarg->type = COMLL_REAL;
                lastarg->dval = 0.0;
            }

            break;

        case 'l': /* Long reference rate name (for Swaps, CMT,...) */

            if (!lastarg)
            {
                /* No reference rate is mentionned: set the default one : NONE */
                strncpy(long_ref_rate, "", strlen(""));
                long_ref_rate[strlen("")] = '\0';
            }
            else
            {
                if (lastarg->type != COMLL_STRING)
                    return serror("%s: %s", GRERR_BAD_REF_RATE, lastarg->sval);
                else
                {
                    strncpy(long_ref_rate, lastarg->sval, strlen(lastarg->sval));
                    long_ref_rate[strlen(lastarg->sval)] = '\0';
                }

                lastarg->type = COMLL_REAL;
                lastarg->dval = 0.0;
            }

            break;

        case 'b': /* Basis */

            if ((lastarg->type != COMLL_STRING) || (err = interp_basis(lastarg->sval, b)))
            {
                return serror("%s: %s", GRERR_BAD_BASIS, lastarg->sval);
            }
            lastarg->type = COMLL_REAL;
            lastarg->dval = 0.0;

            break;

        case 'c': /* Compounding */

            if ((lastarg->type != COMLL_STRING) || (err = interp_compounding(lastarg->sval, c)))
                return serror("%s: %s", GRERR_BAD_COMPD, lastarg->sval);

            lastarg->type = COMLL_REAL;
            lastarg->dval = 0.0;

            break;

        case 's': /* Start date */

            *s = (Date)DTOL(lastarg->dval);

            if (*s < gd->event_dates[gd->I] && type != COMLL_F_CVG)
                return serror("%s: %d", GRERR_BAD_START, *s);

            break;

        case 'e': /* End date */

            *e_nfp = (Date)DTOL(lastarg->dval);

            switch (type)
            {
                /* If end is an nfp (3 for instance) it is in month for CVG, FRA, DF */
            case COMLL_F_CVG:
            case COMLL_F_FRA:
            case COMLL_F_DF:

                if (*e_nfp < 500)
                    *c = SRT_MONTHLY;

                break;
            } /* end switch on type (for 'e') */

            break;

        } /* end switch for code[i] */

        /* Move to the next argument (goes from s <- e <- c <- b <- r...) */
        if (lastarg)
            lastarg = lastarg->prev;

    } /* end for (i=0; i < code_len; i++) */

    return NULL;

} /* END check_financial_arguments(...) */

/* ============================================================================== */

/* -----------------------------------------------------------------------------

    FUNCTION        :   parse_financial_args

    CALLED BY       :   grf_interp_fin_func

    DESCRIPTION     :   This function checks the arguments of a GRFN financial
                        function against the contents of a coded string
                        comprising of special characters defined in the
                        GRFN symbol table.  It does not overwrite any information,
                                                unlike check_  financial_args

                        Example: The string "secb" will be decoded as types
                        of 4 args: (s)tart, (e)nd, (c)ompounding, (b)asis.

                        If there is a "start" and an "end", direction is set.

   ----------------------------------------------------------------------------- */

Err parse_financial_args(
    String          name,
    COMLLType       type,
    COMLL_PTR       lastarg,
    GrfnDeal*       gd,
    Date*           s,
    long*           e_nfp,
    BasisCode*      b,
    SrtCompounding* c,
    long*           aux,
    String          short_ref_rate,
    String          long_ref_rate)
{
    /*..Local variables..*/

    int    i, code_len;
    int    n = 0;
    String code; /* Code for the argument types */
    Err    err;

    /*..Main Body..*/

    code = NULL;

    /* Set the expected sequence of arguments for the function */
    /* s = start, e = end, c = compounding, b = basis, d = double, t = string */
    switch (type)
    {
    case COMLL_F_PAYOFF:
    case COMLL_F_CMS:
    case COMLL_F_ACCINT:
        code = 0;
        break;

    case COMLL_F_SWAP:
        code = "secblr";
        break;

    case COMLL_F_LEVEL:
    case COMLL_F_CLEAN:
    case COMLL_F_YIELD:

        code = "secb";
        break;

    case COMLL_F_CVG:
        code = "seb";
        break;

    case COMLL_F_FRA:

        code = "sebr";
        break;

    case COMLL_F_DF:

        code = "se";
        break;

    } /* end switch(type) */

    /* Return if code was not set up above */

    if (code != NULL)
        code_len = strlen(code);
    else
        return "Incorrect code";

    /* Makes the link between the code letter and the COMLL_ */
    for (i = 0; i < code_len; i++)
    {
        switch (code[i])
        {
        case 'r': /* short reference Rate name */

            if (!lastarg)
            {
                /* No reference rate is mentionned: set the default one : CASH */
                strncpy(short_ref_rate, "CASH", strlen("CASH"));
                short_ref_rate[strlen("CASH")] = '\0';
            }
            else
            {
                if (lastarg->type != COMLL_STRING)
                    return serror("%s: %s", GRERR_BAD_REF_RATE, lastarg->sval);
                else
                {
                    strncpy(short_ref_rate, lastarg->sval, strlen(lastarg->sval));
                    short_ref_rate[strlen(lastarg->sval)] = '\0';
                }
            }

            break;

        case 'l': /* Long reference rate name (for Swaps, CMT,...) */

            if (!lastarg)
            {
                /* No reference rate is mentionned: set the default one : NONE */
                strncpy(long_ref_rate, "", strlen(""));
                long_ref_rate[strlen("")] = '\0';
            }
            else
            {
                if (lastarg->type != COMLL_STRING)
                    return serror("%s: %s", GRERR_BAD_REF_RATE, lastarg->sval);
                else
                {
                    strncpy(long_ref_rate, lastarg->sval, strlen(lastarg->sval));
                    long_ref_rate[strlen(lastarg->sval)] = '\0';
                }
            }

            break;

        case 'b': /* Basis */

            if ((lastarg->type != COMLL_STRING) || (err = interp_basis(lastarg->sval, b)))
            {
                return serror("%s: %s", GRERR_BAD_BASIS, lastarg->sval);
            }

            break;

        case 'c': /* Compounding */

            if ((lastarg->type != COMLL_STRING) || (err = interp_compounding(lastarg->sval, c)))
                return serror("%s: %s", GRERR_BAD_COMPD, lastarg->sval);

            break;

        case 's': /* Start date */

            *s = (Date)DTOL(lastarg->dval);

            if (*s < gd->event_dates[gd->I] && type != COMLL_F_CVG)
                return serror("%s: %d", GRERR_BAD_START, *s);

            break;

        case 'e': /* End date */

            *e_nfp = (Date)DTOL(lastarg->dval);

            switch (type)
            {
                /* If end is an nfp (3 for instance) it is in month for CVG, FRA, DF */
            case COMLL_F_CVG:
            case COMLL_F_FRA:
            case COMLL_F_DF:

                if (*e_nfp < 500)
                    *c = SRT_MONTHLY;

                break;
            } /* end switch on type (for 'e') */

            break;

        } /* end switch for code[i] */

        /* Move to the next argument (goes from s <- e <- c <- b <- r...) */
        if (lastarg)
            lastarg = lastarg->prev;

    } /* end for (i=0; i < code_len; i++) */

    return NULL;

} /* END check_financial_arguments(...) */

/* ============================================================================== */

/* ------------------------------------------------------------------------------

    FUNCTION    :   pay_function

    CALLED BY   :   grfn_interp_fin_func

    DESCRIPTION :   Pay requires a special treatment also because it
                                        triggers payment reports.
                                        In case of simulations, the report does not have to be
                                        done (the report should only happen if in the past)

                                        The comll should look like this when it is passed in:

                                           (top)                   (bottom)
                                            Date          -->  [Amount Expression]
                                            REAL                  REAL/COMLL

                                        (only the arguments within the PAY function are past)
                                        lastargs points to the last argument passed to PAY:
                                        the payment date

                                        If we are in the past and DATE is in the past,
                    don't do anything.

                                        If we are in the past and DATE in in the future:
                                        adds TODAY as the start date of the discounting:
                                                AMOUNT   -->  PAYFROMPAST(TODAY,DATE)
                                                (top)               (bottom)
                                        that will correspond to the following operations
                                                    AMOUNT * DF(TODAY,paydate)
                                        and will allow for proper computations once the
                                        COMLL is evaluated without failing on the use of DF

                    If we are in the future, use the COMLL_PAY type,
                                        setting the discretisation date as the start date
                                        of the discounting:
                                                        AMOUNT --> PAY..(NOW,DATE),
                    that will correspond to the following operations
                                                    AMOUNT * DF(NOW,paydate)
                                        once the COMLL is evaluated


   ------------------------------------------------------------------------------ */

#include "swp_h_all.h"

static Err pay_function(String name, GrfnSymbol* gs, GrfnDeal* gd, COMLL_PTR lastarg)
{
    double      pv_of_cashflow = 0.0;
    Err         err            = NULL;
    COMLL_PTR   top, bot, comll;
    int         und_index = 0;
    double      df_start;
    SRT_Boolean require_df = SRT_YES;

    /* Go to the first (top) element of this COMLL: should be the pay date */
    top = comll_gototop(lastarg);

    /* Makes sure that pay is not used anywhere else but in the last column of the tableau */
    if (gd->J != gd->sswidth - 1)
    {
        return serror("%s", GRERR_PAY_LAST_COL);
    }

    /* Make sure the pay date (top element) is a real; if not return an error message */
    test_real(1, name, top);

    /* Checks that the payment date is after the event date */
    if (gd->nowdt > top->dval)
        return serror("%s: PAY(%d)", GRERR_BAD_END_OR_NFP, DTOL(top->dval));

    /* Adds the "Reference to Pay" status to the cell: the cell should only contain PAY */
    GrfnAddStatus(GrfnCGCell(gd).status, GRFNCSPAYREF);

    /* Pay is used in the PAST */
    if (gd->is_history)
    {
        /* Pay is used in the past with a payment date before today (+ end of day)*/
        if ((long)top->dval < (gd->today + gd->end_of_day_payment))
        {
            /* Adds an extra element in the COMLL with a PAYINPAST type ( to keep track of payment)
             */
            bot       = comll_gotobot(top);
            bot       = comll_insert_after(bot);
            bot->type = COMLL_F_PAYINPAST;

            /* There will be no need of using df: the payment is in the past */
            require_df = SRT_NO;
        }
        /* Pay is used in the past, with (complex) payment in the future : make it a [ AMOUNT *
         * PAYFROMPAST ] COMLL */
        else
        {
            /* Adds an extra element in the COMLL with a PAYFROMPAST type */
            bot       = comll_gotobot(top);
            bot       = comll_insert_after(bot);
            bot->type = COMLL_F_PAYFROMPAST;

            /* Set TODAY as the date to which the cash flow will be discounted */
            df_start = gd->today;

        } /* END PAY used in past with complex amount */

    } /* END if (gd->is_history)  */
    else
    /* Pay is used in the FUTURE, with a simple or complex payment expression: make it a [ AMOUNT *
     * PAY ] COMLL */

    {
        /* Adds an extra element in the COMLL with a PAY type */
        bot       = comll_gotobot(top);
        bot       = comll_insert_after(bot);
        bot->type = COMLL_F_PAY;

        /* Set NOW as the date to which the cash flow will be discounted (date in the future) */
        df_start = gd->nowdt;

    } /* END if !(gd->is_history) */

    /* Keeps track of thepayment date in the COMLL_F_PAY... element */
    bot->dval = top->dval;

    /* If necessary, make initialisation for df calculations */
    if (require_df == SRT_YES)
    {
        /* Stores the index of the DOMESTIC underlying (just in case...) */
        bot->ivec[0] = 0;

        /* Stores the discount factor dates (df_start and pay date) in the bottom element (as for a
         * DF type ) */
        bot->dfindlen = 2;
        bot->dfind    = srt_calloc(2, sizeof(Date));
        bot->dfind[0] = DTOL(df_start);
        bot->dfind[1] = DTOL(top->dval);
    }

    /* Removes the top element of the COMLL, where the payment date was stored (keeping address of
     * top the same) */
    comll = top->next;
    memcpy(top, comll, sizeof(COMLL_STR));
    top->next->prev = top;
    top->prev       = NULL;
    comll->prev     = NULL;
    comll->next     = NULL;
    comll_free_node(comll);

    return NULL;

} /* END  static Err pay_function (...) */

static Err payn_function(String name, GrfnSymbol* gs, GrfnDeal* gd, COMLL_PTR lastarg)
{
    double      pv_of_cashflow = 0.0;
    Err         err            = NULL;
    COMLL_PTR   top, bot, comll;
    int         und_index = 0;
    double      df_start;
    SRT_Boolean require_df = SRT_YES;

    /* Go to the first (top) element of this COMLL: should be the pay date */
    top = comll_gototop(lastarg);

    /* Makes sure that pay is not used anywhere else but in the last column of the tableau */
    /*    if (gd->J != gd->sswidth - 1)
            {
            return serror("%s",GRERR_PAY_LAST_COL);
        }
    */

    /* Make sure the pay date (top element) is a real; if not return an error message */
    test_real(1, name, top);

    /* Checks that the payment date is after the event date */
    if (gd->nowdt > top->dval)
        return serror("%s: PAY(%d)", GRERR_BAD_END_OR_NFP, DTOL(top->dval));

    /* Adds the "Reference to Pay" status to the cell: the cell should only contain PAY */
    //	GrfnAddStatus(GrfnCGCell(gd).status, GRFNCSPAYREF);

    /* Pay is used in the PAST */
    if (gd->is_history)
    {
        /* Pay is used in the past with a payment date before today (+ end of day)*/
        if ((long)top->dval < (gd->today + gd->end_of_day_payment))
        {
            /* Adds an extra element in the COMLL with a PAYINPAST type ( to keep track of payment)
             */
            bot       = comll_gotobot(top);
            bot       = comll_insert_after(bot);
            bot->type = COMLL_F_PAYINPAST;

            /* There will be no need of using df: the payment is in the past */
            require_df = SRT_NO;
        }
        /* Pay is used in the past, with (complex) payment in the future : make it a [ AMOUNT *
         * PAYFROMPAST ] COMLL */
        else
        {
            /* Adds an extra element in the COMLL with a PAYFROMPAST type */
            bot       = comll_gotobot(top);
            bot       = comll_insert_after(bot);
            bot->type = COMLL_F_PAYFROMPAST;

            /* Set TODAY as the date to which the cash flow will be discounted */
            df_start = gd->today;

        } /* END PAY used in past with complex amount */

    } /* END if (gd->is_history)  */
    else
    /* Pay is used in the FUTURE, with a simple or complex payment expression: make it a [ AMOUNT *
     * PAY ] COMLL */

    {
        /* Adds an extra element in the COMLL with a PAY type */
        bot       = comll_gotobot(top);
        bot       = comll_insert_after(bot);
        bot->type = COMLL_F_PAY;

        /* Set NOW as the date to which the cash flow will be discounted (date in the future) */
        df_start = gd->nowdt;

    } /* END if !(gd->is_history) */

    /* Keeps track of thepayment date in the COMLL_F_PAY... element */
    bot->dval = top->dval;

    /* If necessary, make initialisation for df calculations */
    if (require_df == SRT_YES)
    {
        /* Stores the index of the DOMESTIC underlying (just in case...) */
        bot->ivec[0] = 0;

        /* Stores the discount factor dates (df_start and pay date) in the bottom element (as for a
         * DF type ) */
        bot->dfindlen = 2;
        bot->dfind    = srt_calloc(2, sizeof(Date));
        bot->dfind[0] = DTOL(df_start);
        bot->dfind[1] = DTOL(top->dval);
    }

    /* Removes the top element of the COMLL, where the payment date was stored (keeping address of
     * top the same) */
    comll = top->next;
    memcpy(top, comll, sizeof(COMLL_STR));
    top->next->prev = top;
    top->prev       = NULL;
    comll->prev     = NULL;
    comll->next     = NULL;
    comll_free_node(comll);

    return NULL;

} /* END  static Err pay_function (...) */

/* ============================================================================ */

/******************************************************************************\


    FUNCTION        :   grfn_interp_fin_func

    CALLED BY       :   grfn_f_lang

    DESCRIPTION     :   Interpret financial function
                                                Function called by the YACC parser (through the
                                                grfn_interp_func function) when encountering a
                                                NAME(...) type of input string that correspond
                                                to a financial name

                                                The COMLL lastarg points to the last argument
                                                of the function, which is the top of the COMLL
                                                (parsed backwards)

                        NB. Arguments to financial functions must be
                        deterministic.

    *                                                                         *
    *                       OPTIONAL ARGUMENTS                                *
    *                       ------------------                                *
    *                                                                         *
    * If we are looking at an event date in the past (prior to calculation    *
    * date), then we can have a maximum of 3 "optional" arguments specified   *
    *                                                                         *
    * Optional argument 1 = Name of Underlying  (string)                      *
    *   "        "      2 = Reference Rate Code (string) OR value (real)      *
    *   "        "      3 = Reference Rate Code (string)                      *
    *                                                                         *
    * When a Reference Rate Code is specified, then we must have a function   *
    * pointer which will return the required rates given the code.            *
    *                                                                         *
    * NB. For an American event (uses AM) no action is required.              *
    *                                                                         *
    * [Note] The following logic will need to be generalised:                 *
    *                                                                         *
    * We assume that for all functions containing optional arguments, the 1st *
    * optional argument is the name of the underlying. This may not be true   *
    * for functions with optional arguments which may have nothing to do with *
    * the underlying!                                                         *
    *                                                                         *
    *                                                                         *


\******************************************************************************/

Err grfn_interp_fin_func(String name, GrfnSymbol* gs, GrfnDeal* gd, COMLL_PTR lastarg)

{
    /*..Declare local variables..*/

    SwapDP           sdp;
    long             start, fixing;
    long             end_nfp;
    SrtCompounding   comp;
    BasisCode        basis;
    SrtDiffusionType difftype;
    Err              err;
    int              i, und_index, opt_args;
    long             aux[3];

    COMLL_PTR top;        /* pointer to the LAST  argument of function */
    COMLL_PTR bot;        /* pointer to the FIRST argument of function */
    COMLL_PTR comll;      /* Temporary local pointer */
    COMLL_PTR comll_temp; /* Temporary local pointer */

    char long_ref_rate[SRTBUFSZ]  = "";
    char short_ref_rate[SRTBUFSZ] = "";
    char szYieldCurve[SRTBUFSZ]   = "";

    long    dstorelen  = 0;
    double* dstore     = NULL;
    long    sstorelen  = 0;
    comllchar(*sstore) = NULL;

    SrtProduct* product = NULL;
    int         n_unds = 0, *n_dfs = NULL;
    long**      dates = NULL;

    /* Points to first and last arguments of the function */
    top = comll_gototop(lastarg);
    bot = comll_gotobot(lastarg);

    /* The pay function overwrites everything : it is treated separately */
    if (gs->type == COMLL_F_PAY)
    {
        if (gd->pass == 2)
        {
            return (pay_function(name, gs, gd, lastarg));
        }
        else
        {
            return NULL;
        }
    }

    /* The pay function overwrites everything : it is treated separately */
    if (gs->type == COMLL_F_PAYN)
    {
        if (gd->pass == 2)
        {
            return (payn_function(name, gs, gd, lastarg));
        }
        else
        {
            return NULL;
        }
    }

    /* Checks function arguments are all deterministic */
    i     = 0;
    comll = bot;

    while (comll && i < gs->ndargs)
    {
        test_deterministic(i, name, comll);
        comll = comll->prev;
        i++;
    }

    /* How many optional arguments were specified for this function in the parsing  */
    opt_args = top->nargs - gs->nrargs - gs->ndargs;

    /* Check the number and type of the optional arguments */
    if (opt_args > gs->no_opt_args)
        return serror("%s: %s", name, GRERR_WRONG_NUM_OPT_ARG);

    if (opt_args == 1 && top->type != COMLL_STRING)
    {
        return serror("%s: %s", name, GRERR_WRONG_TYPE_OPT_ARG);
    }

    /* If we have no optional arguments and we are just extracting
           names of underlyings, then return */
    if (opt_args == 0 && gd->pass == 1)
        return (NULL);

    /* ----------- UNDERLYING NAME EXTRACTION ---------------------- */

    /* Sets comll on what should be the underlying name (first optional argument) */
    comll = top;
    for (i = opt_args; i > 1; i--)
    {
        comll = comll->next;
    }

    /* Set up and update the und_index (on first parse, make the list of names) */
    und_index = 0;
    if (opt_args >= 1 && comll->type == COMLL_STRING)
    {
        /* On first parsing, collect the underlying name */
        if (gd->pass == 1)
        {
            /* Store the name of the current und (irm) in the Grfn list of underlyings */
            i = grfn_store_und_name(comll->sval, (String)gd->domestic_und);

            return (NULL);
        }
        else
            /* On second parse, get the underlying index from the und name */
            if (gd->pass == 2)
        {
            und_index = get_grf_und_index(comll->sval);

            /* Cut underlying name out of the COMLL */
            if (opt_args != 1)
            {
                /* If not at top: just cut COMLL making sure connexions are preserved */
                comll->prev->next = comll->next;
                comll->next->prev = comll->prev;
                comll->next       = NULL;
                comll->prev       = NULL;
                comll_free_list(comll);
            }
            else
            {
                /* Comll and top are the same: copy comll->prev into top
                          (not to modify lastarg adress) and cut comll->prev */
                comll = comll->next;
                memcpy(top, comll, sizeof(COMLL_STR));
                top->next->prev = top;
                top->prev       = NULL;
                comll->prev     = NULL;
                comll->next     = NULL;
                comll_free_list(comll);
            } /* END if opt_args == 1 */

        } /* END if (gd->pass == 2) */

    } /* END  if (opt_args >= 1 && top->type == COMLL_STRING) */

    /* From now on, the underlying name has been stripped out:any optional argument
       left is a short ref rate, a long one or both (opt_args has not been changed) */

    /*   Reset comll on the first optional argument (not the und name anymore) */
    comll = top;
    for (i = opt_args - 1; i > 1; i--)
    {
        comll = comll->next;
    }

    /* -------------------------------------------------------------------------
                             Non-American Historical Events
       ------------------------------------------------------------------------- */

    if (gd->is_history == SRT_YES && grfn_lang_cur_amcell == SRT_NO)
    {
        /* If we have specified coverage (CVG) in the past calculate it! */

        if (gs->type == COMLL_F_CVG)
        {
            err = check_financial_args(
                name,
                gs->type,
                lastarg,
                gd,
                &start,
                &end_nfp,
                &basis,
                &comp,
                aux,
                short_ref_rate,
                long_ref_rate,
                &dstore,
                &dstorelen,
                &sstore,
                &sstorelen);

            if (err)
                return (err);

            /* Make sure the end date is properly defined and set (not only nfp) */
            err = swp_f_setSwapDP(start, end_nfp, comp, basis, &sdp);
            if (err)
                return err;
            end_nfp = add_unit(sdp.end, 0, SRT_BDAY, MODIFIED_SUCCEEDING);

            /* Represent the coverage as a simple REAL COMLL with its value */
            top->dval = coverage(sdp.start, end_nfp, sdp.basis_code);
            top->type = COMLL_REAL;
            comll_cut(top);

            return NULL;
        } /* end if (gs->type == COMLL_F_CVG) */

        /* Return an error message for DF and PVRNG in the PAST */

        if (gs->type == COMLL_F_DF || gs->type == COMLL_F_PVRANGE || gs->type == COMLL_F_CMS ||
            gs->type == COMLL_F_LEVEL || gs->type == COMLL_F_PAYOFF)
        {
            return serror("%s: %s", name, " Historical use of function is illegal!");
        }

        /* Check the number of optional arguments input : there should be at least one ref rate */
        if (((gs->type == COMLL_F_SWAP) && (opt_args < 2)) ||
            ((gs->type != COMLL_F_SWAP) && (opt_args != gs->no_opt_args)))
        {
            return serror("%s: %s", name, GRERR_NO_REFRATE);
        }

        /* If the LAST argument is a string and we have at least 2 optional arguments
           then we MUST have a valid (non-NULL) fixing function */
        if ((top->type == COMLL_STRING) && (gs->no_opt_args >= 2))
        {
            /* get the info */
            if (err = check_financial_args(
                    name,
                    gs->type,
                    lastarg,
                    gd,
                    &start,
                    &end_nfp,
                    &basis,
                    &comp,
                    aux,
                    short_ref_rate,
                    long_ref_rate,
                    &dstore,
                    &dstorelen,
                    &sstore,
                    &sstorelen))
            {
                return err;
            }
            /* Make sure the end date is properly defined and set (not only nfp) */
            err = swp_f_setSwapDP(start, end_nfp, comp, basis, &sdp);
            if (err)
                return err;

            sdp.spot_lag = 2;

            /* calculate the fixing date */
            fixing = add_unit(sdp.start, -sdp.spot_lag, SRT_BDAY, NO_BUSDAY_CONVENTION);
            if (fixing < gd->today + gd->end_of_day_fixing)
            {
                /* Have we got a function to get the reference rates fixings ? */
                if (gd->fixing_fct == NULL)
                {
                    return serror("%s", GRERR_NO_HIST_FNC);
                }

                /* Get the historical fixing for the first ref rate code ( the long one) on that
                 * date */
                err = gd->fixing_fct(DTOL(fixing), comll->sval, &top->dval);

                if (err != NULL)
                    return err;
            }
            else
            {
                /* get the market name from the underlying */
                SrtGetUnderlyingYieldCurveName(gd->domestic_und, szYieldCurve);
                top->dval =
                    swp_f_fra(sdp.start, sdp.end, sdp.basis_code, szYieldCurve, comll->sval);
            }

            /* We have a REAL value returned, so replace COMLL type */
            top->type = COMLL_REAL;

        } /* END if (gs->no_opt_args >= 2) */

        /* Free list after top (after LAST argument) */

        comll_cut(top);

        return NULL;

    } /* END if (gd->is_history == SRT_YES && grfn_lang_cur_amcell == SRT_NO) */

    /* ----------------------------------------------------------------------------
                  We are looking at an event date at calculation date.
                  So check if we have a reference rate and a fixing
                  or just evaluate the financial function - done in the next loop
       ---------------------------------------------------------------------------- */
    if ((gd->nowdt == gd->today) && ((gs->type == COMLL_F_FRA) || (gs->type == COMLL_F_SWAP)))
    {
        /* Check the number of optional arguments */
        if (opt_args == gs->no_opt_args)
        {
            /* If the LAST argument is a string and we have at least 2 optional arguments
                   then we must have a valid (non-NULL) function */
            if ((top->type == COMLL_STRING) && (gs->no_opt_args >= 2))
            {
                /* get the info */
                comll_temp = lastarg;
                if (err = parse_financial_args(
                        name,
                        gs->type,
                        comll_temp,
                        gd,
                        &start,
                        &end_nfp,
                        &basis,
                        &comp,
                        aux,
                        short_ref_rate,
                        long_ref_rate))
                {
                    return err;
                }
                /* Have we got a function to get the reference rates fixings */
                if (gd->fixing_fct != NULL)
                {
                    /* Check to see if the fixing date of the FRA is today, or in the future */
                    /* Make sure the end date is properly defined and set (not only nfp) */
                    err = swp_f_setSwapDP(start, end_nfp, comp, basis, &sdp);
                    if (err)
                        return err;

                    sdp.spot_lag = 2;

                    /* calculate the fixing date */
                    fixing = add_unit(sdp.start, -sdp.spot_lag, SRT_BDAY, NO_BUSDAY_CONVENTION);
                    if (fixing < gd->today + gd->end_of_day_fixing)
                    {
                        /* Get the historical fixing for the ref rate on that date */
                        err = gd->fixing_fct(DTOL(gd->nowdt), comll->sval, &top->dval);
                        if (err == NULL)
                        {
                            /* We have a REAL value returned, so use it (set type) */
                            top->type = COMLL_REAL;

                            /* Free list after top (after LAST argument) */
                            comll_cut(top);
                            return NULL;

                        } /* END if (err == NULL) */
                        /* If error in getting fixing and end of day fixing flag is on : error */
                        else if (gd->end_of_day_fixing == SRT_YES)
                        {
                            return serror("%s %s", err, GRERR_MISSING_FIXING);
                        }
                    }
                    else
                    {
                        /* get the market name from the underlying */
                        SrtGetUnderlyingYieldCurveName(gd->domestic_und, szYieldCurve);
                        top->dval = swp_f_fra(
                            sdp.start, sdp.end, sdp.basis_code, szYieldCurve, comll->sval);
                        /* We have a REAL value returned, so use it (set type) */
                        top->type = COMLL_REAL;

                        /* Free list after top (after LAST argument) */
                        comll_cut(top);
                        return NULL;
                    }
                } /* END gd->hist_fnc != NULL && gd->hist_fnc->fnc != NULL) */
                  /* If no function but end of day fixing flag is on : error */
                else if (gd->end_of_day_fixing == SRT_YES)
                {
                    return serror("%s", GRERR_NO_HIST_FNC);
                }

            } /* END if (gs->no_opt_args >=2) */

        } /* END if (opt_args == gs->no_opt_args) */

    } /* END if (gd->today == gd->nowdt) */

    /* --------------------------------------------------------------------------------
                    Now we are looking at an event date in "future" (i.e. on or after the
                calculation date so discard the optional arguments (if any)) .
       -------------------------------------------------------------------------------- */

    /* --------------------------------------------------------------------------
       comll point now to the first optional argument that is not the und name
       Since we need to be able to collect the short reference rate name, we only
       get rid of the first optionnal arguments for a COMLL_F_SWAP.
       The und name has already been stored in the first parsing of the tableau,
       and removed from the COMLL at the second parsing.
       If there is no optionnal reference rate code, the function
       check_financial_arg() will set a default string
       -------------------------------------------------------------------------- */

    /* Check arguments of financial functions and sets the sdp & ref rates  accordingly */
    if (err = check_financial_args(
            name,
            gs->type,
            lastarg,
            gd,
            &start,
            &end_nfp,
            &basis,
            &comp,
            aux,
            short_ref_rate,
            long_ref_rate,
            &dstore,
            &dstorelen,
            &sstore,
            &sstorelen))
    {
        return err;
    }

    /* Don't need to create a swap for a PVRANGE and PAYOFF */
    if ((gs->type != COMLL_F_PVRANGE) && (gs->type != COMLL_F_PAYOFF))
    {
        /* Make sure the end date is properly defined and set (not only nfp) */
        err = swp_f_setSwapDP(start, end_nfp, comp, basis, &sdp);
        if (err)
            return err;

        sdp.spot_lag = 2;
    }
    /* ----------------------------------------------------------------------------------------
       Prepares the arguments (dates, space for df, cvg...) that will be required for the
       computation of the financial function, and stack them into top
       ---------------------------------------------------------------------------------------- */

    if (gs->type == COMLL_F_YIELD || gs->type == COMLL_F_CLEAN)
    {
        top = top->next->next;
    }
    top->type = gs->type;

    switch (gs->type)
    {
    case COMLL_F_YIELD:
    case COMLL_F_CLEAN:

        top->ptr = (SwapDP*)srt_calloc(1, sizeof(SwapDP));
        memcpy(top->ptr, &sdp, sizeof(SwapDP));
        break;

    case COMLL_F_PAYOFF:

        /* Look up the product and free sstore and dstore */
        product = lookup_product(sstore[0]);

        top->ptr = srt_calloc(1, sizeof(SrtProduct*));
        /* top->ptr will be destructed by GRFN */

        memcpy(top->ptr, &product, sizeof(SrtProduct*));
        top->dval = dstore[0];

        srt_free(dstore);
        srt_free(sstore);

        if (product == NULL)
            return serror("COMLL_F_PAYOFF: Product not found");
        if (und_index < 0)
            return serror("COMLL_F_PAYOFF: Underlying not found");

        /* Request DFs dates if needed */
        if (product->RequestDfDates)
        {
            err = product->RequestDfDates(product, (long)top->dval, &n_unds, &n_dfs, &dates);
            if (err)
                return err;

            /* Save DFs dates for the domestic und */
            top->ivec[0]  = 0;
            top->dfind    = dates[0];
            top->dfindlen = n_dfs[0];

            /* Save DFs dates for the foreign und (insert another COMLL) */
            if (n_unds > 1)
            {
                top           = top->next;
                top->type     = COMLL_F_PAYOFF;
                top->ptr      = NULL;
                top->ivec[0]  = und_index;
                top->dfind    = dates[1];
                top->dfindlen = n_dfs[1];
            }

            /* Clean up */
            free(n_dfs);
            free(dates);

            /* Extra check for quantos */
            if (n_unds > 1 && und_index == 0)
                return serror("Foreign underlying not specified");
        }

        break;

    case COMLL_F_CMS: /* should certainly be merged with COMLL_F_SWAP */

        /* Store und_index in ivec, swaps dates, cvg, sprd ... in top */
        if (und_index >= 0)
        {
            top->ivec[0] = und_index;

            /* Keep the swapdp structure */
            top->ptr = (SwapDP*)srt_calloc(1, sizeof(SwapDP));
            memcpy(top->ptr, &sdp, sizeof(SwapDP));

            /* Create and populate the fixed leg (dates + cvg...): end is properly set */
            sdp.first_full_fixing = sdp.start;
            err =
                grfn_make_swapdates(&sdp, gd->today, SWAP, &top->dfind, &top->cvg, &top->dfindlen);
            if (err)
                return err;

            /* Create and populate the floating leg (backward from same end as fixed leg) */
            sdp.direction = BKWD;
            err           = grfn_make_floating_leg(
                &sdp,
                gd->today,
                SWAP,
                &top->dfindfloat,
                &top->cvgfloat,
                &top->spread,
                &top->dfindfloatlen,
                short_ref_rate);
            if (err)
                return err;

            /* Transform LOGNORMAL/NORMAL into a SrtDiffusionType */
            /* Store it into dstore and free sstore */
            if (err = interp_diffusion_type(sstore[0], &difftype))
                return err;

            dstorelen++;

            dstore = (double*)realloc(dstore, dstorelen * sizeof(double));
            if (dstore == NULL)
                return serror("%s %s", GRERR_MEMORY_ALLOC, name);

            dstore[dstorelen - 1] = (double)difftype;

            if (sstore != NULL)
            {
                srt_free(sstore);
            }

            /* Modify the delay (in dstore[1]) */
            if (dstore[1] < 30000) /* delay is a coverage */
                dstore[1] = (double)sdp.start + dstore[1] * DAYS_IN_YEAR;

            /* Fill top->dstore array - the only one we need - */
            top->dstore    = dstore;
            top->dstorelen = dstorelen;
        }

        break;

    case COMLL_F_SWAP:

        /* Store und_index in ivec, swaps dates, cvg, sprd ... in top */
        if (und_index >= 0)
        {
            top->ivec[0] = und_index;

            /* Create and populate the fixed leg (dates + cvg...): end is properly set */
            sdp.first_full_fixing = sdp.start;
            err =
                grfn_make_swapdates(&sdp, gd->today, SWAP, &top->dfind, &top->cvg, &top->dfindlen);
            if (err)
                return err;

            /* Create and populate the floating leg (backward from same end as fixed leg) */
            sdp.direction = BKWD;
            err           = grfn_make_floating_leg(
                &sdp,
                gd->today,
                SWAP,
                &top->dfindfloat,
                &top->cvgfloat,
                &top->spread,
                &top->dfindfloatlen,
                short_ref_rate);
            if (err)
                return err;

            /* If the long ref rate is mentionned and the short ref rate == CASH, get long spread */
            if (!(strcmp(short_ref_rate, "CASH")) && (strcmp(long_ref_rate, "") != 0))
            {
                top->long_spread = swp_f_spread(sdp.start, sdp.end, long_ref_rate);
                if (top->long_spread == SRT_SPREAD_ERROR)
                    return serror(
                        "Cannot compute long spread for %s(%d,%d)",
                        long_ref_rate,
                        sdp.start,
                        sdp.end);
            }
            else
            /* Otherwise, no long spread is required */
            {
                top->long_spread = 0.0;
            }
        }
        break;

    case COMLL_F_LEVEL:

        /* Store und_index in ivec, swaps dates, cvg,  ... in top */
        if (und_index >= 0)
        {
            top->ivec[0] = und_index;

            sdp.first_full_fixing = sdp.start;
            err =
                grfn_make_swapdates(&sdp, gd->today, SWAP, &top->dfind, &top->cvg, &top->dfindlen);
            if (err)
                return err;
        }
        break;

    case COMLL_F_FRA:

        /* Store und_index */
        if (und_index >= 0)
        {
            top->ivec[0] = und_index;

            /* Fwd Cash */
            top->dfindlen = 2;
            top->cvg      = srt_calloc(1, sizeof(double));
            top->dfind    = srt_calloc(2, sizeof(Date));
            top->dfind[0] = sdp.start;
            top->dfind[1] = bus_date_method(sdp.end, MODIFIED_SUCCEEDING);
            top->cvg[0]   = coverage(top->dfind[0], top->dfind[1], sdp.basis_code);

            /*Spread */
            top->spread        = srt_calloc(1, sizeof(double));
            top->dfindfloatlen = 0;
            top->spread[0]     = swp_f_spread(top->dfind[0], top->dfind[1], short_ref_rate);
            if (top->spread[0] == SRT_SPREAD_ERROR)
                return serror(
                    "Error in getting %s(%d,%d) spread",
                    short_ref_rate,
                    top->dfind[0],
                    top->dfind[1]);
        }

        break;

    case COMLL_F_DF:

        /* Store und_index */
        if (und_index >= 0)
        {
            top->ivec[0] = und_index;

            top->dfindlen = 2;
            top->cvg      = srt_calloc(1, sizeof(double));
            top->dfind    = srt_calloc(2, sizeof(Date));
            top->dfind[0] = sdp.start;
            top->dfind[1] = sdp.end;
            top->cvg[0]   = coverage(sdp.start, sdp.end, sdp.basis_code);
        }

        break;

    case COMLL_F_CVG:

        top->dval = coverage(sdp.start, sdp.end, sdp.basis_code);
        top->type = COMLL_REAL;

        break;

    case COMLL_F_ACCINT:

        /* Find first i, s.t. ith element is less than start */

        for (i = 0; i < gd->auxlen[aux[0]]; i++)
        {
            if (DTOL(floor(gd->aux[aux[0]][i])) <= sdp.start)
                break;
        }

        if (i == gd->auxlen[aux[0]])
        {
            top->dval = 0.0;
        }
        else
        {
            top->dval = coverage(DTOL(floor(gd->aux[aux[0]][i])), sdp.start, sdp.basis_code) *
                        gd->aux[aux[1]][i];
        }

        top->type = COMLL_REAL;

        break;

    case COMLL_F_PVRANGE:

        if (gs->ndargs == 3)
            err = grfn_pvrange_setup(name, gd, top, und_index);
        else
            err = grfn_dirtyprice_setup(name, gd, top, und_index);

        break;

    } /* end switch(gs->type) */

    /* Just make this instruction one COMLL */
    if (top->next != NULL)
        comll_cut(top);

    return err;
}
