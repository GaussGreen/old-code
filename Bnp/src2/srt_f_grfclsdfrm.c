/* -------------------------------------------------------------------------

   FILE NAME	: srt_f_grfclsdfrm.c

   PURPOSE		: functions to price swaptions, caps or bondoptions without
                  the Grfn interface. THese functions take simplified
                                  descriptions of the products to be priced, and construct
                                  arrays of Grfn cells to feed to the Grfn function
                                  If LGM can be used, then switch to closed form
   ------------------------------------------------------------------------- */
/*
   Add srt_f_newgrfclsdfrm() for a full GrfnTab computation
   JLX
*/

#include "grf_h_public.h"
#include "srt_h_all.h"
#include "srt_h_closedform.h"
#include "srt_h_grfclsdfrm.h"

static void grfn_make_bondoption_string(
    String          bndopt_str,
    SwapDP*         sdp,
    double          strike,
    double          bk,
    double          acc_int,
    Date            pay_date,
    SrtReceiverType rec_pay);

static void grfn_mk_resetcap_spr(
    Date*           eventdates,
    GrfnCell**      sprdsht,
    SwapDP*         sdp,
    double          strike,
    SrtReceiverType rec_pay,
    Date*           dfind,
    double*         cvg,
    long            nrows,
    long            ncols,
    char*           basis_str,
    int             num_caplets_skipped);

static void grfn_mk_bdt_bnd_opt_sht(
    GrfnCell***     sprdshtptr,
    SwapDP*         sdp,
    double          strike,
    double          bk,
    double          acc_int,
    SrtReceiverType rec_pay,
    Date*           eventdatesptr,
    double*         cvg,
    int             numeventdates,
    int             i_start);
static void grfn_mk_bdt_cap_flr_sht(
    GrfnCell***     sprdshtptr,
    SwapDP*         sdp,
    double          strike,
    double          bk,
    double          acc_int,
    SrtReceiverType rec_pay,
    Date*           eventdatesptr,
    double*         cvg,
    int             numeventdates,
    int             numcaplets,
    int             istart);

///// added by Albert Wang 08/25/03 - begin
static void grfn_make_swaption_fix_string(
    String  swptn_str,
    SwapDP* sdp,
    long    start,
    long    end,
    double  strike,
    String  und_name,
    String  ref_rate);

static void grfn_make_swaption_floating_string(
    String swptn_str, SwapDP* sdp, long start, long end, String und_name, String ref_rate);
///// added by Albert Wang 08/25/03 - end

#define SHIFT 0.0001

/* --------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------
   FUNCTION	:       srt_f_grfn_clsdfrm
   DESCRIPTION	:   Prices a swaption, a cap or a bond option using Grfn
                    If LGM can be used, a switch send the calculation to a
                                        real closed form
   --------------------------------------------------------------------------- */

Err srt_f_grfn_clsdfrm(
    SrtUndPtr       und,
    SrtGrfnParam*   grfnparam,
    SwapDP*         sdp,
    double          strike,
    double          bk,
    SrtReceiverType rec_pay,
    StructType      type,
    String          ref_rate_code,
    double*         answer)
{
    Err          err = NULL;
    long         ncols, nrows, numeventdates;
    Date*        eventdates = NULL;
    GrfnCell**   sprdsht    = NULL;
    SrtIOStruct* iolist;
    SrtMdlType   mdl_type;
    SrtMdlDim    mdl_dim;
    String       und_name;
    Date         today;

    /* Gets the Model type (LGM, Cheyette,...) */
    err = get_underlying_mdltype(und, &mdl_type);

    /* Gets the Number of Factors in the Model */
    err = get_underlying_mdldim(und, &mdl_dim);

    /* IF LGM or EtaBeta, there is a REAL closed form: go to it */
    if (((mdl_type == ETABETA || mdl_type == LGM || mdl_type == NEWLGM) ||
         (grfnparam->closed_form_type == QUICK_CLSDFRM))
        ///// added by Albert Wang 08/25/03 - begin
        && (type != SIMPLEMIDAT)
        ///// added by Albert Wang 08/25/03 - end
    )
    {
        err = srt_f_closed_form(
            und, mdl_type, mdl_dim, sdp, strike, bk, rec_pay, type, ref_rate_code, answer);

        if (err)
            return err;
    }
    else
    {
        /* There is no closed form available: have to build a full grfn tableau */
        if (mdl_type == BDT)
        {
            err = grfn_SwapDP_to_BDT_GrfnCells(
                &numeventdates,
                &eventdates,
                &nrows,
                &ncols,
                &sprdsht,
                und,
                sdp,
                strike,
                bk,
                rec_pay,
                type);
            if (err)
                return err;
        }

        else
        {
            und_name = get_underlying_name(und);
            today    = get_today_from_underlying(und);
            err      = grfn_SwapDP_to_GrfnCells(
                &numeventdates,
                &eventdates,
                &nrows,
                &ncols,
                &sprdsht,
                today,
                sdp,
                strike,
                bk,
                rec_pay,
                type,
                und_name,
                ref_rate_code);
            if (err)
                return err;
            if (nrows == 0)
                return serror("Product has expired in srt_f_grfn_clsdfrm");
        }

        err = srt_f_IOstructcreate(&iolist, "clsdfrm");

        /* Call Grfn with this tableau */
        if (!err)
            err = srt_f_grfn(
                und,
                grfnparam,
                numeventdates,
                &eventdates,
                &nrows,
                &ncols,
                &sprdsht,
                0,
                0,
                0,
                0,
                0,
                iolist,
                0,
                0);

        if (!err)
            err = srt_f_IOstructgetpremiumval(*iolist, answer);

        if (!err)
            err = srt_f_IOstructfree(&iolist);

        if (eventdates)
            srt_free(eventdates);

        if (sprdsht)
        {
            grfn_free_GrfnCellmatrix(sprdsht, nrows, ncols);
        }
    }

    return err;

} /* END srt_f_grfn_clsdfrm(...) */

/* --------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------
   FUNCTION	:       srt_f_grfn_clsdfrm_vega
   DESCRIPTION	:   Computes Vegas for srt_f_grfn_clsdfrm from TS sigmas
   --------------------------------------------------------------------------- */

Err srt_f_grfn_clsdfrm_vega(
    SrtUndPtr       undptr,
    SrtGrfnParam*   grfnparam,
    SwapDP*         sdp,
    double          strike,
    double          bk,
    SrtReceiverType rec_pay,
    StructType      type,
    String          ref_rate_code,
    double*         initial_price,
    double**        sigma_vega,
    long*           num_sig,
    double**        tau_vega,
    long*           num_tau)
{
    TermStruct *initial_ts, *shifted_ts;
    Date        today;
    SrtMdlType  mdl_type;
    SrtMdlDim   mdl_dim;
    double *    sigma[6], *tau[2];
    double      shifted_price;
    long        i;
    Err         err;

    today = get_today_from_underlying(undptr);

    err = get_underlying_mdltype(undptr, &mdl_type);
    if (err)
    {
        return err;
    }
    err = get_underlying_mdldim(undptr, &mdl_dim);
    if (err)
    {
        return err;
    }
    if ((mdl_type != LGM) && (mdl_type != CHEY))
    {
        return serror(
            "Error in srt_f_grfn_clsdfrm_vega: ts vega works only with LGM and CHEY underlyings ");
    }
    if (mdl_dim != ONE_FAC)
    {
        return serror("Error in srt_f_grfn_clsdfrm_vega: ts vega works only with 1F underlyings ");
    }
    err = get_underlying_ts(undptr, &initial_ts);
    if (err)
    {
        return err;
    }

    err = srt_f_grfn_clsdfrm(
        undptr, grfnparam, sdp, strike, bk, rec_pay, type, ref_rate_code, initial_price);
    if (err)
    {
        return err;
    }

    err = srt_f_display_IRM_OneFac_TermStruct(
        initial_ts,
        &(sigma[0]),
        &(sigma[1]),
        &(sigma[2]),
        &(sigma[3]),
        &(sigma[4]),
        &(sigma[5]),
        num_sig,
        &(tau[0]),
        &(tau[1]),
        num_tau);
    if (err)
    {
        return err;
    }

    *sigma_vega = srt_calloc(*num_sig, sizeof(double));
    if (!*sigma_vega)
    {
        return serror("Error in srt_f_grfn_clsdfrm_vega: memory allocation (1) ");
    }
    *tau_vega = srt_calloc(*num_tau, sizeof(double));
    if (!*tau_vega)
    {
        return serror("Error in srt_f_grfn_clsdfrm_vega: memory allocation (2) ");
    }

    for (i = 0; i < *num_sig; i++)
    {
        sigma[1][i] += SHIFT;
        err = srt_f_init_IRM_TermStruct(
            today,
            sigma,
            2,
            *num_sig,
            tau,
            2,
            *num_tau,
            mdl_type,
            mdl_dim,
            0.0,
            0.0,
            /* BETAETA */
            0.0,
            0.0,
            0.0,
            0.0,
            0.0, /* vasicek parms */
            0,
            0,
            NULL,
            &shifted_ts);
        if (err)
        {
            return err;
        }

        set_irund_ts(undptr, shifted_ts);

        err = srt_f_grfn_clsdfrm(
            undptr, grfnparam, sdp, strike, bk, rec_pay, type, ref_rate_code, &shifted_price);
        if (err)
        {
            return err;
        }

        err = srt_f_free_IRM_TermStruct(&shifted_ts);
        if (err)
        {
            return err;
        }

        sigma[1][i] -= SHIFT;
        (*sigma_vega)[i] = (shifted_price - *initial_price) / SHIFT;
    }

    for (i = 0; i < *num_tau; i++)
    {
        tau[1][i] += SHIFT;
        err = srt_f_init_IRM_TermStruct(
            today,
            sigma,
            2,
            *num_sig,
            tau,
            2,
            *num_tau,
            mdl_type,
            mdl_dim,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0, /* vasicek parms */
            0,
            0,
            NULL,
            &shifted_ts);
        if (err)
        {
            return err;
        }

        set_irund_ts(undptr, shifted_ts);

        err = srt_f_grfn_clsdfrm(
            undptr, grfnparam, sdp, strike, bk, rec_pay, type, ref_rate_code, &shifted_price);
        if (err)
        {
            return err;
        }

        err = srt_f_free_IRM_TermStruct(&shifted_ts);
        if (err)
        {
            return err;
        }

        tau[1][i] -= SHIFT;
        (*tau_vega)[i] = (shifted_price - *initial_price) / SHIFT;
    }

    srt_free(sigma[0]);
    srt_free(tau[0]);
    srt_free(sigma[1]);
    srt_free(tau[1]);
    srt_free(sigma[2]);
    srt_free(sigma[3]);
    srt_free(sigma[4]);
    srt_free(sigma[5]);

    set_irund_ts(undptr, initial_ts);

    return NULL;

} /* END srt_f_grfn_clsdfrm_vega(...) */

/* --------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------
   FUNCTION	:       grfn_SwapDP_to_GrfnCells
   DESCRIPTION	:   Given a SwapDP (swaps date parameters), a strike and a type,
                    generate the GrfnCells and EventDates corresponding to a
                                        cap, a swaption or a bond option.
   --------------------------------------------------------------------------- */

Err grfn_SwapDP_to_GrfnCells(
    long*           numeventdatesptr,
    Date**          eventdatesptr,
    long*           nrowsptr,
    long*           ncolsptr,
    GrfnCell***     sprdshtptr,
    Date            today,
    SwapDP*         sdp,
    double          strike,
    double          bk,
    SrtReceiverType rec_pay,
    StructType      type,
    String          und_name,
    String          ref_rate_code)
{
    String     basis_str;
    long       numeventdates;
    Date*      fixing_dates;
    Date*      start_dates;
    Date*      end_dates;
    Date*      pay_dates;
    int        num_dates;
    int        num_pay_dates;
    Date       first_pay_date;
    long       nrows, ncols;
    GrfnCell** sprdsht;
    int        num_caplets_skipped;
    SwapDP     caplet;
    int        i;
    ///// added by Albert Wang 08/25/03 - begin
    int j;
    ///// added by Albert Wang 08/25/03 - end
    Err     err        = NULL;
    double  accrue_int = 0.0;
    double* cvg;

    switch (type)
    {
        /* Need to calculate correct accrued interest for bond option, which is then
           subtracted from the bond strike */
    case BOND_OPTION:
        nrows = ncols = numeventdates = 1;
        swp_f_make_FixedLegDatesAndCoverages(
            sdp, today, &pay_dates, &num_pay_dates, &start_dates, &end_dates, &cvg, &num_dates);

        fixing_dates    = (Date*)malloc(1 * sizeof(Date));
        fixing_dates[0] = add_unit(start_dates[0], -sdp->spot_lag, SRT_BDAY, MODIFIED_SUCCEEDING);

        sprdsht    = GrfnCellmatrix(nrows, ncols, GRFN_DEF_ARGBUFSZ);
        accrue_int = acc_int_fct(*sdp, strike);
        if (accrue_int != 0)
            first_pay_date = pay_dates[1];

        grfn_make_bondoption_string(
            sprdsht[0][0].sval, sdp, strike, bk, accrue_int, first_pay_date, rec_pay);
        sprdsht[0][0].type = GRFNSCELL;
        srt_free(start_dates);
        srt_free(end_dates);
        srt_free(pay_dates);
        srt_free(cvg);

        break;

    case SWAPTION:
        nrows = ncols = numeventdates = 1;
        swp_f_make_FixedLegDatesAndCoverages(
            sdp, today, &pay_dates, &num_pay_dates, &start_dates, &end_dates, &cvg, &num_dates);

        fixing_dates    = (Date*)malloc(1 * sizeof(Date));
        fixing_dates[0] = add_unit(start_dates[0], -sdp->spot_lag, SRT_BDAY, MODIFIED_SUCCEEDING);

        sprdsht = GrfnCellmatrix(nrows, ncols, GRFN_DEF_ARGBUFSZ);
        err     = grfn_make_swaption_string(
            sprdsht[0][0].sval, sdp, strike, rec_pay, und_name, ref_rate_code);
        if (err)
            return err;

        sprdsht[0][0].type = GRFNSCELL;

        srt_free(start_dates);
        srt_free(end_dates);
        srt_free(pay_dates);
        srt_free(cvg);

        break;

        ///// added by Albert Wang 08/25/03 - begin
    case SIMPLEMIDAT:

        swp_f_make_FixedLegDatesAndCoverages(
            sdp,
            today,
            &pay_dates,
            &num_pay_dates,
            &start_dates,
            &end_dates,
            &cvg,
            &num_dates);  /// NB: num_dates == # of fix leg start dates = # of fix leg end dates

        /// allocate mem for fixing_dates
        numeventdates = num_dates;
        fixing_dates  = (Date*)malloc(numeventdates * sizeof(Date));
        //// propagate fixing dates
        for (i = 0; i < numeventdates; ++i)
        {
            fixing_dates[i] =
                add_unit(start_dates[i], -sdp->spot_lag, SRT_BDAY, MODIFIED_SUCCEEDING);
        }

        nrows = num_dates;
        ncols = 3;  // 3 columns (floating leg, fix leg, PV)

        //// allocate mem for grfn tableau
        sprdsht = GrfnCellmatrix(nrows, ncols, GRFN_DEF_ARGBUFSZ);

        ///// column 1 contains floating leg
        for (i = 0; i < nrows; ++i)
        {
            grfn_make_swaption_floating_string(
                sprdsht[i][0].sval,
                sdp,
                start_dates[i],
                (long)(sdp->direction == FWD ? sdp->nfp : sdp->end),
                und_name,
                ref_rate_code);
        }

        ///// column 2 contains fixed leg
        for (i = 0; i < nrows; ++i)
        {
            grfn_make_swaption_fix_string(
                sprdsht[i][1].sval,
                sdp,
                start_dates[i],
                end_dates[i],
                strike,
                und_name,
                ref_rate_code);
        }

        ///// column 3 contains PV
        for (i = 0; i < nrows - 1; ++i)
        {
            if (rec_pay == SRT_PAYER)
            {
                sprintf(sprdsht[i][2].sval, "max(c[0,i] - c[1,i], PV[2])");
            }

            else
            {
                sprintf(sprdsht[i][2].sval, "max(c[1,i] - c[0,i], PV[2])");
            }
        }
        {
            //// last payoff
            if (rec_pay == SRT_PAYER)
            {
                sprintf(sprdsht[nrows - 1][2].sval, "max(c[0,i] - c[1,i], 0)");
            }

            else
            {
                sprintf(sprdsht[nrows - 1][2].sval, "max(c[1,i] - c[0,i], 0)");
            }
        }

        for (i = 0; i < nrows; ++i)
        {
            for (j = 0; j < ncols; ++j)
            {
                sprdsht[i][j].type = GRFNSCELL;
            }
        }

        srt_free(start_dates);
        srt_free(end_dates);
        srt_free(pay_dates);
        srt_free(cvg);
        break;
        ///// added by Albert Wang 08/25/03 - end

        /* The logic that determines the cap here is the same as elsewhere
          in the code, and not reproduced here */
    case CAPFLOOR:
        caplet = *sdp;
        /* Get the real dates of the cap */
        err = swp_f_make_FloatLegDatesAndCoverages(
            sdp,
            today,
            &pay_dates,
            &num_pay_dates,
            &fixing_dates,
            &start_dates,
            &end_dates,
            &cvg,
            &num_dates);

        /* Do not include rates that have set before or on today's date */
        num_caplets_skipped = 0;
        while (fixing_dates[0] <= today && num_pay_dates > 1)
        {
            num_caplets_skipped++;
            pay_dates++;
            fixing_dates++;
            start_dates++;
            end_dates++;
            cvg++;
            num_pay_dates--;
            num_dates--;
        }

        /* Sets the Grfn Tableau: one row for each caplet */
        caplet.direction = BKWD;
        caplet.nfp       = 1;
        nrows = numeventdates = num_dates;
        ncols                 = 1;
        sprdsht               = GrfnCellmatrix(nrows, ncols, GRFN_DEF_ARGBUFSZ);
        for (i = 0; i < num_dates; i++)
        {
            caplet.start             = start_dates[i];
            caplet.first_full_fixing = start_dates[i];
            caplet.end               = end_dates[i];
            err                      = grfn_make_caplet_string(
                sprdsht[i][0].sval, &caplet, strike, rec_pay, und_name, ref_rate_code);
            if (err)
                return err;
            sprdsht[i][0].type = GRFNSCELL;
        }
        start_dates -= num_caplets_skipped;
        end_dates -= num_caplets_skipped;
        pay_dates -= num_caplets_skipped;
        cvg -= num_caplets_skipped;
        srt_free(start_dates);
        srt_free(end_dates);
        srt_free(pay_dates);
        srt_free(cvg);

        break;

    case RESETCAPFLOOR:

        translate_basis(&basis_str, (Message)sdp->basis_code);

        /* Get the schedule */
        err = swp_f_make_FloatLegDatesAndCoverages(
            sdp,
            today,
            &pay_dates,
            &num_pay_dates,
            &fixing_dates,
            &start_dates,
            &end_dates,
            &cvg,
            &num_dates);

        /* DON'T INCLUDE RATES SET BEFORE TODAY
         BUT INCLUDE THOSE THAT HAVE RATE SET ON TODAY (DIFFERENT FOR CAPS) */
        num_caplets_skipped = 0;
        while (fixing_dates[0] < today && num_pay_dates > 1)
        {
            num_caplets_skipped++;
            pay_dates++;
            fixing_dates++;
            start_dates++;
            end_dates++;
            cvg++;
            num_pay_dates--;
            num_dates--;
        }

        nrows = numeventdates = num_dates;
        ncols                 = 2;
        sprdsht               = GrfnCellmatrix(nrows, ncols, GRFN_DEF_ARGBUFSZ);

        grfn_mk_resetcap_spr(
            fixing_dates,
            sprdsht,
            sdp,
            strike,
            rec_pay,
            pay_dates,
            cvg,
            nrows,
            ncols,
            basis_str,
            num_caplets_skipped);

        start_dates -= num_caplets_skipped;
        end_dates -= num_caplets_skipped;
        pay_dates -= num_caplets_skipped;
        cvg -= num_caplets_skipped;
        srt_free(start_dates);
        srt_free(end_dates);
        srt_free(pay_dates);
        srt_free(cvg);

        break;

    default:
        return serror("SwapDp_to_GrfnDeal unknown type %d", type);
        break;
    }

    *sprdshtptr       = sprdsht;
    *ncolsptr         = ncols;
    *nrowsptr         = nrows;
    *numeventdatesptr = numeventdates;
    *eventdatesptr    = fixing_dates;

    return err;
}

/* --------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------
   FUNCTION	:       grfn_make_swaption_string
   DESCRIPTION	:   Prints to swptn_str the Grfn description of a swaption
                    payoff (taking into account the reference rater name for
                                        the floating leg)
   --------------------------------------------------------------------------- */

Err grfn_make_swaption_string(
    String          swptn_str,
    SwapDP*         sdp,
    double          strike,
    SrtReceiverType rec_pay,
    String          und_name,
    String          ref_rate)
{
    long   s;
    long   e_nfp;
    long   c;
    long   b;
    String compd_str;
    String basis_str;

    s     = (long)sdp->start;
    e_nfp = (long)(sdp->direction == FWD ? sdp->nfp : sdp->end);
    c     = (long)sdp->compd;
    b     = (long)sdp->basis_code;

    /* Convert compounding from enum to string */
    translate_compounding(&compd_str, (Message)c);

    /* Convert basis_code from enum to string */
    translate_basis(&basis_str, (Message)b);

    if (rec_pay == SRT_PAYER)
    {
        sprintf(
            swptn_str,
            /* "max(df(now,%d)-df(now,%d)-%.10lf*lvl(%d,%d,\"%s\",\"%s\"),0.0)",
             */
            "max(swap(%d,%d,\"%s\",\"%s\",\"%s\",\"\",\"%s\") - %.10lf , 0.0) *  "
            "lvl(%d,%d,\"%s\",\"%s\")",
            s,
            e_nfp,
            compd_str,
            basis_str,
            und_name,
            ref_rate,
            strike,
            s,
            e_nfp,
            compd_str,
            basis_str);
    }

    else if (rec_pay == SRT_RECEIVER)
    {
        sprintf(
            swptn_str,
            /* "max(-df(now,%d)+df(now,%d)+%.10lf*lvl(%d,%d,\"%s\",\"%s\"),0.0)",
             */
            "max(-swap(%d,%d,\"%s\",\"%s\",\"%s\",\"\",\"%s\") + %.10lf , 0.0) *  "
            "lvl(%d,%d,\"%s\",\"%s\")",
            s,
            e_nfp,
            compd_str,
            basis_str,
            und_name,
            ref_rate,
            strike,
            s,
            e_nfp,
            compd_str,
            basis_str);
    }

    return NULL;
}

/* --------------------------------------------------------------------------- */
/* ---------------------------------------------------------------------------
   FUNCTION	:       grfn_make_caplet_string
   DESCRIPTION	:   Prints to caplet_str the Grfn description of a caplet
                    payoff (taking into account the reference rater name for
                                        the floating rate)
   --------------------------------------------------------------------------- */

Err grfn_make_caplet_string(
    String          caplet_str,
    SwapDP*         sdp,
    double          strike,
    SrtReceiverType rec_pay,
    String          und_name,
    String          ref_rate)
{
    long   s;
    long   e_nfp;
    long   b;
    String basis_str;

    s     = (long)sdp->start;
    e_nfp = (long)sdp->end;
    b     = (long)sdp->basis_code;

    /* Convert basis_code from enum to string */
    translate_basis(&basis_str, (Message)b);

    if (rec_pay == SRT_PAYER)
    {
        sprintf(
            caplet_str,
            "max( fra(%d, %d, \"%s\", \"%s\", \"%s\") - %.10lf , 0.0) *  cvg(%d, %d, \"%s\") * "
            "df(now, %d)",
            s,
            e_nfp,
            basis_str,
            und_name,
            ref_rate,
            strike,
            s,
            e_nfp,
            basis_str,
            e_nfp);
    }

    else if (rec_pay == SRT_RECEIVER)
    {
        sprintf(
            caplet_str,
            /* "max(-df(now,%d)+df(now,%d)+%.10lf*lvl(%d,%d,\"%s\",\"%s\"),0.0)",
             */
            "max( - fra(%d, %d, \"%s\", \"%s\", \"%s\") + %.10lf , 0.0) *  cvg(%d, %d, \"%s\") * "
            "df(now, %d)",
            s,
            e_nfp,
            basis_str,
            und_name,
            ref_rate,
            strike,
            s,
            e_nfp,
            basis_str,
            e_nfp);
    }

    return NULL;
}

/* --------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------
   FUNCTION	:       grfn_make_bondoption_string
   DESCRIPTION	:   Prints to bndopt_str the Grfn description of a bondoption
                    payoff (without taking into account the reference rater name
                                        as there is no floating leg)
   --------------------------------------------------------------------------- */

static void grfn_make_bondoption_string(
    String          bndopt_str,
    SwapDP*         sdp,
    double          strike,
    double          bk,
    double          acc_int,
    Date            pay_date,
    SrtReceiverType rec_pay)
{
    long   s;
    long   e_nfp;
    long   c;
    long   b;
    String compd_str;
    String basis_str;

    s     = (long)sdp->start;
    e_nfp = (long)(sdp->direction == FWD ? sdp->nfp : sdp->end);
    c     = (long)sdp->compd;
    b     = (long)sdp->basis_code;

    /* convert compounding from enum to string */
    translate_compounding(&compd_str, (Message)c);

    /* convert basis_code from enum to string */
    translate_basis(&basis_str, (Message)b);

    if (sdp->direction == FWD)
    {
        e_nfp = add_unit(s, e_nfp * 12 / (long)c, SRT_MONTH, NO_BUSDAY_CONVENTION);
    }

    if (acc_int == 0.0)
    {
        if (rec_pay == SRT_PAYER)
        {
            sprintf(
                bndopt_str,
                "max(%.10lf*df(now,%d)-df(now,%d)-%.10lf*lvl(%d,%d,\"%s\",\"%s\"),0.0)",
                bk,
                s,
                e_nfp,
                strike,
                s,
                e_nfp,
                compd_str,
                basis_str);
        }

        else
        {
            sprintf(
                bndopt_str,
                "max(-%.10lf*df(now,%d)+df(now,%d)+%.10lf*lvl(%d,%d,\"%s\",\"%s\"),0.0)",
                bk,
                s,
                e_nfp,
                strike,
                s,
                e_nfp,
                compd_str,
                basis_str);
        }
    }

    else
    {
        if (rec_pay == SRT_PAYER)
        {
            sprintf(
                bndopt_str,
                "max(%.10lf*df(now,%d)-df(now,%d) "
                "-%.10lf*df(now,%d)-%.10lf*lvl(%d,%d,\"%s\",\"%s\"),0.0)",
                bk + acc_int,
                s,
                e_nfp,
                acc_int,
                pay_date,
                strike,
                s,
                e_nfp,
                compd_str,
                basis_str);
        }

        else
        {
            sprintf(
                bndopt_str,
                "max(-%.10lf*df(now,%d)+df(now,%d)+%.10lf*df(now,%d)+%.10lf*lvl(%d,%d,\"%s\",\"%"
                "s\"),0.0)",
                bk + acc_int,
                s,
                e_nfp,
                acc_int,
                pay_date,
                strike,
                s,
                e_nfp,
                compd_str,
                basis_str);
        }
    }
}

/* --------------------------------------------------------------------------- */

/*****************************************************************************
   FUNCTION		:	grfn_mk_resetcap_spr
   DESCRIPTION	:	PRINT TO SPRDSHT THE GRFN DESCRIPTION OF
                                        A RESET CAP
******************************************************************************/

static void grfn_mk_resetcap_spr(
    Date*           eventdates,
    GrfnCell**      sprdsht,
    SwapDP*         sdp,
    double          strike,
    SrtReceiverType rec_pay,
    Date*           dfind,
    double*         cvg,
    long            nrows,
    long            ncols,
    char*           basis_str,
    int             num_caplets_skipped)

{
    long i;

    /* first caplet */
    sprintf(sprdsht[0][0].sval, "fra(%d,%d,\"%s\")", dfind[0], dfind[1], basis_str);
    sprdsht[0][0].type = GRFNSCELL;
    if ((strike > EPS) && (num_caplets_skipped == 0))
    {
        srt_free(sprdsht[0][0].sval);
        sprdsht[0][0].dval = strike;
        sprdsht[0][0].type = GRFNDCELL;
    }
    srt_free(sprdsht[0][1].sval);
    sprdsht[0][1].type = GRFNBCELL;

    /* next caplets */
    for (i = 1; i < nrows; i++)
    {
        sprintf(sprdsht[i][0].sval, "fra(%d,%d,\"%s\")", dfind[i], dfind[i + 1], basis_str);
        sprdsht[i][0].type = GRFNSCELL;
        switch (rec_pay)
        {
        case SRT_PAYER:
            sprintf(sprdsht[i][1].sval, "max(c[0,i]-c[0,i-1],0)*df(d[i],%d)*%g", dfind[i], cvg[i]);
            break;
        case SRT_RECEIVER:
            sprintf(sprdsht[i][1].sval, "max(-c[0,i]+c[0,i-1],0)*df(d[i],%d)*%g", dfind[i], cvg[i]);
            break;
        case SRT_STRADDLE:
            sprintf(sprdsht[i][1].sval, "abs(c[0,i]-c[0,i-1])*df(d[i],%d)*%g", dfind[i], cvg[i]);
            break;
        default:
            break;
        }
        sprdsht[i][1].type = GRFNSCELL;
    }
}

/* ======================================================================= */

/*****************************************************************************
   FUNCTION	: grfn_SwapDP_to_BDT_GrfnCells
   DESCRIPTION	: GIVEN A SWAPDP (SWAP DATE PARAM) AND STRIKE AND TYPE,
        GENERATE THE GRFNCELLS AND EVENT DATES
        CORRESPONDING TO A CAP OR SWAPTION



   AMENDMENTS:
        Reference	:
        Author          :
        Date            :
        Description     :

******************************************************************************/

Err grfn_SwapDP_to_BDT_GrfnCells(
    long*           numeventdatesptr,
    Date**          eventdatesptr,
    long*           nrowsptr,
    long*           ncolsptr,
    GrfnCell***     sprdshtptr,
    SrtUndPtr       und,
    SwapDP*         sdp,
    double          strike,
    double          bk,
    SrtReceiverType rec_pay,
    StructType      t)
{
    long       numeventdates;
    Date*      eventdates;
    Date*      fixing_dates;
    Date*      start_dates;
    Date*      end_dates;
    Date*      pay_dates;
    int        num_dates;
    int        num_pay_dates;
    Date       last_fixing_date, today, first_pay_date;
    long       nrows, ncols;
    GrfnCell** sprdsht;
    int        num_caplets_skipped;
    int        i, istart;
    Err        err        = NULL;
    double     accrue_int = 0.0;
    double*    cvg;

    today = get_today_from_underlying(und);

    switch (t)
    {
    case BOND_OPTION:
        ncols = 1;

        /* Generates first the theoretical swap dates (dates of payments...)*/
        /* dfind[0] correponds to sdp->start, dfind[num_pay_dates-1] to sdp->end */
        swp_f_make_FixedLegDatesAndCoverages(
            sdp, today, &pay_dates, &num_pay_dates, &start_dates, &end_dates, &cvg, &num_dates);

        fixing_dates    = (Date*)malloc(1 * sizeof(Date));
        fixing_dates[0] = add_unit(start_dates[0], -sdp->spot_lag, SRT_BDAY, MODIFIED_SUCCEEDING);

        /* Need to calculate correct accrued interest which will be subtracted
           from the bond strike */
        accrue_int = acc_int_fct(*sdp, strike);
        if (accrue_int != 0)
            first_pay_date = pay_dates[1];

        /* If spot_lag = 0: fixing= start => no need to generate two different dates*/
        if (sdp->spot_lag = 0)
        {
            istart = 0;
            nrows = numeventdates = num_pay_dates;
        }
        else /* Generate sdp dates + fixing date */
        {
            istart = 1;
            nrows = numeventdates = num_pay_dates + 1;
        }
        eventdates    = (Date*)srt_calloc(numeventdates, sizeof(Date));
        eventdates[0] = fixing_dates[0];
        for (i = 0; i < num_dates; i++)
            eventdates[i + istart] = fixing_dates[i];

        sprdsht = GrfnCellmatrix(nrows, ncols, GRFN_DEF_ARGBUFSZ);

        grfn_mk_bdt_bnd_opt_sht(
            &sprdsht, sdp, strike, bk, accrue_int, rec_pay, eventdates, cvg, numeventdates, istart);

        srt_free(pay_dates);
        srt_free(fixing_dates);
        srt_free(start_dates);
        srt_free(cvg);
        break;

    case SWAPTION:
        ncols      = 1;
        accrue_int = 0.0;

        /* Generates first the theoretical swap dates (dates of payments...)*/
        /* dfind[0] correponds to sdp->start, dfind[num_pay_dates-1] to sdp->end */
        swp_f_make_FixedLegDatesAndCoverages(
            sdp, today, &pay_dates, &num_pay_dates, &start_dates, &end_dates, &cvg, &num_dates);

        fixing_dates    = (Date*)malloc(1 * sizeof(Date));
        fixing_dates[0] = add_unit(start_dates[0], -sdp->spot_lag, SRT_BDAY, MODIFIED_SUCCEEDING);

        /* If spot_lag = 0: fixing= start => no need to generate two different dates*/
        if (sdp->spot_lag = 0)
        {
            istart = 0;
            nrows = numeventdates = num_pay_dates;
        }
        else /* Generate sdp dates + fixing date */
        {
            istart = 1;
            nrows = numeventdates = num_pay_dates + 1;
        }
        eventdates    = (Date*)srt_calloc(numeventdates, sizeof(Date));
        eventdates[0] = fixing_dates[0];
        for (i = 0; i < num_pay_dates; i++)
            eventdates[i + istart] = fixing_dates[i];

        sprdsht = GrfnCellmatrix(nrows, ncols, GRFN_DEF_ARGBUFSZ);

        grfn_mk_bdt_bnd_opt_sht(
            &sprdsht, sdp, strike, bk, accrue_int, rec_pay, eventdates, cvg, numeventdates, istart);

        srt_free(pay_dates);
        srt_free(fixing_dates);
        srt_free(start_dates);
        srt_free(cvg);
        break;

    case CAPFLOOR:
        /* Generates first the theoretical swap dates (dates of payments...)*/
        /* dfind[0] correponds to sdp->start, dfind[num_pay_dates-1] to sdp->end */
        err = swp_f_make_FloatLegDatesAndCoverages(
            sdp,
            today,
            &pay_dates,
            &num_pay_dates,
            &fixing_dates,
            &start_dates,
            &end_dates,
            &cvg,
            &num_dates);

        /**** DON'T INCLUDE CAPLETS THAT HAVE RATE SETS ON OR BEFORE TODAY ***/
        num_caplets_skipped = 0;
        today               = get_today_from_underlying(und);
        while (fixing_dates[0] <= today && num_pay_dates > 1)
        {
            pay_dates++;
            fixing_dates++;
            start_dates++;
            end_dates++;
            cvg++;
            num_pay_dates--;
            num_dates--;
        }

        /* As many columns as caplets + 1 : number of remining dates - 1 */
        ncols = num_pay_dates;

        /* As many rows as fixings + payment dates */
        if (sdp->spot_lag = 0)
        {
            istart = 0;
            nrows = numeventdates = num_pay_dates;
        }
        else /* Generate sdp dates + fixing dates */
        {
            istart = 1;
            nrows = numeventdates = 2 * num_pay_dates - 1;
        }
        eventdates = (Date*)srt_calloc(numeventdates, sizeof(Date));
        sprdsht    = GrfnCellmatrix(nrows, ncols, GRFN_DEF_ARGBUFSZ);

        for (i = 0; i < num_dates - 1; i++)
        {
            if (start_dates[i] != fixing_dates[i])
            {
                eventdates[2 * i]     = fixing_dates[i];
                eventdates[2 * i + 1] = start_dates[i];
            }
            else
            {
                eventdates[i] = fixing_dates[i];
            }
        }
        last_fixing_date              = fixing_dates[num_dates - 1];
        eventdates[numeventdates - 1] = fixing_dates[num_pay_dates - 1];

        sprdsht = GrfnCellmatrix(nrows, ncols, GRFN_DEF_ARGBUFSZ);

        grfn_mk_bdt_cap_flr_sht(
            &sprdsht,
            sdp,
            strike,
            bk,
            accrue_int,
            rec_pay,
            eventdates,
            cvg,
            numeventdates,
            num_pay_dates - 1,
            istart);

        start_dates -= num_caplets_skipped;
        end_dates -= num_caplets_skipped;
        pay_dates -= num_caplets_skipped;
        cvg -= num_caplets_skipped;
        srt_free(start_dates);
        srt_free(end_dates);
        srt_free(pay_dates);
        srt_free(cvg);

        break;

    case RESETCAPFLOOR:
        return serror("Cannot price Reset Cap within BDT");

    default:
        return serror("SwapDp_to_GrfnDeal unknown type %d", t);
        break;
    }

    *sprdshtptr       = sprdsht;
    *ncolsptr         = ncols;
    *nrowsptr         = nrows;
    *numeventdatesptr = numeventdates;
    *eventdatesptr    = fixing_dates;

    return err;
}

/*****************************************************************************
   FUNCTION	: grfn_mk_bdt_bnd_opt_sht
   DESCRIPTION	: PRINT TO SWPTN_STR THE GRFN - BDT DESCRIPTION OF
                  A BOND OPTION PAYOFF

                  when *cvg is passed:
                      => cvg[0] = 0.0000;
                      => cvg[1] corresponds to cvg between
                                date[i_start]
                            and
                                date[i_start+1]
******************************************************************************/

static void grfn_mk_bdt_bnd_opt_sht(
    GrfnCell***     sprdshtptr,
    SwapDP*         sdp,
    double          strike,
    double          bk,
    double          acc_int,
    SrtReceiverType rec_pay,
    Date*           eventdatesptr,
    double*         cvg,
    int             numeventdates,
    int             i_start)
{
    int        i, rp;
    double     cash_flow;
    GrfnCell** sprdsht = *sprdshtptr;

    rp = (rec_pay == SRT_PAYER) ? 1 : -1;

    /* Generates coupons of the bond */
    for (i = i_start + 1; i < numeventdates; i++)
    {
        cash_flow = cvg[i - i_start] * strike;
        if (i == numeventdates - 1)
            cash_flow += 1.00;
        /* If on first coupon date after exrcise, need to add accrued
           before exercise as cvg only takes accrued after exercise*/
        if ((i == i_start + 1) && (acc_int != 0.0))
            cash_flow += acc_int;

        sprdsht[i][0].dval = rp * cash_flow;
        sprdsht[i][0].type = GRFNDCELL;
    }

    /* Strike of the bond option (option on clean price)*/
    cash_flow = bk + acc_int;
    /* If fixing and start are different: two different dates
            => can hard write the strike*/
    if (i_start)
    {
        sprdsht[1][0].dval = -rp * cash_flow;
        sprdsht[1][0].type = GRFNDCELL;

        sprintf(sprdsht[0][0].sval, "max( PV[0] , 0.0 )");
        sprdsht[0][0].type = GRFNSCELL;
    }
    else
    {
        if (rec_pay == SRT_PAYER)
        {
            sprintf(sprdsht[0][0].sval, "max( PV[0] - %.10lf , 0.0 )", cash_flow);
        }
        else if (rec_pay == SRT_RECEIVER)
        {
            sprintf(sprdsht[0][0].sval, "max( %.10lf - PV[0] , 0.0 )", cash_flow);
        }
        sprdsht[0][0].type = GRFNSCELL;
    }
}

/*****************************************************************************
   FUNCTION	: grfn_mk_bdt_cap_flr_sht
   DESCRIPTION	: PRINT TO SWPTN_STR THE GRFN - BDT DESCRIPTION OF
                  A BOND OPTION PAYOFF

                  when *cvg is passed:
                      => cvg[0] = 0.0000 or anything to do with the past dates
                      => cvg[1] corresponds to cvg between
                                date[first_fixing after today]
                            and
                                date[next fixing]

                  istart correpsonds to the row number of the first fixing
                  date after today
******************************************************************************/

static void grfn_mk_bdt_cap_flr_sht(
    GrfnCell***     sprdshtptr,
    SwapDP*         sdp,
    double          strike,
    double          bk,
    double          acc_int,
    SrtReceiverType rec_pay,
    Date*           eventdatesptr,
    double*         cvg,
    int             numeventdates,
    int             numcaplets,
    int             istart)
{
    int        col, row, j, rp, increment;
    double     cash_flow;
    GrfnCell** sprdsht = *sprdshtptr;

    rp        = (rec_pay == SRT_PAYER) ? -1 : 1;
    increment = istart + 1;

    col = 0;
    row = istart;
    for (j = 0; j < numcaplets; j++)
    {
        /* Initial exchange at fixing date :
           pays 100 if SRT_RECEIVER, receives 100 if SRT_PAYER */
        cash_flow              = 1.00;
        sprdsht[row][col].dval = -rp * cash_flow;
        sprdsht[row][col].type = GRFNDCELL;

        /* Moves down to next fixing date
           (there might be a spot lag date in between if istart == 1,
            but not for the last caplet, as there is no following  caplet*/
        if (j == numcaplets - 1)
            row += 1; /* No interfering value date */
        else
            row += increment;

        /* Final exchange at payment date = next fixing date:
           receives 100+K*cvg if SRT_RECEIVER, pays 100+K*cvg if SRT_PAYER */
        cash_flow              = 1.00 + cvg[j + 1] * strike;
        sprdsht[row][col].dval = rp * cash_flow;
        sprdsht[row][col].type = GRFNDCELL;

        /* Prepares following caplet : same row (fixing date - next column */
        col += 1;
    }

    col = numcaplets;
    row = 0;
    for (j = 0; j < numcaplets; j++)
    {
        if (istart == 0)
        {
            sprintf(sprdsht[row][col].sval, "max( PV[0] + c[%d, i], 0.0 )", j);
        }
        else if (istart == 1)
        {
            sprintf(sprdsht[row][col].sval, "max( PV[0] , 0.0 )");
        }
        sprdsht[row][col].type = GRFNSCELL;
        row += increment;
    }
}

/* --------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------
   FUNCTION	:       srt_f_grfn_newclsdfrm
   DESCRIPTION	:   Prices a swaption, a cap or a bond option using Grfn
                    If LGM can be used, a switch send the calculation to a
                                        real closed form
   --------------------------------------------------------------------------- */

Err srt_f_grfn_newclsdfrm(
    SrtUndPtr        und,
    SrtGrfnParam*    grfnparam,
    int              lNumInstr,
    SwapDP*          sdp,
    double*          strike,
    double*          bk,
    SrtReceiverType* rec_pay,
    StructType*      type,
    String*          ref_rate_code,
    double*          answer)
{
    Err          err = NULL;
    long         ncols, nrows, numeventdates;
    Date*        eventdates = NULL;
    GrfnCell**   sprdsht    = NULL;
    SrtIOStruct* iolist;
    SrtMdlType   mdl_type;
    SrtMdlDim    mdl_dim;
    String       und_name;
    Date         today;
    int          i;

    /* Gets the Model type (LGM, Cheyette,...) */
    err = get_underlying_mdltype(und, &mdl_type);

    /* Gets the Number of Factors in the Model */
    err = get_underlying_mdldim(und, &mdl_dim);

    /* IF LGM or EtaBeta, there is a REAL closed form: go to it */
    if ((mdl_type == ETABETA || mdl_type == LGM) || (grfnparam->closed_form_type == QUICK_CLSDFRM))
    {
        for (i = 0; i < lNumInstr; i++)
        {
            err = srt_f_closed_form(
                und,
                mdl_type,
                mdl_dim,
                &(sdp[i]),
                strike[i],
                bk[i],
                rec_pay[i],
                type[i],
                ref_rate_code[i],
                &(answer[i]));

            if (err)
                return err;
        }
    }
    else
    {
        /* There is no closed form available: have to build a full grfn tableau */
        if (mdl_type == BDT)
        {
            return serror("Calibration failed: BDT model is not avaible");
        }
        else
        {
            und_name = get_underlying_name(und);
            today    = get_today_from_underlying(und);
            err      = grfn_SwapDParray_to_GrfnCells(
                &numeventdates,
                &eventdates,
                &nrows,
                &ncols,
                &sprdsht,
                today,
                lNumInstr,
                sdp,
                strike,
                bk,
                rec_pay,
                type,
                und_name,
                ref_rate_code);
            if (err)
                return err;
            if (nrows == 0)
                return serror("Product has expired in srt_f_grfn_clsdfrm");
        }

        err = srt_f_IOstructcreate(&iolist, "clsdfrm");

        /* Call Grfn with this tableau */
        if (!err)
            err = srt_f_grfn(
                und,
                grfnparam,
                numeventdates,
                &eventdates,
                &nrows,
                &ncols,
                &sprdsht,
                0,
                0,
                0,
                0,
                0,
                iolist,
                0,
                0);

        if (!err)
            err = srt_f_IOstructgetpremiumval(*iolist, answer);

        if (!err)
            err = srt_f_IOstructfree(&iolist);

        if (eventdates)
            srt_free(eventdates);

        if (sprdsht)
        {
            grfn_free_GrfnCellmatrix(sprdsht, nrows, ncols);
        }
    }

    return err;

} /* END srt_f_grfn_newclsdfrm(...) */

/* ---------------------------------------------------------------------------
   FUNCTION	:       grfn_SwapDP_to_GrfnCells
   DESCRIPTION	:   Given a SwapDP (swaps date parameters), a strike and a type,
                    generate the GrfnCells and EventDates corresponding to a
                                        cap, a swaption or a bond option.
   --------------------------------------------------------------------------- */

Err grfn_SwapDParray_to_GrfnCells(
    long*            numeventdatesptr,
    Date**           eventdatesptr,
    long*            nrowsptr,
    long*            ncolsptr,
    GrfnCell***      sprdshtptr,
    Date             today,
    int              lNumInstr,
    SwapDP*          sdp,
    double*          strike,
    double*          bk,
    SrtReceiverType* rec_pay,
    StructType*      type,
    String           und_name,
    String*          ref_rate_code)
{
    Date*      fixing_dates;
    Date*      start_dates;
    Date*      end_dates;
    Date*      pay_dates;
    int        num_dates;
    int        num_pay_dates;
    GrfnCell** sprdsht;
    int        i, j;
    Err        err        = NULL;
    double     accrue_int = 0.0;
    double*    cvg;
    int        SprdshtRows;
    int        YesNo;
    Date*      FixingDatesForSort = NULL;
    Date*      RealFixingDates    = NULL;
    Date*      EventDates;

    /*----------------------*/

    FixingDatesForSort = (Date*)malloc(lNumInstr * sizeof(Date));
    RealFixingDates    = (Date*)malloc(lNumInstr * sizeof(Date));

    for (i = 0; i < lNumInstr; i++)
    {
        if (type[i] == CAPFLOOR)
        {
            err = swp_f_make_FloatLegDatesAndCoverages(
                &(sdp[i]),
                today,
                &pay_dates,
                &num_pay_dates,
                &fixing_dates,
                &start_dates,
                &end_dates,
                &cvg,
                &num_dates);
            srt_free(start_dates);
            srt_free(end_dates);
            srt_free(pay_dates);
            srt_free(cvg);

            FixingDatesForSort[i] = fixing_dates[0];
            srt_free(fixing_dates);
        }
        else
        {
            swp_f_make_FixedLegDatesAndCoverages(
                &(sdp[i]),
                today,
                &pay_dates,
                &num_pay_dates,
                &start_dates,
                &end_dates,
                &cvg,
                &num_dates);

            FixingDatesForSort[i] =
                add_unit(start_dates[0], -sdp[i].spot_lag, SRT_BDAY, MODIFIED_SUCCEEDING);
            srt_free(start_dates);
            srt_free(end_dates);
            srt_free(pay_dates);
            srt_free(cvg);
        }
    }

    SprdshtRows = 0;
    for (i = 0; i < lNumInstr; i++)
    {
        RealFixingDates[SprdshtRows] = FixingDatesForSort[i];
        YesNo                        = 0;
        for (j = i + 1; j < lNumInstr; j++)
        {
            if (FixingDatesForSort[i] == FixingDatesForSort[j])
                YesNo = 1;
        }
        if (!(YesNo))
            SprdshtRows++;
    }

    if (SprdshtRows <= 0)
        return serror("No fixing date in srt_f_newgrfclsdfrm() ");

    EventDates = (Date*)malloc(SprdshtRows * sizeof(Date));
    for (i = 0; i < SprdshtRows; i++)
    {
        EventDates[i] = RealFixingDates[i];
    }
    if (RealFixingDates)
        srt_free(RealFixingDates);

    sprdsht = GrfnCellmatrix(SprdshtRows, lNumInstr, GRFN_DEF_ARGBUFSZ);

    for (j = 0; j < SprdshtRows; j++)
    {
        for (i = 0; i < lNumInstr; i++)
        {
            if (FixingDatesForSort[i] == EventDates[j])
            {
                switch (type[i])
                {
                case SWAPTION:
                    err = grfn_make_swaption_string(
                        sprdsht[j][i].sval,
                        &(sdp[i]),
                        strike[i],
                        rec_pay[i],
                        und_name,
                        ref_rate_code[i]);
                    break;

                    /* The logic that determines the cap here is the same as elsewhere
                      in the code, and not reproduced here */
                case CAPFLOOR:
                    err = grfn_make_caplet_string(
                        sprdsht[j][i].sval,
                        &(sdp[i]),
                        strike[i],
                        rec_pay[i],
                        und_name,
                        ref_rate_code[i]);
                    break;

                default:
                    return serror("SwapDpArray_to_GrfnDeal unknown type %d", type);
                    break;
                }
            }
            else
            {
                sprintf(sprdsht[j][i].sval, "0.0");
            }

            sprdsht[j][i].type = GRFNSCELL;
        }
    }

    *sprdshtptr       = sprdsht;
    *ncolsptr         = lNumInstr;
    *nrowsptr         = SprdshtRows;
    *numeventdatesptr = SprdshtRows;
    *eventdatesptr    = EventDates;

    if (FixingDatesForSort)
        srt_free(FixingDatesForSort);

    return err;
}

/* ---------------------------------------------------------------------------
   FUNCTION	:       new_grfn_SwapDP_to_GrfnCells
   DESCRIPTION	:   Given a SwapDP (swaps date parameters), a strike and a type,
                    generate the GrfnCells and EventDates corresponding to a
                                        cap, a swaption or a bond option - modification to include
                                        the cap.
   --------------------------------------------------------------------------- */

static int longcmp(const void* vp, const void* vq)
{
    const long* p = vp;
    const long* q = vq;

    if (*p > *q)
        return (+1);
    else if (*p < *q)
        return (-1);
    else
        return (0);
}

Err new_grfn_SwapDParray_to_GrfnCells(
    long*            numeventdatesptr,
    Date**           eventdatesptr,
    long*            nrowsptr,
    long*            ncolsptr,
    GrfnCell***      sprdshtptr,
    Date             today,
    int              lNumInstr,
    SwapDP*          sdp,
    double*          strike,
    double*          bk,
    SrtReceiverType* rec_pay,
    StructType*      type,
    String           und_name,
    String*          ref_rate_code)
{
    Date*      fixing_dates;
    Date*      start_dates;
    Date*      end_dates;
    Date*      pay_dates;
    int        num_dates;
    int        num_pay_dates;
    GrfnCell** sprdsht;
    int        i, j, SprdshtRows, nn;
    Err        err        = NULL;
    double     accrue_int = 0.0;
    double*    cvg;
    Date **    FixingDatesForSort, **startDatesForSort, **endDatesForSort;
    Date*      RealFixingDates = NULL;
    Date*      EventDates      = NULL;
    int*       FixingDatesNum  = NULL;
    int        maxNumDates;
    int*       index = NULL;
    /* FILE		*out;*/

    /*----------------------*/

    FixingDatesForSort = (Date**)malloc(lNumInstr * sizeof(Date*));
    startDatesForSort  = (Date**)malloc(lNumInstr * sizeof(Date*));
    endDatesForSort    = (Date**)malloc(lNumInstr * sizeof(Date*));
    FixingDatesNum     = (int*)malloc(lNumInstr * sizeof(int));
    index              = (int*)malloc(lNumInstr * sizeof(int));
    maxNumDates        = 0;

    for (i = 0; i < lNumInstr; i++)
    {
        /* CAPFLOOR */
        if (type[i] == CAPFLOOR)
        {
            err = swp_f_make_FloatLegDatesAndCoverages(
                &(sdp[i]),
                today,
                &pay_dates,
                &num_pay_dates,
                &fixing_dates,
                &start_dates,
                &end_dates,
                &cvg,
                &num_dates);

            FixingDatesForSort[i] = fixing_dates; /* add an array of date */
            startDatesForSort[i]  = start_dates;
            endDatesForSort[i]    = end_dates;
            FixingDatesNum[i]     = num_dates;
        }
        else /* SWAPTION */
        {
            swp_f_make_FixedLegDatesAndCoverages(
                &(sdp[i]),
                today,
                &pay_dates,
                &num_pay_dates,
                &start_dates,
                &end_dates,
                &cvg,
                &num_dates);

            FixingDatesForSort[i] = (Date*)malloc(sizeof(Date));
            FixingDatesForSort[i][0] =
                add_unit(start_dates[0], -sdp[i].spot_lag, SRT_BDAY, MODIFIED_SUCCEEDING);
            FixingDatesNum[i] = 1;
        }

        maxNumDates += FixingDatesNum[i];
        srt_free(pay_dates);
        srt_free(cvg);
    }

    /* fusion and sort of dates for events dates list */
    RealFixingDates = (Date*)malloc(maxNumDates * sizeof(Date));
    nn              = 0;
    for (i = 0; i < lNumInstr; i++)
        for (j = 0; j < FixingDatesNum[i]; j++)
        {
            RealFixingDates[nn] = FixingDatesForSort[i][j];
            nn++;
        }
    /* Sort the dates */
    qsort(RealFixingDates, maxNumDates, sizeof(Date), longcmp);

    /* put (unique) dates in events dates*/
    SprdshtRows = i = 0;
    while (i < maxNumDates)
    {
        RealFixingDates[SprdshtRows] = RealFixingDates[i];
        SprdshtRows++;
        i++;
        while ((i < maxNumDates) && (RealFixingDates[SprdshtRows - 1] == RealFixingDates[i]))
            i++;
    }

    /* transfer in event dates */
    EventDates = (Date*)malloc(SprdshtRows * sizeof(Date));
    for (i = 0; i < SprdshtRows; i++)
    {
        EventDates[i] = RealFixingDates[i];
    }
    if (RealFixingDates)
        srt_free(RealFixingDates);

    /* build the final grfn tableau with the right grfn syntax for caplets and swaption*/
    sprdsht = GrfnCellmatrix(SprdshtRows, lNumInstr, GRFN_DEF_ARGBUFSZ);

    for (i = 0; i < lNumInstr; i++)
        index[i] = 0;

    for (j = 0; j < SprdshtRows; j++)
    {
        for (i = 0; i < lNumInstr; i++)
        {
            if ((index[i] < FixingDatesNum[i]) &&
                (FixingDatesForSort[i][index[i]] == EventDates[j]))
            {
                switch (type[i])
                {
                case SWAPTION:
                    err = grfn_make_swaption_string(
                        sprdsht[j][i].sval,
                        &(sdp[i]),
                        strike[i],
                        rec_pay[i],
                        und_name,
                        ref_rate_code[i]);
                    break;

                    /* The logic that determines the cap here is the same as elsewhere
                      in the code, and not reproduced here */
                case CAPFLOOR:
                    sdp[i].start = startDatesForSort[i][index[i]];
                    sdp[i].end   = endDatesForSort[i][index[i]];
                    err          = grfn_make_caplet_string(
                        sprdsht[j][i].sval,
                        &(sdp[i]),
                        strike[i],
                        rec_pay[i],
                        und_name,
                        ref_rate_code[i]);
                    break;

                default:
                    return serror("SwapDpArray_to_GrfnDeal unknown type %d", type);
                    break;
                }
                index[i]++;
            }
            else
            {
                sprintf(sprdsht[j][i].sval, "0.0");
            }

            sprdsht[j][i].type = GRFNSCELL;
        }
    }

    *sprdshtptr       = sprdsht;
    *ncolsptr         = lNumInstr;
    *nrowsptr         = SprdshtRows;
    *numeventdatesptr = SprdshtRows;
    *eventdatesptr    = EventDates;

    for (i = 0; i < lNumInstr; i++)
    {
        if (FixingDatesForSort[i])
            srt_free(FixingDatesForSort[i]);
        if (type[i] == CAPFLOOR)
        {
            if (startDatesForSort[i])
                srt_free(startDatesForSort[i]);
            if (endDatesForSort[i])
                srt_free(endDatesForSort[i]);
        }
    }
    if (FixingDatesForSort)
        srt_free(FixingDatesForSort);
    if (type[i] == CAPFLOOR)
    {
        if (startDatesForSort)
            srt_free(startDatesForSort);
        if (endDatesForSort)
            srt_free(endDatesForSort);
    }

    if (FixingDatesNum)
        srt_free(FixingDatesNum);
    if (index)
        srt_free(index);

    /*
        out = fopen("testsprsht.res","w");
        for(j = 0; j < SprdshtRows; j++)
        {
                for(i = 0; i < lNumInstr ; i++)
                        fprintf(out, "=== %s ", sprdsht[j][i].sval);
                fprintf(out, "\n");
        }
    fclose(out);
        */

    return err;
}

/* -------------------------------------------------------------------------- */

///// added by Albert Wang 08/25/03 - begin

/* ---------------------------------------------------------------------------
   FUNCTION	:       grfn_make_swaption_floating_string
   DESCRIPTION	:   Prints to swptn_str the Grfn description of a swaption floating leg
                    payoff (taking into account the reference rate name for
                                        the floating leg)
   --------------------------------------------------------------------------- */

void grfn_make_swaption_floating_string(
    String swptn_str, SwapDP* sdp, long start, long end, String und_name, String ref_rate)
{
    long   c;
    long   b;
    String compd_str;
    String basis_str;

    c = (long)sdp->compd;
    b = (long)sdp->basis_code;

    /* Convert compounding from enum to string */
    translate_compounding(&compd_str, (Message)c);

    /* Convert basis_code from enum to string */
    translate_basis(&basis_str, (Message)b);

    sprintf(
        swptn_str,
        "swap(%d,%d,\"%s\",\"%s\",\"%s\",\"\",\"%s\") *  lvl(%d,%d,\"%s\",\"%s\",\"%s\")",
        start,
        end,
        compd_str,
        basis_str,
        und_name,
        ref_rate,
        start,
        end,
        compd_str,
        basis_str,
        und_name);
}

/* ---------------------------------------------------------------------------
   FUNCTION	:       grfn_make_swaption_fix_string
   DESCRIPTION	:   Prints to swptn_str the Grfn description of a swaption fix leg
                    payoff (taking into account the reference rate name for
                                        the floating leg)
   --------------------------------------------------------------------------- */

void grfn_make_swaption_fix_string(
    String  swptn_str,
    SwapDP* sdp,
    long    start,
    long    end,
    double  strike,
    String  und_name,
    String  ref_rate)
{
    long   c;
    long   b;
    String compd_str;
    String basis_str;

    c = (long)sdp->compd;
    b = (long)sdp->basis_code;

    /* Convert compounding from enum to string */
    translate_compounding(&compd_str, (Message)c);

    /* Convert basis_code from enum to string */
    translate_basis(&basis_str, (Message)b);

    sprintf(
        swptn_str,
        "%.10lf  *  lvl(%d,%d,\"%s\",\"%s\") + PV[1]",
        strike,
        start,
        end,
        compd_str,
        basis_str);
}

///// added by Albert Wang 08/25/03 - end
