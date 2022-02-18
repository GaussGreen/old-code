/******************************************************************************
  SYSTEM	:SRT
  SUB-SYSTEM  	:GRFN
  MODULE NAME	:SRT_F_GRFAMESWP.C

  AUTHOR      	:Jasbir Malhi
  CREATED	:august 1994
  VERSION     	:1.00
  DESCRIPTION 	:Functions for GRFN closed form AMERICAN SWAPTION


  AMENDMENTS
        Reference	:
        Author		: K L Chau
        Date		: Dec 15, 1994
        Description	: use new srtmkt structure

******************************************************************************/

#include "grf_h_public.h"
#include "srt_h_all.h"
#include "srt_h_grfclsdfrm.h"

static void grfn_str0(
    String          opt_str,
    int             nfp,
    String          compd_str,
    String          basis_str,
    double          strike,
    SrtReceiverType rec_pay)
{
    if (rec_pay == SRT_RECEIVER)
        sprintf(
            opt_str,
            "AM(max((%f-swap(now,%d,\"%s\",\"%s\"))*lvl(now,%d,\"%s\",\"%s\"),PV[0]))",
            strike,
            nfp,
            compd_str,
            basis_str,
            nfp,
            compd_str,
            basis_str);
    else if (rec_pay == SRT_PAYER)
        sprintf(
            opt_str,
            "AM(max((swap(now,%d,\"%s\",\"%s\")-%f)*lvl(now,%d,\"%s\",\"%s\"),PV[0]))",
            nfp,
            compd_str,
            basis_str,
            strike,
            nfp,
            compd_str,
            basis_str);
}

static void grfn_str1(
    String          opt_str,
    int             nfp,
    String          compd_str,
    String          basis_str,
    double          strike,
    SrtReceiverType rec_pay)
{
    if (rec_pay == SRT_RECEIVER)
        sprintf(
            opt_str,
            "max((%f-swap(now,%d,\"%s\",\"%s\"))*lvl(now,%d,\"%s\",\"%s\"),0)",
            strike,
            nfp,
            compd_str,
            basis_str,
            nfp,
            compd_str,
            basis_str);
    else if (rec_pay == SRT_PAYER)
        sprintf(
            opt_str,
            "max((swap(now,%d,\"%s\",\"%s\")-%f)*lvl(now,%d,\"%s\",\"%s\"),0)",
            nfp,
            compd_str,
            basis_str,
            strike,
            nfp,
            compd_str,
            basis_str);
}

/*****************************************************************************
   FUNCTION		: grfn_mk_ameswp_GrfnCells
   DESCRIPTION		: generate the Grfn Tableau describing
                        an American swaption.

   AMENDMENTS		:
        Reference	:
        Author          :
        Date            :
        Description     :

******************************************************************************/

static Err grfn_mk_ameswp_GrfnCells(
    double          start,
    double          end,
    int             nfp,
    int             delay,
    int             compd,
    int             basis,
    double          strike,
    SrtReceiverType rec_pay,
    Date            clcn_date,
    long*           num_eventdates,
    Date**          eventdates,
    long*           nrows,
    long*           ncols,
    GrfnCell***     sprdsht)
{
    Err    err = NULL;
    String compd_str;
    String basis_str;

    /*** allocate space ***/ /* careful if more dates added*/
    *ncols = 1;
    *nrows = *num_eventdates = 2;
    *eventdates              = (Date*)srt_calloc(2, sizeof(Date));
    *sprdsht                 = GrfnCellmatrix(*nrows, *ncols, GRFN_DEF_ARGBUFSZ);

    /*** set grfn event dates ***/

    (*eventdates)[0] = (long)start;
    (*eventdates)[1] = (long)end;

    /* convert compounding from enum to string */
    translate_compounding(&compd_str, (Message)compd);

    /* convert basis from enum to string */
    translate_basis(&basis_str, (Message)basis);

    /*** enter payoff string0 for each date ***/
    grfn_str0((*sprdsht)[0][0].sval, nfp, compd_str, basis_str, strike, rec_pay);
    (*sprdsht)[0][0].type = GRFNSCELL;
    if (err)
        return err;

    /*** enter payoff string1 for each date ***/

    grfn_str1((*sprdsht)[1][0].sval, nfp, compd_str, basis_str, strike, rec_pay);
    (*sprdsht)[1][0].type = GRFNSCELL;
    if (err)
        return err;

    return err;
}

/*****************************************************************************
   FUNCTION		: grfn_ameswp
   DESCRIPTION		: PRICE AN AMERICAN SWAPTION USING GRFN.

   AMENDMENTS		:
        Reference	:
        Author          :
        Date            :
        Description     :

******************************************************************************/
Err grf_ameswp_clsdfrm(
    double          start,
    double          end,
    int             nfp,
    int             delay,
    int             compd,
    int             basis,
    double*         strikes,
    SrtReceiverType rec_pay,   /* pay or rec */
    SrtUndPtr       und,       /* name of underlying to use */
    SrtGrfnParam*   grfnparam, /* model and implementation details */
    double*         answer     /* value returned */
)
{
    Err        err            = NULL;
    Err        err1           = NULL;
    long       ncols          = 0;
    long       nrows          = 0;
    long       num_eventdates = 0;
    Date*      eventdates     = NULL;
    GrfnCell** sprdsht        = NULL;

    SrtIOStruct* iolist; /* list of requests for prices */
    Date         today;

    /* Create a list with a request for the price */
    err = srt_f_IOstructcreate(&iolist, "");
    if (err)
        return (err);

    today = get_today_from_underlying(und);

    /* CONSTRUCT A GRFN SPREADSHEET THAT DESCRIBES THE AMERICAN SWAPTION */
    err = grfn_mk_ameswp_GrfnCells(
        start,
        end,
        nfp,
        delay,
        compd,
        basis,
        strikes[0],
        rec_pay,
        today,
        &num_eventdates,
        &eventdates,
        &nrows,
        &ncols,
        &sprdsht);

    if (!err)
    {
        err = srt_f_grfn(
            und,
            grfnparam,
            num_eventdates,
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

        /* Take the price from the list */
        if (!err)
            err = srt_f_IOstructgetpremiumval(*iolist, answer);
    }
    if (eventdates)
        srt_free(eventdates);
    if (sprdsht)
    {
        grfn_free_GrfnCellmatrix(sprdsht, nrows, ncols);
    }

    err1 = srt_f_IOstructfree(&iolist);
    if (err1)
        return err1;

    return err;
}
