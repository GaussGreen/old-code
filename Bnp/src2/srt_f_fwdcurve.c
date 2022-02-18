/* =======================================================================================

   FILENAME :       srt_f_fwdcurve.c

   PURPOSE:         functions to access  information from the SrtCurvePtr object stored
                                    in the double linked list, refering to a SrtFwdDesc, containing
   a SrtDvdObj

   ======================================================================================= */

#include "srt_h_fwdcurve.h"
#include "srt_h_fwdobj.h"
#include "srt_h_repo_obj.h"
#include "swp_h_all.h"
#include "swp_h_curve_struct.h"
#include "utallhdr.h"

/* Compute the forward from 0.0 to start_time */
double srt_f_forward_from_fwdcrv(double start_time, SrtCurvePtr dvdcrv, SrtCurvePtr repocrv)
{
    SrtDvdObj*  dvd_obj;
    SrtRepoObj* repo_obj;

    SrtErr err = NULL;
    double dvd, repo, fwd, today;
    double start_date;

    /* Check the type of the curve */
    if (!ISCURVETYPE(dvdcrv, DVD_CURVE))
        return 0.0;

    if (!ISCURVETYPE(repocrv, REPO_CURVE))
        return 0.0;

    /* Get the forward objectd attached to the SrtCurvePtr = SrtCurveDesc */
    dvd_obj  = get_dvdobj_from_dvdcrv(dvdcrv);
    repo_obj = get_repoobj_from_repocrv(repocrv);

    /* Get today from the obj */
    err = srt_f_dvdobj_extracttoday(dvd_obj, &today);
    if (err)
        return 0.0;

    /* Rebuilds the dates from the time and today */
    start_date = today + start_time * DAYS_IN_YEAR;

    /* 	Compute the relevant forward from the Fwd object */
    err = srt_f_dvdobj_dvd(dvd_obj, start_date, &dvd);
    if (err)
        return 0.0;

    err = srt_f_repo_obj_repo(repo_obj, start_date, &repo);
    if (err)
        return 0.0;

    fwd = repo * dvd;

    return fwd;

} /* END double srt_f_forward_from_fwdcrv(...) */

/* ---------------------------------------------------------------------------------- */

Err srt_f_init_DividendCurve(
    Date today, char* ccy, double** dividend_curve, int ncols, int nrows, String dvd_name)
{
    SrtDvdObj*    dvd_obj;
    SrtCurveList* curve_list;
    Err           err = NULL;

    /* Get the curves list and check it has been initialised  */
    curve_list = get_curve_list();
    if (curve_list == NULL)
        return serror("No Curve list defined: call SrtInit before");

    /* Builds the full SrtDvdObj from the input Values */
    err = srt_f_dvdobj_init(today, dividend_curve, ncols, nrows, dvd_name, &dvd_obj);
    if (err)
        return err;

    /*  Puts the SrtDvdObj as Curve in the Curve List with fwd_name as the reference */
    err = swp_f_addcurvetolist(
        curve_list, dvd_name, "DVD_CURVE", NULL, ccy, today, 0, NULL, (void*)(dvd_obj));
    if (err)
    {
        return serror("Fatal: (SrtInitEQUnd) Failed to add EQ underlying");
    }

    /* Return a success message */
    return NULL;

} /* END of srt_f_init_DividendCurve(...) */

/* -------------------------------------------------------------------------- */
