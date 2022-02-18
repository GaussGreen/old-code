#include "QuadGauss.h"

Err qg_free_struct(qg_str* qgstr)
{
    Err err = NULL;

    if (qgstr->dates)
    {
        free(qgstr->dates);
        qgstr->dates = NULL;
    }
    if (qgstr->times)
    {
        free(qgstr->times);
        qgstr->times = NULL;
    }
    if (qgstr->sigma1)
    {
        free(qgstr->sigma1);
        qgstr->sigma1 = NULL;
    }
    if (qgstr->e1_1)
    {
        free(qgstr->e1_1);
        qgstr->e1_1 = NULL;
    }
    if (qgstr->e1_2)
    {
        free(qgstr->e1_2);
        qgstr->e1_2 = NULL;
    }
    if (qgstr->sigma2)
    {
        free(qgstr->sigma2);
        qgstr->sigma2 = NULL;
    }
    if (qgstr->e2_1)
    {
        free(qgstr->e2_1);
        qgstr->e2_1 = NULL;
    }
    if (qgstr->e2_2)
    {
        free(qgstr->e2_2);
        qgstr->e2_2 = NULL;
    }
    if (qgstr->rho)
    {
        free(qgstr->rho);
        qgstr->rho = NULL;
    }

    return NULL;
}

Err qg_free_und_struct(SrtUndPtr pUndDesc)
{
    SrtIrDesc* pSrtIrPtr;

    pSrtIrPtr = (SrtIrDesc*)(pUndDesc->spec_desc);
    qg_free_struct(pSrtIrPtr->spec);
    free(pSrtIrPtr->spec);
    free(pUndDesc);
    pUndDesc = NULL;
    return NULL;
}

Err qg_initialize_str(qg_str* qgstr)
{
    qgstr->dates  = NULL;
    qgstr->times  = NULL;
    qgstr->sigma1 = NULL;
    qgstr->e1_1   = NULL;
    qgstr->e1_2   = NULL;
    qgstr->sigma2 = NULL;
    qgstr->e2_1   = NULL;
    qgstr->e2_2   = NULL;
    qgstr->rho    = NULL;

    return NULL;
}

Err qg_init_struct(
    long today,

    /*	TS Dates	*/
    int   num_dates,
    long* dates,

    int num_factor,

    /*	X1 TS	*/
    double* sigma1,
    double* e1_1,
    double* e1_2,

    /*	X2 TS	*/
    double* sigma2,
    double* e2_1,
    double* e2_2,

    /*	Correl TS	*/
    double* rho,

    /*	Structure	*/
    qg_str* qgstr)
{
    Err err = NULL;
    int i;

    err = qg_initialize_str(qgstr);

    qgstr->today      = today;
    qgstr->num_factor = num_factor;

    qgstr->num_dates = num_dates;
    qgstr->dates     = (long*)calloc(num_dates, sizeof(long));
    qgstr->times     = (double*)calloc(num_dates, sizeof(double));

    qgstr->sigma1 = (double*)calloc(num_dates, sizeof(double));
    qgstr->e1_1   = (double*)calloc(num_dates, sizeof(double));
    qgstr->e1_2   = (double*)calloc(num_dates, sizeof(double));

    qgstr->sigma2 = (double*)calloc(num_dates, sizeof(double));
    qgstr->e2_1   = (double*)calloc(num_dates, sizeof(double));
    qgstr->e2_2   = (double*)calloc(num_dates, sizeof(double));

    qgstr->rho = (double*)calloc(num_dates, sizeof(double));

    if ((!qgstr->dates) || (!qgstr->times)

        || (!qgstr->sigma1) || (!qgstr->e1_1) || (!qgstr->e1_2)

    )
    {
        err = "memory allocation failed in qg_init_str";
        goto FREE_RETURN;
    }

    if ((qgstr->num_factor == 2) &&
        ((!qgstr->sigma2) || (!qgstr->e2_1) || (!qgstr->e2_2) || (!qgstr->rho)))
    {
        err = "memory allocation failed in qg_init_str";
        goto FREE_RETURN;
    }

    memcpy(qgstr->dates, dates, num_dates * sizeof(long));
    for (i = 0; i < num_dates; ++i)
    {
        qgstr->times[i] = (dates[i] - today) / 365.0;
    }

    memcpy(qgstr->sigma1, sigma1, num_dates * sizeof(double));
    memcpy(qgstr->e1_1, e1_1, num_dates * sizeof(double));
    memcpy(qgstr->e1_2, e1_2, num_dates * sizeof(double));

    memcpy(qgstr->sigma2, sigma2, num_dates * sizeof(double));
    memcpy(qgstr->e2_1, e2_1, num_dates * sizeof(double));
    memcpy(qgstr->e2_2, e2_2, num_dates * sizeof(double));

FREE_RETURN:

    if (err)
    {
        qg_free_struct(qgstr);
    }

    return err;
}

Err SrtInitQGUnd(/* und name */
                 char* undName,

                 /* mkt name */
                 char* ycname,

                 /*	TS Dates	*/
                 int   num_dates,
                 long* dates,

                 int num_factor,

                 /*	X1 TS	*/
                 double* sigma1,
                 double* e1_1,
                 double* e1_2,

                 /*	X2 TS	*/
                 double* sigma2,
                 double* e2_1,
                 double* e2_2,

                 /*	Correl TS	*/
                 double* rho)
{
    int bCleanUpUndFlag = 1;

    SrtUndPtr     pUndPtr;
    SrtIrDesc*    pIrPtr;
    SrtUndListPtr und_list;
    SrtCurvePtr   pYieldCurve;

    qg_str* qgstr = NULL;

    char* ccy;
    long  today;

    Err err = NULL;

    /* Allocate the qg_str structure */
    qgstr = calloc(1, sizeof(qg_str));

    /* Get the yield curve and the spot date from the yield curve name */
    pYieldCurve = lookup_curve(ycname);

    if (!pYieldCurve)
    {
        err = "Cannot find Yield Curve";
        goto FREE_RETURN;
    }

    today = get_today_from_curve(pYieldCurve);
    ccy   = get_curve_ccy(pYieldCurve);

    /* Create the new General underlying */

    pUndPtr = (SrtUndPtr)calloc(1, sizeof(SrtUndDesc));

    pUndPtr->underl_type = INTEREST_RATE_UND;
    strcpy(pUndPtr->underl_name, undName);
    strupper(pUndPtr->underl_name);
    strip_white_space(pUndPtr->underl_name);
    strcpy(pUndPtr->underl_lbl, "QG_UND");
    pUndPtr->underl_ccy = ccy;

    /* Create the IR underlying */

    pIrPtr = calloc(1, sizeof(SrtIrDesc));

    strcpy(pIrPtr->mdl_lbl, "QG_UND");
    pIrPtr->mdl_type = QUAD_GAUSS;

    if (num_factor == 1)
    {
        pIrPtr->mdl_dim = ONE_FAC;
    }
    else if (num_factor == 2)
    {
        pIrPtr->mdl_dim = TWO_FAC;
    }
    else
    {
        err = "QUAD_GAUSS cannot accept more than two factors";
        goto FREE_RETURN;
    }

    strcpy(pIrPtr->yc_name, pYieldCurve->curve_name);
    pIrPtr->spec = qgstr;

    /* Init Structure */
    err = qg_init_struct(
        today,

        /*	TS Dates	*/
        num_dates,
        dates,

        num_factor,

        /*	X1 TS	*/
        sigma1,
        e1_1,
        e1_2,

        /*	X2 TS	*/
        sigma2,
        e2_1,
        e2_2,

        /*	Correl TS	*/
        rho,

        /*	Structure	*/
        qgstr);

    pUndPtr->spec_desc = pIrPtr;

    /* Put the underlying into the depot */

    und_list = get_underlying_list();

    err = srt_f_lstins(
        und_list,
        pUndPtr->underl_name,
        0.0,
        OBJ_PTR_UND,
        (void*)pUndPtr,
        &qg_free_und_struct,
        &(pUndPtr->underl_ticker));

FREE_RETURN:

    return err;
}

Err qg_get_struct_from_und(char* und, qg_str** qgstr)
{
    Err        err = NULL;
    SrtUndPtr  pUndPtr;
    SrtIrDesc* pIrPtr;

    // Get the underlying through its name and check it exists
    // Check on the underlying type
    pUndPtr = lookup_und(und);
    if (!pUndPtr)
    {
        err = "Undefined underlying";
        return err;
    }

    if (get_underlying_type(pUndPtr) != INTEREST_RATE_UND)
    {
        return serror("Underlying %s is not of type IR", und);
        goto FREE_RETURN;
    }

    if (get_mdltype_from_irund(pUndPtr) != QUAD_GAUSS)
    {
        err = serror("Underlying %s is not of type QUAD_GAUSS", und);
        goto FREE_RETURN;
    }

    // Extract the information from the underlying
    pIrPtr = (SrtIrDesc*)(pUndPtr->spec_desc);
    *qgstr = (qg_str*)(pIrPtr->spec);

FREE_RETURN:

    return err;
}

Err qg_get_ts_from_und(/* qg und name */
                       char* undName,

                       /* TS Dates and Times */
                       int*     num_dates,
                       long**   dates,
                       double** times,

                       int* num_factor,

                       /*	X1 TS	*/
                       double** sigma1,
                       double** e1_1,
                       double** e1_2,

                       /*	X2 TS	*/
                       double** sigma2,
                       double** e2_1,
                       double** e2_2,

                       /*	Correl TS	*/
                       double** rho)
{
    Err     err   = NULL;
    qg_str* qgstr = NULL;

    err = qg_get_struct_from_und(undName, &qgstr);

    *num_dates = qgstr->num_dates;

    *dates = (long*)calloc(qgstr->num_dates, sizeof(long));
    *times = (double*)calloc(qgstr->num_dates, sizeof(double));
    *rho   = (double*)calloc(qgstr->num_dates, sizeof(double));

    *sigma1 = (double*)calloc(qgstr->num_dates, sizeof(double));
    *e1_1   = (double*)calloc(qgstr->num_dates, sizeof(double));
    *e1_2   = (double*)calloc(qgstr->num_dates, sizeof(double));

    *sigma2 = (double*)calloc(qgstr->num_dates, sizeof(double));
    *e2_1   = (double*)calloc(qgstr->num_dates, sizeof(double));
    *e2_2   = (double*)calloc(qgstr->num_dates, sizeof(double));
    if ((!(*dates)) || (!(*times)) || (!(*rho)) || (!(*sigma1)) || (!(*e1_1)) || (!(*e1_2)) ||
        (!(*sigma2)) || (!(*e2_1)) || (!(*e2_2)))
    {
        err = "Allocation failed in qg_get_ts_from_und";
        goto FREE_RETURN;
    }

    *num_factor = qgstr->num_factor;

    memcpy(*dates, qgstr->dates, qgstr->num_dates * sizeof(long));
    memcpy(*times, qgstr->times, qgstr->num_dates * sizeof(long));

    memcpy(*sigma1, qgstr->sigma1, qgstr->num_dates * sizeof(double));
    memcpy(*e1_1, qgstr->e1_1, qgstr->num_dates * sizeof(double));
    memcpy(*e1_2, qgstr->e1_2, qgstr->num_dates * sizeof(double));

    memcpy(*sigma2, qgstr->sigma2, qgstr->num_dates * sizeof(double));
    memcpy(*e2_1, qgstr->e2_1, qgstr->num_dates * sizeof(double));
    memcpy(*e2_2, qgstr->e2_2, qgstr->num_dates * sizeof(double));

    memcpy(*rho, qgstr->rho, qgstr->num_dates * sizeof(double));

FREE_RETURN:

    if (err)
    {
        if (*dates)
        {
            free(*dates);
            *dates = NULL;
        }
        if (*times)
        {
            free(*times);
            *times = NULL;
        }
        if (*rho)
        {
            free(*rho);
            *rho = NULL;
        }
        if (*sigma1)
        {
            free(*sigma1);
            *sigma1 = NULL;
        }
        if (*e1_1)
        {
            free(*e1_1);
            *e1_1 = NULL;
        }
        if (*e1_2)
        {
            free(*e1_2);
            *e1_2 = NULL;
        }
        if (*sigma2)
        {
            free(*sigma2);
            *sigma2 = NULL;
        }
        if (*e2_1)
        {
            free(*e2_1);
            *e2_1 = NULL;
        }
        if (*e2_2)
        {
            free(*e2_2);
            *e2_2 = NULL;
        }
    }

    return err;
}
