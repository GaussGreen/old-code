/* ============================================================================

   FILENAME:	SRT_F_TS_IRM_2F_INIT.C

   PURPOSE:		1.-	Initialise a volatility term structure for a two factor model,
                                        straight from raw data (dates + values).
                                2.-	Retrieve all useful quantities at any date from
                                        a given Term Struct.

   TWO FACTOR EQUIVALENT OF SRT_F_TS_IRM_INIT.C

   ============================================================================   */
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_ts_irm_2f.h"

static Err make_2f_J_vectors(TermStruct* l);

/* ----------------------------------------------------------------------------------
                                                                PART I
                                MAIN FUNCTION THAT DEFINES A TERM STRUCT
   ---------------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------- */

/* Function required to free the TermStruct once attached into a linked list */

Err srt_f_irm2ftsvalfree(void* tsvalptr)
{
    TwoFacIrmTermStructVal* tsval = (TwoFacIrmTermStructVal*)tsvalptr;
    Err                     err   = NULL;
    srt_free(tsval);

    return err;
}

/* ---------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------
   The MAIN function to call when a Interest Rate Model TermStructure has
   to be initialised with a TWO FACTOR MODEL
   ------------------------------------------------------------------------- */
Err srt_f_init_IRM_TwoFac_TermStruct(
    TermStruct** ts, /* Pointer to the Term Struct (Filled on output) */

    Date today, /* Today's date */

    double** sig_data, /* Sigma Curve: [0]Date [1]sig1 [2]sig2 [3]ppdRho */
    int      sig_cols, /* Number of Columns (4) */
    int      num_sig,  /* Number of dates for discretisation */
    double** tau_data, /* Tau Curve: [0]Date [1]ppdTau1 [2]ppdTau2 */
    int      tau_cols, /* Number of Coulmns (3) */
    int      num_tau,  /* Number of dates for discretisation */

    SrtMdlType mdl_type,

    /* SMILE */
    double beta, /* For Cheyette Beta */

    /* CORRELATION */
    double alpha,
    double gamma,
    double ppdRho,

    /* MIXED SMILES */
    double omega)
{
    int cur_sig,                 /* Index of the current sigma strarting from 0 */
        cur_tau;                 /* Index of the current tau strarting from 0 */
    TwoFacIrmTermStructVal* ts2; /* One particular atom of the Term Struct */
    int over_sig,                /* SRT_Boolean: whether the current sigma is above the last one */
        over_tau;                /* SRT_Boolean: whether the current tau is above the last one */
    Err      err = NULL;
    int      i;
    double** ppdSigmaValues;
    double** ppdTauValues;
    long     ticker;

    /* Checks on the number of columns in Sigma : one value, or two/four(corr)/five(corr+smile)/six
     * columns */
    if ((sig_cols != 2) && (sig_cols != 4) && (sig_cols != 5) && (sig_cols != 6))
    {
        if ((sig_cols != 1) || (num_sig > 1))
            return serror("Need TWO/FOUR/FIVE/SIX columns or ONE SINGLE VALUE for the Sigma Curve");
    }

    /* Checks on the number of columns in Tau : one value or two/three columns */
    if ((tau_cols != 2) && (tau_cols != 3))
    {
        if ((tau_cols != 1) || (num_tau > 1))
            return serror("Need TWO/THREE columns or ONE SINGLE VALUE for the Tau Curve");
    }

    /* Allocate memory for a new array of data for the full term structure*/
    ppdSigmaValues = dmatrix(0, 5, 0, num_sig - 1);
    ppdTauValues   = dmatrix(0, 2, 0, num_tau - 1);

    /* If only one sigma value input, make a fake term structure */
    if (sig_cols == 1)
    {
        /* Date */
        ppdSigmaValues[0][0] = today + 365.0;

        /* Sigma 1 */
        ppdSigmaValues[1][0] = sig_data[0][0];

        /* Beta 1 */
        ppdSigmaValues[2][0] = beta;

        /* Sigma 2 */
        ppdSigmaValues[3][0] = ppdSigmaValues[1][0] * alpha;

        /* Beta 2 */
        ppdSigmaValues[4][0] = beta + omega;

        /* Rho */
        ppdSigmaValues[5][0] = ppdRho;
    }
    else
        /* Transfer the input data into a real size matrix */
        for (i = 0; i < num_sig; i++)
        {
            /* Date */
            ppdSigmaValues[0][i] = sig_data[0][i];

            /* Sigma 1 */
            ppdSigmaValues[1][i] = sig_data[1][i];

            /* Sigma 2 and Rho */
            if (sig_cols <= 2)
            {
                /* If correlation parameters input directly to the function (alpha, ppdRho) */
                ppdSigmaValues[3][i] = ppdSigmaValues[1][i] * alpha;
                ppdSigmaValues[5][i] = ppdRho;
            }
            else if (sig_cols == 4)
            {
                /* Four columns input : Date - Sig1 - Sig2 - Rho */
                ppdSigmaValues[3][i] = sig_data[2][i];
                ppdSigmaValues[5][i] = sig_data[3][i];
            }
            else if (sig_cols == 5)
            {
                /* Five columns input : Date - Sig1 - Beta - Sig2 - Rho */
                ppdSigmaValues[3][i] = sig_data[3][i];
                ppdSigmaValues[5][i] = sig_data[4][i];
            }
            else if (sig_cols == 6)
            {
                /* Six columns input : Date - Sig1 - Beta - Sig2 - Beta2 - Rho */
                ppdSigmaValues[3][i] = sig_data[3][i];
                ppdSigmaValues[5][i] = sig_data[5][i];
            }

            /* Beta */
            if (sig_cols <= 4)
            {
                /* If Beta input directly to the function */
                ppdSigmaValues[2][i] = beta;
            }
            else
            {
                /* More than four columns : Date - Sig1 - Beta - Sig 2 - ... */
                ppdSigmaValues[2][i] = sig_data[2][i];
            }

            /* Omega */
            if (sig_cols <= 5)
            {
                /* Omega input directly to the function */
                ppdSigmaValues[4][i] = ppdSigmaValues[2][i] + omega;
            }
            else
            {
                /* Six columns : Date - Sig1 - Beta1 - Sig 2 - Beta2 - Rho  */
                ppdSigmaValues[4][i] = sig_data[4][i];
            }

        } /* END for i loop on number of rows in sigma */

    /* Overwrites the Values of Beta for explicit models */
    if ((mdl_type == LGM) || (mdl_type == LGM_STOCH_VOL))
    {
        for (i = 0; i < num_sig; i++)
            ppdSigmaValues[2][i] = 0.0;
    }
    else if ((mdl_type == CHEY) || (mdl_type == CHEY_STOCH_VOL))
    {
        for (i = 0; i < num_sig; i++)
            ppdSigmaValues[2][i] = 1.0;
    }

    /* Overwrites the Values of Beta 2 for non mixed models */
    if (mdl_type != MIXED_BETA)
    {
        for (i = 0; i < num_sig; i++)
            ppdSigmaValues[4][i] = ppdSigmaValues[2][i];
        0.0;
    }

    /* If only one tau value input, make a fake term structure */
    if (tau_cols == 1)
    {
        ppdTauValues[0][0] = today + 365.0;
        ppdTauValues[1][0] = tau_data[0][0];
        ppdTauValues[2][0] = 1.0 / (1 / ppdTauValues[1][0] + gamma);
    }
    else
        /* Transfer the input data into a real size matrix */
        for (i = 0; i < num_tau; i++)
        {
            ppdTauValues[0][i] = tau_data[0][i];
            ppdTauValues[1][i] = tau_data[1][i];

            /* If correlation parameters input directly to the function ( gamma )*/
            if (tau_cols != 3)
            {
                ppdTauValues[2][i] = 1.0 / (1 / tau_data[1][i] + gamma);
            }
            else
                ppdTauValues[2][i] = tau_data[2][i];
        }

    /* If the data exists but is NULL: consider whatever is after as garbage */
    for (i = 0; i < num_sig; i++)
    {
        if (!ppdSigmaValues[0][i] || !ppdSigmaValues[1][i])
            break;
    }
    num_sig = i;
    for (i = 0; i < num_tau; i++)
    {
        if (!ppdTauValues[0][i] || !ppdTauValues[1][i])
            break;
    }
    num_tau = i;

    /* Initialise indexes */
    cur_tau  = 0;
    cur_sig  = 0;
    over_sig = SRT_FALSE;
    over_tau = SRT_FALSE;

    /* Skip dates <= today */
    while ((cur_sig < num_sig) && ((ppdSigmaValues[0][cur_sig] - today) <= 0))
    {
        cur_sig++;
    }
    while ((cur_tau < num_tau) && ((ppdTauValues[0][cur_tau] - today) <= 0))
    {
        cur_tau++;
    }

    /* Return error if there is no acceptable sigma or tau value */
    if ((cur_sig == num_sig) && (cur_tau == num_tau))
    {
        free_dmatrix(ppdSigmaValues, 0, 6, 0, num_sig - 1);
        free_dmatrix(ppdTauValues, 0, 2, 0, num_tau - 1);
        return serror("init_TWOFAC_TermStruct err: no current data !!!");
    }

    /* Makes sure there is no negative values for sigma */
    for (i = 0; i < num_sig; i++)
    {
        if ((ppdSigmaValues[1][i] < 0.0) || (ppdSigmaValues[2][i] < 0.0))
        {
            free_dmatrix(ppdSigmaValues, 0, 5, 0, num_sig - 1);
            free_dmatrix(ppdTauValues, 0, 2, 0, num_tau - 1);
            return serror("Cannot input negative volatility");
        }
    }

    /* Makes sure the correlation is properly bounded */
    for (i = 0; i < num_sig; i++)
    {
        if ((ppdSigmaValues[5][i] < -1.0) || (ppdSigmaValues[5][i] > 1.0))
        {
            free_dmatrix(ppdSigmaValues, 0, 5, 0, num_sig - 1);
            free_dmatrix(ppdTauValues, 0, 2, 0, num_tau - 1);
            return serror("Correlation must be between -1.0 and 1.0");
        }
    }

    /* Create the Term Struct ts */
    srt_f_lstcreate(ts, "TermStruct");

    /* Main loop to fill the Term Struct as long as there are data in */
    while ((over_tau == SRT_FALSE) || (over_sig == SRT_FALSE))
    {
        /* Creation of a new TermStruct Atom */
        ts2 = (TwoFacIrmTermStructVal*)srt_calloc(1, sizeof(TwoFacIrmTermStructVal));

        /* Fill the Term Struct Atom with initial data */
        ts2->exp_ts[0].sig = ppdSigmaValues[1][cur_sig];
        ts2->exp_ts[1].sig = ppdSigmaValues[3][cur_sig];

        ts2->exp_ts[0].tau = ppdTauValues[1][cur_tau];
        ts2->exp_ts[1].tau = ppdTauValues[2][cur_tau];

        ts2->exp_ts[0].beta = ppdSigmaValues[2][cur_sig];
        ts2->exp_ts[1].beta = ppdSigmaValues[4][cur_sig];

        ts2->alpha = ppdSigmaValues[3][cur_sig] / ppdSigmaValues[1][cur_sig];
        ts2->gamma = 1.0 / ppdTauValues[2][cur_tau] - 1.0 / ppdTauValues[1][cur_tau];
        ts2->rho   = ppdSigmaValues[5][cur_sig];

        ts2->omega = ppdSigmaValues[4][cur_sig] - ppdSigmaValues[2][cur_sig];

        /* We create a new date. Where does it come from? Sigma or Tau Struct? */
        /* 1st case: Date comes from both structs, date_sigma = date_tau */
        if (((over_sig == SRT_FALSE) && (over_tau == SRT_FALSE)) &&
            (ppdSigmaValues[0][cur_sig] == ppdTauValues[0][cur_tau]))
        {
            /* Set up specific data */
            ts2->val_origin = BOTH_DATE;
            ts2->date       = DTOL(ppdSigmaValues[0][cur_sig]);
            ts2->time       = (ts2->date - today) * YEARS_IN_DAY;

            /* Update indexes */
            if (cur_sig < num_sig - 1)
            {
                cur_sig++;
            }
            else
            {
                over_sig = SRT_TRUE;
            }
            if (cur_tau < num_tau - 1)
            {
                cur_tau++;
            }
            else
            {
                over_tau = SRT_TRUE;
            }
        }
        /* 2nd case: Date comes from sigma struct, date_sigma < date_tau */
        else if (
            (over_tau == SRT_TRUE) || (((over_sig == SRT_FALSE) && (over_tau == SRT_FALSE)) &&
                                       (ppdSigmaValues[0][cur_sig] < ppdTauValues[0][cur_tau])))
        {
            /* Set up specific data */
            ts2->val_origin = SIGMA_DATE;
            ts2->date       = DTOL(ppdSigmaValues[0][cur_sig]);
            ts2->time       = (ts2->date - today) * YEARS_IN_DAY;

            /* Update indexes */
            if (cur_sig < num_sig - 1)
            {
                cur_sig++;
            }
            else
            {
                over_sig = SRT_TRUE;
            }
        }
        /* 3rd case: date_tau < date_sigma */
        else if (
            (over_sig == SRT_TRUE) || (((over_sig == SRT_FALSE) && (over_tau == SRT_FALSE)) &&
                                       (ppdSigmaValues[0][cur_sig] > ppdTauValues[0][cur_tau])))
        {
            /* Set up specific data */
            ts2->val_origin = TAU_DATE;
            ts2->date       = DTOL(ppdTauValues[0][cur_tau]);
            ts2->time       = (ts2->date - today) * YEARS_IN_DAY;

            /* Update indexes */
            if (cur_tau < num_tau - 1)
            {
                cur_tau++;
            }
            else
            {
                over_tau = SRT_TRUE;
            }
        }

        /* Insert the TermStructAtom in the TS (==linked ts sorted by key = date) */
        srt_f_lstins(
            *ts,
            "TwoFacTsAtom",
            ts2->date,
            OBJ_PTR_IRM2f_TermStruct,
            (void*)ts2,
            &srt_f_irm2ftsvalfree,
            &ticker);
    }

    /* Free memory */
    free_dmatrix(ppdSigmaValues, 0, 5, 0, num_sig - 1);
    free_dmatrix(ppdTauValues, 0, 2, 0, num_tau - 1);

    /* Computes LGM Stuff */
    err = make_2f_F_Psi_vectors(*ts);
    err = make_2f_G_H_vectors(*ts);

    /* Return the error (if any) */
    return err;

} /* END srt_f_init_IRM_TwoFac_TermStruct(...) */

/* ------------------------------------------------------------------------- */

/* From Term Struct structure to arrays */
Err srt_f_display_IRM_TwoFac_TermStruct(
    TermStruct* ts,
    double**    ppdSigmaDates,
    double**    ppdSigma1,
    double**    ppdBeta1,
    double**    ppdSigma2,
    double**    ppdBeta2,
    double**    ppdRho,
    long*       plNumSigmas,
    double**    ppdTauDates,
    double**    ppdTau1,
    double**    ppdTau2,
    long*       plNumTaus)
{
    SrtLst* l = ts->head;

    *plNumSigmas = *plNumTaus = 0;

    *ppdSigmaDates = (double*)malloc(sizeof(double));
    *ppdTauDates   = (double*)malloc(sizeof(double));

    *ppdSigma1 = (double*)malloc(sizeof(double));
    *ppdBeta1  = (double*)malloc(sizeof(double));
    *ppdSigma2 = (double*)malloc(sizeof(double));
    *ppdBeta2  = (double*)malloc(sizeof(double));
    *ppdRho    = (double*)malloc(sizeof(double));

    *ppdTau1 = (double*)malloc(sizeof(double));
    *ppdTau2 = (double*)malloc(sizeof(double));

    while (l != NULL)
    {
        if ((((TwoFacIrmTermStructVal*)l->element->val.pval)->val_origin == TAU_DATE) ||
            (((TwoFacIrmTermStructVal*)l->element->val.pval)->val_origin == BOTH_DATE))
        {
            (*plNumTaus)++;
            *ppdTauDates = realloc(*ppdTauDates, (*plNumTaus) * sizeof(double));
            *ppdTau1     = realloc(*ppdTau1, (*plNumTaus) * sizeof(double));
            *ppdTau2     = realloc(*ppdTau2, (*plNumTaus) * sizeof(double));
            (*ppdTauDates)[*plNumTaus - 1] =
                (long)(((TwoFacIrmTermStructVal*)l->element->val.pval)->date);
            (*ppdTau1)[*plNumTaus - 1] =
                ((TwoFacIrmTermStructVal*)l->element->val.pval)->exp_ts[0].tau;
            (*ppdTau2)[*plNumTaus - 1] =
                ((TwoFacIrmTermStructVal*)l->element->val.pval)->exp_ts[1].tau;
        }

        if ((((TwoFacIrmTermStructVal*)l->element->val.pval)->val_origin == SIGMA_DATE) ||
            (((TwoFacIrmTermStructVal*)l->element->val.pval)->val_origin == BOTH_DATE))
        {
            (*plNumSigmas)++;
            *ppdSigmaDates = realloc(*ppdSigmaDates, (*plNumSigmas) * sizeof(double));
            *ppdSigma1     = realloc(*ppdSigma1, (*plNumSigmas) * sizeof(double));
            *ppdBeta1      = realloc(*ppdBeta1, (*plNumSigmas) * sizeof(double));
            *ppdSigma2     = realloc(*ppdSigma2, (*plNumSigmas) * sizeof(double));
            *ppdBeta2      = realloc(*ppdBeta2, (*plNumSigmas) * sizeof(double));
            *ppdRho        = realloc(*ppdRho, (*plNumSigmas) * sizeof(double));
            (*ppdSigmaDates)[*plNumSigmas - 1] =
                (long)(((TwoFacIrmTermStructVal*)l->element->val.pval)->date);
            (*ppdSigma1)[*plNumSigmas - 1] =
                ((TwoFacIrmTermStructVal*)l->element->val.pval)->exp_ts[0].sig;
            (*ppdBeta1)[*plNumSigmas - 1] =
                ((TwoFacIrmTermStructVal*)l->element->val.pval)->exp_ts[0].beta;
            (*ppdSigma2)[*plNumSigmas - 1] =
                ((TwoFacIrmTermStructVal*)l->element->val.pval)->exp_ts[1].sig;
            (*ppdBeta2)[*plNumSigmas - 1] =
                ((TwoFacIrmTermStructVal*)l->element->val.pval)->exp_ts[1].beta;
            (*ppdRho)[*plNumSigmas - 1] = ((TwoFacIrmTermStructVal*)l->element->val.pval)->rho;
        }

        l = l->next;
    }

    return NULL;

} /* END of Err srt_f_display_2f_ts	(...) */

/* ----------------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------

        STATIC FUNCTIONS TO PRECOMPUTE SOME INTEGRALS IN THE TERM STRUCTURE

   ---------------------------------------------------------------------------- */

/* Computes F and Psi vectors */
/* Stores values corresponding to the time steps of the Term Struct */
/* CAREFUL: F(ti) WILL BE STORED AT TIME STEP ti+1 AND THE SAME FOR PSI */
Err make_2f_F_Psi_vectors(TermStruct* l)
{
    TwoFacIrmTermStructVal *p, /* Previous Atom */
        *c;                    /* Current Atom */
    SrtLst *lp,                /* Previous Element in chain ts */
        *lc;                   /* Current element in chain ts */
    double dt,                 /* delta_t between current and previous */
        lambda,                /* Mean Reversion */
        temp0,                 /* exp (-lambda * dt) */
        temp1;                 /* (1.0 - exp (-lambda * dt)) / lambda */
    int i;
    Err err = NULL;

    /* Initialisation */
    /* First element of the ts */
    lc = l->head;
    c  = (TwoFacIrmTermStructVal*)lc->element->val.pval;
    /* F and Psi are initialised at 1 and 0 */
    for (i = 0; i <= 1; i++)
    {
        c->exp_ts[i].F   = 1.0;
        c->exp_ts[i].Psi = 0.0;
    }

    /* Main loop on all elements of the Term Struct */
    while (lc->next != NULL)
    {
        /* Compute dt */
        /* ---------- */
        /* ABOVE First Element */
        if (lc->previous != NULL)
        {
            /* Set previous */
            lp = lc->previous;
            p  = (TwoFacIrmTermStructVal*)lp->element->val.pval;
            /* Compute dt */
            dt = c->time - p->time;
        }
        else
        /* If firsty element dt = t */
        {
            dt = c->time;
        }

        /* Increment previous to current and current to next */
        p  = c;
        lc = lc->next;
        c  = (TwoFacIrmTermStructVal*)lc->element->val.pval;

        /* Loop on factors */
        for (i = 0; i <= 1; i++)
        {
            lambda = 1.0 / p->exp_ts[i].tau;
            temp0  = exp(-lambda * dt);

            /* Lambda is significative */
            if (fabs(lambda) > EPS)
            {
                temp1 = (1.0 - temp0) / lambda;
            }
            else
            /* Lambda is 0: temp converges to dt */
            {
                temp1 = dt;
            }

            /* Chain rule on F and Psi (see paper by S. Rady and OVE) */
            c->exp_ts[i].F   = p->exp_ts[i].F * temp0;
            c->exp_ts[i].Psi = p->exp_ts[i].Psi + p->exp_ts[i].F * temp1;
        }
    }
    return err;

} /* END static Err make_2f_F_Psi_vectors(...) */

/* ------------------------------------------------------------------------------ */

/* Computes F and Psi matrices (why vectors?) */
/* Stores values corresponding to the time steps of the Term Struct */
/* CAREFUL: G(ti) WILL BE STORED AT TIME STEP ti+1 AND THE SAME FOR H */
/* ------------------------------------------------------------------ */
/* F AND PSI VECTORS ARE TO BE MADE FIRST */
Err make_2f_G_H_vectors(TermStruct* l)
{
    TwoFacIrmTermStructVal *p, *c;
    SrtLst *                lp, *lc;
    double                  var[2], covar, lambda[2], dt, temp0[2], /* exp (-lambda * dt) */
        temp1[2],     /* (1.0 - exp (-lambda * dt) / lambda */
        temp2[2],     /* (exp (2.0 * lambda * dt) - 1.0) / (2 * lambda) */
        tempc1,       /* (exp (lambda1 + lambda2) * dt - 1.0 )/ (lambda1 + lambda2) */
        tempc2[2][2]; /* (exp (-lambda1 - lambda2) * dt - 1.0 )/ (lambda1 + lambda2) */
    int i, j;
    Err err = NULL;

    /* Initialise G and H at 0.0 */
    lc = l->head;
    c  = (TwoFacIrmTermStructVal*)lc->element->val.pval;
    for (i = 0; i <= 1; i++)
    {
        for (j = 0; j <= 1; j++)
        {
            c->reb_ts[i][j].G = 0.0;
            c->reb_ts[i][j].H = 0.0;
        }
    }

    /* Main loop on chain ts elements */
    while (lc->next != NULL)
    {
        /* Compute dt */
        if (lc->previous != NULL)
        {
            lp = lc->previous;
            p  = (TwoFacIrmTermStructVal*)lp->element->val.pval;
            dt = c->time - p->time;
        }
        else
        {
            dt = c->time;
        }

        p  = c;
        lc = lc->next;
        c  = (TwoFacIrmTermStructVal*)lc->element->val.pval;

        /* Compute vars, lambdas and temps */
        for (i = 0; i <= 1; i++)
        {
            var[i]    = pow(p->exp_ts[i].sig, 2.0);
            lambda[i] = 1.0 / p->exp_ts[i].tau;
            temp0[i]  = exp(-lambda[i] * dt);
            if (fabs(lambda[i]) > EPS)
            {
                temp1[i] = (1.0 - temp0[i]) / lambda[i];
                temp2[i] = (exp(2 * lambda[i] * dt) - 1.0) / (2 * lambda[i]);
            }
            else
            /* Lambda is close to 0: quantities converge to dt */
            {
                temp1[i] = temp2[i] = dt;
            }
        }

        /* Compute covar and crossed temps */
        covar = p->rho * p->exp_ts[0].sig * p->exp_ts[1].sig;
        if (fabs(lambda[0] + lambda[1]) > EPS)
        {
            tempc1 = (exp((lambda[0] + lambda[1]) * dt) - 1.0) / (lambda[0] + lambda[1]);
        }
        else
        /* Lambda1 + lambda2 is approximately 0: quantities converge */
        {
            tempc1 = dt;
        }
        /* Tempc2 */
        for (i = 0; i <= 1; i++)
        {
            for (j = 0; j <= 1; j++)
            {
                if (fabs(lambda[i] + lambda[j]) > EPS)
                {
                    tempc2[i][j] =
                        (exp(-(lambda[i] + lambda[j]) * dt) - 1.0) / (lambda[i] + lambda[j]);
                }
                else
                /* Lambda1 + lambda2 is approximately 0: quantities converge */
                {
                    tempc2[i][j] = -dt;
                }
            }
        }

        /* Compute G and H with chain rule (paper by OVE) */
        /* G */
        for (i = 0; i <= 1; i++)
        {
            c->reb_ts[i][i].G = p->reb_ts[i][i].G + temp2[i] * var[i] / pow(p->exp_ts[i].F, 2.0);
        }
        c->reb_ts[0][1].G = c->reb_ts[1][0].G =
            p->reb_ts[0][1].G + tempc1 * covar / p->exp_ts[0].F / p->exp_ts[1].F;

        /* H */
        for (i = 0; i <= 1; i++)
        {
            for (j = 0; j <= 1; j++)
            {
                c->reb_ts[i][j].H =
                    temp0[i] * p->reb_ts[i][j].H +
                    temp0[i] * p->exp_ts[i].F * p->exp_ts[j].F * p->reb_ts[i][j].G * temp1[j] +
                    (i == j ? var[i] : covar) * (temp1[i] + tempc2[i][j]) / lambda[j];
            }
        }
    } /* END of main loop */
    return err;

} /* END Err make_2f_G_H_vectors (...) */

/* ------------------------------------------------------------------------------ */

/* For Jumping Numeraire */
/* Not checked: in God we trust */
static Err make_2f_J_vectors(TermStruct* l)
{
    TwoFacIrmTermStructVal *p, *c;
    SrtLst *                lp, *lc;
    double                  dt, temp, temp2;
    Err                     err = NULL;

    lc             = l->head;
    c              = (TwoFacIrmTermStructVal*)lc->element->val.pval;
    c->exp_ts[0].J = 0.0;
    c->exp_ts[1].J = 0.0;

    while (lc->next != NULL)
    {
        if (lc->previous != NULL)
        {
            lp = lc->previous;
            p  = (TwoFacIrmTermStructVal*)lp->element->val.pval;
            dt = c->time - p->time;
        }
        else
        {
            dt = c->time;
        }
        p  = c;
        lc = lc->next;
        c  = (TwoFacIrmTermStructVal*)lc->element->val.pval;

        temp = exp(dt / p->exp_ts[0].tau);
        if (fabs(1.0 / p->exp_ts[0].tau) > EPS)
        {
            temp2 = (temp - 1.0) * p->exp_ts[0].tau;
        }
        else
        {
            temp2 = dt;
        }

        c->exp_ts[0].J = p->exp_ts[0].J + p->exp_ts[0].sig * temp2 / p->exp_ts[0].F;

        temp = exp(dt / p->exp_ts[1].tau);
        if (fabs(1.0 / p->exp_ts[1].tau) > EPS)
        {
            temp2 = (temp - 1) * p->exp_ts[1].tau;
        }
        else
        {
            temp2 = dt;
        }

        c->exp_ts[1].J = p->exp_ts[1].J + p->exp_ts[1].sig * temp2 / p->exp_ts[1].F;
    }
    return err;

} /* END of Err make_2f_J_vectors (...) */

/* ------------------------------------------------------------------------------ */
