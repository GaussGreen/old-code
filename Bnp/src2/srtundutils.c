/* ===============================================================================

   FILENAME	:	SrtUndUtils.c

   PURPOSE:     Utility Functions to deal with underlying for the outside world

   =============================================================================== */

#include "SrtUndUtils.h"

#include "srt_h_all.h"
#include "srt_h_und_fct.h"

/* ----------------------------------------------------------------------------------- */

SRT_Boolean SrtIsUnderlyingDefined(char* und_name)
{
    SrtUndPtr und;

    /* Get the underlying*/
    und = lookup_und(und_name);

    /* Check if it is empty or not */
    if (!und)
        return SRT_NO;
    else
        return SRT_YES;

} /* END SRT_Boolean SrtIsUnderlyingDefined(...) */

/* ----------------------------------------------------------------------------------- */

SRT_Boolean SrtIsUnderlyingInterestRate(char* und_name)
{
    SrtUndPtr und;

    /* Get the underlying*/
    und = lookup_und(und_name);

    /* Check if it is defined */
    if (!und)
        return SRT_NO;

    /* Check the type of the underlying and return an answer accordingly */
    if ISUNDTYPE (und, INTEREST_RATE_UND)
        return SRT_YES;
    else
        return SRT_NO;

} /* END SRT_Boolean SrtIsUnderlyingInterestRate(...) */

/* ----------------------------------------------------------------------------------- */

SRT_Boolean SrtIsUnderlyingIrOneFactor(char* und_name)
{
    SrtUndPtr und;
    SrtMdlDim mdl_dim;
    Err       err;

    /* Get the underlying*/
    und = lookup_und(und_name);

    /* Check if it is defined */
    if (!und)
        return SRT_NO;

    /* Check the underlying is an interest rate */
    if (!(ISUNDTYPE(und, INTEREST_RATE_UND)))
        return SRT_NO;
    else
    {
        /* Get the dimension of the model */
        err = get_underlying_mdldim(und, &mdl_dim);

        /* Return an answer accordingly */
        if (mdl_dim == ONE_FAC)
            return SRT_YES;
        else
            return SRT_NO;
    }

} /* END SRT_Boolean SrtIsUnderlyingIrOneFactor(...) */

/* ----------------------------------------------------------------------------------- */

SRT_Boolean SrtIsUnderlyingIrTwoFactor(char* und_name)
{
    SrtUndPtr und;
    SrtMdlDim mdl_dim;
    Err       err;

    /* Get the underlying*/
    und = lookup_und(und_name);

    /* Check if it is defined */
    if (!und)
        return SRT_NO;

    /* Check the underlying is an interest rate */
    if (!(ISUNDTYPE(und, INTEREST_RATE_UND)))
        return SRT_NO;
    else
    {
        /* Get the dimension of the model */
        err = get_underlying_mdldim(und, &mdl_dim);

        /* Return an answer accordingly */
        if (mdl_dim == TWO_FAC)
            return SRT_YES;
        else
            return SRT_NO;
    }

} /* END SRT_Boolean SrtIsUnderlyingIrTwoFactor(...) */

/* ----------------------------------------------------------------------------------- */

void SrtGetUnderlyingModelName(char* und_name, char* mdl_name)
{
    SrtUndPtr  und;
    SrtMdlType mdl_type;
    SrtMdlDim  mdl_dim;
    Err        err;

    /* Get the underlying*/
    und = lookup_und(und_name);

    /* Check if it is defined */
    if (!und)
        strcpy(mdl_name, "NOT DEFINED");
    else
    {
        /* Get the type of the model */
        err = get_underlying_mdltype(und, &mdl_type);

        /* Get the dimension of the model */
        err = get_underlying_mdldim(und, &mdl_dim);

        /* Get the corresponding model */
        err = srt_f_translate_model(mdl_type, mdl_dim, mdl_name);
    }

} /* END void SrtGetUnderlyingModelName(...) */

/* ----------------------------------------------------------------------------------- */

void SrtGetUnderlyingYieldCurveName(char* und_name, char* yc_name)
{
    SrtUndPtr und;

    /* Get the underlying*/
    und = lookup_und(und_name);

    /* Check if it is defined */
    if (!und)
        strcpy(yc_name, "NOT DEFINED");
    else
        /* Get the yc of the underlying */
        strcpy(yc_name, get_ycname_from_underlying(und));

} /* END void SrtGetUnderlyingYieldCurveName(...) */

/* ----------------------------------------------------------------------------------- */

void SrtGetUnderlyingCurrency(char* und_name, char* ccy_name)
{
    SrtUndPtr und;

    /* Get the underlying*/
    und = lookup_und(und_name);

    /* Check if it is defined */
    if (!und)
        strcpy(ccy_name, "NOT DEFINED");
    else
        /* Get the yc of the underlying */
        strcpy(ccy_name, get_underlying_ccy(und));

} /* END void SrtGetUnderlyingCurrency(...) */

/* ----------------------------------------------------------------------------------- */

long SrtGetUnderlyingTicker(char* und_name)
{
    SrtUndPtr und;
    long      ticker;

    /* Get the underlying*/
    und = lookup_und(und_name);

    /* Check if the underlying exists */
    if (und == NULL)
        return 0;

    /* Get the ticker attached to the underlying */
    ticker = get_underlying_ticker(und);

    /* return the ticker */
    return ticker;

} /* END long SrtGetUnderlyingTicker(...) */

/* ----------------------------------------------------------------------------------- */

/* Allocate memory for arrays where the underlying TS is displayed */
Err SrtDisplayUndTermStruct(
    char* UndName,

    double*** SigmaCurve,
    long*     NumSigmaRows,
    long*     NumSigmaCols,

    double*** TauCurve,
    long*     NumTauRows,
    long*     NumTauCols,

    double* Alpha,
    double* Gamma,
    double* Rho,

    double* Beta,
    double* Omega,

    double* VoVol)
{
    Err         err = NULL;
    SrtUndPtr   und;
    TermStruct* ts;
    SrtMdlDim   mdl_dim;
    SrtMdlType  mdl_type;
    double *    temp1, *temp2, *temp3;
    int         i;

    /* Get the underlying through its name and check it exists */
    und = lookup_und(UndName);
    if (!und)
        return serror("Underlying %s not defined", UndName);

    /* Check on the underlying type and model dimension */
    if (!ISUNDTYPE(und, INTEREST_RATE_UND))
        return serror("%s is not an Interest Rate Underlying: cannot display TS", UndName);
    err = get_underlying_mdldim(und, &mdl_dim);
    if (err)
        return err;

    /* Get the TermStructure attached to the underlying */
    err = get_underlying_ts(und, &ts);
    if (err)
        return err;

    /* Aloocate memory and populates the output TermStruct according to the Model dimension */
    if (mdl_dim == ONE_FAC)
    {
        err = get_underlying_mdltype(und, &mdl_type);

        if (mdl_type == LGM_STOCH_VOL)
        {
            *NumSigmaCols = 6;
            (*SigmaCurve) = (double**)srt_malloc((*NumSigmaCols) * sizeof(double*));
            *NumTauCols   = 2;
            (*TauCurve)   = (double**)srt_malloc((*NumTauCols) * sizeof(double*));

            err = srt_f_display_IRM_OneFac_TermStruct(
                ts,
                &((*SigmaCurve)[0]),
                &((*SigmaCurve)[1]),
                &((*SigmaCurve)[2]),
                &((*SigmaCurve)[3]),
                &((*SigmaCurve)[5]),
                &((*SigmaCurve)[4]),
                NumSigmaRows,
                &((*TauCurve)[0]),
                &((*TauCurve)[1]),
                NumTauRows);
            if (err)
            {
                for (i = 0; i < (*NumSigmaCols); i++)
                    srt_free(*SigmaCurve[i]);
                srt_free(*SigmaCurve);
                (*SigmaCurve) = NULL;

                for (i = 0; i < (*NumTauCols); i++)
                    srt_free(*TauCurve[i]);
                srt_free(*TauCurve);
                (*TauCurve) = NULL;
                return err;
            }
        }
        else
        {
            *NumSigmaCols = 3;
            (*SigmaCurve) = (double**)srt_malloc((*NumSigmaCols) * sizeof(double*));
            *NumTauCols   = 2;
            (*TauCurve)   = (double**)srt_malloc((*NumTauCols) * sizeof(double*));

            err = srt_f_display_IRM_OneFac_TermStruct(
                ts,
                &((*SigmaCurve)[0]),
                &((*SigmaCurve)[1]),
                &((*SigmaCurve)[2]),
                &temp1,
                &temp2,
                &temp3,
                NumSigmaRows,
                &((*TauCurve)[0]),
                &((*TauCurve)[1]),
                NumTauRows);

            srt_free(temp1);
            srt_free(temp2);
            srt_free(temp3);

            if (err)
            {
                for (i = 0; i < (*NumSigmaCols); i++)
                    srt_free(*SigmaCurve[i]);
                srt_free(*SigmaCurve);
                (*SigmaCurve) = NULL;

                for (i = 0; i < (*NumTauCols); i++)
                    srt_free(*TauCurve[i]);
                srt_free(*TauCurve);
                (*TauCurve) = NULL;

                return err;
            }
        }

        *Alpha = 0.0;
        *Gamma = 0.0;
        *Rho   = 0.0;

        *Beta = (*SigmaCurve)[1][0];

        *VoVol = find_vovol(0.0, ts);

    } /* END if (mdl_dim == ONE_FAC) */
    else
    {
        *NumSigmaCols = 6;
        (*SigmaCurve) = (double**)malloc((*NumSigmaCols) * sizeof(double*));
        *NumTauCols   = 3;
        (*TauCurve)   = (double**)malloc((*NumTauCols) * sizeof(double*));

        err = srt_f_display_IRM_TwoFac_TermStruct(
            ts,
            &((*SigmaCurve)[0]),
            &((*SigmaCurve)[1]),
            &((*SigmaCurve)[2]),
            &((*SigmaCurve)[3]),
            &((*SigmaCurve)[4]),
            &((*SigmaCurve)[5]),
            NumSigmaRows,
            &((*TauCurve)[0]),
            &((*TauCurve)[1]),
            &((*TauCurve)[2]),
            NumTauRows);
        if (err)
        {
            for (i = 0; i < (*NumSigmaCols); i++)
                srt_free(*SigmaCurve[i]);
            srt_free(*SigmaCurve);
            (*SigmaCurve) = NULL;

            for (i = 0; i < (*NumTauCols); i++)
                srt_free(*TauCurve[i]);
            srt_free(*TauCurve);
            (*TauCurve) = NULL;
            return err;
        }

        *Alpha = (*SigmaCurve)[3][0] / (*SigmaCurve)[1][0];
        *Gamma = 1.0 / (*TauCurve)[2][0] - 1.0 / (*TauCurve)[1][0];
        *Rho   = (*SigmaCurve)[5][0];

        *Beta  = (*SigmaCurve)[2][0];
        *Omega = (*SigmaCurve)[4][0] - (*SigmaCurve)[2][0];

        *VoVol = 0.0;

    } /* END if (mdl_dim == TWO_FAC) */

    /* Return a success message */
    return err;

} /* END Err SrtDisplayUndTermStruct(...) */

/* ---------------------------------------------------------------------------------------- */