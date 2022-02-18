/* ========================================================================

        FILE NAME:  srt_f_calibpsCalibParamss.c

        PURPOSE:    functions to set calibration psCalibParamseters,
                    by default or according to user inputs

        FUNCTION:   Err srt_f_set_CalibpsCalibParamss(...)

   ======================================================================== */

#include "srt_h_all.h"
#ifdef PVMPI
#include "parallel.h"
#endif

static Err srt_f_set_default_calib_params(SrtMdlType eModelType, SrtCalibParam* psCalibParams)
{
    Err err = NULL;

    if (!psCalibParams)
        return serror("Null psCalibParams pointer passed in default_CalibpsCalibParamss");

    memset(psCalibParams, 0, sizeof(SrtCalibParam));

    /* The Default Algorithm is Levenberg - 5 iter */
    psCalibParams->eAlgoType = LEVENBERG_MARQUARDT;
    psCalibParams->lNumIter  = 5;

    /* For calibration using Sobenberg */
    psCalibParams->lNumSobenPoints     = 200;
    psCalibParams->lNumBestSobenPoints = 50;
    psCalibParams->lNumClusters        = 5;
    psCalibParams->szSobenMethod       = "BARC";

    /* For calibration using Simplex or Annealing */
    psCalibParams->dSimplexInitScale = 0.120;
    psCalibParams->dSimplexTol       = 1.0e-8;
    psCalibParams->dMaxTemp          = 1000;
    psCalibParams->dMinTemp          = 1.0e-8;
    psCalibParams->dDecreaseFactor   = 20.00;

    /* For pre-sampling using Sobol */
    psCalibParams->lNumSobolPaths = 0;

    /* Use LGM as a starting point or as a Fixed Point approximation */
    psCalibParams->bLgmStartPoint = SRT_YES;
    psCalibParams->lLgmNumIter    = 5;
    psCalibParams->eLgmAlgoType   = LEVENBERG_MARQUARDT;

    /* Calibration forcing a smooth sigma */
    psCalibParams->bSmoothSigma       = SRT_NO;
    psCalibParams->dSmoothSigmaWeight = 0.10;

    /* By default, the correlation will be fitted as well */
    psCalibParams->eCalibType = GLOBAL_CALIB;

    /* By default, use the old method one */
    psCalibParams->bAggregate = SRT_NO;

    /* The way the Tau, Omega and psCalibParamseters are optimised or kept frozen */
    psCalibParams->bOneTau      = SRT_YES;
    psCalibParams->bFreezeTau   = SRT_NO;
    psCalibParams->bOneBeta     = SRT_YES;
    psCalibParams->bFreezeBeta  = SRT_YES;
    psCalibParams->bOneOmega    = SRT_YES;
    psCalibParams->bFreezeOmega = SRT_YES;

    if (eModelType == LGM)
    {
        psCalibParams->dSigmaMin = 1.00e-04;
        psCalibParams->dSigmaMax = 0.20;
    }

    else
    {
        psCalibParams->dSigmaMin = 1.00e-04;
        psCalibParams->dSigmaMax = 3.00;
    }

    psCalibParams->dLambdaMin = -1.00;
    psCalibParams->dLambdaMax = 1.00;

    psCalibParams->dOptionsWeight = 0.90;

    psCalibParams->dBetaMin = 0.00;
    psCalibParams->dBetaMax = 1.00;

    psCalibParams->dAlphaMin = 0.0001;
    psCalibParams->dAlphaMax = 1000.00;

    psCalibParams->dGammaMin = 0;
    psCalibParams->dGammaMax = 1.0;

    psCalibParams->dRhoMin = -1.00;
    psCalibParams->dRhoMax = 1.00;

    psCalibParams->dOmegaMin = -1.00;
    psCalibParams->dOmegaMax = 0.00;

    return err;
}

/* ------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------
        According to the input taken from the spreadsheet, fill sin the
        SrtCalibpsCalibParams accordingly
   ------------------------------------------------------------------------ */

Err srt_f_set_CalibParams(
    String*        pszParamNames,
    String*        pszParamValues,
    int            iNumCalibParamss,
    SrtCalibParam* psCalibParams,
    SrtMdlType     eModelType)
{
    Err              err = NULL;
    int              i;
    long             tmpLng;
    double           tmpDbl;
    SrtCalibAlgoType calibalgo;
    SrtCalibType     calibtype;

    if (err = srt_f_set_default_calib_params(eModelType, psCalibParams))
        return err;

    if ((!iNumCalibParamss) || (!pszParamNames) || (!pszParamValues))
        return NULL;

    /* Make all the string capital wihout spaces */
    for (i = 0; i < iNumCalibParamss; i++)
    {
        strupper(pszParamNames[i]);
        strip_white_space(pszParamNames[i]);
    }

    /* If no psCalibParamseters entered, numpsCalibParamss has to be 0 */
    for (i = 0; i < iNumCalibParamss; i++)
    {
        /* Algorithm */
        if ((!strcmp(pszParamNames[i], "ALGO")) || (!strcmp(pszParamNames[i], "ALGORITHM")) ||
            (!strcmp(pszParamNames[i], "CALIBALGO")))
        {
            err = srt_f_interp_calibalgo(pszParamValues[i], &calibalgo);
            if (err)
            {
                return err;
            }
            else
            {
                psCalibParams->eAlgoType = calibalgo;
            }
        }
        else
            /* Type */
            if ((!strcmp(pszParamNames[i], "TYPE")) || (!strcmp(pszParamNames[i], "CALIBTYPE")) ||
                (!strcmp(pszParamNames[i], "CORRELATION")))
        {
            err = srt_f_interp_calibtype(pszParamValues[i], &calibtype);
            if (err)
            {
                return err;
            }
            else
            {
                psCalibParams->eCalibType = calibtype;
            }
        }
        else
            /* Number of iterations used in the algorithm */
            if ((!strcmp(pszParamNames[i], "NITER")) || (!strcmp(pszParamNames[i], "ITERNUMBER")) ||
                (!strcmp(pszParamNames[i], "ITERNUMB")))
        {
            if (sscanf(pszParamValues[i], "%ld", &tmpLng) == 1)
            {
                if (tmpLng >= 0)
                {
                    psCalibParams->lNumIter = tmpLng;
                }
                else
                {
                    return serror("NITER must be >= 0");
                }
            }
            else
            {
                return serror("Failed to set NITER");
            }
        }
        else
            /* Number of points in Sobeneberg */
            if ((!strcmp(pszParamNames[i], "NPTS")) || (!strcmp(pszParamNames[i], "PTS")) ||
                (!strcmp(pszParamNames[i], "SOBNBGPTS")))
        {
            if (sscanf(pszParamValues[i], "%ld", &tmpLng) == 1)
            {
                if (tmpLng > 0)
                {
                    psCalibParams->lNumSobenPoints = tmpLng;
                }
                else
                {
                    return serror("NPTS must be > 0");
                }
            }
            else
            {
                return serror("Failed to set NPTS");
            }
        }
        else
            /* Number of best points in Sobeneberg */
            if ((!strcmp(pszParamNames[i], "NBPTS")) || (!strcmp(pszParamNames[i], "BESTPTS")) ||
                (!strcmp(pszParamNames[i], "SOBNBGBPTS")))
        {
            if (sscanf(pszParamValues[i], "%ld", &tmpLng) == 1)
            {
                if (tmpLng > 0)
                {
                    psCalibParams->lNumSobenPoints = tmpLng;
                }
                else
                {
                    return serror("NBPTS must be > 0");
                }
            }
            else
            {
                return serror("Failed to set NBPTS");
            }
        }
        else
            /* Number of clusters in Sobeneberg */
            if ((!strcmp(pszParamNames[i], "NCLUST")) || (!strcmp(pszParamNames[i], "CLUSTERS")) ||
                (!strcmp(pszParamNames[i], "SOBNBGCLUST")))
        {
            if (sscanf(pszParamValues[i], "%ld", &tmpLng) == 1)
            {
                if (tmpLng > 0)
                {
                    psCalibParams->lNumClusters = tmpLng;
                }
                else
                {
                    return serror("NCLUST must be > 0");
                }
            }
            else
            {
                return serror("Failed to set NCLUST");
            }
        }
        else
            /* Research method in Sobenberg */
            if ((!strcmp(pszParamNames[i], "RSCMTH")) || (!strcmp(pszParamNames[i], "RESEARCH")))
        {
            if ((!(strcmp(pszParamValues[i], "BARC"))) &&
                (!(strcmp(pszParamValues[i], "BESTFIT"))) &&
                (!(strcmp(pszParamValues[i], "BESTPOINT"))))
            {
                return serror("Unknown research method: %s", pszParamValues[i]);
            }

            else
            {
                psCalibParams->szSobenMethod = pszParamValues[i];
            }
        }
        else
            /* Initial scaling */
            if ((!strcmp(pszParamNames[i], "SCALING")) ||
                (!strcmp(pszParamNames[i], "INIT_SCALING")) ||
                (!strcmp(pszParamNames[i], "INITSCALING")))
        {
            if (sscanf(pszParamValues[i], "%lf", &tmpDbl) == 1)
            {
                if (tmpDbl > 0)
                {
                    psCalibParams->dSimplexInitScale = tmpDbl;
                }
                else
                {
                    return serror("INITSCALING must be > 0");
                }
            }
            else
            {
                return serror("Failed to set INITSCALING");
            }
        }
        else
            /* Fractional tolerance */
            if ((!strcmp(pszParamNames[i], "FTOL")) || (!strcmp(pszParamNames[i], "TOLERANCE")) ||
                (!strcmp(pszParamNames[i], "FRACT_TOL")))
        {
            if (sscanf(pszParamValues[i], "%lf", &tmpDbl) == 1)
            {
                if (tmpDbl > 0)
                {
                    psCalibParams->dSimplexTol = tmpDbl;
                }
                else
                {
                    return serror("FTOL must be > 0");
                }
            }
            else
            {
                return serror("Failed to set FTOL");
            }
        }
        else
            /* Initial temperature */
            if ((!strcmp(pszParamNames[i], "MAXTEMP")) || (!strcmp(pszParamNames[i], "INITTEMP")))
        {
            if (sscanf(pszParamValues[i], "%lf", &tmpDbl) == 1)
            {
                if (tmpDbl > 0)
                {
                    psCalibParams->dMaxTemp = tmpDbl;
                }
                else
                {
                    return serror("INITTEMP must be > 0");
                }
            }
            else
            {
                return serror("Failed to set INITTEMP");
            }
        }
        else
            /* Final temperature */
            if ((!strcmp(pszParamNames[i], "MINTEMP")) || (!strcmp(pszParamNames[i], "FINALTEMP")))
        {
            if (sscanf(pszParamValues[i], "%lf", &tmpDbl) == 1)
            {
                if (tmpDbl > 0)
                {
                    psCalibParams->dMinTemp = tmpDbl;
                }
                else
                {
                    return serror("FINALTEMP must be > 0");
                }
            }
            else
            {
                return serror("Failed to set FINALTEMP");
            }
        }
        else
            /* Decreasing factor */
            if ((!strcmp(pszParamNames[i], "DECFACT")) ||
                (!strcmp(pszParamNames[i], "DECREASINGFACTOR")) ||
                (!strcmp(pszParamNames[i], "ANNEALINGFACTOR")))
        {
            if (sscanf(pszParamValues[i], "%lf", &tmpDbl) == 1)
            {
                if (tmpDbl > 0)
                {
                    psCalibParams->dDecreaseFactor = tmpDbl;
                }
                else
                {
                    return serror("DECFACT must be > 0");
                }
            }
            else
            {
                return serror("Failed to set DECFACT");
            }
        }
        else
            /* Use only one tau date for calibration */
            if ((!strcmp(pszParamNames[i], "ONETAU")) || (!strcmp(pszParamNames[i], "SINGLETAU")) ||
                (!strcmp(pszParamNames[i], "ONE_TAU")))
        {
            if (!strcmp(pszParamValues[i], "YES"))
            {
                psCalibParams->bOneTau = SRT_YES;
            }
            else if (!strcmp(pszParamValues[i], "NO"))
            {
                psCalibParams->bOneTau = SRT_NO;
            }
            else
            {
                return serror("ONETAU must be YES or NO");
            }
        }
        else
            /* Freeze tau structure */
            if ((!strcmp(pszParamNames[i], "FREEZETAU")) ||
                (!strcmp(pszParamNames[i], "FREEZE_TAU")))
        {
            if (!strcmp(pszParamValues[i], "YES"))
            {
                psCalibParams->bFreezeTau = SRT_YES;
            }
            else if (!strcmp(pszParamValues[i], "NO"))
            {
                psCalibParams->bFreezeTau = SRT_NO;
            }
            else
            {
                return serror("FREEZETAU must be YES or NO");
            }
        }
        else
            /* Number of Sobol paths */
            if ((!strcmp(pszParamNames[i], "SOBOL")) || (!strcmp(pszParamNames[i], "SOBOLPATH")))
        {
            if (sscanf(pszParamValues[i], "%ld", &tmpLng) == 1)
            {
                if (tmpLng >= 0)
                {
                    psCalibParams->lNumSobolPaths = tmpLng;
                }
                else
                {
                    return serror("SOBOL must be >= 0");
                }
            }
            else
            {
                return serror("Failed to set SOBOLPATH");
            }
        }
        else
            /* FRD weight in the calibration criteria */
            if ((!strcmp(pszParamNames[i], "FRD_WEIGHT")) ||
                (!strcmp(pszParamNames[i], "FRDWEIGHT")))
        {
            if (sscanf(pszParamValues[i], "%lf", &tmpDbl) == 1)
            {
                if ((tmpDbl >= 0.0) && (tmpDbl <= 1.0))
                {
                    psCalibParams->dOptionsWeight = tmpDbl;
                }
                else
                {
                    return serror("FRD_WEIGHT must be 0 <= ... <= 1");
                }
            }
            else
            {
                return serror("Failed to set FRDWEIGHT");
            }
        }
        else
            /* Sigma Minimum*/
            if ((!strcmp(pszParamNames[i], "SIGMAMIN")) || (!strcmp(pszParamNames[i], "SIGMA_MIN")))
        {
            if (sscanf(pszParamValues[i], "%lf", &tmpDbl) == 1)
            {
                if (tmpDbl >= 0)
                {
                    psCalibParams->dSigmaMin = tmpDbl;
                }
                else
                {
                    return serror("SIGMAMIN must be >= 0");
                }
            }
            else
            {
                return serror("Failed to set SIGMAMIN");
            }
        }
        else
            /* Sigma Maximum*/
            if ((!strcmp(pszParamNames[i], "SIGMAMAX")) || (!strcmp(pszParamNames[i], "SIGMA_MAX")))
        {
            if (sscanf(pszParamValues[i], "%lf", &tmpDbl) == 1)
            {
                if (tmpDbl >= 0)
                {
                    psCalibParams->dSigmaMax = tmpDbl;
                }
                else
                {
                    return serror("SIGMAMAX must be >= 0");
                }
            }
            else
            {
                return serror("Failed to set SIGMAMIN");
            }
        }
        else
            /* Lambda Minimum*/
            if ((!strcmp(pszParamNames[i], "LAMBDAMIN")) ||
                (!strcmp(pszParamNames[i], "LAMBDA_MIN")))
        {
            if (sscanf(pszParamValues[i], "%lf", &tmpDbl) == 1)
            {
                psCalibParams->dLambdaMin = tmpDbl;
            }
            else
            {
                return serror("Failed to set LAMBDAMIN");
            }
        }
        else
            /* Lambda Maximum*/
            if ((!strcmp(pszParamNames[i], "LAMBDAMAX")) ||
                (!strcmp(pszParamNames[i], "LAMBDA_MAX")))
        {
            if (sscanf(pszParamValues[i], "%lf", &tmpDbl) == 1)
            {
                psCalibParams->dLambdaMax = tmpDbl;
            }
            else
            {
                return serror("Failed to set LAMBDAMAX");
            }
        }
        else
            /* Beta Minimum*/
            if ((!strcmp(pszParamNames[i], "BETAMIN")) || (!strcmp(pszParamNames[i], "BETA_MIN")))
        {
            if (sscanf(pszParamValues[i], "%lf", &tmpDbl) == 1)
            {
                psCalibParams->dBetaMin = tmpDbl;
            }
            else
            {
                return serror("Failed to set BETAMIN");
            }
        }
        else
            /* Beta Maximum*/
            if ((!strcmp(pszParamNames[i], "BETAMAX")) || (!strcmp(pszParamNames[i], "BETA_MAX")))
        {
            if (sscanf(pszParamValues[i], "%lf", &tmpDbl) == 1)
            {
                psCalibParams->dBetaMax = tmpDbl;
            }
            else
            {
                return serror("Failed to set BETAMAX");
            }
        }
        else
            /* Alpha Minimum*/
            if ((!strcmp(pszParamNames[i], "ALPHAMIN")) || (!strcmp(pszParamNames[i], "ALPHA_MIN")))
        {
            if (sscanf(pszParamValues[i], "%lf", &tmpDbl) == 1)
            {
                if (tmpDbl >= 0.0)
                {
                    psCalibParams->dAlphaMin = tmpDbl;
                }
                else
                {
                    return serror("ALPHAMIN must be >= 0.0 ");
                }
            }
            else
            {
                return serror("Failed to set ALPHAMIN");
            }
        }
        else
            /* Alpha Maximum*/
            if ((!strcmp(pszParamNames[i], "ALPHAMAX")) || (!strcmp(pszParamNames[i], "ALPHA_MAX")))
        {
            if (sscanf(pszParamValues[i], "%lf", &tmpDbl) == 1)
            {
                if (tmpDbl >= 0.0)
                {
                    psCalibParams->dAlphaMax = tmpDbl;
                }
                else
                {
                    return serror("ALPHAMAX must be >= 0.0 ");
                }
            }
            else
            {
                return serror("Failed to set ALPHAMAX");
            }
        }
        else
            /* Gamma Minimum*/
            if ((!strcmp(pszParamNames[i], "GAMMAMIN")) || (!strcmp(pszParamNames[i], "GAMMA_MIN")))
        {
            if (sscanf(pszParamValues[i], "%lf", &tmpDbl) == 1)
            {
                psCalibParams->dBetaMin = tmpDbl;
            }
            else
            {
                return serror("Failed to set GAMMAMIN");
            }
        }
        else
            /* Beta Maximum*/
            if ((!strcmp(pszParamNames[i], "GAMMAMAX")) || (!strcmp(pszParamNames[i], "GAMMAMAX")))
        {
            if (sscanf(pszParamValues[i], "%lf", &tmpDbl) == 1)
            {
                psCalibParams->dGammaMax = tmpDbl;
            }
            else
            {
                return serror("Failed to set GAMMAMAX");
            }
        }
        else
            /* Rho Minimum*/
            if ((!strcmp(pszParamNames[i], "RHOMIN")) || (!strcmp(pszParamNames[i], "RHO_MIN")))
        {
            if (sscanf(pszParamValues[i], "%lf", &tmpDbl) == 1)
            {
                if ((tmpDbl >= -1.0) && (tmpDbl <= 1.0))
                {
                    psCalibParams->dRhoMin = tmpDbl;
                }
                else
                {
                    return serror("RHOMIN must be -1.0 <= ... < 1.0");
                }
            }
            else
            {
                return serror("Failed to set RHOMIN");
            }
        }
        else
            /* Rho Maximum*/
            if ((!strcmp(pszParamNames[i], "RHOMAX")) || (!strcmp(pszParamNames[i], "RHO_MAX")))
        {
            if (sscanf(pszParamValues[i], "%lf", &tmpDbl) == 1)
            {
                if ((tmpDbl >= -1.0) && (tmpDbl <= 1.0))
                {
                    psCalibParams->dRhoMax = tmpDbl;
                }
                else
                {
                    return serror("RHOMAX must be -1.0 <= ... < 1.0");
                }
            }
            else
            {
                return serror("Failed to set RHOMAX");
            }
        }
        else
            /* Omega Minimum*/
            if ((!strcmp(pszParamNames[i], "OMEGAMIN")) || (!strcmp(pszParamNames[i], "OMEGA_MIN")))
        {
            if (sscanf(pszParamValues[i], "%lf", &tmpDbl) == 1)
            {
                if ((tmpDbl >= -1.0) && (tmpDbl <= 1.0))
                {
                    psCalibParams->dOmegaMin = tmpDbl;
                }
                else
                {
                    return serror("OMEGAMIN must be -1.0 <= ... < 1.0");
                }
            }
            else
            {
                return serror("Failed to set OMEGAMIN");
            }
        }
        else
            /* Omega Maximum*/
            if ((!strcmp(pszParamNames[i], "OMEGAMAX")) || (!strcmp(pszParamNames[i], "OMEGAMAX")))
        {
            if (sscanf(pszParamValues[i], "%lf", &tmpDbl) == 1)
            {
                if ((tmpDbl >= -1.0) && (tmpDbl <= 1.0))
                {
                    psCalibParams->dOmegaMax = tmpDbl;
                }
                else
                {
                    return serror("OMEGAMAX must be -1.0 <= ... < 1.0");
                }
            }
            else
            {
                return serror("Failed to set OMEGAMAX");
            }
        }
        else
            /* LGM Starting Point */
            if ((!strcmp(pszParamNames[i], "LGMSTARTINGPOINT")) ||
                (!strcmp(pszParamNames[i], "LGM_STARTING_POINT")))
        {
            if (!strcmp(pszParamValues[i], "YES"))
            {
                psCalibParams->bLgmStartPoint = SRT_YES;
            }
            else if (!strcmp(pszParamValues[i], "NO"))
            {
                psCalibParams->bLgmStartPoint = SRT_NO;
            }
            else
            {
                return serror("LGMSTARTINGPOINT must be YES or NO");
            }
        }
        else
            /* LGM iterations */
            if ((!strcmp(pszParamNames[i], "LGMITER")) ||
                (!strcmp(pszParamNames[i], "NUMLGMITER")) ||
                (!strcmp(pszParamNames[i], "LGMNITER")) || (!strcmp(pszParamNames[i], "LGM_ITER")))
        {
            if (sscanf(pszParamValues[i], "%ld", &tmpLng) == 1)
            {
                if (tmpLng > 0)
                {
                    psCalibParams->lLgmNumIter = tmpLng;
                }
                else
                {
                    return serror("LGMITER must be > 0");
                }
            }
            else
            {
                return serror("Failed to set LGMITER");
            }
        }
        else
            /* LGM Algorithm */
            if ((!strcmp(pszParamNames[i], "LGMALGO")) ||
                (!strcmp(pszParamNames[i], "LGMALGORITHM")) ||
                (!strcmp(pszParamNames[i], "ALGOLGM")) ||
                (!strcmp(pszParamNames[i], "ALGORITHM_LGM")))
        {
            err = srt_f_interp_calibalgo(pszParamValues[i], &calibalgo);
            if (err)
            {
                return err;
            }
            else
            {
                psCalibParams->eLgmAlgoType = calibalgo;
            }
        }
        else
            /* Smooth Sigma */
            if ((!strcmp(pszParamNames[i], "SMOOTHSIGMA")) ||
                (!strcmp(pszParamNames[i], "SMOOTH_SIGMA")))
        {
            if (!strcmp(pszParamValues[i], "YES"))
            {
                psCalibParams->bSmoothSigma = SRT_YES;
            }
            else if (!strcmp(pszParamValues[i], "NO"))
            {
                psCalibParams->bSmoothSigma = SRT_NO;
            }
            else
            {
                return serror("SMOOTHSIGMA must be YES or NO");
            }
        }
        else
            /* Smooth sigma weight in the calibration criteria */
            if ((!strcmp(pszParamNames[i], "SMOOTHSIGWEIGHT")) ||
                (!strcmp(pszParamNames[i], "SMOOTH_SIG_WEIGHT")))
        {
            if (sscanf(pszParamValues[i], "%lf", &tmpDbl) == 1)
            {
                if ((tmpDbl >= 0.0) && (tmpDbl <= 1.0))
                {
                    psCalibParams->dSmoothSigmaWeight = tmpDbl;
                }
                else
                {
                    return serror("SMOOTHSIGWEIGHT must be 0 <= ... <= 1");
                }
            }
            else
            {
                return serror("Failed to set SMOOTHSIGWEIGHT");
            }
        }
        else
            /* Aggregate */
            if ((!strcmp(pszParamNames[i], "AGGREGATE")) ||
                (!strcmp(pszParamNames[i], "ALLINONE")) ||
                (!strcmp(pszParamNames[i], "ALL_IN_ONE")) || (!strcmp(pszParamNames[i], "ALLIN1")))
        {
            if (!strcmp(pszParamValues[i], "YES"))
            {
                psCalibParams->bAggregate = SRT_YES;
#ifdef PVMPI
                PRL_CONTEXT.bAggregate = 1;
#endif
            }
            else if (!strcmp(pszParamValues[i], "NO"))
            {
                psCalibParams->bAggregate = SRT_NO;
#ifdef PVMPI
                PRL_CONTEXT.bAggregate = 0;
#endif
            }
            else
            {
                return serror("AGGREGATE must be YES or NO");
            }
        }
        /* Error */
        else
        {
            return serror("Model psCalibParamseter not recognised: %s", pszParamNames[i]);
        }

    } /* END loop over numpsCalibParamss */

    /* Checks that the minimum is below the maximum */
    if (psCalibParams->dSigmaMin > psCalibParams->dSigmaMax)
    {
        return serror(
            "Sigma_min(%f) is bigger than sigma_max(%f)",
            psCalibParams->dSigmaMin,
            psCalibParams->dSigmaMax);
    }
    if (psCalibParams->dLambdaMin > psCalibParams->dLambdaMax)
    {
        return serror(
            "Lambda_min(%f) is bigger than Lambda_max(%f)",
            psCalibParams->dLambdaMin,
            psCalibParams->dLambdaMax);
    }
    if (psCalibParams->dBetaMin > psCalibParams->dBetaMax)
    {
        return serror(
            "Beta_min(%f) is bigger than Beta_max(%f)",
            psCalibParams->dBetaMin,
            psCalibParams->dBetaMax);
    }
    if (psCalibParams->dAlphaMin > psCalibParams->dAlphaMax)
    {
        return serror(
            "Alpha_min(%f) is bigger than Alpha_max(%f)",
            psCalibParams->dAlphaMin,
            psCalibParams->dAlphaMax);
    }
    if (psCalibParams->dGammaMin > psCalibParams->dGammaMax)
    {
        return serror(
            "Gamma_min(%f) is bigger than Gamma_max(%f)",
            psCalibParams->dGammaMin,
            psCalibParams->dGammaMax);
    }
    if (psCalibParams->dRhoMin > psCalibParams->dRhoMax)
    {
        return serror(
            "rho_min(%f) is bigger than rho_max(%f)",
            psCalibParams->dRhoMin,
            psCalibParams->dRhoMax);
    }
    if (psCalibParams->dOmegaMin > psCalibParams->dOmegaMax)
    {
        return serror(
            "Omega_min(%f) is bigger than Omega_max(%f)",
            psCalibParams->dOmegaMin,
            psCalibParams->dOmegaMax);
    }

    /* Returns a success message */
    return NULL;
}

static Err srt_f_set_default_FXCalibParams(SrtFXCalibParam* psFXCalibParams)
{
    Err err = NULL;

    if (!psFXCalibParams)
        return serror("Null psFXCalibParams pointer passed in srt_f_set_default_FXCalibParams");

    memset(psFXCalibParams, 0, sizeof(SrtFXCalibParam));

    /* The Default Algorithm is Levenberg - 5 iter */
    psFXCalibParams->eAlgoType    = LEVENBERG_MARQUARDT;
    psFXCalibParams->lNumIter     = 5;
    psFXCalibParams->bCalibCorrel = SRT_NO;

    return err;
}

Err srt_f_set_FXCalibParams(
    String*          pszFXCalibStrings,
    String*          pszFXCalibValues,
    int              iNumCalibParamss,
    SrtFXCalibParam* psFXCalibParams)

{
    Err              err = NULL;
    int              i;
    long             tmpLng;
    SrtCalibAlgoType calibalgo;

    if (err = srt_f_set_default_FXCalibParams(psFXCalibParams))
        return err;

    if ((!iNumCalibParamss) || (!pszFXCalibStrings) || (!pszFXCalibValues))
        return NULL;

    /* Make all the string capital wihout spaces */
    for (i = 0; i < iNumCalibParamss; i++)
    {
        strupper(pszFXCalibStrings[i]);
        strip_white_space(pszFXCalibStrings[i]);
    }

    for (i = 0; i < iNumCalibParamss; i++)
    {
        /* Algorithm */
        if ((!strcmp(pszFXCalibStrings[i], "ALGO")) ||
            (!strcmp(pszFXCalibStrings[i], "ALGORITHM")) ||
            (!strcmp(pszFXCalibStrings[i], "CALIBALGO")))
        {
            err = srt_f_interp_calibalgo(pszFXCalibValues[i], &calibalgo);
            if (err)
            {
                return err;
            }
            else
            {
                psFXCalibParams->eAlgoType = calibalgo;
            }
        }
        else
            /* Number of iterations used in the algorithm */
            if ((!strcmp(pszFXCalibStrings[i], "NITER")) ||
                (!strcmp(pszFXCalibStrings[i], "ITERNUMBER")) ||
                (!strcmp(pszFXCalibStrings[i], "ITERNUMB")))
        {
            if (sscanf(pszFXCalibValues[i], "%ld", &tmpLng) == 1)
            {
                if (tmpLng >= 0)
                {
                    psFXCalibParams->lNumIter = tmpLng;
                }
                else
                {
                    return serror("NITER must be >= 0");
                }
            }
            else
            {
                return serror("Failed to set NITER");
            }
        }
        else
            /* Calibrate to the correlation structure  */
            if ((!strcmp(pszFXCalibStrings[i], "CALIBCORREL")) ||
                (!strcmp(pszFXCalibStrings[i], "CALIBCORR")) ||
                (!strcmp(pszFXCalibStrings[i], "CORRELCALIB")))
        {
            if (!strcmp(pszFXCalibValues[i], "YES"))
            {
                psFXCalibParams->bCalibCorrel = SRT_YES;
            }
            else if (!strcmp(pszFXCalibValues[i], "NO"))
            {
                psFXCalibParams->bCalibCorrel = SRT_NO;
            }
            else
            {
                return serror("CALIBCORREL must be YES or NO");
            }
        }
        else
        {
            return serror("FXCalib params must be CALIBCORREL, CALIBALGO or NITER");
        }
    }
    return NULL;
}
