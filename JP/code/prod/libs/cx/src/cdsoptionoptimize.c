/*
***************************************************************************
** FILENAME: cdsoptionoptimize.c
**
** Multi-Q smile optimization routines for CDS options.
***************************************************************************
*/

#include "cdsoptionoptimize.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <alib/dtivlo.h>
#include <common/include/drutils.h>
#include <crxmultiq/include/crmultiq.h>
#include <crxflow/include/crcrv.h>

#include <alib/lintrp.h>
#include <alib/rtbrent.h>

#include <cxutils/include/alibconv.h>
#include "cds.h"
#include "cdsbootstrap.h"
#include <cxutils/include/mdmin.h>
#include <cxutils/include/hash.h>
#include <cxutils/include/cxmacros.h>

/*#include "private.h"*/

/*
***************************************************************************
** Multi-Q parameter optimization initialization.
**
** Initializes the state object with the initial guess of the ATM vol and
** Q-distribution.
***************************************************************************
*/
CrxTCdsOptionQOptimization* CrxCdsOptionQOptimizationInit(
    CrxTQDist*      qdist,
    double          vol)
{
    static char routine[] = "CrxCdsOptionQOptimizationInit";

    CrxTCdsOptionQOptimization* state = CrxCdsOptionQOptimizationMake(
        qdist,
        vol,
        NULL,
        0,
        0.0,
        0.0);
    
    if (state == NULL)
        GtoErrMsgFailure (routine);

    return state;
}

typedef struct
{
    int  numVars;     /* Number of independent variables */
    int  volVarIdx;   /* -1 = use input vol, otherwise value within vars */
    int  qVarIdx[6];  /* -N = use existing Q-value of that number
                         >= 0 = index with vars */
} VAR_INDEX;

typedef struct
{
    CrxTQDist*                  qdist;
    double                      vol;
    VAR_INDEX*                  varIndex;
    CrxTCdsOption*              baseDeal;
    CrxTCdsOptionQOptModel*     model;
    TDate                       today;
    int                         numOptions;
    double*                     strikes;
    long*                       optionTypes;
    double*                     weights;
    double*                     strikeVols;
    TCurve*                     discCurve;
    CxTCreditCurve*             sprdCurve;
    double                      recoveryRate;
} QOPT_DATA;

static int rmsCalc (int n, double *x, QOPT_DATA *data, double *rms);

static VAR_INDEX* varIndexMake(
    CrxTCdsOptionQOptModel *model,
    int                     nQs);

static void setVars (VAR_INDEX  *varIndex,
                     double     *vars,
                     double      vol,
                     CrxTQDist  *qdist);

static void getVars (VAR_INDEX  *varIndex,
                     double     *vars,
                     double     *vol,
                     CrxTQDist  *qdist);

/*
***************************************************************************
** Multi-Q parameter optimization.
**
** Designed so that we can pass the results back as an input to allow the
** function to be called in a loop for gradual improvements.
***************************************************************************
*/
CrxTCdsOptionQOptimization* CrxCdsOptionQOptimization(
    CrxTCdsOptionQOptimization* initState,
    CrxTCdsOption*              baseDeal,
    CrxTCdsOptionQOptModel*     model,
    TDate                       today,
    int                         numOptions,
    double*                     strikes,
    long*                       optionTypes,
    double*                     weights,
    double*                     strikeVols,
    TCurve*                     discCurve,
    CxTCreditCurve*             sprdCurve,
    double                      recoveryRate,
    char*                       logfilename
)
{
    static char routine[] = "CrxCdsOptionQOptimization";
    int         status    = FAILURE;

    CrxTCdsOptionQOptimization *state = NULL;

    CxTMultiDimMinState *initMinState = NULL;
    CxTMultiDimMinState *minState = NULL;

    int i;
    int numVars; /* number of variables to be optimized */
    double *vars = NULL;
    VAR_INDEX *varIndex = NULL;
    int     numNonZeroOptions = 0;

    double     vol;          /* working copy of vol */
    CrxTQDist *qdist = NULL; /* working copy of qdist */

    QOPT_DATA data;

    REQUIRE (initState != NULL);
    REQUIRE (baseDeal != NULL);
    REQUIRE (model != NULL);
    REQUIRE (model->numQsNoOpt == 0 || model->qsNoOpt != NULL);
    REQUIRE (today > 0);
    REQUIRE (numOptions > 0);
    REQUIRE (strikes != NULL);
    REQUIRE (optionTypes != NULL);
    REQUIRE (weights != NULL);
    REQUIRE (strikeVols != NULL);
    REQUIRE (discCurve != NULL);
    REQUIRE (sprdCurve != NULL);

    vol   = initState->vol;
    qdist = CrxQDistCopy (initState->qdist);
    if (qdist == NULL)
        goto done; /* failure */

    REQUIRE (vol > 0.0);
    REQUIRE (model->numQsNoOpt < qdist->nQs);

    /*
    ** Validation of model vis-a-vis input values.
    */
    varIndex = varIndexMake (model, qdist->nQs);
    if (varIndex == NULL)
        goto done; /* failure */

    numVars = varIndex->numVars;
    for (i = 0; i < numOptions; ++i)
    {
        if (weights[i] > 0.0)
            ++numNonZeroOptions;
    }

    if (numVars > numNonZeroOptions)
    {
        GtoErrMsg ("%s: Number of independent variables (%d) > number of "
                   "non-zero options (%d)\n", routine, numVars,
                   numNonZeroOptions);
        goto done;
    }

    if (initState->direction != NULL)
    {
        if (initState->direction->numDim1 != numVars ||
            initState->direction->numDim2 != numVars)
        {
            GtoErrMsg ("%s: Non-NULL direction matrix should be size [%d,%d] "
                       "and not [%d,%d]\n", routine, numVars, numVars,
                       initState->direction->numDim1,
                       initState->direction->numDim2);
            goto done;
        }
    }

    vars = NEW_ARRAY(double, numVars);
    if (vars == NULL)
        goto done; /* failure */

    /* prepare the optimizer */
    data.qdist          = qdist;
    data.vol            = vol;
    data.varIndex       = varIndex;
    data.baseDeal       = baseDeal;
    data.model          = model;
    data.today          = today;
    data.numOptions     = numOptions;
    data.strikes        = strikes;
    data.optionTypes    = optionTypes;
    data.weights        = weights;
    data.strikeVols     = strikeVols;
    data.discCurve      = discCurve;
    data.sprdCurve      = sprdCurve;
    data.recoveryRate   = recoveryRate;

	setVars (varIndex, vars, vol, qdist);

    initMinState = CxMultiDimMinStateMake (numVars,
                                           vars,
                                           initState->direction,
                                           initState->iter,
                                           initState->value,
                                           initState->vdiff);
    if (initMinState == NULL)
        goto done;

    /* call the optimizer */
    minState = CxMultiDimMinimization (initMinState,
                                       (CxTMultiDimObjFunc)rmsCalc,
                                       NULL, /* dfunc */
                                       &data,
                                       model->tolerance,
                                       model->maxIter,
                                       model->maxTime,
                                       logfilename);
    if (minState == NULL)
    {
        GtoErrMsg ("%s: Minimization for multi-Q failed inside optimizer\n",
                   routine);
        goto done;
    }

    /* set the vol and q-values from minState */
    getVars (varIndex, minState->x, &vol, qdist);

    state = CrxCdsOptionQOptimizationMake (qdist,
                                           vol,
                                           minState->direction,
                                           minState->iter,
                                           minState->value,
                                           minState->vdiff);

    if (state == NULL)
        goto done; /* failure */

    status = SUCCESS;

 done:

    FREE (vars);
    FREE (varIndex);
    CrxQDistFree (qdist);
    CxMultiDimMinStateFree (initMinState);
    CxMultiDimMinStateFree (minState);
 
    if (status != SUCCESS)
    {
        CrxCdsOptionQOptimizationFree (state);
        state = NULL;
        GtoErrMsgFailure (routine);
    }
   
    return state;
}

/*
** Root mean square difference calculator for Q-optimization.
*/
static int rmsCalc (int n, double *x, QOPT_DATA *data, double *rms)
{
    /*static char routine[] = "rmsCalc";*/
    int         status    = FAILURE;
    
    double      err = 0.0;
    double      totalWeight;

    int         i;

    CrxTQDist  *qdist = data->qdist;
    double      vol   = data->vol;

    CrxTCdsOption     *deal = NULL;
    CrxTCdsOptionCalc *priceCalc = NULL;
    CrxTCdsOptionCalc *volCalc = NULL;

    /* populate qdist and vol from x[] */
    getVars (data->varIndex, x, &vol, qdist);
    
    deal = CrxCdsOptionCopy (data->baseDeal);

    totalWeight = 0.0;

    for (i = 0; i < data->numOptions; ++i)
    {
        double weight       = data->weights[i];
        double strikeVol    = data->strikeVols[i];
        double strikeVolDiff;

        if (weight <= 0.0)
            continue;

        totalWeight += weight;

        deal->strike        = data->strikes[i];
        deal->optionType    = data->optionTypes[i];

        CrxCdsOptionCalcFree(priceCalc);
        priceCalc = CrxCdsOptionCalc (deal,
                                      CRX_DIST_TYPE_Q,
                                      data->today,
                                      vol,
                                      qdist,
                                      data->discCurve,
                                      data->sprdCurve,
                                      data->recoveryRate);
        if (priceCalc == NULL)
            goto done;

        CrxCdsOptionCalcFree(volCalc);
        volCalc = CrxCdsOptionVolCalc (deal,
                                       CRX_DIST_TYPE_LOGNORMAL,
                                       data->today,
                                       priceCalc->price,
                                       NULL,
                                       data->discCurve,
                                       data->sprdCurve,
                                       data->recoveryRate);
        if (volCalc == NULL)
            goto done;

        strikeVolDiff = strikeVol - volCalc->vol;
        err += weight * strikeVolDiff * strikeVolDiff;
    }

    *rms = sqrt(err/totalWeight);
    status = SUCCESS;

 done:

    CrxCdsOptionFree(deal);
    CrxCdsOptionCalcFree(priceCalc);
    CrxCdsOptionCalcFree(volCalc);

    if (status != SUCCESS)
        *rms = 1e100; /* i.e. not good - but do not terminate */

    return SUCCESS;
}

/*
** Construct the VAR_INDEX structure which defines a mapping between the
** independent variables array and the economic values of vol and qdist.
*/
static VAR_INDEX* varIndexMake(
    CrxTCdsOptionQOptModel *model,
    int                     nQs)
{
    static char routine[] = "varIndexMake";
    int         status    = FAILURE;

    static char invalidFormat[] = "%s: Invalid format for noOpt '%s'\n";

    int         numVars;
    int         i;
    
    VAR_INDEX *varIndex = NEW(VAR_INDEX);
    if (varIndex == NULL)
        goto done;

    REQUIRE (nQs <= 6);

    for (i = 0; i < model->numQsNoOpt; ++i)
    {
        char *noOpt = model->qsNoOpt[i];
        long  n1;
        long  n2;
        char *ep;
        
        if (noOpt == NULL || *noOpt == '\0')
            continue; /* ignore NULL or empty rules */

        /* two valid forms for noOpt are an integer or integer=integer */
        /* n1 and n2 in this case are both 1-based */
        n1 = strtol (noOpt, &ep, 10);
        if (n1 <= 0)
        {
            GtoErrMsg (invalidFormat, routine, noOpt);
            goto done;
        }
        if (*ep == '\0')
        {
            n2 = 0;
        }
        else if (*ep == '=')
        {
            char *ep2;
            n2 = strtol (ep+1, &ep2, 10);
            if (*ep2 != '\0' || n2 <= 0)
            {
                GtoErrMsg (invalidFormat, routine, noOpt);
                goto done;
            }
        }
        else
        {
            GtoErrMsg (invalidFormat, routine, noOpt);
            goto done;
        }

        if (n1 < 1 || n1 > nQs)
        {
            GtoErrMsg ("%s: noOpt %ld is not in range [1,%d]\n",
                       routine, n1, nQs);
            goto done;
        }

        if (varIndex->qVarIdx[n1-1] != 0)
        {
            GtoErrMsg ("%s: noOpt rule for %ld defined more than once\n",
                       routine, n1);
            goto done;
        }
        
        if (n2 == 0)
        {
            varIndex->qVarIdx[n1-1] = -n1;
        }
        else if (n2 < 0 || n2 > nQs)
        {
            GtoErrMsg ("%s: noOpt %ld for rule %s is not in range [1,%d]\n",
                       routine, n2, noOpt, nQs);
            goto done;
        }
        else
        {
            varIndex->qVarIdx[n1-1] = -n2;
        }
    }

    if (model->optimizeVol)
    {
        varIndex->volVarIdx = 0;
        numVars = 1;
    }
    else
    {
        varIndex->volVarIdx = -1;
        numVars = 0;
    }
    
    for (i = 0; i < nQs; ++i)
    {
        if (varIndex->qVarIdx[i] == 0)
        {
            varIndex->qVarIdx[i] = numVars;
            ++numVars;
        }
    }
    if (numVars <= 0)
    {
        GtoErrMsg ("%s: Invalid number of independent variables (%d)\n",
                   routine, numVars);
        goto done;
    }

    varIndex->numVars = numVars;
    status = SUCCESS;

 done:

    if (status != SUCCESS)
    {
        GtoErrMsgFailure (routine);
        FREE (varIndex);
        varIndex = NULL;
    }
    return varIndex;
}

/*
** Populates vars given varIndex, vol and qdist.
*/
static void setVars (VAR_INDEX  *varIndex,
                     double     *vars,
                     double      vol,
                     CrxTQDist  *qdist)
{
    int i;
    int j;

    if (varIndex->volVarIdx >= 0)
        vars[varIndex->volVarIdx] = vol;

    for (i = 0; i < qdist->nQs; ++i)
    {
        j = varIndex->qVarIdx[i];
        if (j >= 0)
            vars[j] = qdist->Qs[i];
    }
}

/*
** Populates vol, qdist given varIndex and vars.
*/
static void getVars (VAR_INDEX  *varIndex,
                     double     *vars,
                     double     *vol,
                     CrxTQDist  *qdist)
{
    int i;
    int j;

    if (varIndex->volVarIdx >= 0)
        *vol = vars[varIndex->volVarIdx];

    for (i = 0; i < qdist->nQs; ++i)
    {
        j = varIndex->qVarIdx[i];
        if (j >= 0)
            qdist->Qs[i] = vars[j];
    }
    for (i = 0; i < qdist->nQs; ++i)
    {
        j = varIndex->qVarIdx[i];
        if (j < 0)
            /* j=-1 => use q[0], j=-2 => use q[1] etc */
            qdist->Qs[i] = qdist->Qs[-j-1];
    }
}
