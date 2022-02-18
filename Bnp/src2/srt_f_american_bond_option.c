
#include "math.h"
#include "srt_h_all.h"

#define ONE_HALF 0.5
#define TREE_MESH_SPACING 1.73205080756888
#define TREE_MESH_SPACING_SQUARE 3.0
#define TREE_LIM_IN_STDEV 5

static Err srt_f_norm_call(double InOutMoneyness, double Stdev, double* PV)
{
    (*PV) = InOutMoneyness * norm_accurate(InOutMoneyness / Stdev) +
            Stdev * norm(InOutMoneyness / Stdev);

    return NULL;
}

static Err american_bondoption_payofffunc(double InOutMoneyness, double Stdev, double* payoff)
{
    double PV;
    Err    err = NULL;

    err = srt_f_norm_call(-InOutMoneyness, Stdev, &PV);
    if (err)
        return err;

    (*payoff) = InOutMoneyness + Stdev - PV;

    return NULL;
}

Err srt_f_find_american_bond_option_exercise_frontier_index(
    double* InOutMoneyness,
    double  Stdev,
    long    max_index,
    double* PrevExPrePtr,
    double* Payoff,

    long*   imin,
    double* iEx)
{
    Err     err  = NULL;
    long    imax = max_index, IndexMin = -max_index, IndexMax = max_index, j, k;
    double *Array, RefVal = 0.0;

    Array = dvector(IndexMin, IndexMax);

    for (j = IndexMin; j <= IndexMax; j++)
        Array[j] = (PrevExPrePtr[j] - Payoff[j]);

    (*imin) = IndexMin;
    while ((imax - (*imin)) > 1)
    {
        k = (imax + (*imin)) >> 1;
        if (RefVal > Array[k])
            (*imin) = k;
        else if (RefVal < Array[k])
            imax = k;
        else
            (*imin) = k;
    }

    (*iEx) = (*imin) - 1 + Array[(*imin) - 1] / (Array[(*imin) - 1] + Array[(*imin)]);

    if (Array)
        free_dvector(Array, IndexMin, IndexMax);
    Array = NULL;

    return err;
}

Err srt_f_american_bondoption_exercise_premuim(
    Date EvalDate,  /* EVALUATION DATE */
    Date AmExpDate, /* AMERICAN EXPIRY DATE */
    Err (*payofffunc)(),

    /*OUTPUTS */
    long*    MaxIndexOutput,
    double** xExBdryOutput,
    double** EarlyExPr)
{
    Err     err = NULL;
    double *xExBdry, *CurExPrePtr, *PrevExPrePtr, prob[3];
    long    MinNode = 350, MaxNode = 1000, MaxIndex, NumTimeStep, i, j, imin;
    double  TimeStep = 1.0 / 24, dx, Stdev, *InOutMoneyness, *Payoff, iEx, YrtoAmExpDate;

    /*DEFINITION OF THE DISCRETISATION SCHEME */
    YrtoAmExpDate = (double)(AmExpDate - EvalDate) * YEARS_IN_DAY;
    TimeStep      = min(TimeStep, YrtoAmExpDate / MinNode);
    NumTimeStep   = ((long)(YrtoAmExpDate / TimeStep)) + 1;

    dx       = TREE_MESH_SPACING * sqrt(TimeStep); /* TO VERIFY THE STABILITY CONDITION */
    MaxIndex = (long)DTOL(TREE_LIM_IN_STDEV * YrtoAmExpDate / dx) + 1;

    /* ALLOCATIONS OF THE POINTERS */
    xExBdry        = dvector(1, NumTimeStep);
    PrevExPrePtr   = dvector(-MaxIndex, MaxIndex);
    CurExPrePtr    = dvector(-MaxIndex, MaxIndex);
    Payoff         = dvector(-MaxIndex, MaxIndex);
    InOutMoneyness = dvector(-MaxIndex, MaxIndex);

    for (j = -MaxIndex; j <= MaxIndex; j++)
        InOutMoneyness[j] = j * dx;

    /* SET THE TRANSITION PROBABILITIES */
    prob[0] = ONE_HALF * ONE_HALF * TREE_MESH_SPACING;
    prob[2] = ONE_HALF * ONE_HALF * TREE_MESH_SPACING;
    prob[1] = 1.0 - prob[0] - prob[2];

    /* START OF THE FORWARD ALGORITHM */
    for (i = 1; i <= NumTimeStep; i++) /* LOOP ON THE TIME STEP */
    {
        Stdev = sqrt(TimeStep * i);
        for (j = -MaxIndex + 1; j <= (MaxIndex - 1); j++) /* LOOP ON THE SPACE POINTS */
        {
            CurExPrePtr[j] = prob[0] * PrevExPrePtr[j - 1] + prob[1] * PrevExPrePtr[j] +
                             prob[2] * PrevExPrePtr[j + 1];

            err = (*payofffunc)(InOutMoneyness[j], Stdev, &Payoff[j]);
            if (err)
                return err;

            CurExPrePtr[j]  = max(PrevExPrePtr[j], Payoff[j]);
            PrevExPrePtr[j] = CurExPrePtr[j];

        } /*END OF LOOP ON THE SPACE POINTS */

        PrevExPrePtr[-MaxIndex] = PrevExPrePtr[-MaxIndex + 1];
        PrevExPrePtr[MaxIndex]  = 0.0;

        if (i > 1)
        {
            err = srt_f_find_american_bond_option_exercise_frontier_index(
                InOutMoneyness, Stdev, MaxIndex, PrevExPrePtr, Payoff, &imin, &iEx);
            if (err)
                return (err);

            xExBdry[i] = iEx * dx;
        }
    }

    if (PrevExPrePtr)
        free_dvector(PrevExPrePtr, -MaxIndex, MaxIndex);
    PrevExPrePtr = NULL;
    if (Payoff)
        free_dvector(Payoff, -MaxIndex, MaxIndex);
    Payoff = NULL;
    if (InOutMoneyness)
        free_dvector(InOutMoneyness, -MaxIndex, MaxIndex);
    InOutMoneyness = NULL;

    (*MaxIndexOutput) = MaxIndex;
    (*xExBdryOutput)  = xExBdry;
    (*EarlyExPr)      = CurExPrePtr;

    return (err);
}

Err XLAmericanBondOptionExPremuim(
    long EvalDate,
    long AmExpDate,

    long*    MaxIndex,
    double** xExBdryOutput,
    double** EarlyExPr)
{
    Err err = NULL;

    err = srt_f_american_bondoption_exercise_premuim(
        EvalDate,
        AmExpDate,
        &american_bondoption_payofffunc,

        MaxIndex,
        xExBdryOutput,
        EarlyExPr);
    if (err)
        return err;

    return err;
}
