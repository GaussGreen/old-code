#include "MBSPrepay.h"

#include "MBSPPFuncs.h"
#include "MBSPTUtil.h"
#include "math.h"

#define ALSMM1FUNC(x, y, N) ((y) * (x) - (1. - (x)) * (1. - pow((1. - (x)), (N - 1))))
#define ALSMM2FUNC(p, N) ((1. - (p)-pow(1. - (p), (N) + 1)) / p - (N)*pow(1. - (p), (N)))
#define small 0.00001
///

// explicit initialization for one one S curve and NBPOINTS knot-points for S-curve
// call MBS_Prepay_Engine_Initialize
char* MBS_Prepay_Engine_Params_ReaderExplicit(
    char*              data_dir,
    MBS_Prepay_Engine* prepay_engine,
    MTGPTACYPROG       AgencyProg,
    long               oterm,
    double             gwac,
    double             Speedup,
    double             rational,
    double             FebruaryEffect,
    double*            IndexLevel,
    double*            AmortLevel,
    double*            SeasoningLevel,
    double*            SeasoningMat,
    int                SFSNb,
    double*            Amorts,
    double             Seasonality[12],
    int                Delay,
    int                NbFixedPrepay,
    double*            FixedPrepay,
    int                nuniform,
    int                IndexMat,
    int                EndIndexMat,
    int                IndexFreq)
{
    char*   err;
    int     i;
    double* triggerRates;
    ////
    for (i = 0; i < NbFixedPrepay; i++)
        FixedPrepay[i] /= 100.0;

    if (err = history_reader(data_dir, &(prepay_engine->history)))
        return (err);
    if (err = density_reader(data_dir, &(prepay_engine->density)))
        return (err);
    triggerRates = mbspt_dvector(0, SFSNb - 1);
    MBS_Prepay_Engine_Initialize(
        prepay_engine,
        1,
        NBPOINTS,
        IndexLevel,
        &(AmortLevel),
        Seasonality,
        FebruaryEffect,
        NUMRAMPS,
        SeasoningLevel,
        SeasoningMat,
        Speedup,
        rational,
        SFSNb,
        triggerRates,
        Amorts,
        Delay,
        NbFixedPrepay,
        FixedPrepay,
        AgencyProg,
        oterm,
        gwac,
        0,
        0,
        nuniform,
        IndexMat,
        EndIndexMat,
        IndexFreq);

    return 0;
}  // MBS_Prepay_Engine_Params_ReaderExplicit

// maximal use of file data
// call MBS_Prepay_Engine_Params_ReaderExplicit
char* MBS_Prepay_Engine_Params_Reader(
    char*              data_dir,
    MBS_Prepay_Engine* prepay_engine,
    MTGPTACYPROG       AgencyProg,
    long               oterm,
    double             gwac,
    double             Speedup,
    double             rational,
    // double FebruaryEffect,
    int IndexMat)
{
    FILE*   stream;
    double  slowAmort, fastAmort;
    long    nuniform;
    char*   err;
    char    string[MAXSTRBUFSZ];
    double  IndexLevel[NBPOINTS], *AmortLevel, SeasoningLevel[NUMRAMPS], Seasonality[12];
    double *FixedPrepay, *Amorts = 0;
    double  SeasoningMat[NUMRAMPS];
    int     SFSNb, i, Delay, NbFixedPrepay;
    int     EndIndexMat, IndexFreq;
    double  FebEffect;

    err = 0;
    //	if( !(stream = MBS_Sequential_File_Manager( data_dir, "PrepayFunc.dat" )) )
    if (!(stream = _OpenMBSFile(data_dir, FilePPFUNCS)))
        return (nrerror("Can't open PrepayFunc.dat"));
    fgets(string, 80, stream);
    fscanf(stream, "%ld\n", &SFSNb);
    fgets(string, 80, stream);
    fscanf(stream, "%lf \n", &(slowAmort));
    fgets(string, 80, stream);
    fscanf(stream, "%lf \n", &(fastAmort));
    fgets(string, 80, stream);
    fscanf(stream, "%ld\n", &nuniform);
    fgets(string, 80, stream);
    fscanf(stream, "%ld\n", &EndIndexMat);
    if (EndIndexMat <= 0)
        EndIndexMat = IndexMat;
    fgets(string, 80, stream);
    fscanf(stream, "%ld\n", &IndexFreq);
    fgets(string, 80, stream);
    fscanf(stream, "%lf \n", &FebEffect);
    fgets(string, 80, stream);
    fgets(string, 80, stream);
    AmortLevel = mbspt_dvector(0, NBPOINTS - 1);
    for (i = 0; i < NBPOINTS; i++)
        fscanf(stream, "%lf %lf \n", IndexLevel + i, AmortLevel + i);
    fgets(string, 80, stream);
    for (i = 0; i < NUMRAMPS; i++)
        fscanf(stream, "%lf %lf \n", SeasoningLevel + i, SeasoningMat + i);
    fgets(string, 80, stream);
    for (i = 0; i < 12; i++)
        fscanf(stream, "%lf\n", Seasonality + i);
    fgets(string, 80, stream);
    if (fscanf(stream, "%d \n", &(Delay)) != 1)
        return (nrerror("Can't read prepay Delay"));
    fgets(string, 80, stream);
    if (fscanf(stream, "%d \n", &(NbFixedPrepay)) != 1)
        return (nrerror("Can't read prepay NbFixedPrepay"));
    fgets(string, 80, stream);
    FixedPrepay = mbspt_dvector(0, NbFixedPrepay - 1);
    for (i = 0; i < NbFixedPrepay; i++)
    {
        if (fscanf(stream, "%lf \n", &(FixedPrepay[i])) != 1)
        {
            err = nrerror("Can't read FixedPrepay");
            goto Free_Mem;
        }
    }

    Amorts = mbspt_dvector(0, 2 * SFSNb - 1);
    for (i = 0; i < SFSNb; ++i)
    {
        Amorts[i]         = slowAmort;
        Amorts[i + SFSNb] = fastAmort;
    }

    if (err = MBS_Prepay_Engine_Params_ReaderExplicit(
            data_dir,
            prepay_engine,
            AgencyProg,
            oterm,
            gwac,
            Speedup,
            rational,
            FebEffect,
            IndexLevel,
            AmortLevel,
            SeasoningLevel,
            SeasoningMat,
            2 * SFSNb,
            Amorts,
            Seasonality,
            Delay,
            NbFixedPrepay,
            FixedPrepay,
            nuniform,
            IndexMat,
            EndIndexMat,
            IndexFreq))
        goto Free_Mem;

Free_Mem:
{
    mbspt_free_dvector(AmortLevel, 0, NBPOINTS - 1);
    mbspt_free_dvector(Amorts, 0, 2 * SFSNb - 1);
    mbspt_free_dvector(FixedPrepay, 0, NbFixedPrepay - 1);
    return err;
}
    return err;
}  // MBS_Prepay_Engine_Params_Reader

// various adjustments from mtg rate to 10Y rate and from ratio form to diff form of refi-incentive
// change from CPR to smm
char* MBS_Prepay_Engine_SetParams(MBS_Prepay_Engine* prepay_engine, double refiRate, double Rate10Y)
{
    int    i, j, k;
    double ratio1, ratio2;
    int    half_SFSNb;

    // FIX: may has to be reconciled with agency_adj in history_adj
    // adjust_tweaks for different programs:
    if (prepay_engine->AgencyProg == _FN15)
        prepay_engine->Speedup *= 1.15;
    if (prepay_engine->AgencyProg == _GN)
        prepay_engine->Speedup *= 0.8;
    if (prepay_engine->AgencyProg == _GN15)
        prepay_engine->Speedup *= 0.7;

    // density and history adjustments
    adjust_density(
        &(prepay_engine->density),
        (prepay_engine->Speedup) * (prepay_engine->rational),
        prepay_engine->gwac,
        refiRate,
        Rate10Y,
        prepay_engine->nuniform);
    adjust_history(&(prepay_engine->history), prepay_engine->oterm, refiRate, Rate10Y);

    // convert from CPR to smm
    for (i = 0; i < prepay_engine->NbSCurve; ++i)
        for (j = 0; j < prepay_engine->NbIndexLevel; ++j)
        {
            prepay_engine->AmortLevel[i][j] =
                1.0 - pow(1.0 - prepay_engine->AmortLevel[i][j] / 100.0, 1.0 / 12.0);
        }
    // tweaks
    for (i = 0; i < prepay_engine->NbSCurve; ++i)
        spdbase(
            prepay_engine->AmortLevel[i],
            prepay_engine->Speedup,
            prepay_engine->rational,
            prepay_engine->NbIndexLevel);
    // IndexLevel
    for (j = 0; j < prepay_engine->NbIndexLevel; ++j)
        prepay_engine->IndexLevel[j] = transformIncentive(
            prepay_engine->IndexLevel[j], prepay_engine->gwac, refiRate, Rate10Y);

    // seasoning ramp
    for (j = 0; j < prepay_engine->NbSeasoningRamps; ++j)
        prepay_engine->SeasoningLevel[j] = transformIncentive(
            prepay_engine->SeasoningLevel[j], prepay_engine->gwac, refiRate, Rate10Y);
    // triggers
    ratio1     = 0.96 * prepay_engine->gwac / refiRate;
    ratio2     = 0.20 * ratio1 / 0.95;
    j          = 0;
    k          = 0;
    half_SFSNb = prepay_engine->SFSNb / 2;
    if ((0.96 * prepay_engine->gwac + REFI_THRESHOLD - refiRate) <= 0.01)
    {
        for (i = 0; i < half_SFSNb; i++)
        {
            prepay_engine->TriggerRate[i] =
                ratio1 +  // use linear spacing
                (ratio2 - ratio1) * (double)i / (double)(prepay_engine->SFSNb - i);
            prepay_engine->TriggerRate[i] *= refiRate;
        }
    }  // if
    else
    {  // otherwise cubic spacing
        ratio1 = 0.96 * prepay_engine->gwac + REFI_THRESHOLD - refiRate;
        ratio2 = pow((Rate10Y - 1.0) / ratio1, 1. / 3.) + 1.0;
        for (i = 0; i < half_SFSNb; i++)
            prepay_engine->TriggerRate[i] =
                refiRate - 1.0 + ratio1 * pow(1.0 - (double)i / (double)half_SFSNb * ratio2, 3);
    }  // else
    for (i = 0; i < half_SFSNb; ++i)
    {
        prepay_engine->TriggerRate[i] -= refiRate;
        prepay_engine->TriggerRate[i] += Rate10Y;
        prepay_engine->TriggerRate[i + half_SFSNb] = prepay_engine->TriggerRate[i];
    }  // for i
    /// FixedPrepay: convert from cpr to smm
    for (i = 0; i < prepay_engine->NbFixedPrepay; ++i)
    {
        if (prepay_engine->FixedPrepay[i] > 1.0)
            return (nrerror("Some FixedPrepay > 1.0!"));
        prepay_engine->FixedPrepay[i] = 1.0 - pow(1.0 - prepay_engine->FixedPrepay[i], 1.0 / 12.0);
    };
    return (0);
}  // MBS_Prepay_Engine_SetParams

char* MBS_Prepay_Engine_Check(MBS_Prepay_Engine* prepay_engine)
{
    int    i, j;
    double prev;
    char*  err;
    ///
    if (prepay_engine->NbSCurve < 1)
        return (nrerror("NbSCurve >= 1!"));
    if (prepay_engine->NbIndexLevel < 1)
        return (nrerror("NbSCurve >= 1!"));
    prev = 0.0;
    for (i = 0; i < prepay_engine->NbIndexLevel; ++i)
    {
        if (prepay_engine->IndexLevel[i] <= prev)
            return (nrerror("IndexLevel must be increasing!"));
        prev = prepay_engine->IndexLevel[i];
        for (j = 0; j < prepay_engine->NbSCurve; ++j)
        {
            if (prepay_engine->AmortLevel[j][i] < 0.0)
                return (nrerror("AmortLevel >= 0.0!"));
        }
    }
    for (i = 0; i < 12; ++i)
        if (prepay_engine->Seasonality[i] <= 0.0)
            return (nrerror("Seasonality > 0.0!"));
    if (prepay_engine->FebruaryEffect <= 0.0 || prepay_engine->FebruaryEffect >= 1.0)
        return (nrerror("FebrurayEffect in (0,1)!"));
    if (prepay_engine->NbSeasoningRamps < 1)
        return (nrerror("NbSeasoningRamps >= 1 !"));
    prev = 0.0;
    for (i = 0; i < prepay_engine->NbSeasoningRamps; ++i)
    {
        if (prepay_engine->SeasoningLevel[i] <= prev)
            return (nrerror("SeasoningLevel increasing!"));
        prev = prepay_engine->SeasoningLevel[i];
        if (prepay_engine->SeasoningMat[i] <= 0.0)
            return (nrerror("SeasoningMat > 0!"));
    }
    if (prepay_engine->Speedup <= 0.0 || prepay_engine->rational <= 0.0)
        return (nrerror("Speedup and rational > 0!"));
    if (prepay_engine->SFSNb < 1)
        return (nrerror("SFSNb >= 1!"));

    for (i = 0; i < prepay_engine->SFSNb; ++i)
    {
        if (prepay_engine->Amort[i] < 0.0 || prepay_engine->Amort[i] > 1.0)
            return (nrerror("Amort must be in [0,1]!"));
    }
    if (prepay_engine->Delay < 0)
        return (nrerror("Delay >= 0!"));
    if (prepay_engine->NbFixedPrepay < 0)
        return (nrerror("NbFixedPrepay >= 0"));
    if (prepay_engine->NbFixedPrepay < prepay_engine->Delay)
        return ("NbFixedPrepay < PrepayDelay");
    for (i = 0; i < prepay_engine->NbFixedPrepay; ++i)
        if (prepay_engine->FixedPrepay[i] < 0.0 || prepay_engine->FixedPrepay[i] > 1.0)
            return (nrerror("FixedPrepay must be in [0,1]!"));

    //		if( prepay_engine->agency_char != 'F' && prepay_engine->agency_char != 'G' )
    //return(nrerror( "agency_char: F or G" ));
    if (prepay_engine->oterm < 1 || prepay_engine->oterm > MAXOTERM)
        return (nrerror("oterm must be legal!"));
    if (prepay_engine->gwac < 0.0)
        return (nrerror("gwac > 0!"));
    if (prepay_engine->nuniform < 1 || prepay_engine->nuniform > MAXNBWEIGHTS)
        return (nrerror("nuniform must be legal!"));
    if (prepay_engine->IndexFreq <= 0 || prepay_engine->IndexMat <= 0 ||
        prepay_engine->EndIndexMat <= 0 || prepay_engine->IndexMat < prepay_engine->EndIndexMat)
        return ("Bad Index spec for prepay_engine");

    if (err = HISTORY_Check(&(prepay_engine->history)))
        return (err);
    if (err = MBSDENSITY_Check(&(prepay_engine->density)))
        return (err);
    return (0);
}  // MBS_Preppay_Engine_Check

// the ultimate initialization
void MBS_Prepay_Engine_Initialize(
    MBS_Prepay_Engine* prepay_engine,
    int                NbSCurve,
    int                NbIndexLevel,
    double*            IndexLevel,
    double**           AmortLevel,
    double             Seasonality[],
    double             FebruaryEffect,
    int                NbSeasoningRamps,
    double*            SeasoningLevel,
    double*            SeasoningMat,
    double             Speedup,
    double             rational,
    int                SFSNb,
    double*            TriggerRate,
    double*            Amort,
    int                Delay,
    int                NbFixedPrepay,
    double*            FixedPrepay,
    MTGPTACYPROG       AgencyProg,
    long               oterm,
    double             gwac,
    HISTORY*           history,
    MBSDENSITY*        density,
    long               nuniform,
    int                IndexMat,
    int                EndIndexMat,
    int                IndexFreq)
{
    int i, j;
    MBS_Prepay_Engine_Cleanup(prepay_engine);

    ////
    prepay_engine->NbSCurve     = NbSCurve;
    prepay_engine->NbIndexLevel = NbIndexLevel;
    if (NbIndexLevel > 0 && NbSCurve > 0)
    {
        prepay_engine->IndexLevel = mbspt_dvector(0, NbIndexLevel - 1);
        prepay_engine->AmortLevel = mbspt_dmatrix(0, NbSCurve - 1, 0, NbIndexLevel - 1);
        for (i = 0; i < prepay_engine->NbIndexLevel; ++i)
        {
            prepay_engine->IndexLevel[i] = IndexLevel[i];
            for (j = 0; j < prepay_engine->NbSCurve; ++j)
                prepay_engine->AmortLevel[j][i] = AmortLevel[j][i];
        }
    }
    /////
    for (i = 0; i < 12; ++i)
        prepay_engine->Seasonality[i] = Seasonality[i];
    prepay_engine->FebruaryEffect = FebruaryEffect;

    prepay_engine->NbSeasoningRamps = NbSeasoningRamps;
    if (NbSeasoningRamps > 0)
    {
        prepay_engine->SeasoningLevel = mbspt_dvector(0, prepay_engine->NbSeasoningRamps - 1);
        prepay_engine->SeasoningMat   = mbspt_dvector(0, prepay_engine->NbSeasoningRamps - 1);
        for (i = 0; i < prepay_engine->NbSeasoningRamps; ++i)
        {
            prepay_engine->SeasoningLevel[i] = SeasoningLevel[i];
            prepay_engine->SeasoningMat[i]   = SeasoningMat[i];
        }
    }

    prepay_engine->Speedup  = Speedup;
    prepay_engine->rational = rational;

    prepay_engine->SFSNb = SFSNb;
    if (SFSNb > 0)
    {
        prepay_engine->TriggerRate = mbspt_dvector(0, prepay_engine->SFSNb - 1);
        prepay_engine->Amort       = mbspt_dvector(0, prepay_engine->SFSNb - 1);
        for (i = 0; i < prepay_engine->SFSNb; ++i)
        {
            prepay_engine->TriggerRate[i] = TriggerRate[i];
            prepay_engine->Amort[i]       = Amort[i];
        }
    }

    prepay_engine->Delay         = Delay;
    prepay_engine->NbFixedPrepay = NbFixedPrepay;
    if (prepay_engine->NbFixedPrepay > 0)
    {
        prepay_engine->FixedPrepay = mbspt_dvector(0, prepay_engine->NbFixedPrepay - 1);
        for (i = 0; i < prepay_engine->NbFixedPrepay; ++i)
            prepay_engine->FixedPrepay[i] = FixedPrepay[i];
    }
    prepay_engine->AgencyProg = AgencyProg;
    prepay_engine->oterm      = oterm;
    prepay_engine->gwac       = gwac;

    history_copy(history, &(prepay_engine->history));
    density_copy(density, &(prepay_engine->density));
    prepay_engine->nuniform = nuniform;

    prepay_engine->IndexMat    = IndexMat;
    prepay_engine->EndIndexMat = EndIndexMat;
    prepay_engine->IndexFreq   = IndexFreq;
}  // MBS_Prepay_Engine_Initialize

// copy everything else from prepay_engine to new, except trigger info
void MBS_Prepay_Engine_ReplaceTriggers(
    MBS_Prepay_Engine* new_engine,
    MBS_Prepay_Engine* prepay_engine,
    int                num_buckets,
    double*            triggerRate,
    double*            amort)
{
    MBS_Prepay_Engine_Initialize(
        new_engine,
        prepay_engine->NbSCurve,
        prepay_engine->NbIndexLevel,
        prepay_engine->IndexLevel,
        prepay_engine->AmortLevel,
        prepay_engine->Seasonality,
        prepay_engine->FebruaryEffect,
        prepay_engine->NbSeasoningRamps,
        prepay_engine->SeasoningLevel,
        prepay_engine->SeasoningMat,
        prepay_engine->Speedup,
        prepay_engine->rational,
        num_buckets,
        triggerRate,
        amort,
        prepay_engine->Delay,
        prepay_engine->NbFixedPrepay,
        prepay_engine->FixedPrepay,
        prepay_engine->AgencyProg,
        prepay_engine->oterm,
        prepay_engine->gwac,
        &(prepay_engine->history),
        &(prepay_engine->density),
        prepay_engine->nuniform,
        prepay_engine->IndexMat,
        prepay_engine->EndIndexMat,
        prepay_engine->IndexFreq);
}  // MBS_Prepay_Engine_ReplaceTriggers

void MBS_Prepay_Engine_Cleanup(MBS_Prepay_Engine* prepay_engine)
{
    if (prepay_engine->NbIndexLevel > 0)
        mbspt_free_dvector(prepay_engine->IndexLevel, 0, prepay_engine->NbIndexLevel - 1);
    if (prepay_engine->NbSCurve > 0 && prepay_engine->NbIndexLevel > 0)
        mbspt_free_dmatrix(
            prepay_engine->AmortLevel,
            0,
            prepay_engine->NbSCurve - 1,
            0,
            prepay_engine->NbIndexLevel - 1);
    if (prepay_engine->NbSeasoningRamps > 0)
        mbspt_free_dvector(prepay_engine->SeasoningLevel, 0, prepay_engine->NbSeasoningRamps - 1);
    if (prepay_engine->NbSeasoningRamps > 0)
        mbspt_free_dvector(prepay_engine->SeasoningMat, 0, prepay_engine->NbSeasoningRamps - 1);
    if (prepay_engine->SFSNb > 0)
    {
        mbspt_free_dvector(prepay_engine->TriggerRate, 0, prepay_engine->SFSNb - 1);
        mbspt_free_dvector(prepay_engine->Amort, 0, prepay_engine->SFSNb - 1);
    }
    if (prepay_engine->NbFixedPrepay > 0)
        mbspt_free_dvector(prepay_engine->FixedPrepay, 0, prepay_engine->NbFixedPrepay - 1);
}
/// BS_Prepay_Engine_Cleanup

//////////////////

// init_age is wala, for prepay_delay==0, should be 0
// init_month = month of factor of first paydown
// return history of factors
char* MBS_Prepay_FactorHistory(
    MBS_Prepay_Engine* prepay_engine,
    int                init_age,
    int                init_month,
    int                num_rates,
    double*            input_rates,
    int                num_buckets,
    double*            triggerRate,
    double*            amort,
    int                useInitialSMM,
    double**           factors)
// num_buckets <=0 means uses prepay_engine's SFSNb, TriggerRate, and Amort; assumes there is no
// delay
{
    int               i, j;
    int               newSFSNb       = prepay_engine->SFSNb;
    double*           newTriggerRate = prepay_engine->TriggerRate;
    double*           newAmort       = prepay_engine->Amort;
    int               age, monthOfYear, monthfromsettle;
    int               NbFixedPrepay = prepay_engine->NbFixedPrepay;
    double*           FixedAmort;
    MBS_Prepay_Engine new_engine = {0};
    //
    if (num_buckets > 0)
    {
        newSFSNb       = num_buckets;
        newTriggerRate = triggerRate;
        newAmort       = amort;
    }
    MBS_Prepay_Engine_ReplaceTriggers(
        &new_engine, prepay_engine, newSFSNb, newTriggerRate, newAmort);
    FixedAmort = mbspt_dvector(0, NbFixedPrepay - 1);
    for (i = 0; i < NbFixedPrepay; ++i)
    {
        FixedAmort[i] = 1.0 - prepay_engine->FixedPrepay[i];
    }
    // age matters
    age             = init_age;
    monthOfYear     = init_month;
    monthfromsettle = 0;
    // real work, special handling for initial month
    // MBS_Prepay_Engine_SMM(factors[0], &new_engine, age, input_rates[0], monthOfYear,
    // monthfromsettle);
    for (j = 0; j < new_engine.SFSNb; ++j)
        factors[0][j] = 1.0;
    for (i = 1; i <= num_rates; ++i)
    {
        //#ifdef MBSPTCORRECTPPHISTORY
        MBS_Prepay_Engine_SMM(
            factors[i], &new_engine, age + 1, input_rates[i - 1], monthOfYear, monthfromsettle);
        //#else

        //		MBS_Prepay_Engine_SMM(factors[i], &new_engine, age, input_rates[i-1],
        //monthOfYear, monthfromsettle);
        //#endif//MBSPTCORRECTPPHISTORY

        if (useInitialSMM)
            for (j = 0; j < new_engine.SFSNb; ++j)
            {
                if (i <= NbFixedPrepay)
                    factors[i][j] = factors[i - 1][j] * FixedAmort[i - 1];
                else
                    factors[i][j] = factors[i - 1][j] * (1.0 - factors[i][j]);
            }
        else
            for (j = 0; j < new_engine.SFSNb; ++j)
                factors[i][j] = factors[i - 1][j] * (1.0 - factors[i][j]);
        age++;
        monthfromsettle++;
        monthOfYear++;
        if (monthOfYear > 12)
            monthOfYear -= 12;
    }
    // cleanup
    MBS_Prepay_Engine_Cleanup(&new_engine);
    mbspt_free_dvector(FixedAmort, 0, NbFixedPrepay - 1);
    //
    return 0;
}  // MBS_Prepay_FactorHistory

// return a path of prepays in CPR from a path of rates
// rates[i] < 0.0 if useFixedPrepay > 0 and i < NbFixedPrepay
char* MBS_Prepay_Path(
    MBS_Prepay_Engine* prepay_engine,
    long               TradeYYYYMMDD,
    long               spotYYYYMMDD,
    long               wam,
    int                subgroup,
    int                num_future_rates,
    double*            futureRates,
    int                useFixedPrepay,
    int*               numProjections,
    double             rates[MAXOTERM],
    double             cpr[MAXOTERM],
    double             principal_payments[MAXOTERM])
// subgroup == 0 <=> aggregated, slow group comes first starting with 1 then fast group
// principal_payments assume notional of $1 as of spotYYYYMMDD
{
    long     wala, spot_position, o_position;
    long     dd, spotYYYY, firstHistYYYYMM;
    double   sprd;
    double** factors;
    double*  weights;
    double*  fsmm;
    double   tmp;
    int      spot_month, mo;
    int      i, j, numHRates;
    int      Delay;
    char*    err;

    if (subgroup > prepay_engine->SFSNb || subgroup < 0)
        return (nrerror("Illegal subgroup number"));

    Delay           = prepay_engine->Delay;
    *numProjections = min(wam, num_future_rates);
    // position accounting
    wala = prepay_engine->oterm - wam;

    firstHistYYYYMM = prepay_engine->history.dates[0];
    // spotYYYYMM = spotYYYYMMDD / 100;
    // spotYYYY = spotYYYYMM / 100;
    // spot_month = spotYYYYMM - 100 * spotYYYY;
    Dsplit(spotYYYYMMDD, &spotYYYY, &spot_month, &dd);

    spot_position = month_diff(YYYYMM2FirstOfMon(firstHistYYYYMM), spotYYYYMMDD);
    o_position    = spot_position - wala;

    numHRates = Delay;
    if ((spot_position > prepay_engine->history.Number) || (spot_position - Delay < 0))
        return (nrerror("Incomplete history"));
    // fill in future rates
    // recover mtg sprd and transform into 10YRate space
    sprd =
        (prepay_engine->density.IncentKnots[0] * prepay_engine->gwac -
         prepay_engine->density.irates[0]);
    i = 0;
    if (useFixedPrepay > 0)
    {
        for (i = 0; i < *numProjections; ++i)
        {
            if (i < prepay_engine->NbFixedPrepay)
                rates[i] = 0.0;
            else
            {
                rates[i] = futureRates[i - Delay];
                rates[i] = transformIncentive(rates[i], 1.0, sprd, 0.0);
                if (rates[i] < 0.0)
                    return (nrerror("Some future rates < 0 after transformation"));
            }
        }
    }
    else
    {
        for (i = 0; i < *numProjections; ++i)
        {
            if (i < numHRates)
                rates[i] = prepay_engine->history.hrates[spot_position - Delay + i];
            else
            {
                rates[i] = futureRates[i - Delay];
                rates[i] = transformIncentive(rates[i], 1.0, sprd, 0.0);
                if (rates[i] < 0.0)
                    return (nrerror("Some future rates < 0 after transformation"));
            }
        }
    }
    // get factors
    if (*numProjections <= 0)
        return (0);
    // future factors:
    factors = mbspt_dmatrix(0, *numProjections, 0, prepay_engine->SFSNb - 1);
    mo      = spot_month - Delay;
    while (mo < 1)
        mo += 12;
    MBS_Prepay_FactorHistory(
        prepay_engine, wala - Delay, mo, *numProjections, rates, -1, 0, 0, useFixedPrepay, factors);
    // get distribution as of spot
    weights = mbspt_dvector(0, prepay_engine->SFSNb - 1);
    if (!subgroup)
    {
        if (err = get_Trigger_Distribution(prepay_engine, TradeYYYYMMDD, wam, weights))
            return (err);
    }
    else
    {
        for (i = 0; i < prepay_engine->SFSNb; ++i)
            weights[i] = ((i == subgroup - 1) ? 1.0 : 0.0);
    }
    // now aggregrate
    fsmm = mbspt_dvector(0, *numProjections - 1);
    for (i = 1; i <= *numProjections; ++i)
    {
        tmp = 0.0;
        for (j = 0; j < prepay_engine->SFSNb; ++j)
            tmp += weights[j] * factors[i][j];
        fsmm[i - 1] = tmp;
    }
    for (i = 1; i < *numProjections; ++i)
        cpr[i] = 1.0 - fsmm[i] / fsmm[i - 1];
    cpr[0] = 1.0 - fsmm[0];
    if (*numProjections == wam)
        cpr[wam - 1] = 0.0;
    // override the first few smm's
    // for(i=0; i < prepay_engine->NbFixedPrepay; ++i) cpr[i] = prepay_engine->FixedPrepay[i];
    // before turning smm to cpr, calc factor history:
    if (err = get_initial_principal_payments(
            *numProjections, principal_payments, wam, prepay_engine->gwac, *numProjections, cpr))
        return (err);
    //
    for (i = 0; i < *numProjections; ++i)
        cpr[i] = 1.0 - pow((1.0 - cpr[i]), 12.0);
    if (useFixedPrepay > 0)
        for (i = 0; i < prepay_engine->NbFixedPrepay; ++i)
            rates[i] = -100.0;
    // cleanup
    mbspt_free_dvector(fsmm, 0, *numProjections - 1);
    mbspt_free_dmatrix(factors, 0, *numProjections - 1, 0, prepay_engine->SFSNb - 1);
    mbspt_free_dvector(weights, 0, prepay_engine->SFSNb - 1);

    return 0;
}  // MBS_Prepay_Path

// given dense representation of initial_trigger_distribution
// get sparse representation
void get_smm_weights(
    double* weights, MBS_Prepay_Engine* prepay_engine, double* usweights, double* ufweights)
{
    int     nuniform = prepay_engine->density.Nbweights;
    double* vec;
    int     i;
    int     subgroupSFSNb = prepay_engine->SFSNb / 2;
    double  swt, fwt;
    double* reversedTrigger;
    reversedTrigger = mbspt_dvector(0, subgroupSFSNb - 1);
    for (i = 0; i < subgroupSFSNb; ++i)
        reversedTrigger[subgroupSFSNb - 1 - i] = prepay_engine->TriggerRate[i];
    vec = mbspt_dvector(0, subgroupSFSNb - 1);
    for (i = 0; i < subgroupSFSNb; ++i)
    {
        // set kronecker delta vector
        if (i)
            vec[i - 1] = 0.0;
        vec[i] = 1.0;
        swt    = qinterpp(
            subgroupSFSNb,
            reversedTrigger,
            vec,
            nuniform,
            prepay_engine->density.urates,
            usweights);
        fwt = qinterpp(
            subgroupSFSNb,
            reversedTrigger,
            vec,
            nuniform,
            prepay_engine->density.urates,
            ufweights);
        weights[subgroupSFSNb - 1 - i]                 = 0.5 * swt;
        weights[subgroupSFSNb - 1 - i + subgroupSFSNb] = 0.5 * fwt;
    }
    //
    mbspt_free_dvector(vec, 0, subgroupSFSNb - 1);
}  // get_smm_weights

// ultimate call of prepay projection
// smm is ana arry carrying the SMM of various trigger groups
// output is smm, with 1% <->0.01
// indexvalue: 1% <-> 1
// first payment after settle <-> paymentNb=1
// monthOfYear = 2 for February prepay, though the reset date for getting IndexValue may be March 1
int MBS_Prepay_Engine_SMM(
    double*            output,
    MBS_Prepay_Engine* prepay_engine,
    int                paymentNb,
    double             IndexValue,
    int                monthOfYear,
    int monthFromSettle)  // this assumed there is NO prepay-delay i.e., if there is delay, the
                          // input paymentNb, IndexValue, monthOfYear and monthFromSetlle, etc., have
                          // to adjusted properly
// paymentNb < 0 => base smm = 0.0, only asmm
{
    int    SCurveIndex;
    double smm;
    double seasoningFactor;
    double asmm;
    int    i;
    int    mo;
    // pick appropriate S-curve:
    // FIX: check NbSCurve > 0
    if (prepay_engine->NbSCurve == 1)
        SCurveIndex = prepay_engine->NbSCurve - 1;
    else
    {
        if (monthFromSettle > prepay_engine->NbSCurve + 6)
            SCurveIndex = prepay_engine->NbSCurve - 1;
        else if (
            monthFromSettle <= prepay_engine->NbSCurve + 6 &&
            monthFromSettle >= prepay_engine->NbSCurve - 2)
            SCurveIndex = prepay_engine->NbSCurve - 2;
        else if (monthFromSettle >= 0)
            SCurveIndex = monthFromSettle;
        else
            SCurveIndex = 0;  // FIX: should not have to come to this case
    }
    smm = linterpFlat(
        prepay_engine->NbIndexLevel,
        prepay_engine->IndexLevel,
        prepay_engine->AmortLevel[SCurveIndex],
        IndexValue);

    // Seasonality
    mo = monthOfYear;
    smm *= prepay_engine->Seasonality[mo - 1];
    // Seasoning:
    seasoningFactor = ((double)paymentNb) / linterpFlat(
                                                prepay_engine->NbSeasoningRamps,
                                                prepay_engine->SeasoningLevel,
                                                prepay_engine->SeasoningMat,
                                                IndexValue);
    smm *= min(1.0, seasoningFactor);
    smm = max(0.0, smm);

    if (monthOfYear == 2)
        smm *= prepay_engine->FebruaryEffect;

    // trigger
    for (i = 0; i < prepay_engine->SFSNb; ++i)
    {
        asmm = 0.0;
        if (IndexValue < prepay_engine->TriggerRate[i])
            asmm = prepay_engine->Amort[i];
        output[i] = 1.0 - (1.0 - smm) * (1.0 - asmm);

        //#ifndef CORRECTFEBEFFECT
        //		if (monthOfYear == 2) output[i] *= prepay_engine->FebruaryEffect;
        //#endif//#ifndef CORRECTFEBEFFECT
    }
    return 0;
}  // MBS_Prepay_Engine_SMM

// seedup and steepen S-curves
void spdbase(double* AmortLevel, double Speedup, double rational, int numpts)
{
    int i;

    double S0;

    S0 = AmortLevel[4]; /* speed up and steepen base prepayment */
    for (i = 0; i < numpts; i++)
    {
        if (AmortLevel[i] > S0)
            AmortLevel[i] += (rational - 1.0) * (AmortLevel[i] - S0);
        else
            AmortLevel[i] *= pow((AmortLevel[i] / S0), (rational - 1.0));
        AmortLevel[i] *= Speedup;
    }

}  //*spdbase*/

/*Convert average life to smm*/
double alsmm1(double al, int N) /* al is average life in years, N is original term in months */
{
    double x1, x2, x3, y3;
    int    i;

    x1 = 1. / (al); /* I know f(x1,al,N) >0. */
    x2 = small;     /* ALSMM1FUNC(x2, al, N) <0. if N>al	*/
    x3 = (x1 + x2) / 2.0;
    for (i = 0; i < 30; i++)
    {
        x3 = (x1 + x2) / 2.;
        y3 = ALSMM1FUNC(x3, al - 1., N);
        if (y3 > 0.)
            x1 = x3;
        else
            x2 = x3;
    }

    return x3;
}  //*alsmm1*/

/*convert average life to smm*/
double alsmm2(double al, double wac, int N, int wam, int pt) /* al -avg life, wac -- gross coupon, N
                                                                -- original term wam -- remaining
                                                                term, pt -- N-wam*/
{
    double x1, x2, x3, y3, factor1, factor2, factor3;
    int    i;

    factor1 = pow(1. + wac / 1200., N);
    factor3 = pow(1. + wac / 1200., pt);
    factor2 = factor1 - 1.;
    x1      = 1.;    /* intial point, which gives too short an average life */
    x2      = small; /* which gives too long an average life */
    x3      = (x1 + x2) / 2.0;
    for (i = 0; i < 30; i++)
    {
        x3 = (x1 + x2) / 2.;
        y3 = al - 1. - (factor1 / factor2) * ALSMM2FUNC(x3, wam) -
             ALSMM2FUNC(x3 * (1. + wac / 1200.) - wac / 1200., wam) * factor3 / factor1;
        if (y3 > 0.)
            x1 = x3;
        else
            x2 = x3;
    }

    return x3;
}  //*alsmm2*/

// given historical rates
// get dense representation of the density function of various trigger groups at the start (settle
// date or trade date?)

char* get_Trigger_Distribution(
    MBS_Prepay_Engine* prepay_engine, long spotYYYYMMDD, long wam, double* weights)
{
    int    i;
    double sum, sum1;

    double*  slowAmort;
    double*  fastAmort;
    double** fastFactors;
    double** slowFactors;
    double*  input_rates;
    double*  usweights;
    double*  ufweights;

    double slowASMM;
    double fastASMM;
    int    num_rates;
    ///
    double* urates;
    int     spot_position;
    int     age;
    int     num_buckets;
    int     init_month;
    int     Delay;
    long    yy, dd;

    Delay         = prepay_engine->Delay;
    urates        = prepay_engine->density.urates;
    spot_position = month_diff(YYYYMM2FirstOfMon(prepay_engine->history.dates[0]), spotYYYYMMDD);
    age           = prepay_engine->oterm - wam;
    num_buckets   = prepay_engine->density.Nbweights;
    slowASMM      = prepay_engine->Amort[0];
    fastASMM      = prepay_engine->Amort[prepay_engine->SFSNb - 1];

    num_rates = age;

    if ((spot_position - age - Delay < 0) ||
        (spot_position - age - Delay + num_rates - 1 >= prepay_engine->history.Number))
        return (nrerror("Incomplete History!"));

    slowAmort   = mbspt_dvector(0, num_buckets - 1);
    fastAmort   = mbspt_dvector(0, num_buckets - 1);
    fastFactors = mbspt_dmatrix(0, num_rates, 0, num_buckets - 1);
    slowFactors = mbspt_dmatrix(0, num_rates, 0, num_buckets - 1);
    input_rates = mbspt_dvector(0, num_rates - 1);
    usweights   = mbspt_dvector(0, num_buckets - 1);
    ufweights   = mbspt_dvector(0, num_buckets - 1);

    for (i = 0; i < num_buckets; ++i)
    {
        fastAmort[i] = fastASMM;
        slowAmort[i] = slowASMM;
        usweights[i] = prepay_engine->density.usweights[i];
        ufweights[i] = prepay_engine->density.ufweights[i];
    }

    for (i = 0; i < num_rates; ++i)
    {
        input_rates[i] = prepay_engine->history.hrates[spot_position - age + i - Delay];
    }

    Dsplit(
        Date_Add(
            YYYYMM2FirstOfMon(prepay_engine->history.dates[spot_position - age - Delay]),
            MONTH,
            -Delay,
            0,
            1),
        &yy,
        &init_month,
        &dd);
    MBS_Prepay_FactorHistory(
        prepay_engine,
        -Delay,
        init_month,
        num_rates,
        input_rates,
        num_buckets,
        urates,
        fastAmort,
        0,
        fastFactors);
    MBS_Prepay_FactorHistory(
        prepay_engine,
        -Delay,
        init_month,
        num_rates,
        input_rates,
        num_buckets,
        urates,
        slowAmort,
        0,
        slowFactors);

    sum  = 0.0;
    sum1 = 0.0;
    for (i = 0; i < num_buckets; i++)
    {
        ufweights[i] *= fastFactors[num_rates][i];
        usweights[i] *= slowFactors[num_rates][i];

        sum += ufweights[i];
        sum1 += usweights[i];
    }
    for (i = 0; i < num_buckets; i++)
    {
        ufweights[i] /= sum;
        usweights[i] /= sum1;
    }
    //
    get_smm_weights(weights, prepay_engine, usweights, ufweights);
    //
    mbspt_free_dvector(slowAmort, 0, num_buckets - 1);
    mbspt_free_dvector(fastAmort, 0, num_buckets - 1);
    mbspt_free_dmatrix(fastFactors, 0, num_rates, 0, num_buckets - 1);
    mbspt_free_dmatrix(slowFactors, 0, num_rates, 0, num_buckets - 1);
    mbspt_free_dvector(input_rates, 0, num_rates - 1);
    mbspt_free_dvector(usweights, 0, num_buckets - 1);
    mbspt_free_dvector(ufweights, 0, num_buckets - 1);

    return (0);
}  // getTriggerDistribution

// per $ notional as of the start
// gives actual dollar amount of principal payment
// wam = oterm - payment_nb + 1
// where payment_nb = 1 <-> homeowner makes 1st payment
// and in that case out[0] is 1.0 - first_principal_payment
char* get_initial_principal_payments(
    int num_payments, double* out, int wam, double gwac, int NbFixedPrepay, double* FixedPrepay)
{
    int i;
    //	double u;
    int moving_wam;

    if (num_payments > NbFixedPrepay)
        return (nrerror("Ask for too many principal paydowns"));
    if (num_payments <= 0)
        return (0);  // do nothing
    // u = 1.0/(1.0 + gwac /100.0 /12.0 );
    moving_wam = wam;
    out[0]     = 1.0 - TOTAL_AMORT_RATE(moving_wam, gwac, FixedPrepay[0] * 100.0);
    // out[0] = (1.0 - pow(u, moving_wam-1))/(1.0 - pow(u,moving_wam)) * ( 1.0 - FixedPrepay[0] );

    for (i = 1; (i < num_payments) && (moving_wam > 1); ++i)
    {
        moving_wam--;
        out[i] = out[i - 1] * (1.0 - TOTAL_AMORT_RATE(moving_wam, gwac, FixedPrepay[i] * 100.0));
        // out[i] = out[i-1] * (1.0 - pow(u,moving_wam-1))/(1.0 - pow( u, moving_wam)) * (1.0 -
        // FixedPrepay[i]);
    }
    if (i < num_payments)
    {
        while (i < num_payments)
            out[i++] = 0.0;
    }

    for (i = num_payments - 1; i > 0; --i)
        out[i] = out[i - 1] - out[i];
    out[0] = 1.0 - out[0];

    return 0;
}

/***
void test_MBS_Prepay_Path( MBS_Prepay_Engine * prepay_engine, long spotYYYYMMDD, long  wam, int
num_future_rates)
{
        FILE * fp;
        long firstHistYYYYMM;
        double last_rate, * future_rates;
        int i;
        double rates[MAXOTERM], cpr[MAXOTERM], principalPayDowns[MAXOTERM];
        int numProjections, spot_position;

        firstHistYYYYMM = prepay_engine->history.dates[0];
        spot_position = month_diff( firstHistYYYYMM*100+1, spotYYYYMMDD);
        last_rate = prepay_engine->history.hrates[spot_position] +  (
prepay_engine->density.IncentKnots[0] * prepay_engine->gwac - prepay_engine->density.irates[0] );
        future_rates = mbspt_dvector(0, num_future_rates-1);
        for(i=0; i < num_future_rates; ++i) future_rates[i] = last_rate;

        MBS_Prepay_Path( prepay_engine, spotYYYYMMDD, wam, 0, num_future_rates, future_rates,
                                                  &numProjections, rates, cpr, principalPayDowns);
        fp = fopen("C:/workarea/mortgage/jymodel1/fixedpp.txt","w");
        for(i=0; i < numProjections; ++i)
        {
                fprintf(fp, "rates[%d] = %f, cpr = %f, paydown=%f\n", i, rates[i], cpr[i]*100.0,
principalPayDowns[i]*100.0);
        }

        mbspt_free_dvector( future_rates, 0, num_future_rates-1);
        fclose(fp);
}
***/
