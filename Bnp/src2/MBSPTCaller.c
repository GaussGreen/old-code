#include "MBSPTCaller.h"

#include "MBSPTProdStruct.h"
#include "MBSPrepay.h"
#include "utdates.h"

// tree is both input and output
// tree is assumed to have basic information, such as RefDate and PPY ready
// alloc memory for tree_times and assume times, always tree_times[0] = 0.0 <-> RefDate
// add last tree_time = penutilmate_tree_time + penutimate dt
// return the index idx such that tree->tree_times[idx] <-> last DuePeriodEndDates
// Assume tree times coincide with mtg reset dates, with RefDate (whether settle or spot date does
// not matter as long as they differ by a small amount) and some dates at the end; this means that
// PPY == 11 Assume security already seasoned (i.e.. wala >=0) on settle date
int MBS_BKTree_create_times(
    MBS_BKTree* tree, MBSPT_DealStruct* deal_struct, MBS_Prepay_Engine* prepay_engine)
{
    int numTimeSlices, numDuePrdTimes;
    // int refIsDuePrdEnd;
    int                   i, j, k;
    double*               times;
    long                  date, ref_date, last_date;
    int                   interval;
    int*                  numSubIntervals;
    int                   totalNumSubIntervals;
    long                  Y = 0, MM = 0, DD = 0;
    double                dt;
    MBSPT_CashFlowStruct* cashflow_struct;
    //
    cashflow_struct = &(deal_struct->cashflow_struct);
    ref_date        = tree->RefDate;
    // number of cashflow time slices + settle_date
    // if( ref_date == cashflow_struct->duePeriodEndDates[0] ) refIsDuePrdEnd = 1;
    // else refIsDuePrdEnd = 0;
    // refIsDuePrdEnd=0;//FIX: to be removed!
    // actually, given the way we structure dates, refisDuePrdEnd is always 0!!

    numDuePrdTimes = cashflow_struct->numCashFlows + 1;  // - refIsDuePrdEnd;//add 0.0 at beginning
    // add those needed for for index calculation
    numTimeSlices =
        numDuePrdTimes + prepay_engine->EndIndexMat * 12;  // strictly, need 1 less than this
    // extended cashflow times:
    times    = mbspt_dvector(0, numTimeSlices - 1);
    times[0] = 0.0;
    for (i = 1; i < numDuePrdTimes; ++i)
        times[i] = cashflow_struct->duePeriodEndTimes[i - 1];  //+ refIsDuePrdEnd];

    if (cashflow_struct->numCashFlows > 0)
    {
        last_date = cashflow_struct->duePeriodEndDates[cashflow_struct->numCashFlows - 1];
    }
    else
    {
        Dsplit(ref_date, &Y, &MM, &DD);
        last_date = DateYMMDD(Y, MM, 1);
    }

    interval = 12 / cashflow_struct->cpnFreq;
    for (i = numDuePrdTimes, j = 1; i < numTimeSlices; ++i, ++j)
    {
        date     = Date_Add(last_date, MONTH, j * interval, 0, 1);
        times[i] = DateDiff(ref_date, date, _DACT) / 365.0;
    }

    // calc num of subintervals;
    numSubIntervals      = mbspt_ivector(0, numTimeSlices - 2);
    totalNumSubIntervals = 0;
    for (i = 0; i < numTimeSlices - 1; ++i)
    {
        numSubIntervals[i] = round(((double)tree->PPY) * (times[i + 1] - times[i]));
        if (numSubIntervals[i] == 0)
            numSubIntervals[i] = 1;
        totalNumSubIntervals += numSubIntervals[i];
    }

    // now fill in tree times
    tree->NbTimeSlices  = totalNumSubIntervals + 1 + 1;  // include 0.0 AND one more at end
    tree->NbSpaceNodes  = 2 * tree->NbTimeSlices - 1;
    tree->EndTimeSlices = tree->NbTimeSlices - 1;
    tree->tree_times    = mbspt_dvector(0, tree->NbTimeSlices - 1);
    tree->tree_times[0] = 0.0;
    dt                  = 0.0;
    for (i = 0, j = 1; i < numTimeSlices - 1; ++i)
    {
        dt = (times[i + 1] - times[i]) / ((double)numSubIntervals[i]);
        for (k = 0; k < numSubIntervals[i]; ++k)
        {
            tree->tree_times[j] = tree->tree_times[j - 1] + dt;
            j++;
        }
    }
    tree->tree_times[j] = tree->tree_times[j - 1] + dt;  // tag one more at end
    // clean up
    mbspt_free_dvector(times, 0, numTimeSlices - 1);
    mbspt_free_ivector(numSubIntervals, 0, numTimeSlices - 2);
    return (numDuePrdTimes - 1);
}  // MBS_BKTree_create_times

// get smm's for all speed groups at all space nodes for a particular time slice
// smm[0..prepay_engine->SFSNb-1][0, tree->NbSpaceNodes-1], mem to be alloc by user
// 1% <-> 0.01
char* MBS_PricesGetSMM(
    double**              smm,
    int                   NbPer,
    MBSPT_CashFlowStruct* cashflowstruct,
    MBS_Prepay_Engine*    prepay_engine,
    MBS_BKTree*           tree)
{
    double* indices;
    double* smm_holder;
    char*   err = 0;
    int     i, j;
    int     o_wala;
    int     start, end;
    long    yy, mm, dd;
    //
    indices    = mbspt_dvector(0, tree->NbSpaceNodes - 1);
    smm_holder = mbspt_dvector(0, prepay_engine->SFSNb - 1);
    //
    start  = tree->Bottom[NbPer];
    end    = tree->Top[NbPer];
    o_wala = prepay_engine->oterm - cashflowstruct->numCashFlows;
    // see comments in MBS_PricesBackwardInduct:
    if ((NbPer > prepay_engine->NbFixedPrepay - prepay_engine->Delay) &&
        (NbPer < cashflowstruct->numCashFlows - prepay_engine->Delay))
    {
        if (err = MBS_BKTree_Yield_Approx(
                indices, prepay_engine->IndexFreq, prepay_engine->IndexMat, NbPer, tree))
        {
            goto CLEANUP;
        }
        Dsplit(cashflowstruct->duePeriodEndDates[NbPer - 1], &yy, &mm, &dd);  // only mm is needed
        mm--;  // rate on reset date 20030801, e.g., signifies prepay of 200307!
        if (mm <= 0)
            mm += 12;
        for (i = start; i <= end; ++i)
        {
            MBS_Prepay_Engine_SMM(smm_holder, prepay_engine, o_wala + NbPer, indices[i], mm, NbPer);
            for (j = 0; j < prepay_engine->SFSNb; ++j)
                smm[j][i] = smm_holder[j];
        }
    }
    else
    {
        for (i = start; i <= end; ++i)
        {
            for (j = 0; j < prepay_engine->SFSNb; ++j)
            {
                smm[j][i] = 0.0;
            }
        }
    }
    // cleanup
CLEANUP:
    mbspt_free_dvector(indices, 0, tree->NbSpaceNodes - 1);
    mbspt_free_dvector(smm_holder, 0, prepay_engine->SFSNb - 1);
    return (err);
}  // MBS_PricesGetSMM

///////////////

// a basic function, not to be used directly
// one single speed group at a time
// IO PO prices are of notional $100
// coupon: 5% <-> 5
// smm (1% <-> 0.01) is pre-calculated in MBS_PricesGetSMM
char* MBS_PricesBackwardInduct(
    double*               io,
    double*               po,
    int                   NbPer,
    double                coupon,  // may be replaced later
    double*               smm,
    MBSPT_CashFlowStruct* cashflowstruct,
    MBS_Prepay_Engine*    prepay_engine,
    MBS_BKTree*           tree)
{
    char*    err = 0;
    int      start, end;
    double   sch_amort_rate;
    double **d_dptr, **zeros;
    int      i, j, pos, NbZeros, NbZerosCached;
    double   t, t0, t1;
    int      o_wala;  // i.e. number of homeowner payments made as of settle
    int      NbZerosToSkip;
    // backward induct io and po prices:
    if (err = MBS_BKTree_BackwardInduct(io, 1, NbPer, tree))
        return (err);
    if (err = MBS_BKTree_BackwardInduct(po, 1, NbPer, tree))
        return (err);
    if (NbPer == 0)
        return (err);
    //
    o_wala = prepay_engine->oterm - cashflowstruct->numCashFlows;
    // calc NPV of additional cashflows
    // get basic info
    // sch_amort_rate: special handling at the beginning
    // behavior is always such that smm denotes (possible) cash adjustment due to PrepayDelay
    // sch_amort_rate = foreknown amortization rate, without any delay
    // there is no smm for NbPer in [1..NbFixedPrepay-PrepayDelay] and [last-PrepayDelay..last]
    // sch_amort_rate has to incorporate in the FixedPrepay for [1..NbFixedPrepay]

    if (NbPer <= prepay_engine->NbFixedPrepay)
        sch_amort_rate = 0.0;
    else
        sch_amort_rate = SCHEDULED_AMORT_RATE(
            (prepay_engine->oterm) - (NbPer + o_wala) + 1, (prepay_engine->gwac));

    // get zeros
    if (err = MBS_BKTree_ZeroPrices(&NbZerosCached, &d_dptr, 1, NbPer, tree))
    {
        goto CLEANUP;
    }
    if (NbPer <= cashflowstruct->NbInitialCashFlowToSkip - prepay_engine->Delay)
        return (err);  // has to return only after zeros have been calc and cached!
    NbZerosCached++;
    // arrange for the correct zeros with pp delay and pay delay
    if (NbZerosCached <= 1)
        return ("Not enough zeros for inter/extra-polate");
    NbZeros = prepay_engine->Delay + 1;
    if (NbZeros > cashflowstruct->numCashFlows - NbPer + 1)
        NbZeros = cashflowstruct->numCashFlows - NbPer + 1;
    zeros = mbspt_dmatrix(0, NbZeros - 1, 0, tree->NbSpaceNodes - 1);
    start = tree->Bottom[NbPer];
    end   = tree->Top[NbPer];
    for (i = 0; i < NbZeros; ++i)
    {
        t   = cashflowstruct->CashFlowTimes[NbPer - 1 + i];
        pos = interp_search(
            NbZerosCached, tree->tree_times + NbPer, t);  // position relative to NbPer
        t0 = tree->tree_times[NbPer + pos - 1];
        t1 = tree->tree_times[NbPer + pos];
        for (j = start; j <= end; ++j)
        {
            zeros[i][j] = LogInterp(t, t0, t1, d_dptr[pos - 1][j], d_dptr[pos][j]);
        }
    }
    // add cashflow pv's
    NbZerosToSkip = 0;
    if (cashflowstruct->NbInitialCashFlowToSkip > 0 && prepay_engine->Delay > 0)
        NbZerosToSkip = cashflowstruct->NbInitialCashFlowToSkip - NbPer + 1;

    if (err = MBSPT_AddCashFlowPV(
            io, po, start, end, coupon / 100.0, sch_amort_rate, smm, NbZeros, zeros, NbZerosToSkip))
    {
        goto CLEANUPWITHZEROS;
    }
CLEANUPWITHZEROS:
    mbspt_free_dmatrix(zeros, 0, NbZeros - 1, 0, tree->NbSpaceNodes - 1);
    goto CLEANUP;
    ////
CLEANUP:
    return (err);
}  // MBS_PricesBackwardInduct

// filled with zeros
// then start with backward induct at the last slice, i.e., at NbPer == cashflowstruct->numCashFlows
// user alloc mem
char* mbs_iopo_initialize_prices(double* io, double* po, int start, int end)
{
    set_vec_to_val(io, 0.0, start, end);
    set_vec_to_val(po, 0.0, start, end);
    return (0);
}  // mbs_iopo_initialize_prices

// initialize zeros and we quit once tree->cach->currIdx_zeros == UntilNbPer
// UntilNbPer should be cashflowstruct->numCashFlows
char* mbs_iopo_initialize_zeros(MBS_BKTree* tree, int UntilNbPer)
{
    int      lastFwdIndex;
    double** zeros;
    int      i;
    char*    err = 0;

    for (i = tree->EndTimeSlices - 1; i >= UntilNbPer; --i)
    {
        if (err = MBS_BKTree_ZeroPrices(&lastFwdIndex, &zeros, 0, i, tree))
            return (err);
    }
    return (err);
}  // mbs_iopo_initialize_zeros

// io and po are pointers to double
// containing the unadjusted prices of io and po
// which means:
// set scheduled and unscheduled principal payments = 0 for up to NbFixedPrepay - Prepay-Delay
// and the coresponding io and po prices
// the following make the necessary adjustments due to user provided SMM
// the result are the dirty prices

char* MBS_IOPO_Price_Adjust(
    double*            io,
    double*            po,
    MBSPT_DealStruct*  deal_struct,
    MBS_Prepay_Engine* prepay_engine,
    MBS_BKTree*        tree)

{
    double                io_adjust = 0.0, po_adjust = 0.0;
    double*               bals = 0;
    double *              dfs  = 0, **d_dptr;
    int                   NbFixedPrepay;
    int                   i   = 0, pos;
    char*                 err = 0;
    int                   o_wala, wam;
    MBSPT_CashFlowStruct* cashflowstruct;
    int                   NbZerosCached;
    double                t, t0, t1;
    int                   NbCashFlow;

    cashflowstruct = &(deal_struct->cashflow_struct);
    NbFixedPrepay  = prepay_engine->NbFixedPrepay;
    NbCashFlow     = min((NbFixedPrepay), (cashflowstruct->numCashFlows));

    if (NbCashFlow == 0)
    {
        return (err);
    }

    bals = mbspt_dvector(0, NbCashFlow);  // 1 more so that we have 1.0 at 0
    dfs  = mbspt_dvector(0, NbCashFlow - 1);

    // fill bals
    o_wala  = prepay_engine->oterm - cashflowstruct->numCashFlows;
    wam     = prepay_engine->oterm - o_wala;
    bals[0] = 1.0;
    if (err = get_initial_principal_payments(
            NbCashFlow,
            bals + 1,
            wam,
            prepay_engine->gwac,
            prepay_engine->NbFixedPrepay,
            prepay_engine->FixedPrepay))
    {
        goto CLEANUP;
    }

    if (NbCashFlow ==
        cashflowstruct
            ->numCashFlows)  // if it is the case that after NbCashFlow exhaust all principals
    {
        for (i = 1; i < NbCashFlow; ++i)
        {
            if (bals[i] > 0.0)
                bals[i] = bals[i - 1] - bals[i];
            else
                bals[i] = 0.0;
        }
        bals[NbCashFlow] = 0.0;
    }
    else
    {
        for (i = 1; i <= NbFixedPrepay; ++i)
        {
            if (bals[i] > 0.0)
                bals[i] = bals[i - 1] - bals[i];
            else
                bals[i] = 0.0;
        }
    }

    // prepare OAS adjusted discount factors
    // get zeros
    if (err = MBS_BKTree_ZeroPrices(&NbZerosCached, &d_dptr, 1, 0, tree))
    {
        goto CLEANUP;
    }
    if (NbZerosCached < NbCashFlow)
    {
        err = "Not enough zeros cached";
        goto CLEANUP;
    }
    for (i = 0; i < NbCashFlow; ++i)
    {
        t      = cashflowstruct->CashFlowTimes[i];
        pos    = interp_search(NbZerosCached, tree->tree_times + i + 1, t);
        t0     = tree->tree_times[i + pos];
        t1     = tree->tree_times[i + 1 + pos];
        dfs[i] = LogInterp(t, t0, t1, d_dptr[pos + i][0], d_dptr[pos + i + 1][0]);
    }
    //////
    for (i = cashflowstruct->NbInitialCashFlowToSkip; (i < NbCashFlow) && (i < wam); ++i)
    {
        io_adjust += dfs[i] * (bals[i] - bals[NbCashFlow]);
        po_adjust += dfs[i] * (bals[i] - bals[i + 1]);
    }
    io_adjust *= deal_struct->cpn / ((double)deal_struct->cpnFreq);
    po_adjust *= 100.0;
    (*io) *= bals[NbCashFlow];
    (*po) *= bals[NbCashFlow];
    (*io) += io_adjust;
    (*po) += po_adjust;
    (*io) /= bals[cashflowstruct->NbInitialCashFlowToSkip];
    (*po) /= bals[cashflowstruct->NbInitialCashFlowToSkip];
CLEANUP:
    mbspt_free_dvector(bals, 0, NbCashFlow);
    mbspt_free_dvector(dfs, 0, NbCashFlow - 1);
    return (err);
}  // MBS_IOPO_Price_Adjust
///

// outputs io and po are pointers to single numbers
char* MBSPT_IOPO_Price(
    double*            io,
    double*            po,
    MBSPT_DealStruct*  deal_struct,
    MBS_Prepay_Engine* prepay_engine,
    MBS_BKTree*        tree)
{
    char*                 err = 0;
    double **             smm, *smm_at1;
    double **             io_prices, **po_prices;
    double *              weights, sum_of_weights;
    double                coupon, df;
    int                   i, j;
    int                   wam;  // as we use it, wam is as of trade date
    MBSPT_CashFlowStruct* cashflowstruct;
    //
    cashflowstruct = &(deal_struct->cashflow_struct);
    coupon         = deal_struct->cpn;
    // alloc mem:
    smm       = mbspt_dmatrix(0, prepay_engine->SFSNb - 1, 0, tree->NbSpaceNodes - 1);
    smm_at1   = mbspt_dvector(0, prepay_engine->SFSNb - 1);
    io_prices = mbspt_dmatrix(0, prepay_engine->SFSNb - 1, 0, tree->NbSpaceNodes - 1);
    po_prices = mbspt_dmatrix(0, prepay_engine->SFSNb - 1, 0, tree->NbSpaceNodes - 1);
    weights   = mbspt_dvector(0, prepay_engine->SFSNb - 1);
    // first initialize zeros
    if (err = mbs_iopo_initialize_zeros(tree, cashflowstruct->numCashFlows))
        return (err);
    // initialize IO and PO prices
    wam = cashflowstruct->numCashFlows;
    wam += month_diff(tree->RefDate, cashflowstruct->SettleDate);

    for (j = 0; j < prepay_engine->SFSNb; ++j)
    {
        if (err = mbs_iopo_initialize_prices(
                io_prices[j],
                po_prices[j],
                tree->Bottom[cashflowstruct->numCashFlows],
                tree->Top[cashflowstruct->numCashFlows]))
            goto CLEANUP;
    }
    // compute io po prices backward
    for (i = cashflowstruct->numCashFlows; i >= 0; --i)
    {
        // pre-calc smm
        if (err = MBS_PricesGetSMM(smm, i, cashflowstruct, prepay_engine, tree))
            goto CLEANUP;
        if (i == 1)
        // we will use deterministic smm to calculate new trigger weighting
        {
            for (j = 0; j < prepay_engine->SFSNb; ++j)
            {
                smm_at1[j] = smm[j][1];
                smm[j][0]  = smm_at1[j];
                smm[j][2]  = smm_at1[j];
            }
        }
        //
        for (j = 0; j < prepay_engine->SFSNb; ++j)
        {
            if (err = MBS_PricesBackwardInduct(
                    io_prices[j],
                    po_prices[j],
                    i,
                    coupon,
                    smm[j],
                    cashflowstruct,
                    prepay_engine,
                    tree))
                goto CLEANUP;
        }
    }
    // do weigted average
    (*io) = (*po) = 0.0;
    if (err = get_Trigger_Distribution(prepay_engine, deal_struct->TradeDate, wam, weights))
        goto CLEANUP;
    if (cashflowstruct->NbInitialCashFlowToSkip > 0 && prepay_engine->Delay == 0)
    {
        sum_of_weights = 0.0;
        for (j = 0; j < prepay_engine->SFSNb; ++j)
        {
            weights[j] *= (1.0 - smm_at1[j]);
            sum_of_weights += weights[j];
        }
        for (j = 0; j < prepay_engine->SFSNb; ++j)
        {
            weights[j] /= sum_of_weights;
        }
    }

    //
    for (i = 0; i < prepay_engine->SFSNb; ++i)
    {
        (*io) += weights[i] * io_prices[i][0];
        (*po) += weights[i] * po_prices[i][0];
    }

    if (err = MBS_IOPO_Price_Adjust(io, po, deal_struct, prepay_engine, tree))
        goto CLEANUP;
    //
    df = MBSDF(&(tree->df_data), deal_struct->cashflow_struct.SettleDate, deal_struct->SettleDate);
    df *=
        exp(-tree->OAS / 10000.0 *
            (double)day_count_date(
                deal_struct->cashflow_struct.SettleDate, deal_struct->SettleDate, BASIS_ACT_USD) /
            365.0);
    (*io) /= df;
    (*po) /= df;
    // accrued interest:
    (*io) -= cashflowstruct->accrualFactor * deal_struct->cpn;
// cleanup
CLEANUP:
    mbspt_free_dmatrix(smm, 0, prepay_engine->SFSNb - 1, 0, tree->NbSpaceNodes - 1);
    mbspt_free_dvector(smm_at1, 0, prepay_engine->SFSNb - 1);
    mbspt_free_dmatrix(io_prices, 0, prepay_engine->SFSNb - 1, 0, tree->NbSpaceNodes - 1);
    mbspt_free_dmatrix(po_prices, 0, prepay_engine->SFSNb - 1, 0, tree->NbSpaceNodes - 1);
    mbspt_free_dvector(weights, 0, prepay_engine->SFSNb - 1);
    return (err);
}  // MBSPT_IOPO_Price

char* MBSPT_MBSDF_Prepare_Helper(
    MBSDF_Data* df_data,
    int         NbRates,
    double*     rates,
    int         NbMMTerms,
    DATE_UNIT*  mm_date_units,
    int*        MMTerms,
    DATE_UNIT*  yield_date_units,
    int*        yield_terms,
    // mm:
    YBASIS    mm_ybase,
    DIFFBASIS mm_diffBase,
    int       mm_delay,
    // swap
    YBASIS    yield_ybase,
    DIFFBASIS yield_diffBase,
    int       yield_freq,
    int       yield_delay,
    // zero spec
    int  zero_freq,
    int  interp_freq,
    long RefDate)
{
    char*     err        = 0;
    Zero_Data zero_data  = _MBS_ALL_ZEROS;
    MM_Data   mm_data    = _MBS_ALL_ZEROS;
    YBASIS    zero_ybase = _365;

    double *MMYield, *yields;

    int        i;
    long       ValueDate  = RefDate;
    Yield_Data yield_data = _MBS_ALL_ZEROS;

    // checks:
    if (Dateok(RefDate))
        return ("Bad RefDate in Zero_Prepare");
    if ((yield_freq < 0) || (zero_freq < 0) || (interp_freq <= 0))
        return ("Illegal freq in Zero_Prepare");
    if ((mm_delay < 0) || (yield_delay < 0))
        return ("Illegal delays in Zero_Prepare");
    ////////////////

    MMYield = mbspt_dvector(0, NbMMTerms - 1);
    yields  = mbspt_dvector(0, NbRates - NbMMTerms - 1);

    // form mm_data
    for (i = 0; i < NbMMTerms; ++i)
        MMYield[i] = rates[i];
    if (err = MM_Data_Construct(
            &mm_data,
            mm_ybase,
            mm_diffBase,
            mm_delay,
            NbMMTerms,
            mm_date_units,
            MMTerms,
            MMYield,
            ValueDate))
        goto CLEANUPMM;

    // form yield_data
    for (i = NbMMTerms; i < NbRates; ++i)
        yields[i - NbMMTerms] = rates[i];
    if (err = Yield_Data_Construct(
            &yield_data,
            yield_ybase,
            yield_diffBase,
            yield_delay,
            yield_freq,
            NbRates - NbMMTerms,
            yield_date_units,
            yield_terms,
            yields,
            ValueDate))
        goto CLEANUPSWAP;
    // strip:
    zero_freq = yield_freq;
    if (err = MBS_Strip(
            (&zero_data), &mm_data, &yield_data, interp_freq, zero_freq, zero_ybase, ValueDate))
        goto CLEANUPSWAP;
    if (err = MBSDF_Data_Construct_FromZeroData(df_data, (&zero_data)))
        goto CLEANUPSWAP;
    // clean up
CLEANUPSWAP:
    Yield_Data_Destruct(&yield_data);
    goto CLEANUPMM;
CLEANUPMM:
    MM_Data_Destruct(&mm_data);
    Zero_Data_Destruct((&zero_data));

    mbspt_free_dvector(MMYield, 0, NbMMTerms - 1);
    mbspt_free_dvector(yields, 0, NbRates - NbMMTerms - 1);

    return (err);
}  // MBSPT_MBSDF_Prepare_Helper

char* MBSPT_MBSDF_Prepare_FromRatesAndFile(
    MBSDF_Data* df_data,
    char*       data_dir,
    int         NbRates,
    double*     rates,
    int         NbMMTerms,
    DATE_UNIT*  mm_date_units,
    int*        MMTerms,
    DATE_UNIT*  yield_date_units,
    int*        yield_terms,
    long        RefDate)
{
    char* err = 0;

    FILE*     stream;
    char      string[MAXSTRBUFSZ];
    YBASIS    mm_ybasis, yield_ybasis;
    DIFFBASIS mm_diffbasis, yield_diffbasis;
    int       mm_delay, yield_delay, yield_freq, zero_freq, interp_freq;

    if (!(stream = _OpenMBSFile(data_dir, FileTS)))
        return (nrerror("Can't open TermStruct.dat"));

    fgets(string, 80, stream);
    fscanf(stream, "%s\n", string);
    if (err = readBasis(&mm_ybasis, &mm_diffbasis, string))
        return (err);

    fgets(string, 80, stream);
    fscanf(stream, "%ld\n", &mm_delay);

    fgets(string, 80, stream);
    fscanf(stream, "%s\n", string);
    if (err = readBasis(&yield_ybasis, &yield_diffbasis, string))
        return (err);

    fgets(string, 80, stream);
    fscanf(stream, "%ld\n", &yield_delay);

    fgets(string, 80, stream);
    fscanf(stream, "%ld\n", &yield_freq);

    fgets(string, 80, stream);
    fscanf(stream, "%ld\n", &zero_freq);

    fgets(string, 80, stream);
    fscanf(stream, "%ld\n", &interp_freq);

    //

    if (err = MBSPT_MBSDF_Prepare_Helper(
            df_data,
            NbRates,
            rates,
            NbMMTerms,
            mm_date_units,
            MMTerms,
            yield_date_units,
            yield_terms,
            mm_ybasis,
            mm_diffbasis,
            mm_delay,
            yield_ybasis,
            yield_diffbasis,
            yield_freq,
            yield_delay,
            zero_freq,
            interp_freq,
            RefDate))
        goto CLEANUP;

CLEANUP:

    return (err);
}  // MBSPT_MBSDF_Prepare_FromRatesAndFile

// the following is a special zero_data construction
// for the MBS_BKTree used in the first generation mortgage pricer (as of September 2003)
char* MBSPT_MBSDF_Prepare_FromRatesAndFile_Special(
    MBSDF_Data* df_data, char* data_dir, double* tsycurve, double* swpSprd, long RefDate, long spot)
{
#ifdef MBSPTCORRECTDAYCT
    int       NbMMTerms       = 5;
    DATE_UNIT mm_date_units[] = {DAY, MONTH, MONTH, MONTH, MONTH};
    int       MMTerms[]       = {1, 1, 2, 3, 6};
#else
    int       NbMMTerms          = 6;
    DATE_UNIT mm_date_units[]    = {DAY, MONTH, MONTH, MONTH, MONTH, MONTH};
    int       MMTerms[]          = {1, 1, 2, 3, 6, 12};
#endif  // MBSPTCORRECTDAYCT

#ifdef MBSPTCORRECTDAYCT
    int       NbYields           = 9;
    DATE_UNIT yield_date_units[] = {YEAR, YEAR, YEAR, YEAR, YEAR, YEAR, YEAR, YEAR, YEAR};
    int       yield_terms[]      = {1, 2, 3, 4, 5, 7, 10, 20, 30};
#else
    int       NbYields           = 8;
    DATE_UNIT yield_date_units[] = {YEAR, YEAR, YEAR, YEAR, YEAR, YEAR, YEAR, YEAR};
    int       yield_terms[]      = {2, 3, 4, 5, 7, 10, 20, 30};
#endif  // MBSPTCORRECTDAYCT

    double rates[14];
    int    i;
    long   CrvDate;

#ifdef TREEREFDATEISSETTLEDATE
    CrvDate = Date_Add(RefDate, DAY, (int)spot, 1, 0);
#else
    CrvDate                      = RefDate;
#endif  // TREEREFDATEISSETTLEDATE

    for (i = 0; i < 6; ++i)
        rates[i] = tsycurve[i];
    for (i = 0; i < 8; ++i)
        rates[i + 6] = tsycurve[i + 6] + swpSprd[i] / 100.0;

    return (MBSPT_MBSDF_Prepare_FromRatesAndFile(
        df_data,
        data_dir,
        NbMMTerms + NbYields,
        rates,
        NbMMTerms,
        mm_date_units,
        MMTerms,
        yield_date_units,
        yield_terms,
        CrvDate));
}  // MBSPT_MBSDF_Prepare_FromRatesAndFile_Special

char* MBSPT_Structures_Prepare_FromFile(
    MBS_BKTree*        tree,
    MBSPT_DealStruct*  deal_struct,
    MBS_Prepay_Engine* prepay_engine,
    char*              data_dir,
    long               valuation,
    long               FwdSettle,
    MBSDF_Data*        df_data,
    double             refisprd,
    double Rate10Y,  // needed in MBS_Prepay_Engine_SetParams for converting ratio form of
                     // refi-incent to diff form
    MTGPTACYPROG AgcyProg,
    int          NbSwaptions,
    double*      volcurve,
    long*        volexp,
    double       meanrev,
    long         swapmat,
    long*        dealmat,
    double*      coupon,
    double       inoas,
    long         paydelays,
    double*      pptweaks)
{
    int i;

    // double swaptions_exp[5];
    double* swaptions_exp;
    long    TradeDate, TreeRefDate;

    char* err = 0;

    //	Rate10Y = rates[11];

    TradeDate = valuation;
    //	SpotSettleDate = Date_Add( valuation, DAY, sptdays, 1, 0);

    //#ifdef TREEREFDATEISSETTLEDATE
    //	TreeRefDate = SpotSettleDate;
    //#else
    TreeRefDate = TradeDate;
    //#endif//TREEREFDATEISSETTLEDATE

    swaptions_exp = mbspt_dvector(0, NbSwaptions - 1);
    for (i = 0; i < NbSwaptions; ++i)
    {
        swaptions_exp[i] = (double)volexp[i];
    }

    if (err = MBS_BKTree_Construct_FromFile(
            tree,
            data_dir,
            df_data,
            NbSwaptions,
            swapmat,
            swaptions_exp,
            volcurve,
            meanrev,
            inoas,
            Rate10Y,
            refisprd,
            TreeRefDate))
        goto CLEANUP;

    if (err = MBSPT_DealStruct_Construct_FromFile(
            deal_struct, data_dir, dealmat[0], dealmat[1], coupon[0], paydelays))
        goto CLEANUP;

    if (err = MBSPT_DealStruct_SetTimes(deal_struct, TradeDate, FwdSettle))
        goto CLEANUP;

    // prepay

    if (err = MBS_Prepay_Engine_Params_Reader(
            data_dir,
            prepay_engine,
            AgcyProg,
            dealmat[1],
            coupon[1],
            pptweaks[0],
            pptweaks[1],
            (int)swapmat))
        goto CLEANUP;
    if (err = MBS_Prepay_Engine_SetParams(prepay_engine, tree->RefiRate, tree->RefiBaseIndex))
        goto CLEANUP;
    if (err = MBS_Prepay_Engine_Check(prepay_engine))
        goto CLEANUP;

    // calibartion
    MBS_BKTree_create_times(tree, deal_struct, prepay_engine);

    if (err = MBS_BKTree_Calibrate(tree))
        goto CLEANUP;
    // set cach

    if (err = BKTreeCachStruct_Reset(
            (tree->cach),
            1 + prepay_engine->IndexMat * deal_struct->cpnFreq,  // prepay_engine->IndexFreq,
            prepay_engine->NbFixedPrepay + 2 + (deal_struct->PayDelay / 30)))
        goto CLEANUP;

CLEANUP:
    mbspt_free_dvector(swaptions_exp, 0, NbSwaptions - 1);
    return (err);
}  // MBSPT_Structures_Prepare_FromFile

char* MBSPT_Pricer(
    double             outputs[3],
    MBS_BKTree*        tree,
    MBSPT_DealStruct*  deal_struct,
    MBS_Prepay_Engine* prepay_engine)
{
    double io, po;
    char*  err = 0;
    if (err = MBSPT_IOPO_Price(&io, &po, deal_struct, prepay_engine, tree))
        return (err);
    //
    outputs[0] = io;
    outputs[1] = po;
    outputs[2] = io + po;
    //
    return (err);
}  // MBSPT_Pricer

char* MBSPT_OAS_Solver(
    double* out,
    double  mktPrice,
    char*   data_dir,
    long    valuation,
    long    FwdSettle,
    char*   mkt,
    double  refisprd,
    double  Rate10Y,  // needed in MBS_Prepay_Engine_SetParams for converting ratio form of
                     // refi-incent to diff form
    MTGPTACYPROG AgcyProg,
    int          NbSwaptions,
    double*      volcurve,
    long*        volexp,
    double       meanrev,
    long         swapmat,
    long*        dealmat,
    double*      coupon,
    double       guessed_oas,
    long         paydelays,
    double*      pptweaks,
    double       tolerance,
    int          max_trial,
    double       oas_inc  // in bp
)
{
    MBS_BKTree        tree          = _MBS_ALL_ZEROS;
    MBSPT_DealStruct  deal_struct   = _MBS_ALL_ZEROS;
    MBS_Prepay_Engine prepay_engine = _MBS_ALL_ZEROS;
    MBSDF_Data        df_data       = _MBS_ALL_ZEROS;

    ////
    int    it;
    double curr_oas, prev_oas, px;
    double pxs[3], sprdDur;
    char*  err = 0;

    curr_oas = guessed_oas;
    prev_oas = curr_oas + oas_inc;
    it       = 0;

    while (it <= max_trial)
    {
        if (fabs(curr_oas - prev_oas) < tolerance)
        {
            *out = curr_oas;
            return (err);
        }
        // prepare and price for base price
        if (err = MBSDF_Data_Construct_FromName((&df_data), mkt))
        {
            MBSPT_Structures_Cleanup(&df_data, &tree, &deal_struct, &prepay_engine);
            return (err);
        }
        if (err = MBSPT_Structures_Prepare_FromFile(
                &tree,
                &deal_struct,
                &prepay_engine,
                data_dir,
                valuation,
                FwdSettle,
                &df_data,
                refisprd,
                Rate10Y,
                AgcyProg,
                NbSwaptions,
                volcurve,
                volexp,
                meanrev,
                swapmat,
                dealmat,
                coupon,
                curr_oas,
                paydelays,
                pptweaks))
        {
            MBSPT_Structures_Cleanup(&df_data, &tree, &deal_struct, &prepay_engine);
            return (err);
        }
        if (err = MBSPT_Pricer(pxs, &tree, &deal_struct, &prepay_engine))
        {
            MBSPT_Structures_Cleanup(&df_data, &tree, &deal_struct, &prepay_engine);
            return (err);
        }
        MBSPT_Structures_Cleanup(&df_data, &tree, &deal_struct, &prepay_engine);
        if (fabs(pxs[2] - mktPrice) < tolerance)
        {
            *out = curr_oas;
            return (err);
        }

        MBS_INIT_ZERO(df_data);
        MBS_INIT_ZERO(tree);
        MBS_INIT_ZERO(deal_struct);
        MBS_INIT_ZERO(prepay_engine);

        px = pxs[2];

        // prepare and price for base price for shocked price
        if (err = MBSDF_Data_Construct_FromName((&df_data), mkt))
        {
            MBSPT_Structures_Cleanup(&df_data, &tree, &deal_struct, &prepay_engine);
            return (err);
        }
        if (err = MBSPT_Structures_Prepare_FromFile(
                &tree,
                &deal_struct,
                &prepay_engine,
                data_dir,
                valuation,
                FwdSettle,
                &df_data,
                refisprd,
                Rate10Y,
                AgcyProg,
                NbSwaptions,
                volcurve,
                volexp,
                meanrev,
                swapmat,
                dealmat,
                coupon,
                curr_oas + oas_inc,
                paydelays,
                pptweaks))
        {
            MBSPT_Structures_Cleanup(&df_data, &tree, &deal_struct, &prepay_engine);
            return (err);
        }
        if (err = MBSPT_Pricer(pxs, &tree, &deal_struct, &prepay_engine))
        {
            MBSPT_Structures_Cleanup(&df_data, &tree, &deal_struct, &prepay_engine);
            return (err);
        }

        MBSPT_Structures_Cleanup(&df_data, &tree, &deal_struct, &prepay_engine);

        MBS_INIT_ZERO(df_data);
        MBS_INIT_ZERO(tree);
        MBS_INIT_ZERO(deal_struct);
        MBS_INIT_ZERO(prepay_engine);
        //
        sprdDur = (pxs[2] - px) / oas_inc;
        if (fabs(sprdDur) < tolerance / 100.0)
            return ("Solver fail");
        prev_oas = curr_oas;
        curr_oas = prev_oas + (mktPrice - px) / sprdDur;
        if (fabs(prev_oas - curr_oas) < tolerance)
        {
            *out = curr_oas;
            return (err);
        }
        it++;
    }
    if (it == max_trial)
        return ("Solver fails");
    return (err);

}  // MBSPT_OAS_Solver

// indices[i][j] means the value of index at time slice i and space node j
// indices is pre-alloc by user should be [0..endIndex][0..tree->NbSpaceNodex-1]
// endIndex should be cashflowstruct->numCashFlows
char* MBS_Prepay_GetIndexValue(
    double** indices, int endIndex, MBS_BKTree* tree, MBS_Prepay_Engine* prepay_engine)
{
    char* err = 0;
    int   i;
    // alloc mem
    // first initialize zeros
    if (err = mbs_iopo_initialize_zeros(tree, endIndex))
        return (err);
    for (i = endIndex; i >= 0; --i)
    {
        if (err = MBS_BKTree_Yield_Approx(
                indices[i], prepay_engine->IndexFreq, prepay_engine->IndexMat, i, tree))
            return (err);
    }
    return (err);
}  // MBS_Prepay_GetIndexValue

void MBSPT_Structures_Cleanup(
    MBSDF_Data*        df_data,
    MBS_BKTree*        tree,
    MBSPT_DealStruct*  deal_struct,
    MBS_Prepay_Engine* prepay_engine)
{
    MBSDF_Data_Destruct(df_data);
    MBS_BKTree_Destruct(tree);
    MBSPT_DealStruct_Destruct(deal_struct);
    MBS_Prepay_Engine_Cleanup(prepay_engine);
}  // MBSPT_StructuresCleanup

char* MBS_Prepay_Projection_ParamsFromFile(
    char*        data_dir,
    MTGPTACYPROG AgencyProg,
    long         oterm,
    double       gwac,
    double       Speedup,
    double       rational,
    double       RefiRate,
    double       Rate10Y,
    long         TradeYYYYMMDD,
    long         spotYYYYMMDD,
    long         wam,
    int          subgroup,
    int          num_future_rates,
    double*      futureRates,
    int          useFixedPrepay,
    int*         numProjections,
    double       rates[MAXOTERM],
    double       cpr[MAXOTERM],
    double       principal_payments[MAXOTERM],
    int          IndexMat)
{
    MBS_Prepay_Engine prepay_engine = {0};
    char*             err;

    if (err = MBS_Prepay_Engine_Params_Reader(
            data_dir, &prepay_engine, AgencyProg, oterm, gwac, Speedup, rational, IndexMat))
    {
        MBS_Prepay_Engine_Cleanup(&prepay_engine);
        return (err);
    }
    if (err = MBS_Prepay_Engine_SetParams(&prepay_engine, RefiRate, Rate10Y))
    {
        MBS_Prepay_Engine_Cleanup(&prepay_engine);
        return (err);
    }
    if (err = MBS_Prepay_Engine_Check(&prepay_engine))
    {
        MBS_Prepay_Engine_Cleanup(&prepay_engine);
        return (err);
    }
    err = MBS_Prepay_Path(
        &prepay_engine,
        TradeYYYYMMDD,
        spotYYYYMMDD,
        wam,
        subgroup,
        num_future_rates,
        futureRates,
        useFixedPrepay,
        numProjections,
        rates,
        cpr,
        principal_payments);
    MBS_Prepay_Engine_Cleanup(&prepay_engine);
    return (err);
}  // MBS_Prepay_Projection_ParamsFromFile

char* MBS_Prepay_Projection_ParamsExplicit(
    char*        data_dir,
    MTGPTACYPROG AgencyProg,
    long         oterm,
    double       gwac,
    double       Speedup,
    double       rational,
    double       FebruaryEffect,
    double*      IndexLevel,
    double*      AmortLevel,
    double*      SeasoningLevel,
    double*      SeasoiningMat,
    int          SFSNb,
    double*      Amorts,
    double       Seasonality[12],
    int          Delay,
    int          NbFixedPrepay,
    double*      FixedPrepay,
    int          nuniform,
    double       RefiRateSpread,
    double       Rate10Y,
    long         TradeYYYYMMDD,
    long         spotYYYYMMDD,
    long         wam,
    int          subgroup,
    int          num_future_rates,
    double*      futureRates,
    int          useFixedPrepay,
    int*         numProjections,
    double       rates[MAXOTERM],
    double       cpr[MAXOTERM],
    double       principal_payments[MAXOTERM],
    int          IndexMat,
    int          EndIndexMat,
    int          IndexFreq)
{
    MBS_Prepay_Engine prepay_engine = {0};
    char*             err;

    if (err = MBS_Prepay_Engine_Params_ReaderExplicit(
            data_dir,
            &prepay_engine,
            AgencyProg,
            oterm,
            gwac,
            Speedup,
            rational,
            FebruaryEffect,
            IndexLevel,
            AmortLevel,
            SeasoningLevel,
            SeasoiningMat,
            SFSNb,
            Amorts,
            Seasonality,
            Delay,
            NbFixedPrepay,
            FixedPrepay,
            nuniform,
            IndexMat,
            EndIndexMat,
            IndexFreq))
    {
        MBS_Prepay_Engine_Cleanup(&prepay_engine);
        return (err);
    }
    if (err = MBS_Prepay_Engine_SetParams(&prepay_engine, RefiRateSpread, Rate10Y))
    {
        MBS_Prepay_Engine_Cleanup(&prepay_engine);
        return (err);
    }
    if (err = MBS_Prepay_Engine_Check(&prepay_engine))
    {
        MBS_Prepay_Engine_Cleanup(&prepay_engine);
        return (err);
    }
    err = MBS_Prepay_Path(
        &prepay_engine,
        TradeYYYYMMDD,
        spotYYYYMMDD,
        wam,
        subgroup,
        num_future_rates,
        futureRates,
        useFixedPrepay,
        numProjections,
        rates,
        cpr,
        principal_payments);
    MBS_Prepay_Engine_Cleanup(&prepay_engine);
    return (err);
}  // MBS_Prepay_Projection_ParamsExplicit