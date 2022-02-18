#ifndef MBSPREPAY_H
#define MBSPREPAY_H

#include "MBSPPFuncs.h"
#include "MBSPTProdStruct.h"
#include "MBSPTUtil.h"

typedef struct MBS_Prepay_Engine_Tag
{
    int      NbSCurve;      // Number of basic S-curves input, usually = NBPOINTS
    int      NbIndexLevel;  // Number of Index Values in S-curve, usually = NBPOINTS
    double*  IndexLevel;    // Index level that define the prepayment function
    double** AmortLevel;    // Corresponding annual CPR level or Avg Life
                          // AmortLevel[0..NbSCurve-1][0..NbIndexLevel-1]
    double  Seasonality[12];  // Seasonality level for each month
    double  FebruaryEffect;
    int     NbSeasoningRamps;  // Number of seasoning ramps, usually 3
    double* SeasoningLevel;    // Levels of Refi/Coupon that define seasoning (cf SeasoningMat)
    double* SeasoningMat;      // Maturities that define seasoning
    double  Speedup,           // Speed up the PPF
        rational;              // rationality speed up
    int          SFSNb;        // Number of SFS per group
    double*      TriggerRate;  // Trigger rate for each SFS
    double*      Amort;        // Additional amortization once the trigger is reached for each SFS
    int          Delay;        // prepay_delay
    int          NbFixedPrepay;
    double*      FixedPrepay;
    MTGPTACYPROG AgencyProg;
    long         oterm;
    double       gwac;
    HISTORY      history;
    MBSDENSITY   density;
    long         nuniform;
    // index:
    int IndexMat;     // in years
    int EndIndexMat;  // if < IndexMat use it as index at end and progressively length as we
                      // backward induct, this is an approximation technique
    int IndexFreq;
} MBS_Prepay_Engine;

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
    double*            SeasoiningMat,
    int                SFSNb,
    double*            Amorts,
    double             Seasonality[12],
    int                Delay,
    int                NbFixedPrepay,
    double*            FixedPrepay,
    int                nuniform,
    int                IndexMat,
    int                EndIndexMat,
    int                IndexFreq);

char* MBS_Prepay_Engine_Params_Reader(
    char*              data_dir,
    MBS_Prepay_Engine* prepay_engine,
    MTGPTACYPROG       AgencyProg,
    long               oterm,
    double             gwac,
    double             Speedup,
    double             rational,
    // double FebruaryEffect,
    int IndexMat);

char* MBS_Prepay_Engine_SetParams(
    MBS_Prepay_Engine* prepay_engine, double refiRate, double Rate10Y);

char* MBS_Prepay_Engine_Check(MBS_Prepay_Engine* prepay_engine);

void MBS_Prepay_Engine_Initialize(
    MBS_Prepay_Engine* prepay_engine,
    int                NbSCurve,
    int                NbIndexLevel,
    double*            IndexLevel,
    double**           AmortLevel,
    double             Seasonality[],
    double             FebruaryEffect,
    int                NbSeasonaingRamps,
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
    int                IndexFreq);

void MBS_Prepay_Engine_ReplaceTriggers(
    MBS_Prepay_Engine* new_engine,
    MBS_Prepay_Engine* prepay_engine,
    int                num_buckets,
    double*            triggerRate,
    double*            amort);

void MBS_Prepay_Engine_Cleanup(MBS_Prepay_Engine* prepay_engine);
///
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
    double**           factors);
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
    double             principal_payments[MAXOTERM]);

void get_smm_weights(
    double* weights, MBS_Prepay_Engine* prepay_engine, double* usweights, double* ufweights);
int MBS_Prepay_Engine_SMM(
    double*            output,
    MBS_Prepay_Engine* prepay_engine,
    int                age,
    double             IndexValue,
    int                monthOfYear,
    int                monthFromSettle);
/////////////////////////////
/***********Functions:*********************/

void spdbase(double* AmortLevel, double Speedup, double rational, int numpts);

double alsmm1(double al, int N);
double alsmm2(double al, double wac, int N, int wam, int pt);

char* get_Trigger_Distribution(
    MBS_Prepay_Engine* prepay_engine, long spotYYYYMMDD, long wam, double* weights);

char* get_initial_principal_payments(
    int num_payments, double* out, int wam, double gwac, int NbFixedPrepay, double* FixedPrepay);

// void test_MBS_Prepay_Path( MBS_Prepay_Engine * prepay_engine, long spotYYYYMMDD, long  wam, int
// num_future_rates);

#endif  // MBSPREPAY_H