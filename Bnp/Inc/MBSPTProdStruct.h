#ifndef MBSPTPRODSTRUCT_H
#define MBSPTPRODSTRUCT_H

typedef enum MTGPTACYPROG
{
    _ILLEGAL = -1,
    _FN      = 0,  // FN Conventional
    _FN15,
    _GN,  // GN Conventional
    _GN15
} MTGPTACYPROG;

MTGPTACYPROG MBSPT_ProgramLookup(char* code);

struct tagMBSPT_DealStruct;
typedef struct tagMBSPT_DealStruct MBSPT_DealStruct;
struct tagMBSPT_CashFlowStruct;
typedef struct tagMBSPT_CashFlowStruct MBSPT_CashFlowStruct;

struct tagMBSPT_CashFlowStruct
{
    long DealMaturityDate;  // for FNM, this is always the first of the month, this plus 24 days =
                            // last payment day
    long   SettleDate;  // spotsettle, i.e., trade date
    int    oterm;       // in months, for 30Y it is 360
    int    cpnFreq;
    double cpn;            // 6 for 6%
    int    PayDelay;       // in days, for FNM it is 24
    double accrualFactor;  // this multiplied with face amount (100 for $100) and cpn (5 for 5%)
                           // gives the accrued interest
    int   numCashFlows;  // number of coupons that the new owner on SettleDate will own
    int   NbInitialCashFlowToSkip;
    long* duePeriodEndDates;  // always the first of month for FNM, last duePeriodEndDate =
                              // DealMaturityDate
    double* duePeriodEndTimes;  // (duePeriodEndDates - SettleDate)/365
    long*   CashFlowDates;      // duePeriodEndDate + PayDelay
    double* CashFlowTimes;
};

char* MBSPT_CashFlowStruct_Construct(
    MBSPT_CashFlowStruct* cashflow_struct, long spotSettleDate, MBSPT_DealStruct* deal_struct);

void MBSPT_CashFlowStruct_Destruct(MBSPT_CashFlowStruct* cashflow_struct);

struct tagMBSPT_DealStruct
{
    long DealMaturityDate;  // for FNM, this is always the first of the month, this plus 24 days =
                            // last payment day
    MTGPTACYPROG         AgencyProgram;  // for FNM conventional it is _FNCONVENTIONAL
    int                  oterm;          // in months, for 30Y it is 360
    int                  cpnFreq;        // 12
    double               cpn;            // 6 for 6%
    int                  PayDelay;       // in days, for FNM it is 24, for FH Arm it is 44
    long                 TradeDate;
    long                 SettleDate;  // fwdsettle
    MBSPT_CashFlowStruct cashflow_struct;
};

char* MBSPT_DealStruct_Construct(
    MBSPT_DealStruct* deal_struct,
    long              DealMaturityDate,
    MTGPTACYPROG      AgencyProgram,
    int               oterm,
    int               cpnFreq,
    double            cpn,
    int               PayDelay);

char* MBSPT_DealStruct_Construct_FromFile(
    MBSPT_DealStruct* deal_struct,
    char*             data_dir,
    long              DealMaturityDate,
    int               oterm,
    double            cpn,
    int               PayDelay);

void MBSPT_DealStruct_Destruct(MBSPT_DealStruct* deal_struct);

char* MBSPT_DealStruct_SetTimes(MBSPT_DealStruct* deal_struct, long TradeDate, long settleDate);

char* MBSPT_AddCashFlowPV(
    double*  io,
    double*  po,  // in an dout
    int      start,
    int      end,
    double   coupon,
    double   sch_amort_rate,
    double*  smm,
    int      nb_zeros,
    double** zeros,
    int      NbZerosToSkip);

#endif  // MBSPTPRODSTRUCT_H