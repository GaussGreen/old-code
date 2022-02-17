#ifndef TBACALIB_H
#define TBACALIB_H

#include "MBSPPFuncs.h"
#include "MBSPTCalib.h"
#include "MBSPTCaller.h"
#include "TBAUtil.h"

/*******************************
**Functions:
********************************/
double ExpDecay(double a, double t);

void Zero_6m(TERM_DATA *term_data);
// char* Zero_6m_new (TERM_DATA *term_data  , Zero_Data * zero_data);

void yldconv(
    double *Zero,    /* output: zeros at money market points */
    double *ZeroS,   /* output: zeros at mm points */
    long *Period,    /* output: length of money market periods */
    double *MMYield, /* Money Market rates (O/N  , 1  ,2  ,3  ,6  ,12 m) */
    char AoS,        /* 'A'nnual or 'S'emi-Annual frequency */
    long ValueDate,  /* Value date for the benchmark instruments */
    int MMB          /* Money Market basis (360 or 365) */
);

void coupint(
    double *SwapYield); /* Array of interpolated and benchmark swap yields */

void spotrate(
    double *Zero,      /* output: Zero rates at 6m interval */
    double *ZeroS,     /* output: ZeroS rates at 6m interval */
    long *Period,      /* output: length of periods at 6m interval */
    double *SwapYield, /* array of interpolated and benchmark swap yields */
    double *YieldS,    /* array of interpolated and benchmark yieldS */
    char AoS,          /* 'A'nnual or 'S'emi-Annual frequency */
    long ValueDate,    /* Value date for the benchmark instruments */
    char *YBase        /* Year basis for the benchmark swaps ("ACT" or "365") */
);

/***
char * Time_new(
                          MBS_BKTree * tree  ,
                          MBSPT_CashFlowStruct * cashflow_struct  ,
                          MBS_Prepay_Engine * prepay_engine  ,
                  TREE_DATA * tree_data  ,
                  TERM_DATA * term_data  ,
                  DEAL_DATA * deal_data
                  );
***/

void Time(long ValueDate,      /* Value date */
          int Ppy,             /* Number of period per year in the tree */
          long Maturity,       /* Maturity of the deal */
          int NbExer,          /* Number of exercise dates */
          long *Exer,          /* Array of exercise dates */
          int NbFwd,           /* Number of forward prices */
          long *FwdMat,        /* Array of forward maturities */
          double *Strike,      /* Array of the corresponding strikes */
          char Freq,           /* Frequency of the underlying cash-flows */
          char EoA,            /* European or American option */
          TREE_DATA *tree_data /* Output: Structure of tree data */
);

void Tree_Alloc(TREE_DATA *tree_data, int NbNodes);

/***
char *    Build_Tree_New (
                        MBS_BKTree * tree  ,
                        Zero_Data * zero_data  ,
                        Vol_Surface * vol_surface  ,
                        double OAS  ,
                                                TREE_DATA       *tree_data);
**/

char *Build_Tree(TERM_DATA *term_data, TREE_DATA *tree_data,
                 PREPAY_DATA *prepay_data);

/***
char * Forward_new( TREE_DATA * old_tree  ,
                                 MBS_BKTree * new_tree  ,
                                 Zero_Data * zero_data );
**/

void Forward(
    double *ZeroCoupon, /* Output: zero coupon price at each node in the tree */
    double *FwdRate,    /* Output: forward rate at each node in the tree */
    double *Zero,       /* Zero rates at money market points and 6m interval */
    long *Period,       /* Length of money market periods and 6m intervals */
    double *Length,     /* Length of each time step */
    int NbNodes,        /* Total number of nodes (or time steps) */
    char AoS);

void Forward_Prepay(
    double *FwdPrep, /* Output: prepayment at each node in the tree */
    double Spot,     /* Spot prepayment */
    int NbNodes);    /* Total number of nodes (or time steps) */

/**
char * Forward_Vol_new( TREE_DATA * old_tree  ,
                                           MBS_BKTree * new_tree  ,
                                           Vol_Surface * vol_surface  ,
                                           Zero_Data   * zero_data );
**/

char *Forward_Vol(
    double
        *FwdVol,  /* Output: forward rate volatility at each node in the tree */
    double *Jump, /* Output: Size of the jump (in log space) */
    double *FwdRate,     /* Forward rate at each node */
    double *Length,      /* Length of each time step */
    double *Zero,        /* Zero rates at money market points and 6m interval */
    long *Period,        /* Length of money market periods and 6m intervals */
    double *SwaptionVol, /* Volatility of the benchmark swaptions */
    int *SwaptionExp,    /* Expiration (in years) of the benchmark swaptions */
    double Beta,         /* Mean reversion coefficient */
    int IndexMat,        /* Index maturity in years */
    char IndexFreq,      /* Index frequency */
    char AoS,            /* Yield curve type (annual or semi-annual) */
    int Ppy,             /* Number of periods per year in the tree */
    int NbNodes);        /* Total number of nodes (or time steps) */

void Prepay_Vol(double *PrepayVol, /* Output: instantaneous volatility of
                                      forward prepayment */
                double SpotVol,    /* Spot volatility of prepayment */
                int NbNodes);      /* Total number of nodes (or time steps) */

void Conv_Adj(
    double *FwdRate, /* Calculate convexity adjustment for discounting rates */
    double *FwdRateS, double *FwdVol, double *Length, double *Delta,
    int NbNodes, double Beta);

char *Drift_new(TREE_DATA *old_tree, MBS_BKTree *new_tree);

char *Drift(
    double *Drift, /* Output: drift of the 1 period interest rate in the tree */
    int *NodeMax,  /* Output: maximum node reached in the tree (both top and
                      bottom) */
    double Jump,   /* Size of the jump (in log space) */
    double *ZeroCoupon, /* Zero coupon price at each node in the tree */
    double *FwdRate,    /* Forward rate at each node in the tree */
    double *FwdVol,     /* Forward rate volatility at each node in the tree */
    double *Length,     /* Length of each time step */
    double Beta,        /* Mean reversion coefficient (constant) */
    int NbNodes);       /* Total number of nodes (or time steps) */

void NodeMax2(
    double *Jump2,     /* Output: size of the prepayment jump (in log space) */
    int *NodeMax2,     /* Output: maximum node reached in the tree (both top and
                          bottom) */
    double SpotVol,    /* Spot volatility of prepayment */
    double *PrepayVol, /* Instantaneous volatility of forward prepayment */
    double Beta,       /* Mean reversion of prepayment */
    double *Length,    /* Length of each time step */
    int NbNodes);      /* Total number of nodes (or time steps) */

void Print_Term(TERM_DATA *term_data, TREE_DATA *tree_data);

void Lattice_new(double *Discount, /* Output: one period discount factor at each
                                      node at period NbPer */
                 double *DiscountS, /* Output: one period discount factor
                                       including Option Adjusted spread */
                 double *pu, /* Output: array of probabilities in the up state
                                in the interest rate tree */
                 double *pd, /* Output: array of probabilities in the down state
                                in the interest rate tree */
                 double *p0, /* Output: array of probabilities in the flat state
                                in the interest rate tree */
                 int NbPer,  /* Current time period */
                 MBS_BKTree *tree);

void Lattice(
    double *Discount,  /* Output: one period discount factor at each node at
                          period NbPer */
    double *DiscountS, /* Output: one period discount factor including Option
                          Adjusted spread */
    double *pu,        /* Output: array of probabilities in the up state in the
                          interest rate tree */
    double *pd, /* Output: array of probabilities in the down state in the
                   interest rate tree */
    double *p0, /* Output: array of probabilities in the flat state in the
                   interest rate tree */
    int NbPer,  /* Current time period */
    int N,      /* Number of additional periods at the beginning of the tree */
    TERM_DATA *term_data,  /* Structure of term structure data */
    TREE_DATA *tree_data); /* Structure of tree data */

void Zero_Calc_new(
    double ***Zero, /* Set of zeros maturing on each coupon date */
    //	int             *NbZero  ,                                /* Current
    //number of zeros being priced */ 	int             TotNbZero  , /* Total
    //number of zeros needed to calculate the par yield index */ 	int Reset  , /*
    //Reset the zero price to 1 if TRUE */
    int NbPer, /* Current time period */
    //	int             TotPer  ,                                 /* Total
    //number of period in the rate tree including the N additional periods */
    int useOAS, MBS_BKTree *tree);

void Zero_Calc(double **Zero, /* Set of zeros maturing on each coupon date */
               int *NbZero,   /* Current number of zeros being priced */
               int TotNbZero, /* Total number of zeros needed to calculate the
                                 par yield index */
               int Reset,     /* Reset the zero price to 1 if TRUE */
               double *Discount, /* One period discount factor at each node at
                                    period NbPer */
               double *pu, /* Array of probabilities in the up state in the
                              interest rate tree */
               double *pd, /* Array of probabilities in the down state in the
                              interest rate tree */
               double *p0, /* Array of probabilities in the flat state in the
                              interest rate tree */
               int NbPer,  /* Current time period */
               int TotPer, /* Total number of period in the rate tree including
                              the N additional periods */
               TREE_DATA *tree_data); /* Structure of tree data */

void Zero_Price_new(
    double *Zero, /* Array of zero prices at the current period */
    int Reset,    /* Reset the zero prices to 1 if TRUE */
    int NbPer,    /* Current time period */
    MBS_BKTree *tree);

void Zero_Price(double *Zero, /* Array of zero prices at the current period */
                int Reset,    /* Reset the zero prices to 1 if TRUE */
                double *Discount, /* One period discount factor at each node at
                                     period NbPer */
                double *pu, /* Array of probabilities in the up state in the
                               interest rate tree */
                double *pd, /* Array of probabilities in the down state in the
                               interest rate tree */
                double *p0, /* Array of probabilities in the flat state in the
                               interest rate tree */
                int NbPer,  /* Current time period */
                TREE_DATA *tree_data); /* Structure of tree data */

void Dev(double *Price,    /* Array of prices that has to be discounted */
         double *Discount, /* One period discount factor at each node at period
                              NbPer */
         double *pu, /* Array of probabilities in the up state in the interest
                        rate tree */
         double *pd, /* Array of probabilities in the down state in the interest
                        rate tree */
         double *p0, /* Array of probabilities in the flat state in the interest
                        rate tree */
         int NbPer,  /* Current time period */
         TREE_DATA *tree_data); /* Structure of tree data */

void Par_Yield(
    double
        *ParYield, /* Output: par yield at every node at the current period */
    double **Zero, /* Set of zeros maturing on each coupon date */
    int TotNbZero, /* Total number of zeros needed to calculate the par yield
                      index */
    int F,         /* Frequency of underlying payment as a integer */
    int IndexF,    /* Frequency of index payment as a integer */
    int NodeMax,   /* Maximum node reached in the tree (both top and bottom) */
    int NbPer);    /* Current time period */

/****
char * Calibrate(
        MBS_BKTree * tree  ,
        MBSPT_CashFlowStruct * cashflow_struct  ,
        double OAS  ,
        MBS_Prepay_Engine * prepay_engine  ,
        Zero_Data * zero_data  ,
        Vol_Surface * vol_surface  ,
/////
        TREE_DATA * tree_data  ,
        TERM_DATA * term_data  ,
        DEAL_DATA * deal_data);
***/

char *copyTermStruct(MBS_BKTree *tree, MBSPT_DealStruct *deal_struct,
                     TREE_DATA *tree_data);

#endif