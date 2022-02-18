#ifndef TBAPRODSTRUCT_H
#define TBAPRODSTRUCT_H

#include "MBSPPFuncs.h"
#include "MBSPTCalib.h"
#include "MBSPTProdStruct.h"
#include "MBSPrepay.h"
#include "TBAUtil.h"

//****Functions:****//

/***
void    IOPO_Price_new (      double          *IO,                                    // IO price in
the lattice at the current period // double             *PO,                                   // PO
price in the lattice at the current period // double          **Zero, double          *Discount, //
One period discount factor at each node at period NbPer // double          *Discount1, // One period
discount factor without OAS // double          *Amort,                                 // Array of
main amortization factor in the interest rate lattice // double          Coupon, // Underlying
coupon (annualized) // int             CouponNb,                               // Coupon Number: =
deal_data->Term at the end of the tree // int             F,                                      //
Frequency of underlying payment as an integer // int             CpPaid, // Is the current period a
coupon payment date // int             NbPer,                                  // Current time
period // int		Lockout,                                // Lockout date of the SFS expressed
in number of nodes // int             ExEnd,                                  // Maturity of the IO
expressed in number of nodes // TREE_DATA       *tree_data,                             // Structure
of tree data // DEAL_DATA       *deal_data,                             // Structure of deal data
                        MBSPT_CashFlowStruct * cashflowstruct,
                        MBS_Prepay_Engine    *prepay_engine,
                        MBS_BKTree      *tree
);
**/

void IOPO_Price(
    double*    IO, /* IO price in the lattice at the current period */
    double*    PO, /* PO price in the lattice at the current period */
    double**   Zero,
    double*    Discount,   /* One period discount factor at each node at period NbPer */
    double*    Discount1,  /* One period discount factor without OAS */
    double*    pu,         /* Array of probabilities in the up state in the interest rate tree */
    double*    pd,         /* Array of probabilities in the down state in the interest rate tree */
    double*    p0,         /* Array of probabilities in the flat state in the interest rate tree */
    double*    Amort,      /* Array of main amortization factor in the interest rate lattice */
    double     Coupon,     /* Underlying coupon (annualized) */
    int        CouponNb,   /* Coupon Number: = deal_data->Term at the end of the tree */
    int        F,          /* Frequency of underlying payment as an integer */
    int        CpPaid,     /* Is the current period a coupon payment date */
    int        NbPer,      /* Current time period */
    int        Lockout,    /* Lockout date of the SFS expressed in number of nodes */
    int        ExEnd,      /* Maturity of the IO expressed in number of nodes */
    TREE_DATA* tree_data,  /* Structure of tree data */
    DEAL_DATA* deal_data); /* Structure of deal data */
#endif