/*******************************************************************************
*                Copyright (c) 1995 PARIBAS Capital Markets Group              *
********************************************************************************

    SYSTEM          :   GRF

    MODULE NAME     :   GRF_F_LNGSYMTAB

    PURPOSE         :   Symbol table for GRFN reserved words

    AUTHOR          :

    DATE            :

    VERSION         :

    DESCRIPTION     :   This module conatins a static variable which contains
                        detailed information on all the valid GRFN functions
                        and varibales including their data type  , number of
                        mandatory and optional arguments.


                        Financial function names:
                            swap  , fra  , level  , df...

                        Arithmetic function names:
                            exp  , log  , ...

                        Cell function names:
                            colmax  , ...

                        Pre-declared variables names:
                            va  , vb  , cpn  , bal  , ...

                        Deterministic variables:
                            i  , j

     NOTE           :   Strings like "ACT/365"  , are not represented explicitly
                        in the symbol table. They are checked for before
                        searching the symbol table. In interp_comll_name()
                        these strings  , are translated into enumerated types
                        and passed around in the "ival" field of a COMLL.
                        We know they have this particular type (and are not
                        integers) by the type of the comll (e.g. comll_basis)

    FUNCTIONS USED  :   grfn_interp_comll_name  ,  grfn_symbol_helpstr


    NOTE            :   The symbol "S" can take EITHER no arguments inwhich
                        case a domestic stock underlying is assumed OR one
                        string argument (similar to a function). This allows
                        us to use S BOTH as a single reserved symbol AND as
                        a function of the underlying (Nice!)

                        The reason why the parser CAN differentiate between
                        the two cases is due to S("name") appearing BEFORE S.
                        Consequently we detect the function first. This is
                        therefore a "side-effect".


********************************************************************************
*                           Amendment History                                  *
********************************************************************************

    AMENDED BY      :

    DATE            :

    VERSION         :

    REASON          :

    REQUEST/BUG NO  :

    DESCRIPTION

*******************************************************************************/

/*******************************************************************************
 *                               Import Include Files *
 *******************************************************************************/

#include "grf_h_all.h"

/*******************************************************************************
 *                           Macro Constants *
 *******************************************************************************/

#define GRFNSYMTABSZ (sizeof(grfn_symtab) / sizeof(GrfnSymbol))

/****************************************************************************
 *                                                                           *
 *   Key to arguments:                                                       *
 *                                                                           *
 *   d:  Date                            (deterministic)                     *
 *   s:  Start date                                                          *
 *   e:  End after start                        "                            *
 *   f:  End or number of full periods          "                            *
 *   E:  End or number of months                "                            *
 *   c:  Compounding                            "                            *
 *   b:  Basis                                  "                            *
 *   r:  Rec/Pay                                "                            *
 *   i:  Row reference                          "                            *
 *   j:  Column reference                       "                            *
 *   I:  Auxiliary row reference                "                            *
 *   J:  Auxiliary Col Reference                "                            *
 *   n:  Number                          (Non-deterministic)                 *
 *                                                                           *
 ****************************************************************************/

/****************************************************************************
 *                                                                           *
 *   The structure GrfnSymbol contains the following information:            *
 *                                                                           *
 *   (1) String      Name of symbol                                          *
 *   (2) CLLNameType Type of symbol                                          *
 *   (3) COMLLType   Enumerated type of symbol                               *
 *   (4) int         Number of deterministic (known) arguments               *
 *   (5) int         Number of non-deterministic (unknown) arguments         *
 *   (6) int         Number of optional arguments                            *
 *   (7) int         Extra information field                                 *
 *   (8) String      Help information the symbol                             *
 *                                                                           *
 *   The following GRFN functions take [optional] arguments:                 *
 *                                                                           *
 *       SWAP  , FRA  , LVL  , DF  , PAY  , PVRNG *
 *                                                                           *
 ****************************************************************************/

static GrfnSymbol grfn_symtab[] = {

    "Help on GRFN symbols",
    0,
    0,
    0,
    0,
    0,
    0,
    "To see arguments  , type  =grfn_doc(\"grfn\")",
    "Operators",
    0,
    0,
    0,
    0,
    0,
    0,
    "(Assign) :=  , (Arithmetic) +  ,-  ,*  ,/  ,^  , (Logical) <  ,<=  ,>  "
    ",>=  ,=  ,&&  ,||",
    "Indexing information",
    0,
    0,
    0,
    0,
    0,
    0,
    "All indexing in GRFN starts at 0 (NOT 1)",
    "Not a real symbol",
    0,
    0,
    0,
    0,
    0,
    0,
    "expr1 | expr2 : Evaluate & discard expr1; Evaluate & return expr2",
    "Repeat expresssion",
    0,
    0,
    0,
    0,
    0,
    0,
    "AM (expr) or AM: expr - Repeat expr from d[i] to d[i+1].",

    "IF",
    CLLNArithFunc,
    COMLL_A_IF,
    0,
    3,
    0,
    0,
    "IF(expr1  , expr2  , expr3): if expr1 > .5 then expr2 else expr3",

    "FLOOR",
    CLLNArithFunc,
    COMLL_A_FLOOR,
    0,
    2,
    0,
    0,
    "FLOOR(x  ,y): ((int)x/y) * y",
    "CEIL",
    CLLNArithFunc,
    COMLL_A_CEILING,
    0,
    2,
    0,
    0,
    "CEIL(x  ,y): ((int)(x/y+1/2))*y",
    "EXP",
    CLLNArithFunc,
    COMLL_A_EXP,
    0,
    1,
    0,
    0,
    "EXP(x)",
    "LOG",
    CLLNArithFunc,
    COMLL_A_LOG,
    0,
    1,
    0,
    0,
    "LOG(x)",
    "MAX",
    CLLNArithFunc,
    COMLL_A_MAX,
    0,
    2,
    0,
    0,
    "MAX(x  ,y)",
    "MIN",
    CLLNArithFunc,
    COMLL_A_MIN,
    0,
    2,
    0,
    0,
    "MIN(x  ,y)",
    "ABS",
    CLLNArithFunc,
    COMLL_A_ABS,
    0,
    1,
    0,
    0,
    "ABS(x)",
    "NORM",
    CLLNArithFunc,
    COMLL_A_NORM,
    0,
    1,
    0,
    0,
    "Cumulative Normal : norm(x)",
    "BS",
    CLLNArithFunc,
    COMLL_A_BLKSCHLS,
    0,
    5,
    0,
    0,
    "Black_Scholes : BS(f  ,k  ,vol  ,mat  ,c/p)",
    "BSSABR",
    CLLNArithFunc,
    COMLL_A_BSSABR,
    0,
    9,
    0,
    0,
    "Black_Scholes Sabr : BSSABR(f  ,k  ,betavol  ,alpha  ,beta  ,rho  ,mat  "
    ",VolType:0=Log  ,1=Norm  ,2=Beta  ,c/p)",
    "NUMDAYS",
    CLLNArithFunc,
    COMLL_A_NUMDAYS,
    0,
    11,
    0,
    0,
    "Number of days in a corridor : NUMDAYS(f  ,spot  ,bup  ,bdown  ,vol  "
    ",delay  ,mat  ,nstep  ,0->Brwn Brdg  , 1->No Smile  , 2->Linear Smile  , "
    "Par0  , Par1)",
    "BSN",
    CLLNArithFunc,
    COMLL_A_BLKSCHLSNRM,
    0,
    5,
    0,
    0,
    "Black_Scholes Normal : BSN(f  ,k  ,vol  ,mat  ,c/p)",

    "SPREADNUMER",
    CLLNArithFunc,
    COMLL_A_SPREADNUM,
    0,
    10,
    0,
    0,
    "Spread_numer  : SPREADNUMER(fwdx  ,fwdy  ,sigx  ,sigy  ,rho  ,a  ,b  ,k  "
    ",t  ,c/p)",

    "COLMAX",
    CLLNCellFunc,
    COMLL_C_COLMAX,
    3,
    0,
    0,
    0,
    "Column Maximum : COLMAX[column  , start_row  , end_row]",
    "COLMIN",
    CLLNCellFunc,
    COMLL_C_COLMIN,
    3,
    0,
    0,
    0,
    "Column Minimum : COLMIN[column  , start_row  , end_row]",
    "COLSUM",
    CLLNCellFunc,
    COMLL_C_COLSUM,
    3,
    0,
    0,
    0,
    "Column Sum : colsum[column  , start_row  , end_row]",
    "COLAVG",
    CLLNCellFunc,
    COMLL_C_COLAVG,
    3,
    0,
    0,
    0,
    "Column Average : COLAVG[column  , start_row  , end_row]",
    "COLSORT",
    CLLNCellFunc,
    COMLL_C_COLSORT,
    4,
    0,
    0,
    0,
    "Col Sort       : COLSORT[column  , start_row  , end_row  , index]",
    "ROWMAX",
    CLLNCellFunc,
    COMLL_C_ROWMAX,
    3,
    0,
    0,
    0,
    "Row Maximum : ROWMAX[row  , start_col  , end_col]",
    "ROWMIN",
    CLLNCellFunc,
    COMLL_C_ROWMIN,
    3,
    0,
    0,
    0,
    "Row Minimum : ROWMIN[row  , start_col  , end_col]",
    "ROWMAXIDX",
    CLLNCellFunc,
    COMLL_C_ROWMAXIDX,
    3,
    0,
    0,
    0,
    "Index of row Maximum : ROWMAXIDX[row  , start_col  , end_col]",
    "ROWMINIDX",
    CLLNCellFunc,
    COMLL_C_ROWMINIDX,
    3,
    0,
    0,
    0,
    "Index of row Minimum : ROWMINIDX[row  , start_col  , end_col]",
    "ROWSUM",
    CLLNCellFunc,
    COMLL_C_ROWSUM,
    3,
    0,
    0,
    0,
    "Row Sum       : ROWSUM[row  , start_col  , end_col]",
    "ROWAVG",
    CLLNCellFunc,
    COMLL_C_ROWAVG,
    3,
    0,
    0,
    0,
    "Row Average : ROWAVG[row  , start_col  , end_col]",
    "ROWSORT",
    CLLNCellFunc,
    COMLL_C_ROWSORT,
    4,
    0,
    0,
    0,
    "Row Sort       : ROWSORT[row  , start_col  , end_col  , index]: Index th "
    "largest value (index>=0) of a row",
    "ROWINTERP",
    CLLNCellFunc,
    COMLL_C_ROWINTERP,
    4,
    1,
    0,
    0,
    "ROWINTERP[row  , start_col  , end_col  , x_range (a[...])  , x] : "
    "Approximate y=f(x) by interpolation on a row ",

    "C",
    CLLNCellRef,
    COMLL_CONST,
    2,
    0,
    0,
    0,
    "C[col  ,row] : Main GRFN tableau cell contents",

    "CF",
    CLLNArithFunc,
    COMLL_CONSTF,
    0,
    2,
    0,
    0,
    "CF(col  ,row) : Main GRFN tableau cell contents with args evaluated",

    "PV",
    CLLNPvRef,
    COMLL_T_REF,
    1,
    0,
    0,
    0,
    "PV[j] : Discounted expected value of column j",

    // IMPLEMENTATION OF PVR IN GRFN
    "PVR",
    CLLNPvRRef,
    COMLL_T_RISKY_REF,
    1,
    0,
    0,
    0,
    "PVR[j] : Discounted expected value of column j with default risk (GRFNCR "
    "only)",

    // IMPLEMENTATION OF PVR IN GRFN
    "SOURCE",
    CLLNSourceAssign,
    COMLL_T_PVSOURCE_ASSIGN,
    1,
    1,
    0,
    0,
    "SOURCE(opt  , expr) : Assign expr as the new PVSource for this column; "
    "opt=0(cst)  ,1(disc) (GRFNCR only)",

    "CUM",
    CLLNCumRef,
    COMLL_T_CUMREF,
    1,
    0,
    0,
    0,
    "CUM[j] : cumulative capit. value of column j  , included current row",

    "CUR",
    CLLNCurRef,
    COMLL_T_CURREF,
    1,
    0,
    0,
    0,
    "CUR[j] : current value of column j  , current row",

    "INTERP",
    CLLNAuxFunc,
    COMLL_X_INTERP,
    2,
    1,
    0,
    0,
    "INTERP(xrange (a[...])  , yrange (a[...])  , x) : Approximate x by "
    "interpolation using Auxiliary Ranges",

    "PVINTERP",
    CLLNAuxFunc,
    COMLL_X_PVINTERP,
    3,
    1,
    0,
    0,
    "PVINTERP(AuxIndex  ,FirstCol  ,ColSkip  ,x) : Interp y(x) from A[AuxIndex "
    " ,i] and PV[FirsCol+j*ColSkip]",

    "A",
    CLLNAuxRef,
    COMLL_REAL,
    2,
    0,
    0,
    0,
    "A[col  ,row] : Auxiliary cell contents",

    "D",
    CLLNDateRef,
    COMLL_REAL,
    1,
    0,
    0,
    0,
    "D[i] : Date of ith financial event",

    "YIELD",
    CLLNFinFunc,
    COMLL_F_YIELD,
    6,
    0,
    0,
    0,
    "Bond : yield(start  ,end/nfp  ,compd  ,basis  ,cpn  ,clean)",

    "CLEAN",
    CLLNFinFunc,
    COMLL_F_CLEAN,
    6,
    0,
    0,
    0,
    "Bond : clean(start  ,end/nfp  ,compd  ,basis  ,cpn  ,yield)",

    "CVG",
    CLLNFinFunc,
    COMLL_F_CVG,
    3,
    0,
    0,
    0,
    "CVG(date1  , date2  , basiscode) : Coverage",

    /*  "ACCINT"  ,   CLLNFinFunc  ,    COMLL_F_ACCINT  , 4  ,0  ,0  ,0  ,
                    "ACCINT(daterng  , cfrng  , accdate  , basiscode)"  ,  */

    "SWAP",
    CLLNFinFunc,
    COMLL_F_SWAP,
    4,
    0,
    3,
    0,
    "SWAP(start  , end/nfp  , comp  , basis  , [und]  , [ratecode]  , "
    "[floatcode]) : Swap Rate",
    "LVL",
    CLLNFinFunc,
    COMLL_F_LEVEL,
    4,
    0,
    1,
    0,
    "LVL(start  ,end/nfp  , comp  , basis  ,[und]) : Level Payment",
    "FRA",
    CLLNFinFunc,
    COMLL_F_FRA,
    3,
    0,
    2,
    0,
    "FRA(start  , end/nummonth  , basis  , [und]  , [floatcode]) : Moneymarket "
    "Short Rate",
    "DF",
    CLLNFinFunc,
    COMLL_F_DF,
    2,
    0,
    2,
    0,
    "DF(start  , end  , [und])",
    "PAY",
    CLLNFinFunc,
    COMLL_F_PAY,
    1,
    1,
    0,
    0,
    "PAY(amount  ,date) : Future Cashflow",
    "PAYN",
    CLLNFinFunc,
    COMLL_F_PAYN,
    1,
    1,
    0,
    0,
    "PAYN(amount  ,date) : Future Cashflow (multiple calls allowed)",
    "PVRNG",
    CLLNFinFunc,
    COMLL_F_PVRANGE,
    3,
    0,
    1,
    0,
    "PVRNG(daterange  , cfrange  , startdate <after this>  ,[und])",

    "CMS",
    CLLNFinFunc,
    COMLL_F_CMS,
    7,
    0,
    2,
    0,
    "CMS(start  , end/nfp  , comp  , basis  , vol  , delay/date  ,logn/nrm  , "
    "[und]  , [floatcode]) : CMS Rate",

    "PAYOFF",
    CLLNFinFunc,
    COMLL_F_PAYOFF,
    2,
    0,
    1,
    0,
    "PAYOFF(date  , product  , [und]) : payoff of a pre-initialized product",

    /*  "DRTY"  ,     CLLNFinFunc  ,    COMLL_F_PVRANGE  ,4  ,0  ,0  ,0  ,
            "DRTY(daterange  , cfrange  , startdate(>)  , enddate(<=))"  , */

    "PR",
    CLLNArithFunc,
    COMLL_PRINT,
    1,
    0,
    0,
    0,
    "PR(x) : Print x to grf_sam.dat",
    "HIST",
    CLLNUtil,
    COMLL_HIST,
    1,
    0,
    0,
    0,
    "HIST(\"label\") - Create histogram for labelled cell (expr | "
    "hist\"label\")",

    "I",
    CLLNI,
    COMLL_REAL,
    0,
    0,
    0,
    0,
    "I: Index of current row",
    "J",
    CLLNJ,
    COMLL_REAL,
    0,
    0,
    0,
    0,
    "J: Index of current col",
    "NOW",
    CLLNNow,
    COMLL_REAL,
    0,
    0,
    0,
    0,
    "NOW: Date of this event",
    "NEXT",
    CLLNNxt,
    COMLL_REAL,
    0,
    0,
    0,
    0,
    "NEXT: Date of next event",
    "PREV",
    CLLNPrv,
    COMLL_REAL,
    0,
    0,
    0,
    0,
    "PREV: Date of previous event",

    "R",
    CLLNSample,
    COMLL_S_SHORT_RATE,
    0,
    0,
    0,
    0,
    "R: Instantaneous short rate",
    "PHI",
    CLLNSample,
    COMLL_S_PHI,
    0,
    0,
    0,
    0,
    "PHI: Cheyette Cumulative Volatility",
    "SIGMA",
    CLLNSample,
    COMLL_S_SIGMA,
    0,
    0,
    0,
    0,
    "SIGMA: The volatility (stochastic ?) used in the diffusion",
    "STATEVAR",
    CLLNSample,
    COMLL_S_STATEVAR,
    0,
    0,
    0,
    0,
    "STATEVAR: The value of the statevar discretised in the model",
    "NUMERAIRE",
    CLLNSample,
    COMLL_S_NUMERAIRE,
    0,
    0,
    0,
    0,
    "NUMERAIRE: The value of the numeraire discretised in the model",
    "BT",
    CLLNSample,
    COMLL_S_MMACCOUNT,
    1,
    0,
    0,
    0,
    "BT(\"label\"): The money market account identified by 'label'",
    "S",
    CLLNSample,
    COMLL_S_UNDERLYING,
    1,
    0,
    0,
    0,
    "S(\"label\"): Equity/FX Underlying identified by 'label'",
    "S",
    CLLNSample,
    COMLL_S_STOCK,
    0,
    0,
    0,
    0,
    "S : Equity (Domestic Stock)",
    "CPN",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_CPN,
    "CPN: variable",
    "BAL",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_BAL,
    "BAL: variable",
    "VA",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_VA,
    "VA:  variable",
    "VB",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_VB,
    "VB:  variable",
    "VC",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_VC,
    "VC:  variable",
    "VD",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_VD,
    "VD:  variable",
    "VE",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_VE,
    "VE:  variable",
    "VF",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_VF,
    "VF:  variable",
    "VG",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_VG,
    "VG:  variable",
    "VH",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_VH,
    "VH:  variable",
    "V0",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V0,
    "V0:  variable",
    "V1",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V1,
    "V1:  variable",
    "V2",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V2,
    "V2:  variable",
    "V3",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V3,
    "V3:  variable",
    "V4",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V4,
    "V4:  variable",
    "V5",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V5,
    "V5:  variable",
    "V6",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V6,
    "V6:  variable",
    "V7",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V7,
    "V7:  variable",
    "V8",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V8,
    "V8:  variable",
    "V9",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V9,
    "V9:  variable",
    "V10",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V10,
    "V10:  variable",
    "V11",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V11,
    "V11:  variable",
    "V12",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V12,
    "V12:  variable",
    "V13",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V13,
    "V13:  variable",
    "V14",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V14,
    "V14:  variable",
    "V15",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V15,
    "V15:  variable",
    "V16",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V16,
    "V16:  variable",
    "V17",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V17,
    "V17:  variable",
    "V18",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V18,
    "V18:  variable",
    "V19",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V19,
    "V19:  variable",
    "V20",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V20,
    "V20:  variable",
    "V21",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V21,
    "V21:  variable",
    "V22",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V22,
    "V22:  variable",
    "V23",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V23,
    "V23:  variable",
    "V24",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V24,
    "V24:  variable",
    "V25",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V25,
    "V25:  variable",
    "V26",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V26,
    "V26:  variable",
    "V27",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V27,
    "V27:  variable",
    "V28",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V28,
    "V28:  variable",
    "V29",
    CLLNVar,
    COMLL_VAR,
    0,
    0,
    0,
    GRFN_SS_V29,
    "V29:  variable",
};

/****************************************************************************
 *                                                                           *
 *   Function    :   grfn_symbol_helpstr                                     *
 *                                                                           *
 *   Purpose     :   Returns a string containing (on-line) help information  *
 *                   from the structure grfn_symtab.                         *
 *                                                                           * *
 ****************************************************************************/

String grfn_symbol_helpstr(int i) {
  if (i < 0 || i >= GRFNSYMTABSZ)
    return ("Index is out of range.");
  else
    return (grfn_symtab[i].helpstr);
}

long grfn_symbol_list_size() {
  return (long)(sizeof(grfn_symtab) / sizeof(GrfnSymbol));
}

/*******************************************************************************
 *                                                                              *
 *   Function:   grfn_interp_comll_name *
 *                                                                              *
 *   Purpose:    Checks name against the GRFN symbol table * If name is found
 *then the corresponding contents of the        * symbol table is copied to gs
 *, otherwise an error message is    * returned. *
 *                                                                              *
 *   Note:       A further check against the declaration of user-defined *
 *               variables may be required *
 *                                                                              *
 *   Called by:  grfn_interp_cell_func * grfn_interp_func * grfn_interp_name *
 *                                                                              *
 *******************************************************************************/

Err grfn_interp_comll_name(String name, GrfnSymbol *gs)

{
  int i;

  for (i = 0; i < GRFNSYMTABSZ; i++) {

    /*  If name is a valid symbol  , then copy contents of symbol table
        to gs  , else return an error message
    */

    if (!strcmp(grfn_symtab[i].name, name)) {
      memcpy(gs, &grfn_symtab[i], sizeof(GrfnSymbol));
      return NULL;
    }
  }

  return serror("%s %s", GRERR_UNKNOWN_SYMBOL, name);
}

/****************************
 *   NOT YET IMPLEMENTED     *
 ****************************/

/*  Check that the gd->grng does not contain any duplicated
    or otherwise illegal names.  */

Err grfn_check_ranges(GrfnDeal *gd) { return NULL; }
