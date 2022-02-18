/*--------------------------------------------------------------
        FILE: Product.h
        PURPOSE: Generic product description
        AUTHOR: Dimitri Mayevski
        DATE: 21/11/2002
  --------------------------------------------------------------*/

#ifndef __PRODUCT_H__
#define __PRODUCT_H__

#include "swp_h_utils.h"

// ====================== Structures ========================

enum EMktDataType
{
    MKT_IRSTANDARD,  // dfs only
    MKT_GRFN,        // spec_desc points to SrtSample
    MKT_SABR         // spec_desc points to SSabrMktDesc
};

typedef struct _SMktData
{
    double** dfs;   // dfs by currency
    int      type;  // EMktDataType
    void*    spec_desc;
} SMktData;

typedef struct _SSabrMktDesc
{
    double** sigbeta;  // SABR sigma-beta by ccy
    double** alpha;    // SABR alpha by ccy
    double** beta;     // SABR beta by ccy
    double** rho;      // SABR rho by ccy
    int      type;     // EMktDataType
    void*    spec_desc;
} SSabrMktDesc;

typedef struct _SProductDesc SProductDesc;

typedef Err (*PayoffFuncPtr)(
    SProductDesc* product,   // this
    int           idx,       // ex date index
    double        time,      // may be different from ex[idx] for am[idx] = 1
    SMktData*     mkt_data,  // market data (dfs)
    double*       pv         // ninst sized vector with PVs (for backward calculation)
);                           // resulting PVs are stored in pv as well

enum EProductType
{
    PRODUCT_GRFN,           // GRFN table
    PRODUCT_SWAPTIONS,      // Collection of swaptions, e.g. used for calibration
    PRODUCT_SMM_SWAPTIONS,  // Collection of swaptions for SMM
    PRODUCT_CIF             // Callable Inverse Floater
};

struct _SProductDesc
{
    int       nex, nccy, ninst;
    double*   ex;     // exercise times
    long*     ex_d;   // exercise dates
    int*      am;     // american exercise flags
    int**     nmat;   // num of maturities per ccy per exercise
    double*** mat;    // maturities per ccy per exercise
    long***   mat_d;  // maturities dates per ccy per exercise

    int**   nvol;       // num of vols per ccy per exercise
    long*** vol_start;  // vol start dates per ccy per exercise
    long*** vol_end;    // vol end dates per ccy per exercise

    PayoffFuncPtr Payoff;  // Payoff function

    int   type;  // EProductType
    void* spec_desc;
};

// ====================== Functions ========================

Err ProductDesc_InitGRFN(
    SProductDesc* product,
    char*         und,
    int           nevd,
    long*         evd,
    long          ntabrows,
    long          ntabcols,
    char***       tabstrs,
    int**         tabmask,
    long          auxwidth,
    long*         auxlen,
    double**      aux);

Err ProductDesc_FreeGRFN(SProductDesc* product);

Err ProductDesc_InitSwaptions(SProductDesc* product, SCashFlows* g);

#endif  // #ifndef __PRODUCT_H__