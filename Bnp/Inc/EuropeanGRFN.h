/*******************************************************************************
 *
 * FUNCTION      : EuropeanGRFN.h
 *
 * PURPOSE       :
 *
 * DESCRIPTION   :
 *
 * CALLS         :
 *
 *******************************************************************************/
#ifndef __EUROPEANGRFN_H
#define __EUROPEANGRFN_H

#include "CopulaGenericPayoff.h"
#include "grf_h_mdlcomm.h"
#include "srt_h_all.h"

/* ---------------------------------------------------------------------------

        Definition of the GRFNEuropeanCorrMat

  ---------------------------------------------------------------------------- */
typedef struct
{
    long lNbPair;     /* Nb of pair Name1/Name2 */
    long lNbDate;     /* Nb of date */
    long lStringSize; /* ize max of the string */

    char** cName1;
    char** cName2;

    double*  dDate;
    double** dCorrMat;

} GRFNEuropeanCorrMat;

/* ---------------------------------------------------------------------------

        Definition of the GRFNEuropean_Payoff

  ---------------------------------------------------------------------------- */
typedef struct
{
    GRFNCOMMSTRUCT global;
    FIRSTMktAtT*   local;
    int*           iIndex;

} GRFNEuropean_Payoff;

/* ---------------------------------------------------------------------------

        Definition of the structure SrtSABRUnd

  ---------------------------------------------------------------------------- */
typedef struct
{
    long    lNbDate; /* Nb of Dates in the TermStruct */
    double* dTime;   /* Dates or time vector */
    double* dFwd;    /* forward vector */
    double* dSigma;  /* Sigma sabr vector */
    double* dAlpha;
    double* dBeta;
    double* dRho;

} SrtSABRUnd;

/* ---------------------------------------------------------------------------

        Definition of the SABRUnd

  ---------------------------------------------------------------------------- */
char* SrtInitSABRUnd(
    char*   cUndName,
    char*   cYieldCurveName,
    long    lNbDate,
    double* dTime,
    double* dFwd,
    double* dSigma,
    double* dAlpha,
    double* dBeta,
    double* dRho);

/* ----------------------------------------------------------------------------------------------------------

        free_GRFNEuropeanCorrMat

  -----------------------------------------------------------------------------------------------------------
*/
char* free_GRFNEuropeanCorrMat(GRFNEuropeanCorrMat* ptr);

/* ----------------------------------------------------------------------------------------------------------

        PayOff_GRFNEuropean

  -----------------------------------------------------------------------------------------------------------
*/

Err PayOff_GRFNEuropean(/* Underlying value */
                        int     iNbUnderlying,
                        double* UnderlyingValue,

                        /* For Payoff description */
                        int   iNbProduct,
                        void* PayOffParam,

                        /* OutPut */
                        double* dPV);

/* ----------------------------------------------------------------------------------------------------------

        SrtGRFNEuropean
        Pricing of a european GRFN

  -----------------------------------------------------------------------------------------------------------
*/

Err SrtGRFNEuropean(/* Grfn inputs */
                    int      numeventdates,
                    long*    eventdates,
                    long     tableauRows,
                    long     tableauCols,
                    char***  tableauStrings,
                    int**    tableauMask,
                    long     auxWidth,
                    long*    auxLen,
                    double** aux,

                    /* Copula Pricing Parameters */
                    CopulaType   iCopulaType,
                    long         lNumPaths,
                    long         lNbPoints,
                    SrtMCSamType iMCType,
                    ModelType    iModelType,
                    double       dNStdBMMCalib,

                    /* Dom Market */
                    char* cDomYcName,

                    /* Correlation Matrix */
                    GRFNEuropeanCorrMat* ptrCorrMat,

                    /* OutPuts */
                    double*  dPv,
                    double** grfn_ss /* For GRFN_celll */);

#endif