/****************************************************************
 * Module:    Boots
 * File:    bootsmain
 * Function:links all the element functions    
 * Author:    CA
 * Revision:
 *****************************************************************/
#include "drlstd.h"            /* platform compatibility */
#include <math.h>
#include <string.h>
#ifndef    _WIN32
#include <errno.h>
#endif

#include "drlerr.h"
#include "drlio.h"
#include "drlstr.h"
#include "drltime.h"
#include "drllinsm.h"

#define    _Boots_SOURCE
#include "boots.h"


/*f-------------------------------------------------------------
 * Main Function.
 * Calculates KOV matrix for each specified expiry.                
 * Stores the corresponding values in kovTens. 
 */
int 
BootsMain(
    BootsData* that,  /*(I/O) Model parameters*/
    SwMat* sMat,      /*(I) SwaptionMatrix*/
    Target* target,   /*(I) Target correlations*/
    AddInput* addI,   /*(I) Additional Input Parameters*/
    int nbExpTM,      /*(I) Number of expiries to match*/
    double* swExp,      /*(I) Vector of expiries*/
    int nbMatTM,      /*(I) Number of maturities to match*/
    double* swMat       /*(I) Vector of maturities*/
)
{
    int i, j, iter, iterTest, k, l;
    double minimum;
    int errCode = 1;
    int trial = 0;
    int expiryIndex, e, s, di, dj, dk;
    double sum, start, end, K_inte, res;
    int size;
    double ev;

    size = that->nbMatMax;
    that->nbSwpMat = nbMatTM;
    that->swpMat = swMat;

    /* Test: is first expiry = to first time step ?*/
    /* Induction step*/
    if((swExp[0]>that->dt)||(swExp[0]<that->dt)){
        goto done;
    }
 
    /* Precompute all the fixed matrices */
    if(BootsComputeLambda(that) !=0)
        goto done;

    /* FIRST EXPIRY */
    /* First step at expiry == 0 */
    /* Compute all intermediary matrices at that date */
    if(BootsPrecomputeNewKov(that, 
                            sMat,
                            addI,
                            0,
                            size,
                            swMat)!=0)
        goto done;
    /* Update expIndex*/
    that->expIndex = 0;
    /*Test target correlation parameters for this expiry*/
	if(BootsCheckTarget(that, target)!=SUCCESS)
        goto done;
	/*Build correlation matrix*/
	/* CALL THE OPTIMIZER FUNCTION*/ 
    /* Find the Optimal Thetas set of parameters to build the R matrix*/
    switch(target[0].method){
        case 0:
        /* Exponent: 1 param*/
        if(frprmn1(that,
            target,
            addI->initialThetas,
            1,
            1e-5,
            &iter,
            &minimum)!=0)
            goto done;
        break;

        case 1:
        /* Bochner1: 4 params*/
        if(frprmn1(that,
            target,
            addI->initialThetas,
            4,
            1e-5,
            &iter,
            &minimum)!=0)
            goto done;
        break;

        case 2:
        /* Bochner2: 8 params*/
        if(frprmn1(that,
            target,
            addI->initialThetas,
            8,
            1e-5,
            &iter,
            &minimum)!=0)
        goto done;
        break;

        case 3:
        /* Bochner3: 12 params*/
        if(frprmn1(that,
            target,
            addI->initialThetas,
            12,
            1e-5,
            &iter,
            &minimum)!=0)
        goto done;
        break;
        
        default:
        /* Default not implemented*/
        goto done;

    }

    /* Find new covariance matrix */
    if(BootsGetNewKovMatrix(that, size)!=0)
        goto done;
    /* Update Covariance Tensor, with covariance matrix*/
    if(BootsCopyResults(that, 0)!=0)
        goto done;
    /* END FIRST EXPIRY*/

    j = 1;
    /* LOOP ON EXPIRY*/
    for( i = 1; i < that->nbExpMax; i++){
        /* Update the expiryIndex*/
        that->expIndex = i;
        ev = ((i*that->dt + that->dt) - swExp[j]);

        /* Compute all intermediary matrices at corresponding Expiry*/
        if(ev == 0)
        {
            if(BootsPrecomputeNewKov(
            that,
            sMat,
            addI,
            i,
            size,
            swMat)!=0)
            goto done;

			if(BootsCheckTarget(that, target)!=SUCCESS)
	        goto done;

			/* Call the optimizer function*/
            /* Find the Optimal Thetas set of parameters to build
             * the R matrix and match target correls*/
            /* Do it only when it is required*/
             switch(target[0].method){
             case 0:
             /* Exponent: 1 param*/
            if(frprmn1(that,
                target,
                addI->initialThetas,
                1,
                1e-5,
                &iter,
                &minimum)!=0)
                goto done;
            break;

            case 1:
            /* Bochner1: 4 params*/
            if(frprmn1(that,
                target,
                addI->initialThetas,
                4,
                1e-5,
                &iter,
                &minimum)!=0)
                goto done;
            break;

            case 2:
            /* Bochner2: 8 params*/
            if(frprmn1(that,
                target,
                addI->initialThetas,
                8,
                1e-5,
                &iter,
                &minimum)!=0)
            goto done;
            break;

            case 3:
            /* Bochner3: 12 params*/
            if(frprmn1(that,
                target,
                addI->initialThetas,
                12,
                1e-5,
                &iter,
                &minimum)!=0)
            goto done;
            break;
        
            default:
            /* Default not implemented*/
            goto done;  
        }

            if(BootsGetNewKovMatrix(that, size)!=0)
                goto done;    
            if(BootsCopyResults(that, i)!=0)
                goto done;
            j++;
        }
        else
        {
            BootsPrecomputeFlat(that, i, 1);
        }
    }
 
    errCode = 0;
done:

    if (errCode != 0) {
        DrlErrMsg("BootsMain: failed \n");
    }

    return(errCode);

}


/*f-------------------------------------------------------------
 * Another main function.
 * No calibration required.
 * Instead a R matrix is given, which will be used at every time step.
 * Calculates KOV matrix for each expiry.                
 * Stores the corresponding values in kovTens. 
 */
int 
BootsCalibNoOpt(
    BootsData* that,   /*(I/O) Model parameters */
    SwMat* sMat,       /*(I) Initial SwaptionMatrix*/
    double** corr,     /*(I) The correlation matrix (numMat*numMat)*/
    AddInput* addI,    /*(I) Additional Input Parameters*/
    int numExp,        /*(I) Number of expiries to match*/
    double* swExp,     /*(I) Where the routine "BootsPrecomputeNewKov"
                        *    will be used*/
                       /*    First expiry should be first time step*/
    int numMat,        /*(I) Number of tenors to match*/
    double* swMat      /*(I) The values of the tenors*/
)
{
    int i, j, iter, iterTest, k, l;
    double minimum;
    int errCode = 1;
    int trial = 0;
    int expiryIndex, e, s, di, dj, dk;
    double sum, start, end, K_inte, res;
    int size;
    double ev;

    size = that->nbMatMax;
    that->nbSwpMat = numMat;
    that->swpMat = swMat;

    /* Test: is first expiry = to first time step ?*/
    /* Induction step*/
    if((swExp[0]>that->dt)||(swExp[0]<that->dt)){
        goto done;
    }
 
    /* Precompute all the fixed matrices */
    if(BootsComputeLambda(that) !=0)
        goto done;

    /* Copy correlation matrix into that->CorrMat*/
    DrlMatrixCopySM(&that->corrMat, corr, numMat, numMat);

    /* FIRST EXPIRY */
    /* First step at expiry == 0 */
    /* Compute all intermediary matrices at that date */
    if(BootsPrecomputeNewKov(that, 
                            sMat,
                            addI,
                            0,
                            size,
                            swMat)!=0)
        goto done;
    /* Update expIndex*/
    that->expIndex = 0;
    /* Find new covariance matrix */
    if(BootsGetNewKovMatrix(that, size)!=0)
        goto done;
    /* Update Covariance Tensor, with covariance matrix*/
    if(BootsCopyResults(that, 0)!=0)
        goto done;
    /* END FIRST EXPIRY*/

    j = 1;
    /* LOOP ON EXPIRY*/
    for( i = 1; i < that->nbExpMax; i++){
        /* Update the expiryIndex*/
        that->expIndex = i;
        ev = ((i*that->dt + that->dt) - swExp[j]);

        /* Compute all intermediary matrices at corresponding Expiry*/
        if(ev == 0)
        {
            if(BootsPrecomputeNewKov(
            that,
            sMat,
            addI,
            i,
            size,
            swMat)!=0)
            goto done;

            if(BootsGetNewKovMatrix(that, size)!=0)
                goto done;    
            if(BootsCopyResults(that, i)!=0)
                goto done;
            j++;
        }
        else
        {
            BootsPrecomputeFlat(that, i, 1);
        }
    }
    
    errCode = 0;
done:
   
    if (errCode != 0) {
        DrlErrMsg("BootsMain: failed \n");
    }

    return(errCode);

}
