
/****************************************************************
 * Module:   BOOTS
 * File:     bootsmeth
 * Function: low level calibration routines
 * Author:   CA
 * Revision:    
 *****************************************************************/
#include "drlstd.h"            /* platform compatibility */
#include <math.h>

#include "drlerr.h"
#include "drltime.h"
#include "drlsort.h"
#include "drlio.h"
#include "drllinsm.h"

#define  _Boots_SOURCE
#include "boots.h"

#define    __DEBUG__
#undef    __DEBUG__

/*f-------------------------------------------------------------
 * Build a correlation matrix from a vector Thetas.
 * Vector Thetas[1,n].
 */
int 
BootsGetNewCorrMatrix(
    BootsData* that,      /*(I/O) model parameters */
    Target* targetCorr,   /*(I) Target correlations input*/
    double* Thetas        /*(I) Thetas parameters*/
)
{
    /*Thetas in "vector" format: starts at 1 finish at MAT*/
    int errCode = 1;
    int i,j;
    double sum = 0.0;


    /*Switch according to method */
    switch(targetCorr[that->expIndex].method)
    {
        /*Multiplicative method */
        /*Parameter stored in Thetas[1] */
        case 0:
        for(i = 0; i < that->nbSwpMat; i++){
            for(j = 0; j < that->nbSwpMat; j++){
                CORRM[i][j] = exp(-fabs(Thetas[1])*
                              fabs(that->swpMat[i]-that->swpMat[j]));
            }
        }
        break;

        /*Bochner parametrisation */
        /*Fourier transform of a Lorentzian (for i or log(i) with law +) */
        /*Degrees of freedom: 4*/
        case 1:
            for(i = 0; i < that->nbSwpMat; i++){
                for(j = 0; j < that->nbSwpMat; j++){
                    that->corrMat[i][j] = (fabs(Thetas[3])*
                    exp(-fabs(Thetas[1])*
                    fabs(that->swpMat[i]-that->swpMat[j]))+ 
                    fabs(Thetas[4])*
                    exp(-fabs(Thetas[2])*
                    fabs(log(that->swpMat[i])-log(that->swpMat[j]))))/
                    (double)(fabs(Thetas[3])+fabs(Thetas[4]));
                }
            }
        break;

        /*More Bochner */
        /*Sum of the Fourier transforms of a Lorentzian */
        /*and a Gaussian (for i or log(i) with law +)   */
        /*Degrees of freedom: 8*/
        case 2:
            for(i = 0; i < that->nbSwpMat; i++){
                for(j = 0; j < that->nbSwpMat; j++){
                    that->corrMat[i][j] = 
                    fabs(Thetas[5])*exp(-fabs(Thetas[1])*
                    fabs(that->swpMat[i]-that->swpMat[j]))+
                    fabs(Thetas[6])*exp(-fabs(Thetas[2])*
                    fabs(log(that->swpMat[i])-log(that->swpMat[j])));

                    that->corrMat[i][j] +=
                    fabs(Thetas[7])*exp(-fabs(Thetas[3])*
                    SQR(fabs(that->swpMat[i]-that->swpMat[j])));
                    
                    that->corrMat[i][j] +=
                    fabs(Thetas[8])*exp(-fabs(Thetas[4])*
                    SQR(fabs(log(that->swpMat[i])-log(that->swpMat[j]))));

                    that->corrMat[i][j] /= 
                    (double) (fabs(Thetas[5])+fabs(Thetas[6])+fabs(Thetas[7])+fabs(Thetas[8]));
                }
            }
        break;

        /* More Bochner 
         * Sum of the Fourier transforms of a Lorentzian, a Gaussian,
         * and an exponential (for i or log(i) with law +) 
         * Degrees of freedom: 12 */
        case 3:
            for(i = 0; i < that->nbSwpMat; i++){
                for(j = 0; j < that->nbSwpMat; j++){
                    sum = 0.0;
                    that->corrMat[i][j] = 
                    fabs(Thetas[7])*exp(-fabs(Thetas[1])*
                    fabs(that->swpMat[i]-that->swpMat[j]))+
                    fabs(Thetas[8])*exp(-fabs(Thetas[2])*
                    fabs(log(that->swpMat[i])-log(that->swpMat[j])));
                    sum += (fabs(Thetas[7])+fabs(Thetas[8]));

                    that->corrMat[i][j] +=
                    fabs(Thetas[9])*exp(-fabs(Thetas[3])*
                    SQR(fabs(that->swpMat[i]-that->swpMat[j])));

                    that->corrMat[i][j] +=
                    fabs(Thetas[10])*exp(-fabs(Thetas[4])*
                    SQR(fabs(log(that->swpMat[i])-log(that->swpMat[j]))));
                    sum += (fabs(Thetas[9])+fabs(Thetas[10]));

                    that->corrMat[i][j] +=
                    fabs(Thetas[11])/(1+fabs(Thetas[5])*
                    SQR(fabs(that->swpMat[i]-that->swpMat[j])));
                    that->corrMat[i][j] +=
                    fabs(Thetas[12])/(1+fabs(Thetas[6])*
                    SQR(fabs(log(that->swpMat[i])-log(that->swpMat[j]))));
                    sum += (fabs(Thetas[11])+fabs(Thetas[12]));

                    that->corrMat[i][j] /= sum;
                }
            }
        break;

        default:
        /*Not implemented*/
        goto done;


    }

    errCode = 0;
done:
    if (errCode != 0) {
        DrlErrMsg("BootsGetNewCorrMatrix: failed\n");
    }
    return(errCode);
}

/*f-------------------------------------------------------------
 * 
 * Given a new Correlation Matrix, recalculates the covariance 
 * Matrix.
 *
 */
int 
BootsGetNewKovMatrix(
    BootsData* that,   /*(I/O) Model parameters */
    int numSwap        /*(I) Number of maturities to be matched
                        *    in the swaption line*/
)
{
    int errCode = 0;

    DrlMatrixProductSM(&that->productMat,
                    ALPHAM,
                    numSwap, 
                    numSwap,
                    CORRM,
                    numSwap,
                    numSwap
                    );

    DrlMatrixProductSM(&that->productMat2,
                    that->productMat,
                    numSwap,
                    numSwap,
                    ALPHAM,
                    numSwap,
                    numSwap
                    );

    DrlMatrixProductSM(&that->productMat3,
                    that->productMat2,
                    numSwap,
                    numSwap,
                    INVLM,
                    numSwap,
                    that->nbMatMax
                    );

    DrlMatrixProductSM(&KOVM,
                    INVLMT,
                    that->nbMatMax,
                    numSwap,
                    that->productMat3,
                    numSwap,
                    that->nbMatMax
                    );

    errCode = 0;

    if (errCode != 0) {
        DrlErrMsg("BootsGetNewKovMatrix: failed\n");
    }
    return(errCode);
}

/*f-------------------------------------------------------------
 * Copy the Kij matrix into the covariance tensor.
 * At a given expiry index and with a given strategy.
 * Either K(t+1, i, j) = K(t, i, j)
 * Or K(t+1, i, j) = K(t, i-1, j-1)
 * i and j being absolute indexes.
 */
int 
BootsPrecomputeFlat(
    BootsData* that,    /*(I/O) Model parameters */
    int idxS ,          /*(I) Index of the covariance Tensor*/
    int stationary      /*(I) 0 = relative, 1 = absolute reset*/
)
{
    int i, j;
    int errCode = 1;
    
    
     switch(stationary)
    {   
        /*Relative*/
        case 0:
            for(i = 0; i < that->nbMatMax; i++){
                for(j = 0; j < that->nbMatMax; j++){
                   that->kovTens[i][j][idxS] = that->kovTens[i][j][idxS-1];
                }
            }
        break;

        /*Absolute*/
        case 1:
            for(i = 0; i < that->nbMatMax; i++){
                for(j = 0; j < that->nbMatMax; j++){
                   that->kovTens[i][j][idxS] = that->kovTens[i+1][j+1][idxS-1];
                }
            }
        break;
    
        default:
        goto done;

    }

    errCode = 0;
done:

    if (errCode != 0) {
        DrlErrMsg("BootsUpdateMatrices: failed\n");
    }
    return(errCode);
}

/*f-------------------------------------------------------------
 * Update all Matrices at each expiry index on the expiry timeline.
 * Recalculates the alpha Matrix, Eq.(19), in doc.
 * Refill the lambaMatrix
 * (new Vol Fra/Vol Swap coeffs for the corresponding expiry).
 * Calculate the inverse of the LambdaMatrix
 * and transposed value of the inverse.
 *
 * (!beware shift of one in the indexes with regard to newTimeLine)
 * [newTimeLine[0] = 0]
 */
int 
BootsPrecomputeNewKov(
    BootsData* that,   /*(I/O) Model parameters */
    SwMat* sMat,       /*(I) Input SwaptionMatrix*/
    AddInput* addI,    /*(I) Additional Input Parameters*/
    int expiryIndex,   /*(I) ExpiryIndex for update*/
    int nbSwp,         /*(I) Number of points to use from SwaptionMatrix*/
    double* swpMat     /*(I) Swaption maturities to be used 
                        *    (maturity indexes starting at 0)*/
)
{
    int i, j, di, dj, dk, i1;
    /* start, end */
    int s,e;
    /* start and end of the swaption in years, various indexes */
    double start, end, t1; 
    /* K_inte is an intermediary for the calculation
     * of the second term in Eq.(19)
     */
    double sum = 0.0, K_inte = 0.0;
    int errCode = 1;
    int idx;


    /* Copy nbSwp into the structure variable
     * that->swpMat will now point on swpMat*/
    that->nbSwpMat = nbSwp;
    that->swpMat = swpMat;
    

    /* Find the maturity index*/
#undef  IDS
#define IDS(i, idx) DrlDoubleArrayClosestIdx((*sMat).maturity,that->nbMatMax,swpMat[i],(idx))

    /* Fill alpha matrix in the case expiryIndex = 0 */ 
    if(expiryIndex == 0){
        for(i = 0; i < nbSwp; i++){
            IDS(i, &idx);
            ALPHAM[i][i] = sqrt((sMat->expiry[expiryIndex])*
                           SQR(sMat->swaptionMatrix[expiryIndex][idx]));
        }
    }
    else{
    /* Fill alpha matrix otherwise */
        for(i = 0; i < nbSwp; i++){
            sum = 0.0;
            IDS(i, &idx);
            ALPHAM[i][i] = (sMat->expiry[expiryIndex])*
                            SQR(sMat->swaptionMatrix[expiryIndex][idx]);

            start = sMat->expiry[expiryIndex];
            /* end date of the swaption */
            end = sMat->expiry[expiryIndex] +
                            sMat->maturity[idx];
            /* e is the last swaption index, in timeline ZeroDate */
            IDXZ(end, &e); 
            /* s is the first swaption index, in timeline ZeroDate */
            IDXZ(start, &s); 

            if (s>e){
                DrlErrMsg("Swaption definition error "
                "(start=%d, end=%d, startIndex=%d, endIndex=%d)\n",
                start, end, s,e);
                goto done;
            }
        
            for(di = s; di < e; di++){
                for(dj = s; dj < e; dj++){
                    K_inte = 0.0;
                    /* Should be 1-<=expiryIndex or 0-<expiryIndex
                    *  according to storage convention in KOVT*/
                    for(dk = 0; dk < expiryIndex; dk++){
                        K_inte += KOVT[(di-1-dk)][(dj-1-dk)][dk];
                    }
                    sum += LAMBDA[di][s][e]*LAMBDA[dj][s][e]*K_inte;
                }
            }
            ALPHAM[i][i] -= sum; 

            /* CHECK ALPHAM>0 */
            if(ALPHAM[i][i]<0){
                ALPHAM[i][i] = 0.0001;
            }
            else{
            ALPHAM[i][i] = sqrt(ALPHAM[i][i]);
            }
        }
    }

    /* Fill lambda matrix */
    for(j = 0; j < nbSwp; j++){
        IDS(j, &idx);
        /* start date of the corresponding swaption */
        start = sMat->expiry[expiryIndex];  
        /* end date of the swaption */
        end = sMat->expiry[expiryIndex] + swpMat[j]; 
        /* e is the last swaption index, in timeline newTimeLine */
        IDXZ(end, &e);  
        /* s is the first swaption index, in timeline newTimeLine */
        IDXZ(start, &s); 
        /* can't fill the first line directly cos maturity[0]!=0 */      
        LAMBDAM[0][j] = LAMBDA[s][s][e];  
        for(i = 1; i < that->nbMatMax; i++){
            if(i<=idx){
            t1 = sMat->expiry[expiryIndex] + sMat->maturity[i-1];
            IDXZ(t1, &i1);    
            LAMBDAM[i][j] = LAMBDA[i1][s][e];  
            }
            else{
            LAMBDAM[i][j] = 0.0;
            }
        }
    }
 
    if (DrlPseudoInverseSM(that->nbMatMax,
                        nbSwp,
                        LAMBDAM,
                        &INVLM,
                        1,
                        addI->smoothingCoeffs[(int)(expiryIndex*that->dt)],
                        1)!=0){
        errCode = 1;
        goto done;
    }

    DrlMatrixCopySM(&INVLMT, INVLM, nbSwp, that->nbMatMax);
    DrlMatrixTransposeSM(&INVLMT, nbSwp, that->nbMatMax);

#undef    IDS
    errCode = 0;
done:

    if (errCode != 0) {
        DrlErrMsg("BootsUpdateMatrices: failed\n");
    }
    return(errCode);
}

/*-------------------------------------------------------------
 * Update covariance tensor with new value of covariance matrix.
 * Copy correlations and volatilities to OutputMatrix.
 */
int 
BootsCopyResults(
    BootsData* that,   /*(I/O) Model parameters*/
    int expiryIndex    /*(I) Expiry Index for storage in tensor*/
)
{
    int errCode = 0;
    int i, j;

    /* Copy Covariance matrix at the expiryIndex to Tensor*/
    for(i = 0; i < that->nbMatMax; i++){
        for(j = 0; j < that->nbMatMax; j++){
			if((i==j)&&(KOVM[i][j]<0)){
            KOVT[i][j][expiryIndex] = 0.0;
            }
            else{
            KOVT[i][j][expiryIndex] = KOVM[i][j];
			}
			
        }
    }

    errCode = 0;

    if(errCode != 0){
        DrlErrMsg("BootsCopyResults: failed\n");
    }
    return(errCode);
}

/*f-------------------------------------------------------------
 * Test Function.
 * Test the possible shapes of the correlation matrix for a given target. 
 * This is a direct test for which the R matrix
 * has to match some target correls.
 *                                                    
 */
int 
BootsOptimizeCorrelations(
    BootsData* that,       /*(I) Model parameters */
    Target* target,        /*(I) Target correlations*/
    double* initialThetas  /*(I) Initial thetas*/
)
{

    double  minimum;
    int errCode = 1;
    int iter;


    /* Call the optimizer function*/
    /* Find the Optimal Thetas set of parameters to build
     * the R matrix and match target correls*/
    /* Do it only when it is required*/
    switch(target[0].method){
    case 0:
    /* Exponent: 1 param*/
        if(frprmn1(that,
        target,
        initialThetas,
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
        initialThetas,
        4,
        1e-30,
        &iter,
        &minimum)!=0)
        goto done;
        break;

    case 2:
    /* Bochner2: 8 params*/
        if(frprmn1(that,
        target,
        initialThetas,
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
        initialThetas,
        12,
        1e-5,
        &iter,
        &minimum)!=0)
        goto done;
        break;
        
        default:
        /*Default not implemented*/
        goto done;

    }

    /* Update correlations and covariance and copy results*/
    if(BootsGetNewCorrMatrix(
       that,
       target,
       initialThetas)!=0){
        goto done;
    }

    errCode = 0;
done:

    if (errCode != 0) {
        DrlErrMsg("BootsOptimizeCorrelations: failed \n");
    }
    return(errCode);
}

/*f-------------------------------------------------------------
 *  Test Function.
 *  Compute swaption correlations at a given expiry. 
 *  Swaption expiry index input by user as viewIndex.
 *                                                       
 */
int 
BootsComputeSwaptionCorrelations(
    BootsData* that,              /*(I/O) Model parameters */
    SwMat* sMat,                  /*(I) Initial SwaptionMatrix*/
    double**   outSwMat,          /*(I) Input SwaptionVolMatrix
                                   *    from bootstrapping*/
    double**   outSwCorrMat       /*(O) Output SwaptionCorrMat*/
)
{
    int i, j, s, e1, e2, di, dj, dk;
    double  start, end1, end2, K_inte, sum;
    int errCode = 1;
    outSwCorrMat = DrlMatrixNew(that->nbMatMax, that->nbMatMax);

    /* Test are there enough expiries ?*/
    if(that->viewIndex>=that->nbExpMax){
        DrlErrMsg("More expiry lines to be specified \n");
        goto done;
    }
    start = sMat->expiry[that->viewIndex];
    /* s becomes the first swaption index, in timeline ZeroDate */
    IDXZ(start, &s); 
    for(i = 0; i < (that->nbMatMax/(double)2); i++){
        /* Find the end for the corresponding swaption1*/
        end1 = sMat->expiry[that->viewIndex] + sMat->maturity[i]; 
        /* e1 becomes the last swaption1 index, in timeline ZeroDate */
        IDXZ(end1, &e1); 
        for(j = i; j < (that->nbMatMax/(double)2); j++){
            sum = 0.0;
            /* Find the end for the corresponding swaption */
            end2 = sMat->expiry[that->viewIndex] + sMat->maturity[j]; 
            /* e2 becomes the last swaption2 index, in timeline ZeroDate */
            IDXZ(end2, &e2); 
            
            if ((s>e1)||(s>e2)){
                goto done;
            }
            for(di = s; di < e1; di++){
                for(dj = s; dj < e2; dj++){
                    K_inte = 0;
                    /* dk index should be 1-<=expiryIndex
                     * or 0-<expiryIndex according to conventions */
                    for(dk = 0; dk <= that->viewIndex; dk++){ 
                        K_inte += KOVT[(di-1-dk)][(dj-1-dk)][dk];
                    }
                    sum+= LAMBDA[di][s][e1]*LAMBDA[dj][s][e2]*K_inte;
                }
            }
            /* Ouptut Correlations */
            outSwCorrMat[i][j] = sum/sMat->expiry[that->viewIndex];
            outSwCorrMat[i][j] /= (outSwMat[that->viewIndex][i]*
                               outSwMat[that->viewIndex][j]);
        }
    }

    errCode = 0;
done:
    if (errCode != 0) {
        DrlErrMsg("BootsComputeSwaptionCorrelations: failed \n");
    }
    return(errCode);
}

/*f-------------------------------------------------------------
 * Test Function.
 * Compute implied forward correlations at a given expiry. 
 *                                                       
 */
int 
BootsComputeImpliedForwardCorrelations(
    BootsData* that,         /*(I/O) Model parameters */
    double** outFwdCorrMat   /*(O) Forward CorrelationMatrix*/
)
{
    int i, j, k, i1, i2;
    double inte1, inte2, inte3;
    int errCode = 1;

    outFwdCorrMat = DrlMatrixNew(that->nbMatMax, that->nbMatMax);

    /*Are there enough expiries ?*/
    if(that->viewIndex>=that->nbExpMax){
        DrlErrMsg("More expiry lines to be specified \n");
        goto done;
    }
    
    for(i = 0; i < (that->nbMatMax/(double)2); i++){
        for(j = 0; j <=i; j++){
            inte1 = inte2 = inte3 = 0.0;
            for(k = 0; k <= that->viewIndex; k++){
                i1 = i+that->viewIndex-k;
                i2 = j+that->viewIndex-k;

                inte1 +=that->kovTens[i1][i1][k];
                inte2 +=that->kovTens[i2][i2][k];
                inte3 +=that->kovTens[i1][i2][k];

            }

            outFwdCorrMat[i][j] = inte3;
            outFwdCorrMat[i][j] /= sqrt(inte1);
            outFwdCorrMat[i][j] /= sqrt(inte2);
        }
    }

    errCode = 0;
done:
    if (errCode != 0) {
        DrlErrMsg("BootsComputeImpliedForwardCorrelations: failed \n");
    }
    return(errCode);
}

/*f-------------------------------------------------------------
 * Test Function.
 * Check target, for a given expiryIndex.
 *    
 *    
 *                                                    
 */
int 
BootsCheckTarget(
    BootsData* that, /*(I) Model parameters */
    Target* target   /*(I) Target correlations*/
)
{
    int i, cInd;
    double start, end;
    int errCode = 1;

    cInd = (int)(that->expIndex*that->dt);
    
	/*General*/
    for(i = 0; i < target[cInd].nb; i++)
    {
        if((target[cInd].x[i]<0)||
           (target[cInd].x[i]>that->swpMat[that->nbSwpMat-1])||
           (target[cInd].x[i]<that->swpMat[0])){
           goto done;
		}

        if((target[cInd].y[i]<0)||
           (target[cInd].y[i]>that->swpMat[that->nbSwpMat-1])||
           (target[cInd].x[i]<that->swpMat[0])){
           goto done;
        }

        if((target[cInd].val[i]<-1)||(target[cInd].val[i]>1)){
           goto done;
		}
    }     

	/*If fCorr=1 optimization on Forward correlations*/
	if(that->fCorr ==1){
	    for(i = 0; i < target[cInd].nb; i++)
		{
	        start = (that->expIndex+1)*that->dt;
	        end = (that->expIndex)*that->dt + that->swpMat[that->nbSwpMat-1];
        
		    if((target[cInd].x[i]<start)||
               (target[cInd].x[i]>end)){
               goto done;
			}
		
            if((target[cInd].y[i]<start)||
               (target[cInd].y[i]>end)){
               goto done;
			}
		}
	}
    errCode = 0;
done:
    if (errCode != 0) {
        DrlErrMsg("BootsOptimizeCorrelations: failed \n");
    }
    return(errCode);
}
