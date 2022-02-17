/****************************************************************
 * Module:     Boots
 * File:     bootsout
 * Function: Output routines    
 * Author:     CA
 * Revision:
 *****************************************************************/
#include "drlstd.h"            /* platform compatibility */
#include <math.h>

#include "drlerr.h"
#include "drltime.h"
#include "drlsort.h"
#include "drlio.h"
#include "drllinsm.h"

#define    _Boots_SOURCE
#include "boots.h"

#define    __DEBUG__
#undef    __DEBUG__

/*f-------------------------------------------------------------
 * Output Function.
 * Compare swaption prices from boostrapping 
 * and swaption prices from input data.                                                   
 */
int 
BootsRecomputeSwaptionVol(
    BootsData* that,       /*(I/O) Model parameters */
    SwMat* sMat,           /*(I) Initial SwaptionMatrix*/
    double** outSwpMat,    /*(O) Swaption Matrix from model*/
    double** errSwpMat     /*(O) Errors on output SwaptionMatrix*/
)
{
    int i, j, expiryIndex, s, e, di, dj, dk, i1, i2;
    double iter, minimum, start, end, t1, t2, K_inte, sum;
    int errCode = 1;

#undef  IDX
#define IDX(x, idx) DrlDoubleArrayClosestIdx(that->newTimeLine,that->newNbDates,x,(idx))

    for(expiryIndex = 0; expiryIndex < that->nbExpMax; expiryIndex++){
        for(i = 0; i < (that->nbMatMax/(double)2); i++){
            sum = 0.0;
            start = sMat->expiry[expiryIndex];
            /* find the end for the corresponding swaption */
            end = sMat->expiry[expiryIndex] + sMat->maturity[i]; 
            /* e becomes the last swaption index, in timeline ZeroDate */
            IDX(end, &e); 
            /* s becomes the first swaption index, in timeline ZeroDate */
            IDX(start, &s); 
            if (s>e){
                DrlErrMsg("Swaption definition error "
                "(start=%d, end=%d, startIndex=%d, endIndex=%d)\n",
                start, end, s,e);
                goto done;
            }
            for(di = s; di < e; di++){
                for(dj = s; dj < e; dj++){
                    K_inte = 0;
                    /* dk index should be 1-<=expiryIndex 
                     * or 0-<expiyrIndex according ot conventions */
                    for(dk = 0; dk <= expiryIndex; dk++){ 
                        K_inte += KOVT[(di-1-dk)][(dj-1-dk)][dk];
                    }
                sum+= LAMBDA[di][s][e]*LAMBDA[dj][s][e]*K_inte;
                }
            }
            /* Distance to expected price in bp */
            outSwpMat[expiryIndex][i] = sqrt(sum/sMat->expiry[expiryIndex]);
            errSwpMat[expiryIndex][i] = 10000*
                                        (sqrt(sum/sMat->expiry[expiryIndex])-
                                        sMat->swaptionMatrix[expiryIndex][i]);
        }
    }
    errCode = 0;
done:
#undef    IDX
    if(errCode != 0){
        DrlErrMsg("BootsRecomputeSwaptionPrices: failed \n");
    }
    return(errCode);
}

/*f-------------------------------------------------------------
 * Output Function.
 * The following routine allow to query arbitrary
 * forward/average vol of any rate.
 * Four dates are required.
 * They correspond to:
 * rateObs   : Start of observation period
 * rateReset : End of observation period                                                   
 * rateStart : first Forward constituting "rate"
 * rateLast  : last Forward constituting "rate"
 * (Rate     = Forward if StartIndex == LastIndex,
 *  Rate = SwaptionRate otherwise) 
 */
int 
BootsDataRateVol(
    BootsData* that,        /*(I)*/ 
    double    rateObs,      /*(I)*/
    double    rateReset,    /*(I)*/
    double    rateStart,    /*(I)*/
    double    rateEnd,      /*(I)*/
    int    rateFreq,        /*(I), specific frequency of the rate, 
                             * different from timeline freq*/  
    int    outputType,      /*(I), 0 = %vol, 1 = bp vol*/
    double *outputVol       /*(O)*/
)      
{
    int i, j, k, l, i1, i2, di, dj;
    double inte1, inte2, inte3, K_inte, sum;
    int errCode = 1;
    int rateObsIdx, rateResetIdx, rateStartIdx, rateEndIdx;
	/*For speed*/
	int nbMax;
	nbMax = (that->nbMatMax - 1);

#undef  IDX
#define IDX(x, idx) DrlDoubleArrayClosestIdx(that->newTimeLine,that->newNbDates,x,(idx))

    /* Here convert dates to indexes*/
    IDX(rateObs, &rateObsIdx);   
    IDX(rateReset, &rateResetIdx);   
    IDX(rateStart, &rateStartIdx);   
    IDX(rateEnd, &rateEndIdx);   

    /* Indexes coherence*/
    if(rateObsIdx >= rateResetIdx){
        DrlErrMsg("Observation Index should "
        "be before Reset Index, see prototype \n");
        goto done;
    }
    
    if((rateStartIdx<rateObsIdx)){
        DrlErrMsg("Observation Index should "
                    "be before rate Start Index, see prototype \n");
        goto done;
    }

    if((rateEndIdx<rateStartIdx)){
        DrlErrMsg("startIdx should "
                    "be after endIdx !, see prototype \n");
        goto done;
    }

    /* Double summmation*/
    sum = 0.;
    for(i = rateStartIdx; i <= rateEndIdx; i++){
        for(j = rateStartIdx; j <= rateEndIdx; j++){
            K_inte = 0.;
            /* Time-Integration between T_rateObsIdx and T_rateResetIdx */
            for(k = rateObsIdx; k < rateResetIdx; k++){
                i1 = i-k-1;
                i2 = j-k-1;
				/*Check i1 and i2 within the boundaries */
				if((i1<0)||(i1>nbMax)){
					goto done;
				}
					
				if((i2<0)||(i2>nbMax)){
					goto done;
				}

                K_inte += KOVT[i1][i2][k];
            }
               sum+= LAMBDA[i][rateStartIdx][rateEndIdx+1]*
                     LAMBDA[j][rateStartIdx][rateEndIdx+1]*K_inte;
            }
    }

    /* Normal/LogNormal*/
    if(outputType ==1)
    {
    sum /= (double)(rateReset-rateObs);
    (*outputVol) = sqrt(sum);
    }
    else
    {
    sum = (exp(sum) - 1);
    sum /= (double)(rateReset-rateObs);
    (*outputVol) = sqrt(sum);
    }

    errCode = 0;
 
done:
#undef	IDX
    if (errCode != 0) {
        DrlErrMsg("BootsDataRateVol: failed \n");
    }
    return(errCode);
}

/*f-------------------------------------------------------------
 * Output Function.
 * The following routine allow to query
 * arbitrary forward/average correlation of any rate.
 * Six dates are required.
 * They correspond to:
 * rateObs    : Start of observation period
 * rateReset  : End of observation period                                                   
 * rateStart1 : First Forward constituting "rate1"
 * rateLast1  : Last Forward "rate1"
 * (Rate1     = Forward if Start1 == Last1, Rate1 = SwaptionRate otherwise) 
 * rateStart2 : First Forward constituting "rate2"
 * rateLast2  : Last Forward "rate2" 
 * (Rate2     = Forward if Start2 == Last2, Rate2 = SwaptionRate otherwise)  
 */
int 
BootsDataRateCorrel(
    BootsData* that,        /*(I)*/ 
    double    rateObs,      /*(I)*/
    double    rateReset,    /*(I)*/
    double    rateStart1,   /*(I)*/
    double    rateEnd1,     /*(I)*/
    double    rateStart2,   /*(I)*/
    double    rateEnd2,     /*(I)*/
    int    rateFreq,        /*(I), specific frequency of the rate, 
                             * different from timeline freq*/  
    int    outputType,      /*(I), 0 = %vol, 1 = bp vol*/
    double *outputCorrels   /*(O)*/
)  
{
    int i, j, k, l, i1, i2, di, dj;
    double K_inte, sum, vol1, vol2;
    int errCode = 1;
    int rateObsIdx, rateResetIdx;
    int rateStartIdx1, rateEndIdx1;
    int rateStartIdx2, rateEndIdx2;
	/*For speed*/
	int nbMax;
	nbMax = (that->nbMatMax - 1);
 
#undef  IDX
#define IDX(x, idx) DrlDoubleArrayClosestIdx(that->newTimeLine,that->newNbDates,x,(idx))

    /*Here convert dates to indexes*/
    IDX(rateObs, &rateObsIdx);   
    IDX(rateReset, &rateResetIdx);   
    IDX(rateStart1, &rateStartIdx1);   
    IDX(rateEnd1, &rateEndIdx1);   
    IDX(rateStart2, &rateStartIdx2);   
    IDX(rateEnd2, &rateEndIdx2);   

    /*Indexes coherence*/
    if(rateObsIdx >= rateResetIdx){
        DrlErrMsg("Observation Index should "
                    "be before Reset Index, see prototype \n");
        goto done;
    }
    
    if((rateStartIdx1<rateObsIdx)){
        DrlErrMsg("Observation Index1 should "
                    "be before rate Start Index1, see prototype \n");
        goto done;
    }

    if((rateStartIdx2<rateObsIdx)){
        DrlErrMsg("Observation Index2 should "
                    "be before rate Start Index2, see prototype \n");
        goto done;
    }

    if((rateEndIdx1<rateStartIdx1)){
        DrlErrMsg("startIdx1 should "
                    "be after endIdx1 !, see prototype \n");
        goto done;
    }

    if((rateEndIdx2<rateStartIdx2)){
        DrlErrMsg("startIdx2 should "
                    "be after endIdx2 !, see prototype \n");
        goto done;
    }

    /*Double summmation*/
    sum = 0.;
    for(i = rateStartIdx1; i <= rateEndIdx1; i++){
            for(j = rateStartIdx2; j <= rateEndIdx2; j++){
                K_inte = 0.;
                /*Time-Integration between T_rateObsIdx and T_rateResetIdx*/
                for(k = rateObsIdx; k < rateResetIdx; k++){
                    i1 = i-k-1;
                    i2 = j-k-1;
					/*Check i1 and i2 within the boundaries */
					if((i1<0)||(i1>nbMax)){
						goto done;
					}
					
					if((i2<0)||(i2>nbMax)){
						goto done;
					}

                    K_inte += KOVT[i1][i2][k];
                }
                sum+= LAMBDA[i][rateStartIdx1][rateEndIdx1+1]*
                      LAMBDA[j][rateStartIdx2][rateEndIdx2+1]*
                      K_inte;
            }
    }

    if(outputType ==1)
    {
        sum /= (double)(rateReset-rateObs);
    }
    else
    {
        sum  = (exp(sum)-1);
        sum /= (double)(rateReset-rateObs);
    }
    /*Divide by the two vols*/
    /*They can't be equal to zero*/
    if(BootsDataRateVol(that, 
                 rateObsIdx,        
                 rateResetIdx,    
                 rateStartIdx1,    
                 rateEndIdx1, 
                 rateFreq,
                 outputType,
                 &vol1)!=0){
        errCode = 1;
        goto done;
    }

    if(BootsDataRateVol(that, 
                 rateObsIdx,        
                 rateResetIdx,    
                 rateStartIdx2,    
                 rateEndIdx2,
                 rateFreq,
                 outputType,
                 &vol2)!=0){
        errCode = 1;
        goto done;
    }

    sum /= (vol1*vol2);
    (*outputCorrels) = (sum);        

    errCode = 0;
        
#undef    CONV
done:
#undef IDX
    if (errCode != 0) {
        DrlErrMsg("BootsDataRateCorrel: failed \n");
    }
    return(errCode);
}

/*-------------------------------------------------------------
 * Prints a matrix to fp.
 * nl = number of ligns.
 * nc = number of cols.
 */
int
BootsPrintMat(
    double*** matrix,    /*(I) matrix to be printed */
    int nl,             /*(I) number of cols */
    int nc,             /*(I) number of cols */
    FILE *fp)            /*(O) ouptut file*/
{
    static    char    routine[] = "BootsPrintMat";
    int m, n;
   
    /*Print covariance and lambda tensors*/
    DrlFPrintf(fp, "matrix: \n");
   
    for (m=0; m<nl;m++) {
        for(n=0; n<nc; n++){
            DrlFPrintf(fp, "%7.4f ", (*matrix)[m][n]);
        }
        DrlFPrintf(fp, "\n");
    }

    return(0);
}


/*f-------------------------------------------------------------
 * Test Function.
 * Compute implied forward correlations at a given expiry. 
 *                                                       
 */
int 
BootsImpForCorr(
    BootsData* that,
    int index, 
    double** outMat
)
{
    int i, j, k, i1, i2;
    double inte1, inte2, inte3;
    int errCode = 1;

    /*Are there enough expiries ?*/
    if(index>=that->nbExpMax){
        DrlErrMsg("More expiry lines to be specified \n");
        goto done;
    }
    
    for(i = 0; i < (that->nbMatMax/(double)2); i++){
        for(j = 0; j <=i; j++){
            inte1 = inte2 = inte3 = 0.0;
            for(k = 0; k <= index; k++){
                i1 = (i + index - k);
                i2 = (j + index - k);

                inte1 +=that->kovTens[i1][i1][k];
                inte2 +=that->kovTens[i2][i2][k];
                inte3 +=that->kovTens[i1][i2][k];
            }

            outMat[i][j] = inte3;
            outMat[i][j] /= sqrt(inte1);
            outMat[i][j] /= sqrt(inte2);
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
 * 
 * Transform Forward correlations matrix into Swaption correlations Matrix. 
 *                                                       
 */
int 
BootsGenSwapCorrs(
    BootsData* that,        /*(I) Data*/
    double** inputForCorr,  /*(I) Input forward correlations*/
	double expiry,          /*(I) Start of the swap*/
	int nbFor,              /*(I) Number of maturities to consider*/
    double*** outputForCorr /*(O) Output swap correlation matrix*/
)
{
    int i, j, k, i1, i2, expInd, s,e, Endi, Endj;
    double start, end, t1;
	double tEndi, tEndj;
	double covar, var1, var2;
    int errCode = 1;
	/*Find corresponding expiryIndex in Timeline*/ 
#undef  IDX
#define IDX(x, idx) DrlDoubleArrayClosestIdx(that->newTimeLine,that->newNbDates,x,(idx))
	IDX(expiry, &s);
	
	/*calculate correlations*/
	for(i = 0; i < nbFor; i++){
		for(j = 0; j < nbFor; j++){
			/*Access lambda*/
			tEndi = expiry + (i+1)*that->dt;
			tEndj = expiry + (j+1)*that->dt;

			IDX(tEndi, &Endi);
			IDX(tEndj, &Endj);
			covar = 0;
			var1 = 0;
			var2 = 0;

			/*covariance*/
			for(i1 = s; i1<Endi; i1++){
				for(i2 = s; i2<Endj; i2++){
					covar += LAMBDA[i1][s][Endi]*LAMBDA[i2][s][Endj]*inputForCorr[(i1-s)][(i2-s)];
				}
			}
			/*variance1*/
			for(i1 = s; i1<Endi; i1++){
				for(i2 = s; i2<Endi; i2++){
					var1 += LAMBDA[i1][s][Endi]*LAMBDA[i2][s][Endi]*inputForCorr[(i1-s)][(i2-s)];
				}
			}
			/*variance2*/
			for(i1 = s; i1<Endj; i1++){
				for(i2 = s; i2<Endj; i2++){
					var2 += LAMBDA[i1][s][Endj]*LAMBDA[i2][s][Endj]*inputForCorr[(i1-s)][(i2-s)];
				}
			}

			/*copy correl into new matrix*/
			(*outputForCorr)[i][j] = covar/sqrt(var1*var2);
		}
	}

#undef IDX  
    errCode = 0;
done:
    if (errCode != 0) {
        DrlErrMsg("BootsComputeImpliedForwardCorrelations: failed \n");
    }
    return(errCode);
}
