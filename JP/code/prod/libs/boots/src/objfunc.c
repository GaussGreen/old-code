/****************************************************************
 * Module:    Boots
 * File:    Objective Function
 * Function:Define the objective functin to be minimized     
 * Author:    CA
 * Revision:    
 *****************************************************************/
#include "drlstd.h"            /* platform compatibility */
#include <math.h>
#include "drlmem.h"            /* Drl memory management */
#include "drllinsm.h"
#include "drlio.h"             /* DrlFPrintf */
#include "drlerr.h"
#include "drlsort.h"
#include "drltime.h"

#define NR_END 1

static	double *dvector(long nl, long nh)
{
	double *v;
	v= (double *)MALLOC((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) DrlErrMsg("allocation failure in dvector()\n");
	return v-nl+NR_END;
}
static	void free_dvector(double *v, long nl, long nh)
{
	FREE((void*) (v+nl-NR_END));
}

#define    _Boots_SOURCE
#include "boots.h"

/*f-------------------------------------------------------------
 * Objective function to be minimized.
 * Hard coded target correlations for testing.
 * As usual retuns 0 if OK.
 */
int
BootsObjFunc_initial(
    BootsData* that,   /*(I) Model parameters */
    double* Thetas,    /*(I) Initial parameters*/
    double* res        /*(I) Corresponding value of objective function*/
)
{
    int errCode = 1;
    int i, j, trigger, expInd, nbCorr, x, y;
    double result;
    double inter, inte1, inte2, inte3, targetVal;
    result = 0.0;

    for(i = 0; i < that->nbMatMax; i++){
        result += pow((Thetas[i+1]-(15.256+i)),4);
    }

    *res = result;
    /*ok*/
    errCode = 0;
done:
    if (errCode != 0) {
        DrlErrMsg("BootsObjFunc: failed\n");
    }
    return(errCode); 
}

/*f-------------------------------------------------------------
 * Objective function to be minimized.
 * ! Important: Target.x and Target.y, represent 
 * either absolute maturities, 
 * when swap rates correlations are optimized,
 * or relative (to expIndex) expiry when forward correlations are optimized.
 * 
 * As usual retuns 0 if OK.
 */
int
BootsObjFunc(
    BootsData* that,   /*(I) Model parameters */
    Target* target,    /*(I) Target correlations */
    double* Thetas,    /*(I) Values of Input parameters*/
    double* res        /*(O) Output value of the function*/
)
{
    int errCode = 1;
    int i, j, trigger, expInd, nbCorr, x, y, k;
    double result;
    double inter, inte1, inte2, inte3, targetVal;

    /*CurrentT contains the current target correlation point*/
    Target currentT;
    result = 0.0;

    /*Find the index in that->swpMat, for swap rates correlations*/
#undef     IDXSM
#define    IDXSM(x, idx)    DrlDoubleArrayClosestIdx(that->swpMat, that->nbSwpMat, x, (idx))

    /*Find the relative index for forward rates correlations*/
#undef     IDXE
#define    IDXE(x)    (int)((double)(x - (that->expIndex+1)*that->dt)*that->freq)

    /*Read current expiryIndex */
    expInd = that->expIndex;

    /*Read number of target correlations to match */
    nbCorr = target[(int)(expInd*that->dt)].nb;

    /*Define current target correlation point */
    currentT = target[(int)(expInd*that->dt)];

    if(that->fCorr ==0)
    {
        BootsGetNewCorrMatrix(that, target, Thetas);
        /*Loop on target point */
        for(i = 0; i < nbCorr; i++){
        /*copy the current point coordinates, and point value*/
        IDXSM(currentT.x[i], &x);
        IDXSM(currentT.y[i], &y);

        targetVal = currentT.val[i];
        /*Intermediate for the calculation of correlation */
        inte1 = that->corrMat[x][y];
        /*Increase penalty function */
        result += pow(fabs(inte1 - targetVal),2);        
        }
        *res = result;
    }
    else
    {
    /*Update*/
    BootsGetNewCorrMatrix(that, target, Thetas);
    BootsGetNewKovMatrix(that, that->nbSwpMat);

    /*Loop on target point */
    for(i = 0; i < nbCorr; i++){
        /*copy the current point coordinates, and point value*/
        x = IDXE(currentT.x[i]);
        y = IDXE(currentT.y[i]);
        targetVal = currentT.val[i];

        /*Test: is it a correlation or a vol target ?*/
        /*Target == Volatility*/
        if(x == y){
            /*Volatility to match */
            inte1 = sqrt(fabs(KOVM[x][y]));
            /*Increase of the penalty function */
            result += pow(fabs(inte1-targetVal),2);
        }
        /*Target == Correlation*/
        else{
            /*Intermediaries for the calculation of correlation */
            inte1 = fabs(KOVM[x][x]);
            inte2 = fabs(KOVM[y][y]);
            inte3 = fabs(KOVM[x][y]);
            
            /*Test before division by vol */
            if((inte1 < 1e-20)||(inte2 < 1e-20)){
            if((inte1 < 1e-20)&&(inte2 < 1e-20)&&(inte3 < 1e-10)){
            inter = 0.0;
            goto next;
            }
            DrlErrMsg("Division by zero\n");
            goto done;
            }
            inter = fabs(KOVM[x][y])/(sqrt(fabs(KOVM[x][x]))*sqrt(fabs(KOVM[y][y])));
next:
            /*Increase penalty function */
            result += pow(fabs(inter - targetVal),2);            
        }
    }

#undef IDXSM
#undef IDXE
    *res = result;
    }
    errCode = 0;
done:
    if (errCode != 0) {
        DrlErrMsg("BootsObjFunc: failed\n");
    }
    return(errCode); 
}


/*f-------------------------------------------------------------
 * Derivative of the objective function.
 * 
 */
int
BootsdObjFunc(
    BootsData* that,      /*(I) Model parameters */
    Target* target,       /*(I) Target correlations */
    double* thetas,       /*(I) Values of Input parameters*/
    double* derivative    /*(O) Value of the derivative*/
)
{
    int errCode = 1;
    int i;
    int vectSiz;
    double *newVector;
    double h = 1e-5;
    double inter;
    double right, left;

    if(target[(int)(that->expIndex*that->dt)].method ==0){
        vectSiz = 1;
    }
    else{
        vectSiz = (4*target[(int)(that->expIndex*that->dt)].method +1);
    }

    newVector = dvector(0,vectSiz);
    /*Copy vector */
    DrlVectUnaryOper(newVector, vectSiz, "=", thetas);

    for(i = 1; i < vectSiz; i++){
    inter = thetas[i];
    newVector[i] = inter + h;
    thetas[i] = inter - h;
        
    if((BootsObjFunc(that, target, newVector, &right)!=0)||(BootsObjFunc(that, target, thetas, &left)!=0)){
    DrlErrMsg("Failure in objective function, right or left\n");
    goto done;
    }

    derivative[i] = (right-left)/(2*h);
    newVector[i] = inter;
    thetas[i] = inter;
    }

///////////////////////////////
	/*Revert to initial function*/
	BootsGetNewCorrMatrix(that, target, thetas);



    /*ok*/
    errCode = 0;
done:
    free_dvector(newVector,0, vectSiz);
    if (errCode != 0) {
        DrlErrMsg("BootsdObjFunc: failed\n");
    }
    return(errCode);
}

