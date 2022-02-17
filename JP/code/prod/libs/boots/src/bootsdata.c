/****************************************************************
 * Module:    BOOTS
 * File:    BootsData
 * Function:BootsData management and memory routines.    
 * Author:    CA
 * Revision:
 *****************************************************************/
#include "drlstd.h"            /* platform compatibility */
#include "drlio.h"
#include "drlerr.h"
#include "drltime.h"
#include "drlinter.h"
#include "drllinsm.h"

#include <math.h>

#define   _Boots_SOURCE
#include "boots.h"


static    int
BootsMemoryAlloc(
    BootsData* that,  /*(I/O) Model parameters */
    int nbMaxMat,     /*(I) Total number of maturities*/
    int nbMaxExp,     /*(I) Total number of expiries*/
    int freq          /*(I) Frequency*/
);

static    int
BootsMemoryFree(
    BootsData *that   /*(I/O) Model parameters */ 
);


/*f-------------------------------------------------------------
 * Allocate a new BootsData structure.
 * Returns a pointer to a new Bootsdata structure.
 * Returns NULL if failed.                                                            
 */
BootsData*
BootsNew(
    int nbMatMax,  /*(I) Total number of maturities*/
    int nbExpMax,  /*(I) Total number of expiries*/
    int freq       /*(I) Value of the frequency*/
)
{
static    char        routine[] = "BootsNew";
    BootsData*    that = (BootsData*)NULL;
    int        status = FAILURE;

    if ((that = NEW(BootsData)) == NULL) {
        DrlErrMsg("%s: memory allocation failed\n", routine);
        goto done;
    }

    if(BootsMemoryAlloc(that, nbMatMax, nbExpMax, freq) != 0)
        goto done;

    status = SUCCESS;
done:
    if (status != SUCCESS) {
        if (that != NULL) {
            FREE(that);
            that = NULL;
        }
        DrlErrMsg("%s: failed\n", routine);
    }
    return(that);
}

/*f-------------------------------------------------------------
 * Copies a BootsData structure.
 * ZeroBond, Timeline, lambdaMatrix, lambdaTensor,
 * covariance matrix, covariance tensor.     
 * And zero curve.
 */
BootsData*
BootsNewCopy(
    BootsData* that    /*(I/O) Model parameters */
)
{
static    char routine[] = "BootsNewCopy";
    int    status = FAILURE;
    int i, j, k;
    BootsData *copy = NULL;
    
    /*Allocate a BootsData structure*/
    copy = BootsNew(that->nbMatMax, that->nbExpMax, that->freq);
    if (copy == NULL) goto done;

    /*Copy data*/
    copy->newNbDates = that->newNbDates; 
    copy->nbExpMax = that->nbExpMax;
    copy->nbMatMax = that->nbMatMax;
	copy->fCorr = that->fCorr;
    copy->freq = that->freq;
    copy->dt = that->dt;

    for(i = 0; i < that->newNbDates; i++) {
        copy->newTimeLine[i] = that->newTimeLine[i];
        copy->newZeroBond[i] = that->newZeroBond[i];        
    }
    
    for(i = 0; i <= that->newNbDates; i++){
        for(j = 0; j <= that->newNbDates; j++){
            for(k = 0; k <= that->newNbDates; k++){
                copy->lambda[i][j][k] =  that->lambda[i][j][k];
            }
        }
    }

    for(i = 0; i < that->newNbDates; i++){
       for(j = 0; j < that->newNbDates; j++){
           for(k = 0; k < that->nbExpMax; k++){
               copy->kovTens[i][j][k] =  that->kovTens[i][j][k];
           }
       }
    }
    
    DrlMatrixCopySM(
        &(copy->kovMat),
        that->kovMat,
        that->nbMatMax,
        that->nbMatMax);

    DrlMatrixCopySM(
        &(copy->lambdaMat),
        that->lambdaMat,
        that->nbMatMax,
        that->nbMatMax);

    /* Copy ZcCurve */
    copy->fZcCurve = DrlDCurveNewCopy(that->fZcCurve);
    if (copy->fZcCurve == NULL) goto done;
   
    /* Made it through OK */
    status = SUCCESS;
done:
    if (status != SUCCESS) {
        BootsFree(copy);
        DrlErrMsg("%s: failed\n", routine);
    }
    return(copy);
}

/*------------------------------------------------------
 * Allocate Memory  
 */
static    int
BootsMemoryAlloc(
    BootsData* that,  /*(I/O) Model parameters */
    int nbMaxMat,     /*(I) Total number of maturities*/
    int nbMaxExp,     /*(I) Total number of expiries*/
    int freq          /*(I) Frequency*/
)
{
    int    errCode=1,i,j;
    int nbDatesNewTimeLine;
    int nbMaxMaturity;
    int nbMaxExpiry;
    double* oneDiag;

#undef    CHECK
#define    CHECK(ptr)    if ((ptr) == NULL) {errCode = (-4); goto done;}

    /*Allocate newTimeLine and newZeroBond maximal sizes */
    nbDatesNewTimeLine =  (int)(freq*(nbMaxMat+nbMaxExp) + 1);
    that->newNbDates = nbDatesNewTimeLine;

    /*Maximal number of maturities */
    nbMaxMaturity = (int)(freq*nbMaxMat);
    nbMaxExpiry = (int)(freq*nbMaxExp);
    that->nbMatMax = nbMaxMaturity;                    
    that->nbExpMax = nbMaxExpiry;

    that->newTimeLine = NULL;
    that->newZeroBond = NULL;

    CHECK(that->newTimeLine = NEW_ARRAY(double, nbDatesNewTimeLine))
    CHECK(that->newZeroBond = NEW_ARRAY(double, nbDatesNewTimeLine))

    /*Allocate Lambdas */
    that->lambda=NULL;
    CHECK(that->lambda=NEW_ARRAY(double**,(nbDatesNewTimeLine+1)))
    for (i=0; i<=nbDatesNewTimeLine; i++) {
        that->lambda[i]=NULL;
        CHECK(that->lambda[i]=NEW_ARRAY(double*,(nbDatesNewTimeLine+1)))
        for(j=0; j<=nbDatesNewTimeLine; j++){
            that->lambda[i][j]=NULL;
            CHECK(that->lambda[i][j]=NEW_ARRAY(double,(nbDatesNewTimeLine+1)))
            that->lambda[i][j][0]=1.0;
        }
    }

    /*Fill diagonal of alphaMat with something*/
    oneDiag = NULL;
    CHECK(oneDiag = NEW_ARRAY(double, nbMaxMaturity));
    for(i = 0; i < nbMaxMaturity; i++)
        oneDiag[i] = 1.0;

    /*Allocate intermediate matrices*/
    that->lambdaMat = NULL;
    that->invLambdaMat = NULL;
    that->invLambdaMatT = NULL;
    that->alphaMat = NULL;
    that->corrMat = NULL;
    that->kovMat = NULL;

    CHECK(that->lambdaMat = DrlMatrixNew(nbMaxMaturity, nbMaxMaturity))
    CHECK(that->invLambdaMat = DrlMatrixNew(nbMaxMaturity, nbMaxMaturity))
    CHECK(that->invLambdaMatT = DrlMatrixNew(nbMaxMaturity, nbMaxMaturity))
    CHECK(that->alphaMat = DrlMatrixNewDiag(nbMaxMaturity, oneDiag))
    CHECK(that->corrMat = DrlMatrixNew(nbMaxMaturity, nbMaxMaturity))
    CHECK(that->kovMat = DrlMatrixNew(nbMaxMaturity, nbMaxMaturity))    

    /*    Allocate the covariance tensor. 
    *    kovTens stores the values of the
    *    covariance matrix at each expiry 
    *    date.
    *    Since the covariance matrix is MAT*MAT,
    *    and kovTens is (MAT+EXP)*(MAT+EXP), a
    *    "filling" method is used to provide 
    *    additional data: $BootsCopyResults$.
    */
    that->kovTens = NULL;
    CHECK(that->kovTens = NEW_ARRAY(double**, (nbDatesNewTimeLine)))
    for(i=0; i<(nbDatesNewTimeLine); i++){
        that->kovTens[i] = NULL;
        CHECK(that->kovTens[i] = NEW_ARRAY(double*, (nbDatesNewTimeLine)))
        for(j=0; j<(nbDatesNewTimeLine); j++){
            that->kovTens[i][j] = NULL;
            CHECK(that->kovTens[i][j] = NEW_ARRAY(double, nbMaxExpiry))
        }
    }

    /*Allocate memory for intermediarie smatrices */
    that->productMat = NULL;
    that->productMat2 = NULL;
    that->productMat3 = NULL;

    CHECK(that->productMat = DrlMatrixNew(nbMaxMaturity, nbMaxMaturity))
    CHECK(that->productMat2 = DrlMatrixNew(nbMaxMaturity, nbMaxMaturity))
    CHECK(that->productMat3 = DrlMatrixNew(nbMaxMaturity, nbMaxMaturity))    
    /********************************************************************/
    /*
     *
     */
    /*Ok*/
    errCode = 0;


done:

    /*Free intermediary oneDiag vector */
    if(oneDiag != NULL)FREE_ARRAY(oneDiag);
    if (errCode != 0) {
         DrlErrMsg("BootsMemoryAlloc: malloc failure ");
    }
    return(errCode);

#undef    CHECK
}

/*f-----------------------------------------------------------------------
 * Free a BootsData structure.
 *                                                             
 */
int
BootsFree(
    BootsData *that   /*(I/O) Model parameters */
)
{
    BootsMemoryFree(that);
    if (that != NULL) FREE((char*) that);
    return(0);
}

/*------------------------------------------------------------------------
 * Frees memory allocated by BootsAlloc.
 * Returns 0 iff OK.
 */
static    int
BootsMemoryFree(
    BootsData *that   /*(I/O) Model parameters */ 
)
{
    int    i,j;
    int nbMaxMaturity, nbMaxExpiry, nbDates;

    /*Maximal number of expiries, maturities, dates in newTimeLine */
    nbDates = that->newNbDates;
     nbMaxMaturity = that->nbMatMax;
    nbMaxExpiry = that->nbExpMax; 

    if ((that != NULL)){
    FREE_ARRAY(that->newTimeLine);
    FREE_ARRAY(that->newZeroBond);

    /*Free Lambda*/
    for(i = 0; i<(nbDates+1); i++){
    for(j = 0; j<(nbDates+1); j++)
        if(that->lambda[i][j] != NULL)
        FREE_ARRAY(that->lambda[i][j]);
    }
    
    for(i = 0; i<(nbDates+1); i++){
        if(that->lambda[i] != NULL)
        FREE_ARRAY(that->lambda[i]);
    }

    if(that->lambda != NULL)
    FREE_ARRAY(that->lambda);
    /*End FreeLambda*/

    /*Free kovTens*/
    for(i = 0; i<(that->newNbDates); i++){
        for(j = 0; j<(that->newNbDates); j++){
            if(that->kovTens[i][j]!=NULL)
            FREE_ARRAY(that->kovTens[i][j]);
        }
    }
    
    for(i = 0; i<(that->newNbDates); i++){
        if(that->kovTens[i]!=NULL)
        FREE_ARRAY(that->kovTens[i]);
    }

    if(that->kovTens!=NULL)
    FREE_ARRAY(that->kovTens);
    /*End Free koTens*/
    
    /*Free all other matrices*/
    if(that->kovMat !=NULL)
    DrlMatrixFree(that->kovMat, nbMaxMaturity, nbMaxMaturity);
    
    if(that->alphaMat !=NULL)
    DrlMatrixFree(that->alphaMat, nbMaxMaturity, nbMaxMaturity);
    
    if(that->corrMat !=NULL)
    DrlMatrixFree(that->corrMat, nbMaxMaturity, nbMaxMaturity);
    
    if(that->lambdaMat !=NULL)
    DrlMatrixFree(that->lambdaMat, nbMaxMaturity, nbMaxMaturity);
    
    if(that->invLambdaMat !=NULL)
    DrlMatrixFree(that->invLambdaMat, nbMaxMaturity, nbMaxMaturity);
    
    if(that->invLambdaMatT !=NULL)
    DrlMatrixFree(that->invLambdaMatT, nbMaxMaturity, nbMaxMaturity);
    
    if(that->productMat !=NULL)
    DrlMatrixFree(that->productMat, nbMaxMaturity, nbMaxMaturity);
    
    if(that->productMat2 !=NULL)
    DrlMatrixFree(that->productMat2, nbMaxMaturity, nbMaxMaturity);
    
    if(that->productMat3 !=NULL)
    DrlMatrixFree(that->productMat3, nbMaxMaturity, nbMaxMaturity);
    
    }

    return(SUCCESS);
}

/*-------------------------------------------------------------
 * Prints all cached coefficients to a file pointer.
 * Covariance and Lambda matrix printed at a given expiry.
 * Expiry given as a number of years.
 */
int
BootsPrintCoeff(
    BootsData *that, /*(I) model parameters */
    double expiry,   /*(I) expiry in years */
    FILE *fp)         /*(O) ouptut file*/
{
    static    char    routine[] = "BootsPrintCoeff";
    int m, n;
    int expIndex;
    /*Conversion to index*/
    expIndex = (int)((expiry-that->dt)*that->freq);

    DrlFPrintf(fp, "%s: \n", routine);

    /*Print input data parameters*/
    DrlFPrintf(fp, "Input Data Parameters:\n");
    DrlFPrintf(fp, "Number of dates in timeline:");
    DrlFPrintf(fp, "%d \n", that->newNbDates);
    DrlFPrintf(fp, "Total number of expiries:");
    DrlFPrintf(fp, "%d \n", that->nbExpMax);
    DrlFPrintf(fp, "Total number of maturities:");
    DrlFPrintf(fp, "%d \n", that->nbMatMax);
    DrlFPrintf(fp, "Frequency:");
    DrlFPrintf(fp, "%d \n", that->freq);

    /*Print newtimeline and interpolated values*/
    DrlFPrintf(fp, "NEW Timeline Info:\n");
    for (m=0; m<(that->newNbDates);m++) {
        DrlFPrintf(fp, "%7.4f ", that->newTimeLine[m]);
        DrlFPrintf(fp, "\n");
    }

    DrlFPrintf(fp, "NEW InterpolatedZeroBond Info:\n");
    for (m=0; m<(that->newNbDates);m++) {
        DrlFPrintf(fp, "%7.4f ", that->newZeroBond[m]);
        DrlFPrintf(fp, "\n");
    }

    /*Print covariance and lambda tensors*/
    DrlFPrintf(fp, "covariance at expiry date:");
    DrlFPrintf(fp, "%7.1f \n", expiry);

    for (m=0; m<(that->nbMatMax);m++) {
        for(n=0; n<(that->nbMatMax); n++){
            DrlFPrintf(fp, "%7.4f ", that->kovTens[m][n][expIndex]);
        }
        DrlFPrintf(fp, "\n");
    }

    DrlFPrintf(fp, "lambda matrix at expiry date:");
    DrlFPrintf(fp, "%7.1f \n", expiry);

    for (m=(expIndex+1); m<(that->newNbDates+1);m++) {
        for(n=(expIndex+2); n<(that->newNbDates+1); n++){
            DrlFPrintf(fp, "%7.4f ", that->lambda[m][expIndex+1][n]);
        }
        DrlFPrintf(fp, "\n");
    }

    return(0);
}

/*-------------------------------------------------------------
 * Read a BootsData structure.
 *                                                             
 * 
 * Reads a set of parameters from a wrapper interface
 * Returns 0 iff OK.
 */
int 
BootsWrapRead
(
    BootsData **that,/*(O) Object*/
    DCurve *zcCurve, /*(I) Zero curve*/
    int nbMaxMat,    /*(I) Max Number of Maturities*/
    int nbMaxExp,    /*(I) Max Number of Expiries*/
    int freq         /*(I) Freq to specify the biggest matrices sizes we need*/  
    )    
{
static    char    routine[] = "BootsWrapRead";
    int    status = FAILURE, 
        nZDates = DRLTCURVE_NUMITEMS(zcCurve);

#undef    CHECK
#define    CHECK(cond)    {if (!(cond)) {DrlErrMsg("%s: assertion failed (%s)\n",\
            routine, #cond); goto done;}}

    *that = (BootsData*)NULL;

    /*Check array of arrays */
    CHECK(zcCurve != NULL);
    CHECK(DRLTCURVE_NUMITEMS(zcCurve) >= 1);

    /*Now allocate memory */
    *that = BootsNew(nbMaxMat, nbMaxExp, freq);
    if (*that == NULL) goto done;

    /*Set the zero curve */
    (*that)->fZcCurve = zcCurve;
    (*that)->freq = freq;
    (*that)->dt = (double)1/(double)freq;

    /*Made it through */
    status = SUCCESS;
done:
    if (status != 0) {
        if (*that != (BootsData*)NULL) BootsFree(*that);
        *that = (BootsData*)NULL;
        DrlErrMsg("%s: failed (code %d)\n", routine, status);
    }
    return(status);
#undef    CHECK
}

/*f-------------------------------------------------------------
 * First build a TimeLine according to frequency, nbMaturity, nbFrequency.  
 * Interpolate the zero curve accordingly.
 * As ouput newTimeLine and newZeroBond
 * 
 * Then compute once and for all the lambda coefficients.
 * See doc for available forms.
 * For the moment only Rebonato approximation available.
 */
int
BootsComputeLambda(
    BootsData *that     /* (I/O) model parameters */
)    
{
    static    const char    routine[] = "BootsComputeLambda";
    int        status = FAILURE;
    /* i index for the current forward,
     * s starting date of the swap, e end date of the swap */
    int i,s,e;
    double freq;
    DDate* newTimeLineDate = NULL;
    /*vector used in the spline interpolation*/
    double *y2 = NULL; 
    double *diff = NULL;
    /*values of the derivatives at the beginning 
    * and end of the Zerobond vector */  

    freq = (double)that->freq;
    y2 = NEW_ARRAY(double, that->newNbDates);
    diff = NEW_ARRAY(double, that->newNbDates);
    newTimeLineDate = NEW_ARRAY(DDate, that->newNbDates);

    /* Set time accural */
    that->dt = (double) 1/(double)freq;

    /*Build time line according to swaption frequency */
    /*And corresponding DDate timeline */
    newTimeLineDate[0] = (that->fZcCurve)->ValueDate;
    that->newZeroBond[0] = 1;
    for(i = 1; i < (that->newNbDates); i++){
    (that->newTimeLine)[i] = (double) i/freq;
    DrlDDateAdvanceYears(newTimeLineDate[i-1], that->dt, &newTimeLineDate[i]);
    DrlDCurveDiscFact(that->fZcCurve, 
                        newTimeLineDate[i],
                        &(that->newZeroBond)[i]);
    }

    /*Rebonato approximation
    *lambda[i][s][e] = omega[i][s][e]*li/swap
    *see written doc for notations*/
    for(s = 0; s < (NBD); s++){
        for(e = (s+1); e <= (NBD); e++){
            for(i = s; i<=(e-1); i++){
                /* Time added for homogeneity, see paper */
                LAMBDA[i][s][e] = (ZERO[i]-ZERO[i+1]);
                LAMBDA[i][s][e] /= (ZERO[s]-ZERO[e]);
                LAMBDA[i][s][e] *= sqrt(that->dt);
            }    
        }
    }

    status = SUCCESS;
done:

    if(y2!=NULL){
    FREE_ARRAY(y2);
    }
    if(diff!=NULL){
    FREE_ARRAY(diff);
    }
    if(newTimeLineDate!=NULL){
    FREE_ARRAY(newTimeLineDate);
    }
    if (status!= SUCCESS) {
        DrlErrMsg("%s: failed\n", routine);
    }
    return(status);
}

/*-------------------------------------------------------------
 * Fill a swaption matrix.
 * Maturity and expiry constraints should be given in years.                                                            
 * !!!!!!!!!!Initial swaption matrix should be "annual".
 * (No frequency conversion needed in the reading phase)!!!!!!!!!!!!
 * Input swaption matrix should be 30*30.
 * Output swaption matrix complete at freq level.
 * Output swaption matrix 30*60.
 * Reads a set of parameters from a wrapper interface
 * Returns 0 iff OK.
 */
int
BootsProcessSwaptionMatrix(
    SwMat    *swaptionMatrix,    /*(I/O) Swaption matrix 
                                  * before and after interpolation*/
    int     nbExpiryToMatch,     /*(I) Number of expiries to match*/
    double  *expiryToMatch,      /*(I) Values of the points
                                  * to match (in expiry)*/
    int     nbMaturityToMatch,   /*(I) Number of maturities to match*/
    double     *maturityToMatch, /*(I) Values of the of the points
                                  * to match (in maturity)*/
    double    terminalVol,       /*(I) Total decrease in vol*/        
    double    terminalSlope,     /*(I) Interpolation type*/
    int     freq                 /*(I) Frequency*/
)            
{
static    char    routine[] = "BootsProcessSwaptionMatrix";
    int    status = FAILURE, i, j, k;
    int intpInd = 30;

    /*Variables for Interpolation */
    double * valueKnown=NULL;
    double * valueInit=NULL;
    double point;
    double slope;
    double alpha;
    double epsilon;
    double beta;
    double signBeta;
    double dt;


#undef  MTR
#define MTR(x, y)  (*swaptionMatrix).swaptionMatrix[x][y]

    dt = 1/(double)freq;
    /*Allocation*/
    valueKnown = NEW_ARRAY(double, swaptionMatrix->nbMat);
    /*For spline interpolation*/
    valueInit = NEW_ARRAY(double, swaptionMatrix->nbMat);

    /*Smooth Initial Matrix */
    /*Smooth on maturities*/
    for(j = 0; j < nbExpiryToMatch; j++){
        /*Allocate valueKnown*/
        for(i = 0; i < nbMaturityToMatch; i++){
            valueKnown[i] =
            MTR((int)expiryToMatch[j]-1,(int)maturityToMatch[i]-1);
            valueInit[i] = 1;
        }

        DrlSplineInterp1dInit(
            maturityToMatch,        
            valueKnown,        
            nbMaturityToMatch,            
            1e12,    
            1e12,    
            valueInit);

        for(i = 0; i < freq*(2*intpInd); i++){
            DrlSplineInterp1dInterp(
                maturityToMatch,            /* (I) array of X(i) (i=0,..,n-1) */
                valueKnown,                 /* (I) array of X(i) (i=0,..,n-1) */
                valueInit,                  /* (I) from SplineInit */
                nbMaturityToMatch,          /* (I) # of points */
                (double)(i+1)/(double)freq, /* (I) point to interpolate */
                &point);                    /* (O) interpolated value */

            MTR(((int)expiryToMatch[j]-1),i) = point;
        }
    }


    /*Smooth on expiries*/
    /*Spline on expiry*SQR(vol)*/
    for(k = 0; k < (2*intpInd)*freq; k++){
        for(j = 0; j < nbExpiryToMatch; j++){
            valueKnown[j] = 
            expiryToMatch[j]*
            (MTR((int)expiryToMatch[j]-1,k))*
            (MTR((int)expiryToMatch[j]-1,k));
        }
        
        DrlSplineInterp1dInit(
            expiryToMatch,        
            valueKnown,        
            nbExpiryToMatch,            
            1e12,    
            1e12,    
            valueInit);

        for(i = 0; i < intpInd*freq; i++){
            DrlSplineInterp1dInterp(
                expiryToMatch,              /* (I) array of X(i) (i=0,..,n-1) */
                valueKnown,                 /* (I) array of X(i) (i=0,..,n-1) */
                valueInit,                  /* (I) from SplineInit */
                nbExpiryToMatch,            /* (I) # of points */
                (double)(i+1)/(double)freq, /* (I) point to interpolate */
                &point);                    /* (O) interpolated value */

                MTR(i,k) = point;
                MTR(i,k) /= (double)(i+1)/(double)freq;
                MTR(i,k) = sqrt(MTR(i,k));

            }
        }

    /*Use Rebonato form to interpolate after 30 years*/
    /*Continuous derivative*/
	/*Multiply by 100 to take into account the change in vol format (*100)*/ 
    alpha = terminalVol;
    if(fabs(alpha)<1e-8)
        alpha = 0.01;
    if((alpha>0.5)||(alpha<-0.5))
        goto done;
    epsilon = terminalSlope;
    if((epsilon>100)||(epsilon<-100))
        goto done;
    for(i = 0; i < intpInd*freq; i++){
        slope = (MTR(i,intpInd*freq-1)-MTR(i,intpInd*freq-2));
        beta = 5*pow((double)intpInd,0.8)*slope/fabs(alpha) + epsilon;
        
        if(beta != 0){
        signBeta = beta/fabs(beta);
        }
        else
        signBeta = 1;    
    
        for(j = intpInd*freq; j <(2*intpInd*freq); j++){
            MTR(i,j) = (MTR(i,intpInd*freq-1)
                     +  signBeta*alpha 
                     -  signBeta*alpha*
                        exp(-fabs(beta)*
                        (pow((double)(j+1)/(double)freq,0.2)
                     -  pow((double)intpInd,0.2)))*
                        (1+signBeta*epsilon*
                        (pow((double)(j+1)/(double)freq,0.2)
                     -  pow((double)intpInd/(double)freq,0.2))));
        }
    }

    status = 0;

#undef MTR
done:
    if (status != 0) {
        DrlErrMsg("%s: failed (code %d)\n", routine, status);
    }

    /*Free vectors for interpolation */
    if(valueKnown != NULL){
    FREE_ARRAY(valueKnown);
    }

    if(valueInit != NULL){
    FREE_ARRAY(valueInit);
    }
    return(status);
}

