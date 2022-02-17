/****************************************************************************/
/*      Calculation routines customised for a turbo swap                    */
/****************************************************************************/
/*     FXTOOLS.C                                                         */
/****************************************************************************/


/*
$Header$
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cupslib.h"




/****** Hyb3_Gfunc ****************************************/
/**                                                     
        local volatility mapping function              
        local vol = \sigma_t * Hyb3_Gfunc(M,a1,a2,a3)/M     
                                                       
 *****************************************************/

int Hyb3_Gfunc(
      double *val, /**<(O)                     */
      double M,    /**<(I) S/F, moneyness      */
      double a1,   /**<(I) skew                */
      double a2,   /**<(I) curva               */
      double a3)   /**<(I) scale               */
{
    int     status = FAILURE;
    double  x, z, tanha3x, secha3x;

    if(M < TINY)
    {
        DR_Error("G function not defined for S/F < TINY.\n");
        goto RETURN;
    }    
    
    if(a3 < TINY)
    {
        *val = M;
        return(SUCCESS);
    }


    x        =   log(M);            /* log moneyness */
    z        =   exp(a3 * x);       /* ie, z = M^a3  */
    tanha3x  =   (z - 1. / z) / (z + 1. / z);
    secha3x  =   2. / (z + 1./z);

    *val     =   M * (1. + a1 * tanha3x + a2 * (1. - secha3x));

  
    status = SUCCESS;

RETURN:

    return(status);
}

/****** Hyb3_GDashfunc ************************************/
/**                                                     
        dg/dM                                          
        where g = Hyb3_Gfunc                                
                                                       
  ****************************************************/

int     Hyb3_GDashfunc(double *val, /**<(O)                 */
                  double M,    /**<(I) S/F, moneyness  */
                  double a1,   /**<(I) skew            */
                  double a2,   /**<(I) curva           */
                  double a3)   /**<(I) scale           */
{

    int     status = FAILURE;
    double  x, z, sinha3x, cosha3x;

    if(a3 < TINY)
    {
        *val = 1;
        status = SUCCESS;
        goto RETURN;
    }

    if(M < TINY)
    {
        (*val) = (1 -a1 + a2);
        status = SUCCESS;
        goto RETURN;
    }    

    if(Hyb3_Gfunc(val, 
             M,    
             a1,   
             a2,   
             a3) != SUCCESS){
        goto RETURN;    
    }


    x        =   log(M);            /* log moneyness */
    z        =   exp(a3 * x);       /* ie, z = M^a3  */
    sinha3x  =   0.5 * (z - 1./z);
    cosha3x  =   0.5 * (z + 1./z);

    /* derivative with respect to log(M), multiplied by 1/M */
    (*val)     /=   M;
    (*val)     +=   a3 * (a1 + a2 * sinha3x) / (cosha3x * cosha3x);

    status = SUCCESS;

RETURN:

    return(status);
}


/****** Hyb3_GDriftfunc    ******************************/
/**                                                     
        the functions needed for the calculation of 
        the drift process under FX smile
        dg/dM     where g = Hyb3_Gfunc     
        and K'(x).x (= x/G(x)) as well as both are needed 
        for the drift calculation (and comes for free)
                                                       
  ****************************************************/

int     Hyb3_GDriftfunc(double *val, /**<(O) dg/dM                */
                       double *val0, /**<(O) x/g                 */
                       double M,    /**<(I) S/F, moneyness  */
                       double a1,   /**<(I) skew            */
                       double a2,   /**<(I) curva           */
                       double a3)   /**<(I) scale           */
{


    int     status = FAILURE;
    double  x, z, gFunc,  sinha3x, tanha3x, secha3x;

    
    if(M < TINY)
    {
        (*val) = (1. -a1 + a2);
        (*val0) = 1./(1. - a1 + a2);
        status = SUCCESS;
        goto RETURN;
    }

    
    if(a3 < TINY)
    {
        *val = 1.;
        *val0 = 1.;
        return(SUCCESS);
    }

    x        =   log(M);            /* log moneyness */
    z        =   exp(a3 * x);       /* ie, z = M^a3  */
    sinha3x  =   0.5 * (z - 1./z);
    secha3x  =   2. / (z + 1./z);
    tanha3x  =   sinha3x * secha3x;

    gFunc    =   M * (1. + a1 * tanha3x + a2 * (1. - secha3x));

    /* derivative with respect to log(M), multiplied by 1/M */
    (*val)     =   gFunc/M;
    (*val)    +=   a3 * (a1 + a2 * sinha3x) * secha3x*secha3x; 

    /* produce the m/g */
    (*val0) = M/gFunc;
    status = SUCCESS;

RETURN:

    return(status);
}



/****Hyb3_Kfunc ***********************************************/
/**                                                         
     Calculates the general smile mapping function K(x)    
     ie, K(x) where K'(x) = 1/g(x)                         
     g(x) is the local volatility mapping function         
                                                           
 *********************************************************/

int     Hyb3_Kfunc(
               double  *val, /**<(O)                        */
               double  M,    /**<(I) S/F, moneyness         */
               double  a1,   /**<(I) skew                   */
               double  a2,   /**<(I) curva                  */
               double  a3)   /**<(I) scale                  */
{
    int    status = FAILURE;
    double x, z, alpha, beta, gamma, delta, kappa_up, kappa_low;
    double a, b, c; /*rational fraction decomposition*/
   

    /*? Initially a3 < TINY ?*/
    if(a3 < 0)
    {
        DR_Error("Smile function not defined for negative scale parameters.\n");
        goto RETURN;
    }

    /*special lognormal case*/
    if ( ((fabs(a1)<TINY) && (fabs(a2)<TINY)) || (fabs(a3)<TINY) ) {
        (*val) = log(M);
        status    = SUCCESS;
        goto RETURN;
    }

    alpha  =  1. + a1 + a2;
    beta   =  -2. * a2;
    gamma  =  1. - a1 + a2; 
    delta  =  beta * beta - 4. * alpha * gamma;
        
    x   =   log(M); /* log moneyness */
    z   =   exp(a3 * x);

    /*Special cases alpha = 0, gamma = 0*/
    /*alpha == 0 case*/
    if((fabs(alpha) < TINY))
    {
        /*Two subcases, beta = 0 or gamma = 0*/
        if((fabs(beta)<TINY))
        {   
            *val = z*z/2.0 - 0.5 + log(z);
            *val /= gamma;
            *val /= a3;
            return(SUCCESS);

        }else if((fabs(gamma)<TINY))
        {
            *val = z - 1/z;
            *val /= beta;
            *val /= a3;
            return(SUCCESS);
        }else
        {
            *val = (z - 1)/beta;
            *val -= gamma*(log(beta * z + gamma) - log(beta + gamma))/((beta*beta));
            *val += (log(z) - log(beta*z + gamma) + log(beta + gamma))/gamma;
            *val /= a3;
            return(SUCCESS);
        }
    }

    /*gamma == 0 case*/
    if((fabs(gamma) < TINY))
    {
        if (fabs(beta) < TINY)
        {
            *val = log(z) - 1.0/(2*z*z);
            *val += 0.5;
            *val /= (a3 * alpha);
            return (SUCCESS);
        }
        else
        {
            a    = -alpha/(beta*beta);
            b    = 1/beta;
            c    = (alpha*alpha)/(beta*beta);
            *val = (log(alpha*z + beta)- log(alpha + beta))/alpha;
            *val += a*log(z) - b/z + b + c*log(alpha*z + beta)/alpha - c*log(alpha + beta)/alpha;
            *val /= a3;
            return(SUCCESS);
        }
    }

    
    /* kappa  =  \int_1^{y^c} 1/(alpha * z^2 + beta * z+ gamma) dz */
    kappa_up  =  (delta > 0)?
                 log(  (2. * alpha * z + beta - sqrt(delta) ) / 
                       (2. * alpha * z + beta + sqrt(delta) ) ) / sqrt(delta) :
                 2*atan( (2. * alpha * z + beta) / sqrt(-delta) ) / sqrt(-delta);
    kappa_low =  (delta > 0)?
                 log(  (2. * alpha + beta - sqrt(delta) ) / 
                       (2. * alpha + beta + sqrt(delta) ) ) / sqrt(delta) :
                 2*atan( (2. * alpha + beta) / sqrt(-delta) ) / sqrt(-delta);

    /* deal with the case delta == 0 (we overwrite kappa_up and kappa_low) */
    if(fabs(delta) < TINY)
    {
        kappa_up  = -2/(2.0 * alpha * z + beta);
        kappa_low = -2/(2.0 * alpha     + beta);
    }
    

    if ( a3*x < 50.0 ) 
    {
        *val      =  a3*x / (a3 * gamma) + 
                    0.5 * (1. / alpha - 1. / gamma) * ( log(alpha*z*z+beta*z+gamma) - log(alpha+beta+gamma) ) / a3 -
                     0.5 *  beta * (1. / alpha + 1. /gamma )  * ( kappa_up - kappa_low ) / a3; 
    }
    else /* treat this case spearedly for speed and stability purpose (AlexK) */
    {
        /* approximate log(alpha*z*z) by log(alpha) + 2.0*(a3*x) + neglect */
        /* log( exp( a3*x) ) = a3*x : more efficient                       */
        *val      =  a3*x / (a3 * gamma) + 
                     0.5 * (1. / alpha - 1. / gamma) * ( log(alpha) + 2.0*a3*x - log(alpha+beta+gamma) ) / a3 -
                     0.5 *  beta * (1. / alpha + 1. /gamma )  * ( kappa_up - kappa_low ) / a3; 

    }

    status    = SUCCESS;

RETURN:

    return(status);
}


/****  Hyb3_KDashTimesX ***********************************/
/**                                                     
       calculates K'(x).x (= x/G(x))                   
                                                       
 *****************************************************/

int Hyb3_KDashTimesX(double *val,  /**<(O)                  */
                double  M,    /**<(I) S/F, moneyness   */
                double  a1,   /**<(I) scale            */
                double  a2,   /**<(I) curva            */
                double  a3)   /**<(I) scale            */
{
    
    int     status = FAILURE;


    if(M < TINY) 
    {
        (*val) = 1/(1 - a1 + a2);
        status = SUCCESS;
        goto RETURN;
    }


    if(Hyb3_Gfunc(val,M,a1,a2,a3) == FAILURE) goto RETURN;

    *val   =  M/ *val;    

    status =  SUCCESS;

RETURN:
    
    return(status);
}



/****  Hyb3_UpdateSmileCache ***********************************/
/**                                                     
       updated the tabulated K values needed for inverse K map 
       plus the functions needed for the drift calculation
                                                       
 *****************************************************/
int Hyb3_UpdateSmileCache(HYB3_FXSMILE_CACHE *smileCache , /* (O) contains the tabulated K and g values */
                        int               iCache,
                        long             *tMinTree,
                        long             *tMaxTree,
                        int               NbTP,
                        int               Ppy,
                        double            NbSigmaDBL,
                        double           *ExVol,
                        double           *ExMidNode,
                        int               t,
                        double            a1,
                        double            a2, 
                        double            a3)
{

    int status = FAILURE; 
    
    int i , nPts ;
    double maxKT, minKT;
    double minMony, maxMony, point;
    double inder, finder, K, h;
    long tMin, tMax, l;
    double maxStd = NbSigmaDBL + 5.0;  /* +5.0 arbitrary. To ensure coverage of tabulation */
    int loopCounter;


    /*First find the limits in time;                             */
    /*ie the indexes such that on the time interval [tMin, tMax] */
    /*The same smile index is in use                             */
    
    /*The method in use depends on the forward Index*/
    tMin = tMinTree[t] + 1;
    tMax = MIN( tMaxTree[t] + 1, NbTP);


    /* Find the limits in K dimension using an approximate method */
    /* This method is based on the fact that the Ex dim tree bound is initially  */
    /* computed using a maxStd input. We use this (+ a bit) to estimate the      */
    /* min (and max) K over the timepoints that share the same smile index       */
    minKT = ExMidNode[tMin] - maxStd*ExVol[tMin]*sqrt((tMin)/(double)Ppy);
    maxKT = ExMidNode[tMin] + maxStd*ExVol[tMin]*sqrt((tMin)/(double)Ppy);


    for(l = tMin; l<=tMax; l++)
    {
        minKT = MIN( ExMidNode[l] - maxStd*ExVol[l]*sqrt((l)/(double)Ppy), minKT);
        maxKT = MAX( ExMidNode[l] + maxStd*ExVol[l]*sqrt((l)/(double)Ppy), maxKT);
    }
    
    
    minMony = 1;
    maxMony = 1;
    
    /* Find min Moneyness by iteratively stepping down at sqrt(10) steps in log scale */      
    K = 0;
    loopCounter = 0;
    while(K>minKT)
    {
        minMony /= sqrt(10.0);
        if(Hyb3_Kfunc(&(K), 
                      minMony, 
                      a1,
                      a2,
                      a3)!=SUCCESS)
        {
            goto RETURN;
        }
        /* if Kfunc doesn't converge to < minKT in 10,000 iterations, return error */
        loopCounter++;
        if (loopCounter > 10000)
        {
            DR_Error("Unable to converge on min Moneyness calculation after 10,000 iterations - Smile parameters "
                "are stressing to model too much and should be adjusted\n");
            goto RETURN;
        }
    }
    /*Check K(min) < minKT*/
    if(K>minKT)
    {
        DR_Error("Pb in tabK, k(min) is not smaller than minKT");
        goto RETURN;
    }
    
    /* Find max Moneyness by iteratively stepping up at sqrt(10) steps in log scale */      
    K = 0;
    loopCounter = 0;
    while(K<maxKT)
    {
        maxMony *= sqrt(10.0);
        if(Hyb3_Kfunc(&(K), 
                      maxMony, 
                      a1,
                      a2,
                      a3) != SUCCESS )
        {
            goto RETURN;
        }
        /* if Kfunc doesn't converge to > maxKT in 10,000 iterations, return error */
        loopCounter++;
        if (loopCounter > 10000)
        {
            DR_Error("Unable to converge on max Moneyness calculation after 10,000 iterations - Smile parameters "
                "are stressing to model too much and should be adjusted\n");
            goto RETURN;
        }
    }
    /*Check K(min) < minKT*/
    if(K<maxKT)
    {
        DR_Error("Pb in tabK, k(max) is not bigger than maxKT");
        goto RETURN;
    }
    
    point = minMony;
    h = (log(maxMony) - log(minMony))/(double)NBSPLINE;
    
    maxMony *= exp(h*log(10.0));
    
    i = 0;
    while(point<=maxMony)
    {
        smileCache->X[iCache][i] = point;
        if(Hyb3_Kfunc(&(smileCache->K[iCache][i]), 
                      point, 
                      a1,
                      a2,
                      a3) != SUCCESS)
        {
            goto RETURN;
        }
        
        point *= exp(h);
        i++;
    }
    
    nPts = i;
    smileCache->nbPtSPL[iCache] = nPts;
    

    /*Calculation of initial and final derivative*/
    inder  = (smileCache->X[iCache][1] - smileCache->X[iCache][0])/
                        (smileCache->K[iCache][1]-smileCache->K[iCache][0]);

    finder = (smileCache->X[iCache][i-1] - smileCache->X[iCache][i-2])/
                        (smileCache->K[iCache][i-1]-smileCache->K[iCache][i-2]);
    

    /*Init spline interpolation for the K^(-1) function */
    if(SplineInterp1dInit(smileCache->K[iCache],     
                          smileCache->X[iCache],      
                          smileCache->nbPtSPL[iCache],            
                          inder,     
                          finder,    
                          smileCache->SPL[iCache]) != SUCCESS)
    {
        DR_Error("Pb in Spline Initialisation");
        goto RETURN;
    }


    /*Calculation of initial and final derivative*/
    inder  = (smileCache->K[iCache][1] - smileCache->K[iCache][0])/
                        (smileCache->X[iCache][1]-smileCache->X[iCache][0]);

    finder = (smileCache->K[iCache][i-1] - smileCache->K[iCache][i-2])/
                        (smileCache->X[iCache][i-1]-smileCache->X[iCache][i-2]);

    /* Init spline interpolation for the K function*/
    if(SplineInterp1dInit(smileCache->X[iCache],     
                          smileCache->K[iCache],      
                          smileCache->nbPtSPL[iCache],            
                          inder,     
                          finder,    
                          smileCache->SPL_Inv[iCache]) != SUCCESS)
    {
        DR_Error("Pb in Spline Initialisation");
        goto RETURN;
    }


    /* ------------------------------------------------------------------ */
    /* last step: update the remainder of the smile functions             */
    /* the g-function for the drift calculation as well as the K function */
    /* tabulate the gdash kdashTimesX and K function if necessary         */
    /* ------------------------------------------------------------------ */

    for (i= 0; i < nPts; ++i )
    {
        /* new: interpolateg-function and K function for change of drift as well */
        if(Hyb3_GDriftfunc(&smileCache->gd[iCache][i],
                           &smileCache->kdX[iCache][i],
                            smileCache->X[iCache][i], /* NormalisedEx */
                            a1,
                            a2,
                            a3) == FAILURE) goto RETURN;

    } /* loop over smile points i */
    
    

    /*Calculation of initial and final derivative and initialising spline */

    inder  = (smileCache->gd[iCache][1] - smileCache->gd[iCache][0])/
        (smileCache->X[iCache][1]-smileCache->X[iCache][0]);

    finder = (smileCache->gd[iCache][nPts-1] - smileCache->gd[iCache][nPts-2])/
        (smileCache->X[iCache][nPts-1]-smileCache->X[iCache][nPts-2]);
    

    if(SplineInterp1dInit(smileCache->X[iCache],     
                          smileCache->gd[iCache],      
                          smileCache->nbPtSPL[iCache],            
                          inder,     
                          finder,    
                          smileCache->gd_SPL[iCache] )!=SUCCESS)
    {
        DR_Error("Pb in Spline Initialisation for the g' interpolation");
        goto RETURN;
    }
    

    /*Calculation of initial and final derivative and initialising spline */

    inder  = (smileCache->kdX[iCache][1] - smileCache->kdX[iCache][0])/
        (smileCache->X[iCache][1]-smileCache->X[iCache][0]);

    finder = (smileCache->kdX[iCache][nPts-1] - smileCache->kdX[iCache][nPts-2])/
        (smileCache->X[iCache][nPts-1]-smileCache->X[iCache][nPts-2]);
    

    if(SplineInterp1dInit(smileCache->X[iCache],     
                          smileCache->kdX[iCache],      
                          smileCache->nbPtSPL[iCache],            
                          inder,     
                          finder,    
                          smileCache->kdX_SPL[iCache])!=SUCCESS)
    {
        DR_Error("Pb in Spline Initialisation for x/g function");
        goto RETURN;
    }



    status = SUCCESS;
RETURN:

    return status;
}



/****  Hyb3_BuildFXSmileCache***********************************/
/**                           
        allocate memory for smile and calculate the values                          
                                                       
 *****************************************************/
int Hyb3_BuildFXSmileCache(HYB3_TREE_DATA  *tree_data , /* (I/O) tree data */
                           double          *ExVol,
                           double          *ExMidNode )
{

    int status = FAILURE; 
    HYB3_FXSMILE_CACHE *smileCacheL = &(tree_data->FXsmileCache);

    int t, smileInUse, smileIdx;


    /* nothing to be done here, the caching has been done before */
    if (smileCacheL->isCached)
    {
        return SUCCESS;
    }


    /* "hack": the last smile data gives us number of possible smile sections    */
    smileCacheL->nbCachePts = tree_data->SmileIndex[ tree_data->NbTP] + 1; 


    /* allocate necessary memory for the smile function */
    smileCacheL->X       = (double**)DR_Matrix (DOUBLE, 0, smileCacheL->nbCachePts, 0 , NBSPLINE + 3);
    smileCacheL->K       = (double**)DR_Matrix (DOUBLE, 0, smileCacheL->nbCachePts, 0 , NBSPLINE + 3);
    smileCacheL->SPL     = (double**)DR_Matrix (DOUBLE, 0, smileCacheL->nbCachePts, 0 , NBSPLINE + 3);
    smileCacheL->SPL_Inv = (double**)DR_Matrix (DOUBLE, 0, smileCacheL->nbCachePts, 0 , NBSPLINE + 3);
    smileCacheL->gd      = (double**)DR_Matrix (DOUBLE, 0, smileCacheL->nbCachePts, 0 , NBSPLINE + 3);
    smileCacheL->gd_SPL  = (double**)DR_Matrix (DOUBLE, 0, smileCacheL->nbCachePts, 0 , NBSPLINE + 3);
    smileCacheL->kdX     = (double**)DR_Matrix (DOUBLE, 0, smileCacheL->nbCachePts, 0 , NBSPLINE + 3);
    smileCacheL->kdX_SPL = (double**)DR_Matrix (DOUBLE, 0, smileCacheL->nbCachePts, 0 , NBSPLINE + 3);

    smileCacheL->nbPtSPL = (int*)DR_Array( INT, 0, smileCacheL->nbCachePts );
    

    if ((smileCacheL->X       == NULL) ||
        (smileCacheL->K       == NULL) ||
        (smileCacheL->SPL     == NULL) ||
        (smileCacheL->SPL_Inv == NULL) ||
        (smileCacheL->nbPtSPL == NULL) ||
        (smileCacheL->gd      == NULL) ||
        (smileCacheL->gd_SPL  == NULL) ||
        (smileCacheL->kdX     == NULL) ||
        (smileCacheL->kdX_SPL == NULL) )
    {
        DR_Error("Unable to allocate cache memory for X, K, SPL, gd, gd_SPL, kdX and kd_SPL  arrays.");
        goto RETURN;
    }


    /* go through the tree and cache the necessary information */
    smileInUse = -999;    
    for ( t = 0; t < tree_data->NbTP; ++t)
    {
        /* force recalculation if new parameterisation */ 
        smileIdx = tree_data->SmileIndex[t];
        if(smileIdx > smileCacheL->nbCachePts )
        {
            DR_Error("Hyb3_BuildFXSmileCache: Internal Error: smileIdx exceeds the allowed cached size!");
            goto RETURN;
        }

        if ( smileIdx != smileInUse)
        {

            /* update the smile information */
            smileInUse = tree_data->SmileIndex[t];

            if (Hyb3_UpdateSmileCache( smileCacheL, 
                            smileIdx, 
                            tree_data->tMin, 
                            tree_data->tMax, 
                            tree_data->NbTP, 
                            tree_data->Ppy,
                            tree_data->NbSigmaDBL, 
                            ExVol, 
                            ExMidNode, 
                            t, 
                            tree_data->A1[smileIdx], 
                            tree_data->A2[smileIdx], 
                            tree_data->A3[smileIdx]) != SUCCESS )
            {
                DR_Error("Could not tabulate the K-function");
                goto RETURN;
            }


        } /* a new smile Index needs calculation */
    } /* for t */

    /* set cache flag as it has been calculated */
    smileCacheL->isCached = TRUE;

    status = SUCCESS;
RETURN:

    return status;
}






/*****Hyb3_FillGrid_2d ************************************************************************/
/**
 * This function populates the SpotEx slice with Ex values based on 
 * the smile parameters at SmileIdx[t-1]
 *                                      
 * NOTE: this function should only be called from lattice.c
 *       in particular, when calling it with t, we MUST also have
 *       a1 = tree_data->A1C[t-1]  or equivalently = tree_data->A1[SmileIdx[t-1]] etc
 *
 ***************************************************************************************/
int     Hyb3_FillGrid_2d(HYB3_TREE_DATA *tree_data,
                         TSLICE          SpotEx,
                         TSLICE          gDash,     /* (O) (only under smile) g-function at each node */
                         TSLICE          kdashTimesX,/* (O)(only under smile) x*K' at each node */
                         TSLICE          kVar,      /* (O) (only under smile+change) K_t(FX) for each node */
                         double         *FwdEx,
                         double         *ExMidNode,
                         double         *ExVol,
                         int             t,
                         int             updateGdash)
{
    int status = FAILURE;
    int offset;
    
    double  Zidx; /* first point in third dimension */
    double  X;    /* position in Gaussian space     */
    double  du, LengthJ;
    double  Jump1, Jump2;

    HYB3_FXSMILE_CACHE *smileCacheL = &(tree_data->FXsmileCache);

    int    Top1,      Bottom1;           /* Tree limits (IR domestic (nominally foreign)) */
    int   *Top2,     *Bottom2;           /* Tree limits (equity)                          */

    double  *SpotExL      = NULL;
    double  *gDashL       = NULL;
    double  *kdashTimesXL = NULL;
    double  *kVarL        = NULL;

    const double  FwdExL     = FwdEx[t];
    const double  ExMidNodeL = ExMidNode[t];
    double  InvK;
    double  a1, a2, a3;             /* smile parameters   */
    double  a1Next, a2Next, a3Next; /* for the g-function */
    int isLognormal, changeInK;  


    /*Limit values in X*/
    double minX, maxX;
    int     i,j, nSpline, nSplineN;
    double a, b, h;
    double * xa;
    double * ya;
    double * y2a;
    int jlo, jhi, jspl, jloN, jhiN;
    double deltaEx;

    /* New variables*/
    int smilePt, gFctPt;
    int NbTP;
    double expJump;

    double * gda, * gd2a; 
    double * kdxa, * kdx2a;
    double *Xg, *Yg, *k2Inv; 


    Top1      = tree_data->OutTop1[t];
    Top2      = tree_data->OutTop2[t];
    Bottom1   = tree_data->OutBottom1[t];
    Bottom2   = tree_data->OutBottom2[t];

    NbTP = tree_data->NbTP;

    /* Previous period coefficients [t-1] */
    LengthJ   = tree_data->LengthJ[t-1];
    du        = sqrt (JUMPCOEFF * LengthJ);
        
    Jump1     = tree_data->Aweight[1][t-1] * du;   
    Jump2     = tree_data->Aweight[2][t-1] * du;   
    expJump   = exp(Jump2);

    /* get the smile parameter for t-1 (!) */
    a1 = tree_data->A1[tree_data->SmileIndex[t-1] ];
    a2 = tree_data->A2[tree_data->SmileIndex[t-1] ];
    a3 = tree_data->A3[tree_data->SmileIndex[t-1] ];

    a1Next = tree_data->A1[tree_data->SmileIndex[t] ];
    a2Next = tree_data->A2[tree_data->SmileIndex[t] ];
    a3Next = tree_data->A3[tree_data->SmileIndex[t] ];


    /* K_t(FX) - K_t-1(FX) only needs calculation if index changes */
    changeInK = ( ( fabs(a1Next-a1) > TINY) || (fabs(a2Next-a2) > TINY) 
                  || (fabs(a3Next-a3) > TINY) );


    isLognormal = ((fabs(a1)<TINY) && (fabs(a2)<TINY)) || (fabs(a3)<TINY);


    /* build a smile cache */
    if (!smileCacheL->isCached)
    {
        Hyb3_BuildFXSmileCache(tree_data,  ExVol, ExMidNode);
    }


    /*Inserting the inverse function in fill grid*/
    /*special lognormal case*/
    if( isLognormal ){


        /* constant node width in mapped space */
        for(i = Bottom1 ; i <= Top1; i++)
        {
            /* prepare to address Asset2 */
            offset  = Hyb3_Node_Offset(2,i,0,t,tree_data);
            SpotExL = SpotEx + offset;
            
            Zidx  = ExMidNodeL + i * Jump1 + Jump2 * Bottom2[i];            
            X     = Zidx; 
        
            SpotExL[Bottom2[i]] = FwdExL * exp(X);

            for (j = (Bottom2[i]+1); j <= Top2[i]; j++)
            {
                SpotExL[j] = SpotExL[j-1] * expJump;
            }
        }/*for i */

        /* in case the previous calculation was smiley 
           and we need to know the change of the K-function for the
           drift calculation, then continue and calculate g-drift and K-drift */
        if (!changeInK || !updateGdash) /* "illegal" goto, try to avoid this */
        {
            status    = SUCCESS;
            goto RETURN;
        }
    }

    /*Update Smile index previously in use*/
    smilePt = tree_data->SmileIndex[t-1];

    nSpline   = smileCacheL->nbPtSPL[smilePt];
    minX      = smileCacheL->K[smilePt][0];
    maxX      = smileCacheL->K[smilePt][nSpline-1];

    /* Update SplineVectors  */ 
    xa        = smileCacheL->K[smilePt];
    ya        = smileCacheL->X[smilePt];
    y2a       = smileCacheL->SPL[smilePt];

    /* in order to use index starting from 1 (for quick search algorithm) */
    xa --;
    ya --;
    y2a --;

    /* for the g' and k' * X function we need smile at t, not t-1 (!))  */
    if (updateGdash)
    { 
        gFctPt    = tree_data->SmileIndex[t];
        nSplineN  = smileCacheL->nbPtSPL[gFctPt];

        gda  = smileCacheL->gd[gFctPt];  gda--;
        gd2a = smileCacheL->gd_SPL[gFctPt]; gd2a--;
        kdxa  = smileCacheL->kdX[gFctPt]; kdxa--;
        kdx2a = smileCacheL->kdX_SPL[gFctPt]; kdx2a--;
        
        Xg  = smileCacheL->X[gFctPt]; Xg--; /* use the smile point of t not t-1!! */
        Yg  = smileCacheL->K[gFctPt]; Yg--; /* use the smile point of t not t-1!! */

        /* for the K-function fit */
        k2Inv = smileCacheL->SPL_Inv[gFctPt]; k2Inv--;

    }

    /* constant node width in mapped space */
    for(i = Bottom1 ; i <= Top1; i++)
    {
        /* prepare to address Asset3 */
        offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
        SpotExL  = SpotEx + offset;
        gDashL   = gDash + offset;
        kdashTimesXL  = kdashTimesX + offset;
        kVarL    = kVar + offset;
        
        Zidx = ExMidNodeL + i * Jump1 + Jump2 * Bottom2[i];
        
        X = Zidx; 

        deltaEx = Jump2*(Top2[i] - Bottom2[i]);
        
        if ( X < minX )
        {
            DR_Error("x smaller than tree_data->K[0].\n");
            goto RETURN;
        }

        if ( (X+deltaEx) > maxX)
        {
            DR_Error("x bigger than tree_data->K[n-1].\n");
            goto RETURN;
        }
    
        jlo = 1;
        jloN = 1;
        for (j = Bottom2[i]; j <= Top2[i]; j++)
        {
            jhi = nSpline;
            while (jhi-jlo > 1) 
            {
                jspl=(jhi+jlo) >> 1;
                if (xa[jspl] > X) jhi=jspl;
                else jlo=jspl;
            }
            h=xa[jhi]-xa[jlo];
            if (fabs (h) < TINY) 
            {
                DR_Error("bad xa input.\n");
                return(FAILURE);
            }

            a=(xa[jhi]-X)/h;
            b=(X-xa[jlo])/h;
            InvK=a*ya[jlo]+b*ya[jhi]+((a*a*a-a)*y2a[jlo]+
                   (b*b*b-b)*y2a[jhi])*(h*h)/6.0;


            SpotExL[j] =  FwdExL * InvK; 
            X          += Jump2;

            /* update the GDash and KDashTimesX only when required */
            if (updateGdash)
            {
                /* in case we have different smile indices, we need to interpolate further */
                if (smilePt != gFctPt) 
                {
                    jhiN = nSplineN;
                    while (jhiN-jloN > 1) 
                    {
                        jspl=(jhiN+jloN) >> 1;
                        if (Xg[jspl] > InvK ) jhiN = jspl;
                        else jloN = jspl;
                    }
                }
                else
                {
                    jloN = jlo;
                    jhiN = jhi;
                }
                h=Xg[jhiN]-Xg[jloN];
                if (fabs (h) < TINY*1e-5 ) /* needs smaller error margin in this space ! */
                {
                    /* when FX points are lying to closely together or even 0.0 */
                    /* then simply use the lower point in the tabulated vector  */
                    a =1.0;
                }
                else
                {
                    a=(Xg[jhiN]-InvK)/h;
                }
                b=1.-a;

                /* avoid spline in out of bound areas */
                if ( (a > 1.1) || (a < -0.1) ) 
                {
                    /*DR_Error("FX spline out of bounds: possibly too big change in"
                             "smile parameter between t=%d and %d\n", t-1,t);
                    goto RETURN; */
                    if(Hyb3_GDriftfunc(&(gDashL[j]),
                                   &(kdashTimesXL[j]),
                                   InvK, /* NormalisedEx */
                                   a1Next,
                                   a2Next,
                                   a3Next) == FAILURE) goto RETURN;

                    if (changeInK)
                    {   
                        if(Hyb3_Kfunc(&(kVarL[j]),
                                   InvK, /* NormalisedEx */
                                   a1Next,
                                   a2Next,
                                   a3Next) == FAILURE) goto RETURN;

                    } /* changeInK */
                }
                else
                {
                    gDashL[j] = a*gda[jloN]+b*gda[jhiN]+((a*a*a-a)*gd2a[jloN]
                                +(b*b*b-b)*gd2a[jhiN])*(h*h)/6.0;
                    kdashTimesXL[j] = a*kdxa[jloN]+b*kdxa[jhiN]
                                +((a*a*a-a)*kdx2a[jloN]+(b*b*b-b)*kdx2a[jhiN])*(h*h)/6.0;
                                
                    if (changeInK)
                    {   
                        kVarL[j] = a*Yg[jloN]+b*Yg[jhiN]
                                +((a*a*a-a)*k2Inv[jloN]+(b*b*b-b)*k2Inv[jhiN])*(h*h)/6.0;
                    } /* changeInK */
                }
            } /* updateGDash */        
        }
    }/*for i */

    status = SUCCESS;
    
RETURN:

    return(status);

} /* Hyb3_FillGrid_2d */





/*****Hyb3_FillGrid_3d ************************************************************************/
/**
 * This function populates the SpotEx slice with Ex values based on 
 * the smile parameters at SmileIdx[t-1]
 *                                      
 * NOTE: this function should only be called from lattice.c
 *       in particular, when calling it with t, we MUST also have
 *       a1 = tree_data->A1C[t-1]  or equivalently = tree_data->A1[SmileIdx[t-1]] etc
 *
 ***************************************************************************************/
int     Hyb3_FillGrid_3d(HYB3_TREE_DATA *tree_data,
                         TSLICE          SpotEx,    /* (O) spot FX data for each node point */
                         TSLICE          gDash,     /* (O) (only under smile) g-function at each node */
                         TSLICE          kdashTimesX,/* (O)(only under smile) x*K' at each node */
                         TSLICE          kVar,      /* (O) (only under smile+change) K_t(FX) for each node */
                         double         *FwdEx,
                         double         *ExMidNode,
                         double         *ExVol,
                         int             t,
                         int             updateGdash)
{
    int status = FAILURE;
    int offset;
    
    HYB3_FXSMILE_CACHE *smileCacheL = &(tree_data->FXsmileCache);

    double  Zidx; /* first point in third dimension */
    double  X;    /* position in Gaussian space     */
    double  du, LengthJ;
    double  Jump3, Jump4, Jump5;

    int    Top1,      Bottom1;           /* Tree limits (IR foreign)       */
    int   *Top2,     *Bottom2;           /* Tree limits (IR domestic)      */
    int  **Top3,    **Bottom3;           /* Tree limits (third asset)      */

    double  *SpotExL      = NULL;
    double  *gDashL       = NULL;
    double  *kdashTimesXL = NULL;
    double  *kVarL        = NULL;

    const double  FwdExL     = FwdEx[t];
    const double  ExMidNodeL = ExMidNode[t];
    double  InvK;                   /* normalised FX/Eq */
    double  a1, a2, a3;             /* smile parameters   */
    double  a1Next, a2Next, a3Next; /* for the g-function */
    int isLognormal, changeInK;  

    /*Limit values in X*/
    double minX, maxX;
    int     i,j,k, nSpline, nSplineN;
    double a, b, h;
    double * xa;
    double * ya;
    double * y2a;
    int klo, khi, kspl, kloN, khiN;
    double deltaEx;

    /*New variables*/
    int smilePt, gFctPt;
    int NbTP;
    double expJump;

    double * gda, * gd2a; 
    double * kdxa, * kdx2a;
    double *Xg, *Yg, *k2Inv; 


    
    Top1      = tree_data->OutTop1[t];
    Top2      = tree_data->OutTop2[t];
    Top3      = tree_data->OutTop3[t];
    Bottom1   = tree_data->OutBottom1[t];
    Bottom2   = tree_data->OutBottom2[t];
    Bottom3   = tree_data->OutBottom3[t];

    NbTP = tree_data->NbTP;

    /* Previous period coefficients [t-1] */
    LengthJ   = tree_data->LengthJ[t-1];
    du        = sqrt (JUMPCOEFF * LengthJ);
        
    Jump3     = tree_data->Aweight[3][t-1] * du;   
    Jump4     = tree_data->Aweight[4][t-1] * du;   
    Jump5     = tree_data->Aweight[5][t-1] * du; 
    expJump   = exp(Jump5);


    /* get the smile parameter for t-1 (!) */
    a1 = tree_data->A1[tree_data->SmileIndex[t-1] ];
    a2 = tree_data->A2[tree_data->SmileIndex[t-1] ];
    a3 = tree_data->A3[tree_data->SmileIndex[t-1] ];

    a1Next = tree_data->A1[tree_data->SmileIndex[t] ];
    a2Next = tree_data->A2[tree_data->SmileIndex[t] ];
    a3Next = tree_data->A3[tree_data->SmileIndex[t] ];

    /* K_t(FX) - K_t-1(FX) only needs calculation if index changes */
    changeInK = ( ( fabs(a1Next-a1) > TINY) || (fabs(a2Next-a2) > TINY) 
                  || (fabs(a3Next-a3) > TINY) );


    isLognormal = ((fabs(a1)<TINY) && (fabs(a2)<TINY)) || (fabs(a3)<TINY);

/*    fprintf(fp, "%d %d %d %d %d %d %d %d %d  %d %d %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f\n"
        ,t, updateGdash, changeAtPrevT,changeAtT,changeInK,isLognormal,smileCacheL->isCached,
        smileCacheL->SmlIndPInUse,tree_data->SmileIndex[t-1],smileCacheL->gFuncInUse,tree_data->SmileIndex[t],
        a1,a2,a3,a1Next,a2Next,a3Next);*/
    if (!smileCacheL->isCached)
    {
        Hyb3_BuildFXSmileCache(tree_data,  ExVol, ExMidNode);
    }

    /*Inserting the inverse function in fill grid*/
    /*special lognormal case*/
    if( isLognormal )
    {
        /* constant node width in mapped space */
        for(i = Bottom1 ; i <= Top1; i++)
        {
            for(j = Bottom2[i]; j<= Top2[i]; j++)
            {
                /* prepare to address Asset3 */
                offset  = Hyb3_Node_Offset(3,i,j,t,tree_data);
                SpotExL = SpotEx + offset;
                
                Zidx  = ExMidNodeL + i * Jump3 + j * Jump4;
            
                Zidx += Jump5 * Bottom3[i][j];
            
                X = Zidx; 
            
                SpotExL[Bottom3[i][j]] =  FwdExL * exp(X);

                for (k = (Bottom3[i][j]+1); k <= Top3[i][j]; k++)
                {
                  SpotExL[k] = SpotExL[k-1] * expJump;
                }
            }/* for j */
        }/*for i */

        /* in case the previous calculation was smiley 
           and we need to know the change of the K-function for the
           drift calculation, then continue and calculate g-drift and K-drift */
        if (!changeInK || !updateGdash) /* "illegal" goto, try to avoid this */
        {
            status    = SUCCESS;
            goto RETURN;
        }
    }


    /*Update Smile index previously in use*/
    smilePt = tree_data->SmileIndex[t-1];

    nSpline   = smileCacheL->nbPtSPL[smilePt];
    minX      = smileCacheL->K[smilePt][0];
    maxX      = smileCacheL->K[smilePt][nSpline-1];

    /* Update SplineVectors  */ 
    xa        = smileCacheL->K[smilePt];
    ya        = smileCacheL->X[smilePt];
    y2a       = smileCacheL->SPL[smilePt];

    /* in order to use index starting from 1 (for quick search algorithm) */
    xa --;
    ya --;
    y2a --;


    /* for the g' and k' * X function we need smile at t, not t-1 (!))  */
    if (updateGdash)
    { 
        gFctPt    = tree_data->SmileIndex[t];
        nSplineN  = smileCacheL->nbPtSPL[gFctPt];

        gda  = smileCacheL->gd[gFctPt];  gda--;
        gd2a = smileCacheL->gd_SPL[gFctPt]; gd2a--;
        kdxa  = smileCacheL->kdX[gFctPt]; kdxa--;
        kdx2a = smileCacheL->kdX_SPL[gFctPt]; kdx2a--;
        
        Xg  = smileCacheL->X[gFctPt]; Xg--; /* use the smile point of t not t-1!! */
        Yg  = smileCacheL->K[gFctPt]; Yg--; /* use the smile point of t not t-1!! */

        /* for the K-function fit */
        k2Inv = smileCacheL->SPL_Inv[gFctPt]; k2Inv--;

    }



    /* constant node width in mapped space */
    for(i = Bottom1 ; i <= Top1; i++)
    {
        for(j = Bottom2[i]; j<= Top2[i]; j++)
        {
            /* prepare to address Asset3 */
            offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
            SpotExL  = SpotEx + offset;
            gDashL   = gDash + offset;
            kdashTimesXL  = kdashTimesX + offset;
            kVarL    = kVar + offset;
            
            Zidx = ExMidNodeL + i * Jump3 + j * Jump4;
            
            Zidx += Jump5 * Bottom3[i][j];
            
            X = Zidx; 

            deltaEx = Jump5*(Top3[i][j] - Bottom3[i][j]);
            
            if ( X < minX)
            {
                DR_Error("x smaller than tree_data->K[0].\n");
                goto RETURN;
            }

            if ((X+deltaEx) > maxX)
            {
                DR_Error("x bigger than tree_data->K[n-1].\n");
                goto RETURN;
            }
        
            klo  = 1;
            kloN = 1;
            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
            {
                khi = nSpline;
                while (khi-klo > 1) 
                {
                    kspl=(khi+klo) >> 1;
                    if (xa[kspl] > X) khi=kspl;
                    else klo=kspl;
                }
                h=xa[khi]-xa[klo];
                if (fabs (h) < TINY) 
                {
                    DR_Error("bad xa input.\n");
                    return(FAILURE);
                }

                a=(xa[khi]-X)/h;
                b=(X-xa[klo])/h;
                InvK=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+
                            (b*b*b-b)*y2a[khi])*(h*h)/6.0;
 
                SpotExL[k] = FwdExL * InvK; 
                X         += Jump5;


                /* update the GDash and KDashTimesX only when required */
                if (updateGdash)
                {
                    /* in case we have different smile indices, we need to interpolate further */
                    if (smilePt != gFctPt) 
                    {
                        khiN = nSplineN;
                        while (khiN-kloN > 1) 
                        {
                            kspl=(khiN+kloN) >> 1;
                            if (Xg[kspl] > InvK ) khiN = kspl;
                            else kloN = kspl;
                        }
                    }
                    else
                    {
                        kloN = klo;
                        khiN = khi;
                    }
                    h=Xg[khiN]-Xg[kloN];
                    if (fabs (h) < TINY*1e-5 ) /* needs smaller error margin in this space ! */
                    {
                        /* when FX points are lying to closely together or even 0.0 */
                        /* then simply use the lower point in the tabulated vector  */
                        a =1.0;
                    }
                    else
                    {
                        a=(Xg[khiN]-InvK)/h;
                    }
                    b=1.-a;

                    /* avoid spline in out of bound areas */
                    if ( (a > 1.1) || (a < -0.1) ) 
                    {
                        /*DR_Error("FX spline out of bounds: possibly too big change in"
                                 "smile parameter between t=%d and %d\n", t-1,t);
                        goto RETURN; */
                        if(Hyb3_GDriftfunc(&(gDashL[k]),
                                       &(kdashTimesXL[k]),
                                       InvK, /* NormalisedEx */
                                       a1Next,
                                       a2Next,
                                       a3Next) == FAILURE) goto RETURN;

                        if (changeInK)
                        {   
                            if(Hyb3_Kfunc(&(kVarL[k]),
                                       InvK, /* NormalisedEx */
                                       a1Next,
                                       a2Next,
                                       a3Next) == FAILURE) goto RETURN;
                        } /* changeInK */
                    }
                    else
                    {
                        gDashL[k] = a*gda[kloN]+b*gda[khiN]+((a*a*a-a)*gd2a[kloN]
                                +(b*b*b-b)*gd2a[khiN])*(h*h)/6.0;
                        kdashTimesXL[k] = a*kdxa[kloN]+b*kdxa[khiN]
                                +((a*a*a-a)*kdx2a[kloN]+(b*b*b-b)*kdx2a[khiN])*(h*h)/6.0;
                                                                               
                        if (changeInK)
                        {
                            kVarL[k] = a*Yg[kloN]+b*Yg[khiN]
                                +((a*a*a-a)*k2Inv[kloN]+(b*b*b-b)*k2Inv[khiN])*(h*h)/6.0;
                        } /* changeInK */
                    }
                } /* updateGDash */
            } /* for k */
        }/* for j */
    }/*for i */

    status = SUCCESS;
    
RETURN:

    return(status);
    
} /* Hyb3_FillGrid_3d */



/* ------------------------------------------------------------------------- */
/* FillGridEqFwd_2d: calculate the forward prices using the dividiend stream */
/* ------------------------------------------------------------------------- */
int     Hyb3_FillGridEqFwd_2d(HYB3_TREE_DATA  *tree_data,
                              HYB3_DEV_DATA   *dev_data,
                              int        t)
{
    int status = FAILURE;
    int offset;
    
    register double   FwdFactor=0.0;   /* Ratio of fwd eq to spot eq prices  */

    int    Top1,      Bottom1;           /* Tree limits (IR domestic)       */
    int   *Top2,     *Bottom2;           /* Tree limits (equity      )      */
 
    int     i,j;

    double  *Discount_1D[3] = {NULL, NULL, NULL};
 
    double  *SpotEqL;
    double  *FwdEqL;

    int     CvDiscF;
    
    double  Length = tree_data->Length[t];

    Top1    = tree_data->OutTop1[t];
    Top2    = tree_data->OutTop2[t];
    Bottom1 = tree_data->OutBottom1[t];
    Bottom2 = tree_data->OutBottom2[t];

    CvDiscF = tree_data->CvDisc[0];

    offset = Hyb3_Node_Offset(1,0,0,t,tree_data);

    Discount_1D[0] = dev_data->Discount_1D[0] + offset;
    Discount_1D[1] = dev_data->Discount_1D[1] + offset;
    Discount_1D[2] = dev_data->Discount_1D[2] + offset;

    for(i = Bottom1 ; i <= Top1; i++)
    {
        /*From outside */
        FwdFactor = pow(Discount_1D[CvDiscF][i],
                        -tree_data->NodeSettleTime[t]/Length);

        offset = Hyb3_Node_Offset(2,i,0,t,tree_data);

        SpotEqL = dev_data->EqSpot + offset;
        FwdEqL  = dev_data->EqFwd  + offset;
        
        for(j = Bottom2[i]; j<= Top2[i]; j++)
        {
            FwdEqL[j] = SpotEqL[j]*FwdFactor;
        }
    }/*for i */

    status = SUCCESS;
    return(status);

} /* FillGridEqFwd_2d */



/* ------------------------------------------------------------------------- */
/* FillGridEqFwd_3d: calculate the forward prices using the dividiend stream */
/* ------------------------------------------------------------------------- */
int     Hyb3_FillGridEqFwd_3d(HYB3_TREE_DATA  *tree_data,
                              HYB3_DEV_DATA   *dev_data,
                              int        t)
{
    int status = FAILURE;
    int offset;
    
    register double   FwdFactor=0.0;   /* Ratio of fwd eq to spot eq prices  */

    int    Top1,      Bottom1;           /* Tree limits (IR foreign)       */
    int   *Top2,     *Bottom2;           /* Tree limits (IR domestic)      */
    int  **Top3,    **Bottom3;           /* Tree limits (third asset)      */

    int     i,j,k;

    int   EquityOn  = !(tree_data->TreeType == TTYPE_FX2IR);
    int   EquityDom = EquityOn && !(tree_data->TreeType == TTYPE_EQF2IR);

    double  *Discount_1D[3] = {NULL, NULL, NULL};
    double  *Discount_2D[3] = {NULL, NULL, NULL};
 
    double  *SpotEqL;
    double  *FwdEqL;

    int     CvDiscD;
    int     CvDiscF;
    
    double  Length = tree_data->Length[t];

    Top1    = tree_data->OutTop1[t];
    Top2    = tree_data->OutTop2[t];
    Top3    = tree_data->OutTop3[t];
    Bottom1 = tree_data->OutBottom1[t];
    Bottom2 = tree_data->OutBottom2[t];
    Bottom3 = tree_data->OutBottom3[t];

    CvDiscD = tree_data->CvDisc[1];
    CvDiscF = tree_data->CvDisc[0];

    /*From outside */

    offset = Hyb3_Node_Offset(1,0,0,t,tree_data);

    Discount_1D[0] = dev_data->Discount_1D[0] + offset;
    Discount_1D[1] = dev_data->Discount_1D[1] + offset;
    Discount_1D[2] = dev_data->Discount_1D[2] + offset;

    for(i = Bottom1 ; i <= Top1; i++)
    {
        offset = Hyb3_Node_Offset(2,i,0,t,tree_data);

        Discount_2D[0] = dev_data->Discount_2D[0] + offset;
        Discount_2D[1] = dev_data->Discount_2D[1] + offset;
        Discount_2D[2] = dev_data->Discount_2D[2] + offset;

        for(j = Bottom2[i]; j<= Top2[i]; j++)
        {
            if (EquityDom)
                FwdFactor = pow(Discount_2D[CvDiscD][j], 
                                -tree_data->NodeSettleTime[t]/Length);
            else 
                FwdFactor = pow(Discount_1D[CvDiscF][i], 
                                -tree_data->NodeSettleTime[t]/Length);

            offset  = Hyb3_Node_Offset(3,i,j,t,tree_data);

            SpotEqL = dev_data->EqSpot + offset;
            FwdEqL  = dev_data->EqFwd  + offset;
            
            for (k = (Bottom3[i][j]); k <= Top3[i][j]; k++)
            {
                FwdEqL[k] = SpotEqL[k]*FwdFactor;
            }
        }/* for j */
    }/*for i */

    status = SUCCESS;
    return(status);

} /* FillGridEqFwd_3d */

