/******************************************************************************
 * Module:      Q3
 * Submodule:
 * File:        mcutil.c        
 * Function:    
 * Author:      Interest Rates DR
 * Revision:    $Header: $
 *****************************************************************************/
#include "math.h"
#include "crxq.h"
#include <crxflow/include/crxerror.h>


#if !defined(NULL)
#define NULL 0
#define UNDEF_NULL
#endif

/*f----------------------------------------------------------------------------
 *
 * ran2
 *
 * Generates uniform deviates.  Copied from Numerical Recipes.
 *
 */
double ran2(long *idum)
{

    double temp;

    static long iy = 0;
    static long iv[NTAB];
    static long idum2 = 123456789;
    long k;

    int  j;

    if (*idum <= 0 || !iy)
    {
        if (-(*idum) < 1) 
        {
            *idum = 1;
        }
        else
        {
            *idum = -(*idum);
        }
        idum2 = (*idum);

        for (j = NTAB + 7; j >= 0; j--)
        {
            k = (*idum) / IQ1;
            *idum = (IA1 * (*idum - k * IQ1)) - (IR1 * k);
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum) / IQ1;
    *idum = (IA1 * (*idum - k * IQ1)) - (IR1 * k);
    if (*idum < 0) *idum += IM1;
    k = idum2 / IQ2;
    idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
    if (idum2 < 0) idum2 += IM2;
    j = iy / NDIV;
    iy = iv[j] - idum2;
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp = AM * iy) > RNMX)
    {
        return RNMX;
    }
    else
    {
        return temp;
    }

}


/*f----------------------------------------------------------------------------
 *
 * ran2x2
 *
 * Generates pairs of uniform deviates.
 *
 */
void ran2x2(
    long   *seed,  /* (I/O) */
    double *X      /* (O)   */
    )
{

    X[0] = ran2(seed);
    X[1] = ran2(seed);

}


/*f----------------------------------------------------------------------------
 *
 * sobseq
 *
 * Generate Sobol sequences.  Copied from Numerical Recipes.
 *
 */
void sobseq(
    int    n,  /* (I) n<0 Initialize routine, n>0 generate n dim point */
    double x[] /* (O) */
    )
{

    unsigned long j,k,l; /* Changed from int. */
    unsigned long i,im,ipp;
    static double fac;
    static unsigned long in,ix[MAXDIM+1],*iu[MAXBIT+1];
    static unsigned long mdeg[MAXDIM+1]={0,1,2,3,3,4,4};
    static unsigned long ip[MAXDIM+1]={0,0,1,1,2,1,4};
    static unsigned long iv[MAXDIM*MAXBIT+1]={
        0,1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};

    if (n < 0) 
    {
        for (j=1,k=0;j<=MAXBIT;j++,k+=MAXDIM) iu[j] = &iv[k];
        for (k=1;k<=MAXDIM;k++) 
        {
            for (j=1;j<=mdeg[k];j++) iu[j][k] <<= (MAXBIT-j);
            for (j=mdeg[k]+1;j<=MAXBIT;j++) 
            {
                ipp=ip[k];
                i=iu[j-mdeg[k]][k];
                i ^= (i >> mdeg[k]);
                for (l=mdeg[k]-1;l>=1;l--) 
                {
                    if (ipp & 1) i ^= iu[j-l][k];
                    ipp >>= 1;
                }
                iu[j][k]=i;
            }
        }
        fac=1.0/(1L << MAXBIT);
        in=0;
    } 
    else 
    {
        im=in;
        for (j=1;j<=MAXBIT;j++) 
        {
            if (!(im & 1)) break;
            im >>= 1;
        }
        if (j > MAXBIT) DR_Error("sobseq: MAXBIT too small.");
        im=(j-1)*MAXDIM;
        for (k=1;k<=(unsigned int) MIN(n,MAXDIM);k++) 
        { /* Add cast */
            ix[k] ^= iv[im+k];
            x[k]=ix[k]*fac;
        }
        in++;
    }

}



/*f----------------------------------------------------------------------------
 *
 * BoxMuller
 *
 * Generates independent n(0,1) deviates of dimensionality dimX.  
 * Only dimX = (1,2) is supported.
 *
 */
int BoxMuller(
    long    genTyp, /* (I) Selects generator for uniform deviates. */
    long    dimX,   /* (I) Requested dimensionality of deviates.   */  
    long    numX,   /* (I) Requested number of deviates.           */
    double *X       /* (O) Deviates.                               */
    )
{

    double fac, rsq, v[2];     
    long   seed   = CRXQ_BV_RAN2_SEED;
    int    status = FAILURE;
    int    i,j;

    if (X == NULL) goto RETURN;

    if (numX < 1 || dimX < 1 || dimX > 2) goto RETURN;
     
    switch (genTyp) {

    case CRXQ_BV_RAN2:
            
        for (i = 0, j = 0; i < numX; i++) {
            do {
                ran2x2(&seed, v);
                v[0] = 2.0 * v[0] - 1.0;
                v[1] = 2.0 * v[1] - 1.0;
                rsq = (v[0] * v[0]) + (v[1] * v[1]);
            } while(rsq >= 1.0 || rsq == 0.0);
            fac = sqrt(-2.0 * log(rsq) / rsq);

            X[j] = v[0] * fac; j++;
            if (dimX > 1)
	    {
		X[j] = v[1] * fac; 
		j++;
	    }
        }
            
        break;

    case CRXQ_BV_SOBOL:

        sobseq(-1, v); /* Initialize. */

        for (i = 0, j = 0; i < numX; i++){
            do{
                sobseq(2, v);
                v[0] = 2.0 * v[0] - 1.0;
                v[1] = 2.0 * v[1] - 1.0;
                rsq = (v[0] * v[0]) + (v[1] * v[1]);
            } while (rsq >= 1.0 || rsq == 0.0);
            fac = sqrt(-2.0 * log(rsq) / rsq);
            
            X[j] = v[0] * fac; j++;
            if (dimX > 1) 
	    {
		X[j] = v[1] * fac; 
		j++;
	    }

        }

        break;

    default: goto RETURN;

    }
    
    status = SUCCESS;

 RETURN:

    return status;

}

/*f----------------------------------------------------------------------------
 *
 * Gauss
 *
 * Linear transformation of independent n(0,1) deviates to produce n(mu,sig)
 * deviates with linear correlation rho.  dimX * numX independent n(0,1)
 * deviates are required to produce numX deviates with dimensionality dimX 
 * and correlation structure rhoX.  Routine only supports dimX = (1, 2).
 * Ouput is formatted such that the dimension index varies most rapidly.
 *
 */
int Gauss(
    long    dimX, /* (I)   Requested dimensionality of deviates.                      */
    long    numX, /* (I)   Requested number of deviates of dimension dimX.            */
    double *sigX, /* (I)   Requested std devs for corr. deviates.                     */
    double *muX,  /* (I)   Requested means for corr. deviates.                        */
    double *rhoX, /* (I)   Requested linear correlation structure.                    */
    double *X     /* (I/O) dimX * numX n(0,1) deviates to transform on input.         */
    /*       numX corr. n(mu,sig) deviates of dimension dimX on output. */
    )
{

    double sig0 =0., sig1 = 0., mu0 = 0., mu1 = 0.;
    long   i;
    int    status = FAILURE;

    if (X == NULL || sigX == NULL || muX == NULL) goto RETURN;

    if (numX < 1 || dimX < 1 || dimX > 2) goto RETURN;

    /* NULL allowed for rhoX in 1d case. */
    if (rhoX == NULL && dimX != 1) goto RETURN;
        
    if (fabs(*rhoX) > 1.0) goto RETURN;

    mu0  = muX[0]; sig0 = sigX[0];    

    if (dimX > 1)
    {
        mu1 = muX[1]; sig1 = sigX[1];
    }
        
    for (i = 0; i < numX * dimX; ) {
        if (dimX > 1)
            X[i+1] = sig1*(*rhoX*X[i] + sqrt(1. - *rhoX**rhoX)*X[i+1]) + mu1;
        X[i] = sig0 * X[i] + mu0; 
        i += dimX;
    }

    status = SUCCESS;

 RETURN:

    return status;

}


#ifdef UNDEF_NULL
#undef NULL
#endif
