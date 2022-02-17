/*********************************************************************
 sobol.c
 Original source code from SPRNG Sobol sequence generator.
 See http://sprng.cs.fsu.edu/
 
 Modified for use with Supercube library - M.Huq, Credit DR
 demoted four arrays to local variables that have to be passed through
 the various calls. See soboldat.c. This was done so as to ensure that
 we correctly can have multiple Sobol sequences if needed.
 *********************************************************************/
#include "edginc/coreConfig.hpp"
#include "sobol.h"

#include <stdio.h>
#include <stdlib.h>


CORE_BEGIN_NAMESPACE

extern RNG_DLL int Sobol_edi[DIMENSION];
extern RNG_DLL int Sobol_initdat[MAXORDPOLY][DIMENSION];
extern RNG_DLL int Sobol_deg[DIMENSION];

/***********************************************************
 Initialize the parameters of sobol' sequences
 generating the vc in each dimension
 dimension specifies the number of dimension needed
 0 if successful, -1 if failed
***********************************************************/
int sobinit(int dimension,
            int **Sobol_gen,
            int *Sobol_xn,
            int *Sobol_n,
            int *Sobol_nstream /*pass by reference */
           )
{
    /*extern int Sobol_gen[MAXBIT][DIMENSION];
    extern int Sobol_xn[DIMENSION];
    extern int Sobol_n[DIMENSION];
    extern int Sobol_nstream;*/

    int i, j, k;
    int dirnum; /* direction number */

    if (dimension>DIMENSION) {
        fprintf(stderr,"Too high dimension for initializing the Sobol' sequence\n");
        return(-1);
    }

    /* set the number of quasirandom number streams */
    (*Sobol_nstream)=dimension;

    /* read in all of the initial value
       the initial value m1, m2, ..., md can be selected freely
       provided that each mi is odd and mi<2^i */
    for (i=0; i<MAXORDPOLY; i++)
        for (j=0;j<dimension; j++)
            Sobol_gen[i][j]=Sobol_initdat[i][j];

    /* dimension 0 is generated specially, 2^0 ~ 2^(MAXBIT-1) in each*/
    for (i=0;i<MAXBIT;i++)
        Sobol_gen[i][0]=1<<(MAXBIT-i-1);

    /* generate the other vc in other dimensions */
    for (i=1;i<dimension;i++) {
        /* set up vj, where j from 0 to d-1, vj=2^j*m[i-j]
           the reason the set it as 2^MAX-j-1 is to set the most
           significant bit when translate it to decimal number */
        for (j=0; j<Sobol_deg[i]; j++)
            Sobol_gen[j][i]*=1<<(MAXBIT-j-1);

        /* set up vj for j from d to MAXBIT, using the direction
           number of the primitive polynomial */
        for (j=Sobol_deg[i];j<MAXBIT;j++) {
            /* Compute v[i-d]EOR[v[i-d]/2^d] */
            Sobol_gen[j][i]=Sobol_gen[j-Sobol_deg[i]][i];
            Sobol_gen[j][i]^=Sobol_gen[j][i]>>Sobol_deg[i];

            /* get the direction number */
            dirnum=Sobol_edi[i];

            /* compute the vi by vi=a[1]v[i-1]EORa[2]v[i-2]EOR...EORa[d-1]v[i-d+1] */
            for (k=Sobol_deg[i]-1; k>0; k--) {
                if ((dirnum&1)!=0)
                    Sobol_gen[j][i]^=Sobol_gen[j-k][i];
                dirnum>>=1;
            }
        }
    }

    /* in each dimension the initial value x[0] is 0, n is 0 at each
       dimension in the beginning*/
    for (i=0;i<dimension;i++) {
        Sobol_xn[i]=0;
        Sobol_n[i]=0;
    }

    return(0);
}

#ifndef IRIX
/* __attribute__ ((aligned (16))) */ int expon[4] = {0, MAXBIT, 0, MAXBIT};
/* __attribute__ ((aligned (16))) */
double expo[2] = {1./(double)(1<<MAXBIT), 1./(double)(1<<MAXBIT)};
#endif

/***********************************************************
  Everytime by calling sobvect, quasirandom vector will return
  parameter specifies the dimension of the vector
  parameter x is the vector to be generated
***********************************************************/
int sobvect(int dimension,
            double *x,
            int **Sobol_gen,
            int *Sobol_xn,
            int *Sobol_n,
            int *Sobol_nstream /*pass by reference */
           )
{
    /*  extern int Sobol_gen[MAXBIT][DIMENSION];
    extern int Sobol_xn[DIMENSION];
    extern int Sobol_n[DIMENSION];
    extern int Sobol_nstream; */

    int i, j/*, dim4*/;
    int firstzero; /* the first zero bit in n */

    if (dimension>DIMENSION||dimension>(*Sobol_nstream)) {
        fprintf(stderr,"Too huge dimension or the Sobol sequence hasn't been initialized for such a high dimension yet\n");
        return(-1);
    }

    for (i=0;i<dimension;i++) {
        // find the first zero bit in n stored in Sobol_n
#if !defined(IRIX) && 0
        j = ~Sobol_n[i];
asm ("bsf %1, %%ecx\n\t movl %%ecx, %0\n" : "=r"(firstzero) : "r"(j) : "%ecx");
#else
        // FIXME: double check the output with another generator
        // FIXME: instead of asm, use firstzero = ffs(~Sobol_n[i]) - 1; and check for -1
        firstzero=0;
        for (j=Sobol_n[i];(j&1)!=0;j>>=1)
            firstzero++;
#endif
        // x[n+1]=x[n]EORvc
        Sobol_xn[i]^=Sobol_gen[firstzero][i];

        // divided by 2^MAXBIT
#if 0

        m.d1 = (double)Sobol_xn[i];
        new_exp = (m.i2[1] >> 20) - 30;
        m.i2[1] = (new_exp << 20) | (m.i2[1] & 0xfffff);
        x[i] = m.d1;
#else

        x[i]=(double)(Sobol_xn[i])/((double)(1<<MAXBIT));
#endif
        // increase n by 1
        Sobol_n[i]++;
    }


    return(0);
}

/***********************************************************
 Initialize the parameters of a stream of sobol' sequence
 streamid>=0 if successful, -1 if failed
***********************************************************/
int sobstream(
    int **Sobol_gen,
    int *Sobol_xn,
    int *Sobol_n,
    int *Sobol_nstream /*pass by reference */
)
{
    /*  extern int Sobol_gen[MAXBIT][DIMENSION];
      extern int Sobol_xn[DIMENSION];
      extern int Sobol_n[DIMENSION];
      extern int Sobol_nstream; */

    int j, k;
    int dirnum; /* direction number */
    int streamid; /* the streamid return to the user */

    if ((*Sobol_nstream)+1>DIMENSION) {
        fprintf(stderr,"Cannot obtain more sobol streams\n");
        return(-1);
    }

    /* set the number of quasirandom number streams */
    streamid=(*Sobol_nstream);
    (*Sobol_nstream)++;

    /* read in all of the initial value
       the initial value m1, m2, ..., md can be selected freely
       provided that each mi is odd and mi<2^i */
    for (j=0; j<MAXORDPOLY; j++)
        Sobol_gen[j][streamid]=Sobol_initdat[j][streamid];

    if (streamid==0) {
        /* dimension 0 is generated specially, 2^0 ~ 2^(MAXBIT-1) in each*/
        for (j=0;j<MAXBIT;j++)
            Sobol_gen[j][0]=1<<(MAXBIT-j-1);
    } else
        /* generate the other vc in other dimensions */
    {
        /* set up vj, where j from 0 to d-1, vj=2^j*m[i-j]
           the reason the set it as 2^MAX-j-1 is to set the most
           significant bit when translate it to decimal number */
        for (j=0; j<Sobol_deg[streamid]; j++)
            Sobol_gen[j][streamid]*=1<<(MAXBIT-j-1);

        /* set up vj for j from d to MAXBIT, using the direction
           number of the primitive polynomial */
        for (j=Sobol_deg[streamid];j<MAXBIT;j++) {
            /* Compute v[i-d]EOR[v[i-d]/2^d] */
            Sobol_gen[j][streamid]=Sobol_gen[j-Sobol_deg[streamid]][streamid];
            Sobol_gen[j][streamid]^=Sobol_gen[j][streamid]>>Sobol_deg[streamid];

            /* get the direction number */
            dirnum=Sobol_edi[streamid];

            /* compute the vi by vi=a[1]v[i-1]EORa[2]v[i-2]EOR...EORa[d-1]v[i-d+1] */
            for (k=Sobol_deg[streamid]-1; k>0; k--) {
                if ((dirnum&1)!=0)
                    Sobol_gen[j][streamid]^=Sobol_gen[j-k][streamid];
                dirnum>>=1;
            }
        }
    }

    /* in each dimension the initial value x[0] is 0, n is 0 at each
       dimension in the beginning*/
    Sobol_xn[streamid]=0;
    Sobol_n[streamid]=0;

    return(streamid);
}

/***********************************************************
  Everytime by calling sobolseq, a quasirandom number in the
  specified stream will be generated.
  the parameter streamid specifies the streamid in which
  the next quasirandom number will be generated
  A quasirandom number will return if successful
  -1.0 will return if failed
***********************************************************/
double sobolseq(int streamid,
                int **Sobol_gen,
                int *Sobol_xn,
                int *Sobol_n,
                int *Sobol_nstream /*pass by reference */
               )
{
    int j;
    int firstzero; /* the first zero bit in n */
    double x;

    if (streamid>=(*Sobol_nstream)) {
        fprintf(stderr,"Illegal streamid in function sobolseq(int streamid)\n");
        return(-1.0);
    }

    /* find the first zero bit in n stored in Sobol_n */
#if !defined(IRIX) && 0
    j = ~Sobol_n[streamid];
asm ("bsf %1, %%ecx\n\t movl %%ecx, %0\n" : "=r"(firstzero) : "r"(j) : "%ecx");
#else

    firstzero=0;
    for (j=Sobol_n[streamid];(j&1)!=0;j>>=1)
        firstzero++;
#endif

    /* x[n+1]=x[n]EORvc */
    Sobol_xn[streamid]^=Sobol_gen[firstzero][streamid];

    /* increase n by 1 */
    Sobol_n[streamid]++;

    /* divided by 2^MAXBIT */
    x=(double)(Sobol_xn[streamid])/((double)(1<<MAXBIT));
    return(x);
}

CORE_END_NAMESPACE
