/*--------------------------------------------------------------
        FILE: RandomGen.h
        PURPOSE: Random numbers generation routines
        AUTHOR: Dimitri Mayevski
        DATE: 25/02/2003
  --------------------------------------------------------------*/

#ifndef __RANDOMGEN_H__
#define __RANDOMGEN_H__

#define SHUFFLE_NTAB 32

typedef struct _SRandomGen SRandomGen;

struct _SRandomGen {
  long idum, idum2, iv[SHUFFLE_NTAB], iy;

  Err (*Uniform)(SRandomGen *rg, double *res);
  Err (*Gauss)(SRandomGen *rg, double *res);

  void *spec_desc;
};

typedef struct _SABSDesc {
  double *NormalTable, **Spots;
  long NBuck, NNode, TableSize;
  long **Index;
  int nBrownians, iBrownian, bStrat;
  long nPaths, iPath, iStep;
} SABSDesc;

Err RandomGen_Init(SRandomGen *rg, long seed);
Err ABS_Init(SRandomGen *rg, long seed, long nPaths, int nBrownians,
             int bStrat);
Err ABS_Free(SRandomGen *rg);

#endif // #ifndef __RANDOMGEN_H__