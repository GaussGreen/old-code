/*--------------------------------------------------------------
        FILE: RandomGen.cxx
        PURPOSE: Random numbers generation routines
        AUTHOR: Dimitri Mayevski
        DATE: 25/02/2003
  --------------------------------------------------------------*/

#include "RandomGen.h"
#include "math.h"
#include "srt_h_all.h"

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0 / IM1)
#define IMM1 (IM1 - 1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NDIV (1 + IMM1 / SHUFFLE_NTAB)
#define RNMX_EPS 1.2e-7
#define RNMX (1.0 - RNMX_EPS)

static Err RandomGen_Uniform(SRandomGen *rg, double *res) {
  int j;
  long k;

  k = rg->idum / IQ1; // Compute idum=(IA1*idum) % IM1 without overflows
  rg->idum = IA1 * (rg->idum - k * IQ1) - k * IR1; // by Schrage's method
  if (rg->idum < 0)
    rg->idum += IM1;
  k = rg->idum2 / IQ2;
  rg->idum2 = IA2 * (rg->idum2 - k * IQ2) -
              k * IR2; // Compute idum2=(IA2*idum) % IM2 likewise
  if (rg->idum2 < 0)
    rg->idum2 += IM2;
  j = rg->iy / NDIV; // Will be in the range 0..SHUFFLE_NTAB-1
  rg->iy = rg->iv[j] -
           rg->idum2;   // Here idum is shuffled        , idum and idum2 are
  rg->iv[j] = rg->idum; // combined to generate output
  if (rg->iy < 1)
    rg->iy += IMM1;
  if ((*res = AM * rg->iy) > RNMX)
    *res = RNMX; // Because users don't expect endpoint values
  return NULL;
}

static Err RandomGen_Gauss(SRandomGen *rg, double *res) {
  double temp;
  RandomGen_Uniform(rg, &temp);
  *res = inv_cumnorm_fast(temp);
  return NULL;
}

Err RandomGen_Init(SRandomGen *rg, long seed) {
  int j;
  long k;
  rg->idum =
      (seed > 1 ? seed
                : (seed < -1 ? -seed : 1)); // Be sure to prevent idum = 0
  rg->idum2 = rg->idum;
  for (j = SHUFFLE_NTAB + 7; j >= 0;
       j--) { // Load the shuffle table (after 8 warm-ups)
    k = rg->idum / IQ1;
    rg->idum = IA1 * (rg->idum - k * IQ1) - k * IR1;
    if (rg->idum < 0)
      rg->idum += IM1;
    if (j < SHUFFLE_NTAB)
      rg->iv[j] = rg->idum;
  }
  rg->iy = rg->iv[0];

  rg->Uniform = RandomGen_Uniform;
  rg->Gauss = RandomGen_Gauss;
  rg->spec_desc = NULL;

  return NULL;
}

#define SWAP(a, b)                                                             \
  ltemp = (a);                                                                 \
  (a) = (b);                                                                   \
  (b) = ltemp

static Err FillIndices(SRandomGen *rg) {
  Err err = NULL;
  long *uni = NULL, *idx = NULL;
  int br;
  long i, j, k, ltemp;
  double a;
  SABSDesc *desc = (SABSDesc *)rg->spec_desc;

  uni = (long *)calloc(desc->TableSize, sizeof(long));
  if (!uni) {
    err = serror("Memory failure in FillIndices");
    goto FREE_RETURN;
  }

  if (!desc->bStrat || desc->iStep == 0) {
    for (i = 0; i < desc->TableSize; i++)
      uni[i] = i; // Init uni

    for (br = 0; br < desc->nBrownians; br++) {
      for (i = 0; i < desc->TableSize - 1; i++) // Permute randomly
      {
        rg->Uniform(rg, &a);
        k = (long)(a * (desc->TableSize - i)) + i;
        SWAP(uni[i], uni[k]);
      }
      memcpy(desc->Index[br], uni, desc->TableSize * sizeof(long));
      if (desc->bStrat)
        for (i = 0; i < desc->TableSize; i++)
          desc->Spots[br][i] = desc->NormalTable[uni[i]];
    }
  } else {
    idx = (long *)calloc(desc->TableSize, sizeof(long));
    if (!idx) {
      err = serror("Memory failure in FillIndices");
      goto FREE_RETURN;
    }

    for (br = 0; br < desc->nBrownians; br++) {
      err = indexx_dl(desc->Spots[br], idx, desc->TableSize);
      if (err)
        goto FREE_RETURN;

      // Initialize uni:
      for (j = 0; j < desc->NNode; j++)
        for (i = 0; i < (desc->NBuck >> 1); i++)
          uni[j * desc->NBuck + i] = i * desc->NNode + j;

      for (i = 0; i < (desc->NBuck >> 1); i++) {
        // Permute randomly within each of desc->NBuck/2 new buckets:
        for (j = 0; j < desc->NNode - 1; j++) {
          rg->Uniform(rg, &a);
          k = (long)(a * (desc->NNode - j)) + j;
          SWAP(uni[desc->NBuck * j + i], uni[desc->NBuck * k + i]);
        }
        // Apply antithetic rule:
        for (j = 0; j < desc->NNode; j++)
          uni[desc->NBuck * j + i + (desc->NBuck >> 1)] =
              desc->TableSize - 1 - uni[desc->NBuck * j + i];
      }
      // Permute buckets randomly within each of desc->NNode old nodes:
      for (j = 0; j < desc->NNode; j++)
        for (i = 0; i < desc->NBuck - 1; i++) {
          rg->Uniform(rg, &a);
          k = (long)(a * (desc->NBuck - i)) + i;
          SWAP(uni[desc->NBuck * j + i], uni[desc->NBuck * j + k]);
        }

      // Update indices and spots:
      for (i = 0; i < desc->TableSize; i++) {
        desc->Index[br][idx[i]] = uni[i];
        desc->Spots[br][idx[i]] += desc->NormalTable[uni[i]];
      }
    }
  }

FREE_RETURN:
  free(uni);
  free(idx);
  return err;
}

static Err ABS_Gauss(SRandomGen *rg, double *res) {
  Err err = NULL;
  SABSDesc *desc = (SABSDesc *)rg->spec_desc;

  if (desc->iPath == 0 && desc->iBrownian == 0) {
    err = FillIndices(rg);
    if (err)
      return err;
  }

  *res = desc->NormalTable[desc->Index[desc->iBrownian][desc->iPath]];

  if (++desc->iBrownian == desc->nBrownians) {
    desc->iBrownian = 0;
    if (++desc->iPath == desc->nPaths) {
      desc->iPath = 0;
      ++desc->iStep;
    }
  }
  return NULL;
}

Err ABS_Init(SRandomGen *rg, long seed, long nPaths, int nBrownians,
             int bStrat) {
  SABSDesc *desc = NULL;
  long i;

  RandomGen_Init(rg, seed);

  desc = (SABSDesc *)malloc(sizeof(SABSDesc));
  if (!desc)
    return serror("Memory failure in ABS_Init");
  rg->spec_desc = desc;
  memset(desc, 0, sizeof(SABSDesc));

  desc->nPaths = nPaths;
  desc->nBrownians = nBrownians;
  desc->bStrat = bStrat;

  desc->TableSize = nPaths;

  // Calculate NNode and NBuck (NBuck must be even)
  if (bStrat) {
    if (nPaths < 10)
      return serror("Too few paths for ABS");

    if (desc->TableSize % 2)
      desc->TableSize++;

    desc->NBuck = 2 * (long)(sqrt(desc->TableSize) / 2 + 1e-8);
    while (desc->TableSize % desc->NBuck)
      desc->NBuck += 2;
    desc->NNode = desc->TableSize / desc->NBuck;

    desc->Spots = dmatrix(0, nBrownians - 1, 0, desc->TableSize - 1);
    if (!desc->Spots)
      return serror("Memory failure in ABS_Init");
  }

  desc->NormalTable = (double *)calloc(desc->TableSize, sizeof(double));
  if (!desc->NormalTable)
    return serror("Memory failure in ABS_Init");

  // Fill the table of normal deviates:
  for (i = 0; i < desc->TableSize; i++)
    desc->NormalTable[i] = inv_cumnorm_fast((i + 0.5) / desc->TableSize);

  desc->Index = (long **)calloc(nBrownians, sizeof(long *));
  if (!desc->Index)
    return serror("Memory failure in ABS_Init");
  memset(desc->Index, 0, nBrownians * sizeof(long *));

  for (i = 0; i < nBrownians; i++) {
    desc->Index[i] = (long *)calloc(desc->TableSize, sizeof(long));
    if (!desc->Index[i])
      return serror("Memory failure in ABS_Init");
  }

  rg->Gauss = ABS_Gauss;
  return NULL;
}

Err ABS_Free(SRandomGen *rg) {
  long i;
  SABSDesc *desc = (SABSDesc *)rg->spec_desc;

  if (desc) {
    free(desc->NormalTable);
    if (desc->Spots)
      free_dmatrix(desc->Spots, 0, desc->nBrownians - 1, 0,
                   desc->TableSize - 1);
    if (desc->Index)
      for (i = 0; i < desc->nBrownians; i++)
        free(desc->Index[i]);
    free(desc->Index);
  }
  free(desc);
  memset(rg, 0, sizeof(SRandomGen));

  return NULL;
}
