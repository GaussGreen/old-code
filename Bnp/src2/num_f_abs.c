/* -------------------------------------------------------------------------
   FILENAME		: num_f_utabs.c
   AUTHOR		: A. Savine
   PURPOSE		: generate brownian paths using abs technique
   ------------------------------------------------------------------------- */
#include "math.h"
#include "num_h_allhdr.h"

static void RandomPermutation(double *RandomArray, long Length, long *seed) {
  long Rand, Counter;
  double TempVar, tmp;

  for (Counter = 0; Counter < Length; Counter++) {
    tmp = uniform_fast(seed);
    Rand = (long)(tmp * (Length - Counter - 1)) + Counter;

    /* Swap variables */
    TempVar = RandomArray[Counter];
    RandomArray[Counter] = RandomArray[Rand];
    RandomArray[Rand] = TempVar;
  }
}

#define SWAP(a, b)                                                             \
  itemp = (a);                                                                 \
  (a) = (b);                                                                   \
  (b) = itemp
#define M 7
#define NSTACK 50

static long QuickSort(double *arr, long *indx, long end, long start) {
  long i, indxt, ir, itemp, j, k, l, jstack = 0;
  unsigned long *istack = NULL;
  double a;

  l = start;
  ir = end;

  istack = lvector(1, NSTACK);

  if (istack == NULL)
    return 1;

  for (;;) {
    if (ir - l < M) {
      for (j = l + 1; j <= ir; j++) {
        indxt = indx[j];
        a = arr[indxt];

        for (i = j - 1; i >= l; i--) {
          if (arr[indx[i]] <= a)
            break;

          indx[i + 1] = indx[i];
        }

        indx[i + 1] = indxt;
      }

      if (jstack == 0)
        break;

      ir = istack[jstack--];
      l = istack[jstack--];
    } else {
      k = (l + ir) >> 1;
      SWAP(indx[k], indx[l + 1]);

      if (arr[indx[l]] > arr[indx[ir]]) {
        SWAP(indx[l], indx[ir]);
      }

      if (arr[indx[l + 1]] > arr[indx[ir]]) {
        SWAP(indx[l + 1], indx[ir]);
      }

      if (arr[indx[l]] > arr[indx[l + 1]]) {
        SWAP(indx[l], indx[l + 1]);
      }

      i = l + 1;
      j = ir;
      indxt = indx[l + 1];
      a = arr[indxt];

      for (;;) {
        do
          i++;
        while (arr[indx[i]] < a);

        do
          j--;
        while (arr[indx[j]] > a);

        if (j < i)
          break;

        SWAP(indx[i], indx[j]);
      }

      indx[l + 1] = indx[j];
      indx[j] = indxt;

      jstack += 2;

      if (jstack > NSTACK)
        return 1;

      if (ir - i + 1 >= j - l) {
        istack[jstack] = ir;
        istack[jstack - 1] = i;
        ir = j - 1;
      } else {
        istack[jstack] = j - 1;
        istack[jstack - 1] = l;
        l = i;
      }
    }
  }

  free_lvector(istack, 1, NSTACK);

  return 0;
}

#undef SWAP
#undef M
#undef NSTACK

static unsigned long IsDivisible(long x, long y) {

  double xx, yy;

  xx = (double)x;
  yy = (double)y;

  return (xx / yy == x / y);
}
/* -------------------------------------------------------------------------- */

/* --------------------------------------------------------------------------
   The generation of a Cube of Anthithetic Balanced Stratified set of
   Gaussian random numbers
   The cube is described as follows:
                          rand[path][brownian][step]
   where:
                - first_path     <= path     <= last_path
                - first_step     <= step     <= last_step
                - first_brownian <= brownian <= last_brownian
   Memory allocation for the cube is done OUTSIDE
  -------------------------------------------------------------------------- */

Err ABSCube(double ***rand, long first_path, long last_path,
            long first_brownian, long last_brownian, long first_step,
            long last_step, long *seed) {
  double *uni1, *uni2, *vect, *spots;
  long i, j, k, ibr, b, NBuck, NPB, temp, *index;
  long num_path = last_path - first_path + 1;
  long num_step = last_step - first_step + 1;
  long gen_path;
  SRT_Boolean odd = SRT_NO;

  /* Number of generated paths must be even (antithetic principle) */
  gen_path = num_path;
  if ((gen_path % 2) != 0) {
    odd = SRT_YES;
    gen_path += 1;
  }

  /* Number of paths per bucket  , must be even */
  if (gen_path == 2) {
    NPB = 2;
  } else {
    NPB = 2 * (long)(sqrt(gen_path) / 2);
    /*NPB = ( mod(gen_path  , NPB) );*/
    while (!IsDivisible(gen_path, NPB))
      NPB += 2;
  }

  /* Number of buckets */
  NBuck = gen_path / NPB;

  /* Memory allocation */
  uni1 = dvector(0, gen_path - 1);
  uni2 = dvector(0, gen_path - 1);
  vect = dvector(0, NPB);
  spots = dvector(0, gen_path - 1);
  index = lngvector(0, gen_path - 1);
  if ((uni1 == NULL) || (uni2 == NULL) || (vect == NULL) || (spots == NULL) ||
      (index == NULL))
    return serror("memory error in ABSCube");

  memset(&spots[0], 0, (long)(gen_path * sizeof(double)));

  /* Creates a uniform mesh (for the number of path) */
  for (j = 0; j < gen_path; j++) {
    index[j] = j;
    uni1[j] = (0.50 + j) / gen_path;
  }

  /* Loop on the number of Brownians */
  for (ibr = first_brownian; ibr <= last_brownian; ibr++) {

    /* Initialises uni2 as uni1 (uniform mesh) */
    memcpy(&uni2[0], &uni1[0], gen_path * sizeof(double));

    /* Permutates randomly the uniform mesh (for uncorrelated Brownians on first
     * step) */
    RandomPermutation(uni2, gen_path, seed);

    /* Initialises first step as regular grid (with random ordering...) */
    for (j = 0; j < gen_path - 1; j++) {
      spots[j] = inv_cumnorm_fast(uni2[j]);
      rand[first_path + j][ibr][first_step] = inv_cumnorm_fast(uni2[j]);
    }

    /* Be careful if it is an odd number of path */
    if (odd == SRT_NO) {
      rand[last_path][ibr][first_step] = inv_cumnorm_fast(uni2[gen_path - 1]);
    }

    /* Loops on all the steps (after the first one) */
    for (k = 1; k < num_step; k++) {
      /* Sorts relative to spot value so far */
      if (QuickSort(spots, index, gen_path - 1, 0))
        return serror("Sorting error in ABSCube");

      /* Initialises uni2 as uni1 : uni2 is a set of NPB blocks of NBuck points
       * (dual to paths) */
      memcpy(uni2, uni1, gen_path * sizeof(double));

      /* Bucket number */
      b = 0;

      /* Treatment of all the increments  , one by one */
      for (j = 0; j < gen_path;) {

        /* Treats all the increments v[0..NPB-1] in the curent bucket # b by
         * INCREASING order */
        for (i = 0; i < NPB / 2; i++) {
          /* Chooses randomly a REMAINING increment to associate to the current
           * path (#i) in current bucket (#b) */
          temp = (long)(uniform_fast(seed) * (NBuck - 1 - b) + i * NBuck);
          vect[i] = uni2[temp];

          /* Antithetic principle on the bucket */
          vect[i + NPB / 2] = 1 - uni2[temp];

          /* Copies the last value of the uniform bucket in place of the chosen
           * value to allow to choose it next time */
          uni2[temp] = uni2[NBuck - b - 1 + i * NBuck];
        }

        /* Permutates randomly these increments in the bucket (so far  , first
         * increment is lower  ,...) */
        RandomPermutation(vect, NPB, seed);

        /* Associates these uniform numbers to their randomly selected path in
         * the bucket */
        for (i = 0; i < NPB; i++) {
          spots[index[j + i]] += inv_cumnorm_fast(vect[i]);
          if (first_path + index[j + i] <= last_path) {
            rand[first_path + index[j + i]][ibr][first_step + k] =
                inv_cumnorm_fast(vect[i]);
          } else {
            if (odd == SRT_NO) {
              rand[last_path][ibr][first_step + k] = inv_cumnorm_fast(vect[i]);
            }
          }
        }

        /* Move on to the next bucket (treated NPB paths...) */
        j += NPB;
        b++;

      } /* end of loop on path */

    } /* end of loop on step */

  } /* end of loop on brownian */

  /* Free memory */
  free_dvector(uni1, 0, gen_path - 1);
  free_dvector(uni2, 0, gen_path - 1);
  free_dvector(vect, 0, NPB);
  free_dvector(spots, 0, gen_path - 1);
  free_lngvector(index, 0, gen_path - 1);

  return NULL;
}
