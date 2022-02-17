////   Various routines of sorting arrays

#include "num_h_allhdr.h"

// S. Galluccio: 21 November 2000
// Created Source File

void sort_simple(unsigned long n, double *arr, double *brr) {
  unsigned long i, ir = n, j, k, l = 1;
  int *istack, jstack = 0;
  double a, b, temp;

  istack = ivector(1, NSTACK);
  for (;;) {
    if (ir - l < MPAR) {
      for (j = l + 1; j <= ir; j++) {
        a = arr[j];
        b = brr[j];
        for (i = j - 1; i >= 1; i--) {
          if (arr[i] <= a)
            break;
          arr[i + 1] = arr[i];
          brr[i + 1] = brr[i];
        }
        arr[i + 1] = a;
        brr[i + 1] = b;
      }
      if (!jstack) {
        free_ivector(istack, 1, NSTACK);
        return;
      }
      ir = istack[jstack];
      l = istack[jstack - 1];
      jstack -= 2;
    } else {
      k = (l + ir) >> 1;
      DSWAP(arr[k], arr[l + 1]);
      DSWAP(brr[k], brr[l + 1]);
      if (arr[l + 1] > arr[ir]) {
        DSWAP(arr[l + 1], arr[ir]);
        DSWAP(brr[l + 1], brr[ir]);
      }
      if (arr[l] > arr[ir]) {
        DSWAP(arr[l], arr[ir]);
        DSWAP(brr[l], brr[ir]);
      }
      if (arr[l + 1] > arr[l]) {
        DSWAP(arr[l + 1], arr[l]);
        DSWAP(brr[l + 1], brr[l]);
      }
      i = l + 1;
      j = ir;
      a = arr[l];
      b = brr[l];
      for (;;) {
        do
          i++;
        while (arr[i] < a);
        do
          j--;
        while (arr[j] > a);
        if (j < i)
          break;
        DSWAP(arr[i], arr[j]);
        DSWAP(brr[i], brr[j]);
      }
      arr[l] = arr[j];
      arr[j] = a;
      brr[l] = brr[j];
      brr[j] = b;
      jstack += 2;
      if (jstack > NSTACK)
        goto end; // nrerror("NSTACK too small in sort2.");
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
end:;
}

/////////////////////////////// indexing and ranking of a given array
///////////////////////////

void index_data(unsigned long n, double *arr, unsigned long *indx) {
  unsigned long i, indxt, ir = n, itemp, j, k, l = 1;
  int jstack = 0, *istack;
  double a;

  istack = ivector(1, NSTACK);
  for (j = 1; j <= n; j++)
    indx[j] = j;
  for (;;) {
    if (ir - l < MPAR) {
      for (j = l + 1; j <= ir; j++) {
        indxt = indx[j];
        a = arr[indxt];
        for (i = j - 1; i >= 1; i--) {
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
      LSWAP(indx[k], indx[l + 1]);
      if (arr[indx[l + 1]] > arr[indx[ir]]) {
        LSWAP(indx[l + 1], indx[ir]);
      }
      if (arr[indx[l]] > arr[indx[ir]]) {
        LSWAP(indx[l], indx[ir]);
      }
      if (arr[indx[l + 1]] > arr[indx[l]]) {
        LSWAP(indx[l + 1], indx[l]);
      }
      i = l + 1;
      j = ir;
      indxt = indx[l];
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
        LSWAP(indx[i], indx[j]);
      }
      indx[l] = indx[j];
      indx[j] = indxt;
      jstack += 2;
      if (jstack > NSTACK)
        goto end; // nrerror("NSTACK too small in indexx.");
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
end:

  free_ivector(istack, 1, NSTACK);
}

////////////////////////////// Ranks an indexed table
///////////////////////////////////

void rank_data(unsigned long n,
               unsigned long *indx,  // input: the indexed table
               unsigned long *irank) // output: the rank table
{
  unsigned long j;

  for (j = 1; j <= n; j++)
    irank[indx[j]] = j;
}

///////////////////////////////////////////////////////////////////////////////////////

void SortSeries(long n, double *ra) {
  long i, ir, j, l;
  double rra;

  if (n < 2)
    return;
  l = (n >> 1) + 1;
  ir = n;
  for (;;) {
    if (l > 1) {
      rra = ra[--l];
    } else {
      rra = ra[ir];
      ra[ir] = ra[1];
      if (--ir == 1) {
        ra[1] = rra;
        break;
      }
    }
    i = l;
    j = l + l;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j + 1])
        j++;
      if (rra < ra[j]) {
        ra[i] = ra[j];
        i = j;
        j <<= 1;
      } else
        j = ir + 1;
    }
    ra[i] = rra;
  }
}

///////////////////////////////////////////////////////////////////////////////////////
// Sorts a series in a descending  , rather than ascending order
///////////////////////////////////////////////////////////////////////////////////////
void SortSeriesDescending(long n, double *ra) {
  long i, ir, j, l;
  double rra;

  if (n < 2)
    return;
  l = (n >> 1) + 1;
  ir = n;
  for (;;) {
    if (l > 1) {
      rra = ra[--l];
    } else {
      rra = ra[ir];
      ra[ir] = ra[1];
      if (--ir == 1) {
        ra[1] = rra;
        break;
      }
    }
    i = l;
    j = l + l;
    while (j <= ir) {
      if (j < ir && ra[j] > ra[j + 1])
        j++;
      if (rra > ra[j]) {
        ra[i] = ra[j];
        i = j;
        j <<= 1;
      } else
        j = ir + 1;
    }
    ra[i] = rra;
  }
}

///////////////////////////////////////////////////////////////////////////////////////

void sort_simple2(unsigned long n, double **arr) {
  unsigned long i, ir = n, j, k, l = 1;
  int *istack, jstack = 0;
  double a, b, temp;

  istack = ivector(1, NSTACK);
  for (;;) {
    if (ir - l < MNUM) {
      for (j = l + 1; j <= ir; j++) {
        a = arr[j][0];
        b = arr[j][1];
        for (i = j - 1; i >= 1; i--) {
          if (arr[i][0] <= a)
            break;
          arr[i + 1][0] = arr[i][0];
          arr[i + 1][1] = arr[i][1];
        }
        arr[i + 1][0] = a;
        arr[i + 1][1] = b;
      }
      if (!jstack) {
        free_ivector(istack, 1, NSTACK);
        return;
      }
      ir = istack[jstack];
      l = istack[jstack - 1];
      jstack -= 2;
    } else {
      k = (l + ir) >> 1;
      DSWAP(arr[k][0], arr[l + 1][0]);
      DSWAP(arr[k][1], arr[l + 1][1]);

      if (arr[l + 1][0] > arr[ir][0]) {
        DSWAP(arr[l + 1][0], arr[ir][0]);
        DSWAP(arr[l + 1][1], arr[ir][1]);
      }
      if (arr[l][0] > arr[ir][0]) {
        DSWAP(arr[l][0], arr[ir][0]);
        DSWAP(arr[l][1], arr[ir][1]);
      }
      if (arr[l + 1][0] > arr[l][0]) {
        DSWAP(arr[l + 1][0], arr[l][0]);
        DSWAP(arr[l + 1][1], arr[l][1]);
      }
      i = l + 1;
      j = ir;
      a = arr[l][0];
      b = arr[l][1];
      for (;;) {
        do
          i++;
        while (arr[i][0] < a);
        do
          j--;
        while (arr[j][0] > a);
        if (j < i)
          break;
        DSWAP(arr[i][0], arr[j][0]);
        DSWAP(arr[i][1], arr[j][1]);
      }
      arr[l][0] = arr[j][0];
      arr[j][0] = a;

      arr[l][1] = arr[j][1];
      arr[j][1] = b;

      jstack += 2;
      if (jstack > NSTACK)
        goto end;
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

end:

  free_ivector(istack, 1, NSTACK);
}

///////////////////////////////////////////////////////////////////////////////////////

void shell(unsigned long n, double **a) {
  unsigned long i, j, inc;
  double v, w;

  inc = 1;
  do {
    inc *= 3;
    inc++;
  } while (inc <= n);
  do {
    inc /= 3;
    for (i = inc + 1; i <= n; i++) {
      v = a[i][0];
      w = a[i][1];
      j = i;
      while (a[j - inc][0] > v) {
        a[j][0] = a[j - inc][0];
        a[j][1] = a[j - inc][1];

        j -= inc;
        if (j <= inc)
          break;
      }
      a[j][0] = v;
      a[j][1] = w;
    }
  } while (inc > 1);
}

Err indexx_ll(long *f, long *pI, long n) {
  Err err = NULL;
  long i, indxt, itemp, j, k, l = 0, ir = n - 1;
  unsigned long *istack = NULL;
  int jstack = 0;
  long a;

  istack = lvector(1, NSTACK);
  if (!istack) {
    err = serror("Memory failure in indexx_ll");
    goto FREE_RETURN;
  }

  for (j = 0; j < n; j++)
    pI[j] = j;

  for (;;) {
    if (ir - l < MPAR) { /* Insertion sort when subarray small enough */
      for (j = l + 1; j <= ir; j++) {
        indxt = pI[j];
        a = f[indxt];
        for (i = j - 1; i >= l; i--) {
          if (f[pI[i]] <= a)
            break;
          pI[i + 1] = pI[i];
        }
        pI[i + 1] = indxt;
      }
      if (jstack == 0)
        break;
      ir = istack[jstack--]; /* Pop stack and begin */
      l = istack[jstack--];  /* a new round of partitioning */
    } else {
      k = (unsigned long)(l + ir) >>
          1; /* Choose median of left  , center  , and right */
      LSWAP(pI[k], pI[l + 1]); /* elements as partitioning element a. */
                               /* Rearrange so that a[l] <= a[l+1] <= a[ir] */
      if (f[pI[l]] > f[pI[ir]]) {
        LSWAP(pI[l], pI[ir]);
      }
      if (f[pI[l + 1]] > f[pI[ir]]) {
        LSWAP(pI[l + 1], pI[ir]);
      }
      if (f[pI[l]] > f[pI[l + 1]]) {
        LSWAP(pI[l], pI[l + 1]);
      }
      i = l + 1; /* Initialize pointers for partitioning */
      j = ir;
      indxt = pI[l + 1];
      a = f[indxt]; /* Partitioning element */
      for (;;)      /* Beginning the innermost loop */
      {
        do
          i++;
        while (f[pI[i]] < a); /* Scan up to find element > a */
        do
          j--;
        while (f[pI[j]] > a); /* Scan down to find element < a */
        if (j < i)
          break;             /* Pointers crossed. Partitioning complete. */
        LSWAP(pI[i], pI[j]); /* Exchange elements */
      }                      /* End of innermost loop */
      pI[l + 1] = pI[j];     /* Insert partitioning element */
      pI[j] = indxt;
      jstack += 2;
      /* Push pointers to larger subarray on stack  , process smaller subarray
       * immediately */
      if (jstack > NSTACK) {
        err = serror("Stack overflow in indexx_ll");
        goto FREE_RETURN;
      }
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

FREE_RETURN:
  if (istack)
    free_lvector(istack, 1, NSTACK);
  return err;
}

Err qsort_d(double *arr, long n) {
  Err err = NULL;
  long i, j, k, l = 0, ir = n - 1;
  unsigned long *istack = NULL;
  int jstack = 0;
  double a, temp;

  istack = lvector(1, NSTACK);
  if (!istack) {
    err = serror("Memory failure in qsort_d");
    goto FREE_RETURN;
  }

  for (;;) {
    if (ir - l < MPAR) { /* Insertion sort when subarray small enough */
      for (j = l + 1; j <= ir; j++) {
        a = arr[j];
        for (i = j - 1; i >= l; i--) {
          if (arr[i] <= a)
            break;
          arr[i + 1] = arr[i];
        }
        arr[i + 1] = a;
      }
      if (jstack == 0)
        break;
      ir = istack[jstack--]; /* Pop stack and begin */
      l = istack[jstack--];  /* a new round of partitioning */
    } else {
      k = (unsigned long)(l + ir) >>
          1; /* Choose median of left  , center  , and right */
      DSWAP(arr[k], arr[l + 1]); /* elements as partitioning element a. */
                                 /* Rearrange so that a[l] <= a[l+1] <= a[ir] */
      if (arr[l] > arr[ir]) {
        DSWAP(arr[l], arr[ir]);
      }
      if (arr[l + 1] > arr[ir]) {
        DSWAP(arr[l + 1], arr[ir]);
      }
      if (arr[l] > arr[l + 1]) {
        DSWAP(arr[l], arr[l + 1]);
      }
      i = l + 1; /* Initialize pointers for partitioning */
      j = ir;
      a = arr[l + 1]; /* Partitioning element */
      for (;;)        /* Beginning the innermost loop */
      {
        do
          i++;
        while (arr[i] < a); /* Scan up to find element > a */
        do
          j--;
        while (arr[j] > a); /* Scan down to find element < a */
        if (j < i)
          break;               /* Pointers crossed. Partitioning complete. */
        DSWAP(arr[i], arr[j]); /* Exchange elements */
      }                        /* End of innermost loop */
      arr[l + 1] = arr[j];     /* Insert partitioning element */
      arr[j] = a;
      jstack += 2;
      /* Push pointers to larger subarray on stack  , process smaller subarray
       * immediately */
      if (jstack > NSTACK) {
        err = serror("Stack overflow in qsort_d");
        goto FREE_RETURN;
      }
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

FREE_RETURN:
  if (istack)
    free_lvector(istack, 1, NSTACK);
  return err;
}

Err indexx_dl(double *f, long *pI, long n) {
  Err err = NULL;
  long i, indxt, itemp, j, k, l = 0, ir = n - 1;
  unsigned long *istack = NULL;
  int jstack = 0;
  double a;

  istack = lvector(1, NSTACK);
  if (!istack) {
    err = serror("Memory failure in indexx_dl");
    goto FREE_RETURN;
  }

  for (j = 0; j < n; j++)
    pI[j] = j;

  for (;;) {
    if (ir - l < MPAR) { /* Insertion sort when subarray small enough */
      for (j = l + 1; j <= ir; j++) {
        indxt = pI[j];
        a = f[indxt];
        for (i = j - 1; i >= l; i--) {
          if (f[pI[i]] <= a)
            break;
          pI[i + 1] = pI[i];
        }
        pI[i + 1] = indxt;
      }
      if (jstack == 0)
        break;
      ir = istack[jstack--]; /* Pop stack and begin */
      l = istack[jstack--];  /* a new round of partitioning */
    } else {
      k = (unsigned long)(l + ir) >>
          1; /* Choose median of left  , center  , and right */
      LSWAP(pI[k], pI[l + 1]); /* elements as partitioning element a. */
                               /* Rearrange so that a[l] <= a[l+1] <= a[ir] */
      if (f[pI[l]] > f[pI[ir]]) {
        LSWAP(pI[l], pI[ir]);
      }
      if (f[pI[l + 1]] > f[pI[ir]]) {
        LSWAP(pI[l + 1], pI[ir]);
      }
      if (f[pI[l]] > f[pI[l + 1]]) {
        LSWAP(pI[l], pI[l + 1]);
      }
      i = l + 1; /* Initialize pointers for partitioning */
      j = ir;
      indxt = pI[l + 1];
      a = f[indxt]; /* Partitioning element */
      for (;;)      /* Beginning the innermost loop */
      {
        do
          i++;
        while (f[pI[i]] < a); /* Scan up to find element > a */
        do
          j--;
        while (f[pI[j]] > a); /* Scan down to find element < a */
        if (j < i)
          break;             /* Pointers crossed. Partitioning complete. */
        LSWAP(pI[i], pI[j]); /* Exchange elements */
      }                      /* End of innermost loop */
      pI[l + 1] = pI[j];     /* Insert partitioning element */
      pI[j] = indxt;
      jstack += 2;
      /* Push pointers to larger subarray on stack  , process smaller subarray
       * immediately */
      if (jstack > NSTACK) {
        err = serror("Stack overflow in indexx_dl");
        goto FREE_RETURN;
      }
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

FREE_RETURN:
  if (istack)
    free_lvector(istack, 1, NSTACK);
  return err;
}
