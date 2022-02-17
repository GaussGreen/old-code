/* -------------------------------------------------------------------------
   FILENAME		: num_f_vector_tools.c
   AUTHOR		: Toto 16Dec1999

   ------------------------------------------------------------------------- */
#include "math.h"
#include "num_h_allhdr.h"

#define MY_DBL_EPSILON                                                         \
  (10.0 * DBL_EPSILON) /* a little above the smallest such                     \
                          that 1.0+DBL_EPSILON != 1.0 */

#define LOCAL_DTOL(X) ((long)(X + 1.0e-08))
/*	Fill evenly a vector	*/
void num_f_fill_vector(
    /*	Size of the input vector	*/
    int p,
    /*	Vector byref	*/
    double **v,
    /*	Required size of the output vector	*/
    int n,
    /*	Required size of the first section	*/
    int n0) {
  double *v2 = (double *)calloc(n, sizeof(double));
  double total_lenght = (*v)[p - 1] - (*v)[0];
  int total_number_of_extra_points = n - p;
  int i, free_idx = 0, k, number_of_extra_points;
  double proportional_lenght;

  if (total_number_of_extra_points < 1) {
    free(v2);
    return;
  }

  for (i = 0; i < p - 1; i++)

  {
    v2[free_idx] = (*v)[i];
    free_idx++;

    proportional_lenght = ((*v)[i + 1] - (*v)[i]) / total_lenght;

    number_of_extra_points =
        LOCAL_DTOL(proportional_lenght * total_number_of_extra_points);
    if (i == 0 && number_of_extra_points < n0) {
      number_of_extra_points = n0;
    }

    total_lenght -= (*v)[i + 1] - (*v)[i];

    if (number_of_extra_points > 0) {
      total_number_of_extra_points -= number_of_extra_points;
      for (k = 1; k <= number_of_extra_points; k++) {
        v2[free_idx] = (*v)[i] + k * ((*v)[i + 1] - (*v)[i]) /
                                     (1 + number_of_extra_points);
        free_idx++;
      }
    }
  }

  v2[n - 1] = (*v)[p - 1];
  free(*v);
  *v = v2;
}

/*	Another algorithm... */
void num_f_fill_vector_newalgo(
    /*	Initial size of vector	*/
    int *p,
    /*	Vector byref	*/
    double **v,
    /*	Required size of the output vector	*/
    int n) {
  int i, j, k;
  int n1, n2, n3, nn;
  double *v1, *v2;
  double t;
  double idiff, diff;
  int i0, j0;

#ifdef _DEGUG
  FILE *f = fopen("c:\\toto.txt", "w");

  fprintf(f, "Starting\n");
  fprintf(f, "Existing points: %d\n", *p);
  for (i = 0; i < *p - 1; i++) {
    fprintf(f, "%.2f  , ", (*v)[i]);
  }
  fprintf(f, "%.2f\n", (*v)[*p - 1]);
#endif

  v1 = *v;
  n1 = *p;

  if (n < 2 * n1 - 1) {
    n = 2 * n1 - 1;
  }

#ifdef _DEGUG
  fprintf(f, "Requested points: %d\n", n);
  fprintf(f, "Extra points: %d\n", n - n1);
#endif

  idiff = (v1[n1 - 1] - v1[0]) / (n - 1);
  n -= n1;

  diff = idiff / 1.01;
  do {
    diff *= 1.01;
    n2 = 0;
    for (i = 0; i < n1 - 1; i++) {
      j = (int)((v1[i + 1] - v1[i] - 1.0e-08) / diff);
      if (j < 1) {
        j = 1;
      }
      n2 += j;
    }
  } while (n2 > n);

  nn = n - n2;

  n3 = 0;
  i = 1;
  do {
    j = (int)((v1[i] - v1[i - 1] - 1.0e-08) / idiff);
    if (j < 1) {
      j = 1;
    }

    k = (int)((v1[i] - v1[i - 1] - 1.0e-08) / diff);
    if (k < 1) {
      k = 1;
    }

    n3 += j - k;
    i++;
  } while (n3 < nn);

  i0 = i - 2;

#ifdef _DEGUG
  fprintf(f, "Part 1 points end at pillar: %d\n", (i0 > 0 ? i0 : 0));
#endif

  v2 = (double *)calloc(n, sizeof(double));
  n2 = 0;

  for (i = 0; i < i0; i++) {
    j = (int)((v1[i + 1] - v1[i] - 1.0e-08) / idiff);
    if (j < 1) {
      j = 1;
    }

    t = (v1[i + 1] - v1[i]) / (j + 1);

    for (k = 0; k < j; k++) {
      v2[n2 + k] = v1[i] + (k + 1) * t;
    }

    n2 += j;
  }

#ifdef _DEGUG
  fprintf(f, "Part 1 points: %d\n", n2);
  for (i = 0; i < n2 - 1; i++) {
    fprintf(f, "%.2f  , ", v2[i]);
  }
  if (n2 - 1 >= 0)
    fprintf(f, "%.2f\n", v2[n2 - 1]);
  nn = n2;
#endif

  for (i = i0; i < n1 - 2; i++) {
    j = (int)((v1[i + 1] - v1[i] - 1.0e-08) / diff);
    if (j < 1) {
      j = 1;
    }

    t = (v1[i + 1] - v1[i]) / (j + 1);

    for (k = 0; k < j; k++) {
      v2[n2 + k] = v1[i] + (k + 1) * t;
    }

    n2 += j;
  }

#ifdef _DEGUG
  fprintf(f, "Part 2 points: %d\n", n2 - nn);
  for (i = nn; i < n2 - 1; i++) {
    fprintf(f, "%.2f  , ", v2[i]);
  }
  if (n2 - 1 >= 0)
    fprintf(f, "%.2f\n", v2[n2 - 1]);
  nn = n2;
#endif

  j0 = n - n2;
  t = (v1[n1 - 1] - v1[n1 - 2]) / (j0 + 1);

  for (k = 0; k < j0; k++) {
    v2[n2 + k] = v1[n1 - 2] + (k + 1) * t;
  }

  n2 += j0;

#ifdef _DEGUG
  fprintf(f, "Part 3 points: %d\n", n2 - nn);
  for (i = nn; i < n2 - 1; i++) {
    fprintf(f, "%.2f  , ", v2[i]);
  }
  fprintf(f, "%.2f\n", v2[n2 - 1]);

  fprintf(f, "Total points %d - requested %d", n2, n);
#endif

  if (n2) {
    num_f_concat_vector(p, v, n2, v2);
    num_f_sort_vector(*p, *v);
    num_f_unique_vector(p, *v);

    free(v2);
  }

#ifdef _DEGUG
  fclose(f);
#endif
}

/*	Fill evenly a vector so that the space between two consecutive points is
   less than a specified value */
void num_f_fill_vector_maxtime(
    /*	Initial size of vector	*/
    int *p,
    /*	Vector byref	*/
    double **v,
    /*	Required maxtime	*/
    double maxtime) {
  int i, j, k;
  int n1, n2;
  double *v1, *v2;
  double t;

  v1 = *v;
  n1 = *p;

  v2 = NULL;
  n2 = 0;

  for (i = 0; i < n1 - 1; i++) {
    if (v1[i + 1] - v1[i] > maxtime) {
      j = (int)((v1[i + 1] - v1[i] - 1.0e-08) / maxtime);
      t = (v1[i + 1] - v1[i]) / (j + 1);

      if (n2) {
        v2 = (double *)realloc(v2, (n2 + j) * sizeof(double));
      } else {
        v2 = (double *)calloc(j, sizeof(double));
      }

      for (k = 0; k < j; k++) {
        v2[n2 + k] = v1[i] + (k + 1) * t;
      }

      n2 += j;
    }
  }

  if (n2) {
    num_f_concat_vector(p, v, n2, v2);
    num_f_sort_vector(*p, *v);
    num_f_unique_vector(p, *v);

    free(v2);
  }
}

/*	Add a double	*/
void num_f_add_number(
    /*	Initial size of vector	*/
    int *p,
    /*	Initial vector byref	*/
    double **v1,
    /*	Number to be addded	*/
    double d) {
  *v1 = (double *)realloc(*v1, (*p + 1) * sizeof(double));

  (*v1)[*p] = d;
  (*p) += 1;
}

/*	Concatenate two vectors	*/
void num_f_concat_vector(
    /*	Initial size of vector	*/
    int *p,
    /*	Initial vector byref	*/
    double **v1,
    /*	Size of vector to be concatenated	*/
    int n,
    /*	Vector to be concatenated	*/
    double *v2) {
  int i;
  *v1 = (double *)realloc(*v1, (n + *p) * sizeof(double));

  for (i = 0; i < n; i++) {
    (*v1)[*p + i] = v2[i];
  }

  (*p) += n;
}

/*	Sort vector	*/
#define SWAP(a, b)                                                             \
  temp = (a);                                                                  \
  (a) = (b);                                                                   \
  (b) = temp;
#define M 7
#define NSTACK 50
void num_f_sort_vector(
    /*	Size of vector	*/
    int p,
    /*	Vector	*/
    double *v) {
  int i, ir = p - 1, j, k, l = 0;
  int jstack = -1, *istack;
  double a, temp;

  istack = (int *)calloc(NSTACK, sizeof(int));

  for (;;) {
    if (ir - l < M) {
      for (j = l + 1; j <= ir; j++) {
        a = v[j];
        for (i = j - 1; i >= l; i--) {
          if (v[i] <= a)
            break;
          v[i + 1] = v[i];
        }
        v[i + 1] = a;
      }

      if (jstack == -1)
        break;

      ir = istack[jstack--];
      l = istack[jstack--];
    } else {
      k = (l + ir) >> 1;
      SWAP(v[k], v[l + 1])
      if (v[l] > v[ir]) {
        SWAP(v[l], v[ir])
      }
      if (v[l + 1] > v[ir]) {
        SWAP(v[l + 1], v[ir])
      }
      if (v[l] > v[l + 1]) {
        SWAP(v[l], v[l + 1])
      }
      i = l + 1;
      j = ir;
      a = v[l + 1];
      for (;;) {
        do
          i++;
        while (v[i] < a);
        do
          j--;
        while (v[j] > a);

        if (j < i)
          break;

        SWAP(v[i], v[j]);
      }

      v[l + 1] = v[j];
      v[j] = a;
      jstack += 2;

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
  free(istack);
}

/*	Sort a matrix: sorting is done against the column col of	*/
/*	all columns are swithed
 */
void num_f_sort_matrix(
    /*	Size of vector			*/
    int p,
    /*	Number of row			*/
    int n,
    /*	Column on which we sort	*/
    int row,
    /*	Vector	*/
    double **matrix) {
  int i, ir = p - 1, j, k, l = 0, m;
  int jstack = -1, *istack;
  double *a = NULL, *v = matrix[row];
  double temp;

  istack = (int *)calloc(NSTACK, sizeof(int));
  a = (double *)calloc(n, sizeof(double));

  for (;;) {
    if (ir - l < M) {
      for (j = l + 1; j <= ir; j++) {
        for (m = 0; m < n; m++) {
          a[m] = matrix[m][j];
        }
        for (i = j - 1; i >= l; i--) {
          if (v[i] <= a[row])
            break;

          for (m = 0; m < n; m++) {
            matrix[m][i + 1] = matrix[m][i];
          }
        }

        for (m = 0; m < n; m++) {
          matrix[m][i + 1] = a[m];
        }
      }

      if (jstack == -1)
        break;

      ir = istack[jstack--];
      l = istack[jstack--];
    } else {
      k = (l + ir) >> 1;
      for (m = 0; m < n; m++) {
        SWAP(matrix[m][k], matrix[m][l + 1])
      }
      if (v[l] > v[ir]) {
        for (m = 0; m < n; m++) {
          SWAP(matrix[m][l], matrix[m][ir])
        }
      }
      if (v[l + 1] > v[ir]) {
        for (m = 0; m < n; m++) {
          SWAP(matrix[m][l + 1], matrix[m][ir])
        }
      }
      if (v[l] > v[l + 1]) {
        for (m = 0; m < n; m++) {
          SWAP(matrix[m][l], matrix[m][l + 1])
        }
      }
      i = l + 1;
      j = ir;
      for (m = 0; m < n; m++) {
        a[m] = matrix[m][l + 1];
      }

      for (;;) {
        do
          i++;
        while (v[i] < a[row]);
        do
          j--;
        while (v[j] > a[row]);

        if (j < i)
          break;

        for (m = 0; m < n; m++) {
          SWAP(matrix[m][i], matrix[m][j]);
        }
      }

      for (m = 0; m < n; m++) {
        matrix[m][l + 1] = matrix[m][j];
        matrix[m][j] = a[m];
      }

      jstack += 2;

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
  free(istack);
  free(a);
}

#undef M
#undef NSTACK
#undef SWAP
#undef NRANSI

/*	Remove identical values	*/
void num_f_unique_vector(
    /*	Size of vector byref	*/
    int *p,
    /*	Vector	*/
    double *v) {
  int i = 0, j;

  while (i < *p - 1) {
    if (fabs(v[i + 1] - v[i]) < MY_DBL_EPSILON) {
      for (j = i; j < *p - 1; j++) {
        v[j] = v[j + 1];
      }
      (*p)--;
    } else {
      i++;
    }
  }
}