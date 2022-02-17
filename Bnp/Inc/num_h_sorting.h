#if !defined(_NUMHSORTING_H_)
#define _NUMHSORTING_H_

#define MPAR 7
#define MNUM 7
#define NSTACK 50
#define DSWAP(a, b)                                                            \
  temp = (a);                                                                  \
  (a) = (b);                                                                   \
  (b) = temp
#define LSWAP(a, b)                                                            \
  itemp = (a);                                                                 \
  (a) = (b);                                                                   \
  (b) = itemp

#pragma warning(disable : 4786)

void sort_simple(unsigned long n, double *arr, double *brr);

void index_data(unsigned long n, double *arr, unsigned long *indx);

void rank_data(unsigned long n,
               unsigned long *indx,   // input: the indexed table
               unsigned long *irank); // output: the rank table

void SortSeries(long n, double *ra);

void SortSeriesDescending(long n, double *ra);

void sort_simple2(unsigned long n, double **arr);

void shell(unsigned long n, double **a);

Err indexx_ll(long *f, long *pI, long n);
Err indexx_dl(double *f, long *pI, long n);
Err qsort_d(double *arr, long n);

#endif