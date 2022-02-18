/* ======================================================
   FILENAME:  utmemory.h

   PURPOSE:   memory allocation functions for vectors,
              arrays, f3tensor...
   ====================================================== */

#ifndef UTMEMORY_H
#define UTMEMORY_H

#define SQR(a) (a == 0.0 ? 0.0 : a * a)
#define DSQR(a) (a == 0.0 ? 0.0 : a * a)
#define DMAX(a, b) (a > b ? a : b)
#define DMIN(a, b) (a < b ? a : b)
#define FMAX(a, b) (a > b ? a : b)
#define FMIN(a, b) (a < b ? a : b)
#define LMAX(a, b) (a > b ? a : b)
#define LMIN(a, b) (a < b ? a : b)
#define IMAX(a, b) (a > b ? a : b)
#define IMIN(a, b) (a < b ? a : b)
#ifndef SIGN
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#endif
/* To be find with "find in file":  PVMPI*/
/* The following lines must be included in changes */
#ifdef __cplusplus
extern "C"
{
#endif

    double*        vector(long nl, long nh);
    int*           ivector(long nl, long nh);
    unsigned char* cvector(long nl, long nh);
    unsigned long* lvector(long nl, long nh);
    long*          lngvector(long nl, long nh);
    double*        dvector(long nl, long nh);
    char**         svector(long nl, long nh);
    char**         svector_size(long nl, long nh, long strsize);

    double** matrix(long nrl, long nrh, long ncl, long nch);
    double** dmatrix(long nrl, long nrh, long ncl, long nch);
    char***  smatrix(long nrl, long nrh, long ncl, long nch);
    char***  smatrix_size(long nrl, long nrh, long ncl, long nch, long sz);
    int**    imatrix(long nrl, long nrh, long ncl, long nch);
    long**   lngmatrix(long nrl, long nrh, long ncl, long nch);
    double** submatrix(
        double** a, long oldrl, long oldrh, long oldcl, long oldch, long newrl, long newcl);
    double** convert_matrix(double* a, long nrl, long nrh, long ncl, long nch);

    double*** f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
#define dcube f3tensor

    long*** l3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);

    void free_vector(double* v, long nl, long nh);
    void free_ivector(int* v, long nl, long nh);
    void free_cvector(unsigned char* v, long nl, long nh);
    void free_lvector(unsigned long* v, long nl, long nh);
    void free_lngvector(long* v, long nl, long nh);
    void free_dvector(double* v, long nl, long nh);
    void free_svector(char** v, long nl, long nh);
    void free_svector_size(char** v, long nl, long nh, long size);

    void free_matrix(double** m, long nrl, long nrh, long ncl, long nch);
    void free_dmatrix(double** m, long nrl, long nrh, long ncl, long nch);
    void free_smatrix(char*** m, long nrl, long nrh, long ncl, long nch);
    void free_smatrix_size(char*** m, long nrl, long nrh, long ncl, long nch, long sz);
    void free_imatrix(int** m, long nrl, long nrh, long ncl, long nch);
    void free_lngmatrix(long** m, long nrl, long nrh, long ncl, long nch);
    void free_submatrix(double** b, long nrl, long nrh, long ncl, long nch);
    void free_convert_matrix(double** b, long nrl, long nrh, long ncl, long nch);
    void free_f3tensor(double*** t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
#define free_dcube free_f3tensor

    void free_l3tensor(long*** t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);

    void free_f4tensor(
        double**** t,
        long       m_min,
        long       m_max,
        long       n_min,
        long       n_max,
        long       o_min,
        long       o_max,
        long       p_min,
        long       p_max);

    double**** f4tensor(
        long m_min,
        long m_max,
        long n_min,
        long n_max,
        long o_min,
        long o_max,
        long p_min,
        long p_max);

/* To be find with "find in file":  PVMPI*/
/* The following lines must be included in changes */
#ifdef __cplusplus
}
#endif

#endif /* _NR_UTILS_H_ */
