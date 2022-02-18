/* ===========================================================
   FILENAME :    utmemory.c

   PURPOSE:      memory allocation functions for:
                       - vectors,
                                           - matrixes,
                                           - f3tensor ...
                 of longs, doubles, strings,...
   =========================================================== */

#include "utallhdr.h"

#define NR_END 1
#define FREE_ARG char*

/* -----------------------------------------------------------------
  DESCRIPTION     :Numerical Recipes in 'C' allocation routines, modified for
SORT purposes in the following ways:
(1) we never use floats.
(2) lngvector (vector of longs) added
(3) calloc and free replaced with srt_calloc and srt_free
(4) nrerror no longer calls exit() since this could kill something important.
rather, these functions all return NULL if they fail.  This is done by
redefining nrerror as a call to a function nrerror_print that prints a message
 and returns NULL. Any
function calling these functions should make appropriate error checks.
(5) static arguments to numerical recipes macros
(which were in nrutil.h and hence in every file) have been made extern and
live in this file only.
(6)
   ---------------------------------------------------------------------*/

#undef nrerror
#define nrerror(ERROR_STRING)        \
    {                                \
        nrerror_print(ERROR_STRING); \
        return NULL;                 \
    }

/* modified not to exit */
static void nrerror_print(char error_text[])
/* Numerical Recipes standard error handler */
{
    FILE* f;
    f = fopen("nrerror.out", "a");
    if (f)
    {
        fprintf(f, "Numerical Recipes run-time error...\n");
        fprintf(f, "%s\n", error_text);
        fclose(f);
    }

    /*	fprintf(stderr,"...now exiting to system...\n");
            exit(1); */
}

double* vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
    double* v;

    v = (double*)srt_calloc(nh - nl + 1 + NR_END, sizeof(double));
    if (!v)
        nrerror("allocation failure in vector()");
    return v - nl + NR_END;
}

int* ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
    int* v;

    v = (int*)srt_calloc(nh - nl + 1 + NR_END, sizeof(int));
    if (!v)
        nrerror("allocation failure in ivector()");
    return v - nl + NR_END;
}

unsigned char* cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
    unsigned char* v;

    v = (unsigned char*)srt_calloc(nh - nl + 1 + NR_END, sizeof(unsigned char));
    if (!v)
        nrerror("allocation failure in cvector()");
    return v - nl + NR_END;
}

unsigned long* lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
    unsigned long* v;

    v = (unsigned long*)srt_calloc(nh - nl + 1 + NR_END, sizeof(long));
    if (!v)
        nrerror("allocation failure in lvector()");
    return v - nl + NR_END;
}

long* lngvector(long nl, long nh)
/* KNL: ansi-friendly version for longs, 5/10/94 	 */
/* allocate a long vector with subscript range v[nl..nh] */
{
    long* v;

    v = (long*)srt_calloc(nh - nl + 1 + NR_END, sizeof(long));
    if (!v)
        nrerror("allocation failure in lngvector()");
    return v - nl + NR_END;
}

double* dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
    double* v;

    v = (double*)srt_calloc(nh - nl + 1 + NR_END, sizeof(double));
    if (!v)
        nrerror("allocation failure in dvector()");
    return v - nl + NR_END;
}

char** svector(long nl, long nh)
/* allocate a String vector with subscript range v[nl..nh] */
{
    char** v;

    v = (char**)srt_calloc(nh - nl + 1 + NR_END, sizeof(char*));
    if (!v)
        nrerror("allocation failure in svector()");
    v += NR_END;
    v -= nl;
    return v;
}

char** svector_size(long nl, long nh, long sz)
/* allocate a vector of Strings (length sz) with subscript range v[nl..nh] */
{
    char** v;
    long   i;

    /* allocate pointers to Strings */
    v = (char**)srt_calloc(nh - nl + 1 + NR_END, sizeof(char*));
    if (!v)
        nrerror("allocation failure in svector_size");
    v += NR_END;
    v -= nl;

    /* allocate Strings and set pointers to them */
    for (i = nl; i <= nh; i++)
    {
        v[i] = (char*)srt_calloc(sz, sizeof(char));
        if (!v[i])
        {
            free((FREE_ARG)(v + nl - NR_END));
            v = NULL;
            nrerror("allocation failure (2) in svector_size");
        }
    }

    return v;
}

double** matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long     i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    double** m;

    /* allocate pointers to rows */
    m = (double**)srt_calloc(nrow + NR_END, sizeof(double*));
    if (!m)
        nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl] = (double*)srt_calloc(nrow * ncol + NR_END, sizeof(double));
    if (!m[nrl])
    {
        free((FREE_ARG)(m + nrl - NR_END));
        m = NULL;
        nrerror("allocation failure 2 in matrix()");
    }
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

double** dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long     i, nrow, ncol;
    double** m;

    nrow = nrh - nrl + 1;
    ncol = nch - ncl + 1;
    m    = NULL;

    /* allocate pointers to rows */
    m = (double**)srt_calloc(nrow + NR_END, sizeof(double*));
    if (!m)
        nrerror("allocation failure 1 in dmatrix()");
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl] = (double*)srt_calloc(nrow * ncol + NR_END, sizeof(double));
    if (!m[nrl])
    {
        free((FREE_ARG)(m + nrl - NR_END));
        m = NULL;
        nrerror("allocation failure 2 in dmatrix()");
    }
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

char*** smatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a string matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long    i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    char*** m;

    /* allocate pointers to rows */
    m = (char***)srt_calloc(nrow + NR_END, sizeof(char**));
    if (!m)
        nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl] = (char**)srt_calloc(nrow * ncol + NR_END, sizeof(char*));
    if (!m[nrl])
    {
        free((FREE_ARG)(m + nrl - NR_END));
        m = NULL;
        nrerror("allocation failure 2 in matrix()");
    }
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

char*** smatrix_size(long nrl, long nrh, long ncl, long nch, long sz)
/* allocate a string matrix (size sl) with subscript range m[nrl..nrh][ncl..nch] */
{
    long    i, j;
    char*** m;

    /* allocate memory for the matrix of pointers to char */
    m = smatrix(nrl, nrh, ncl, nch);
    if (!m)
        nrerror("allocation failure 1 in smatrix_size");

    /* allocate memory for the strings themselves*/
    for (i = nrl; i <= nrh; i++)
    {
        for (j = ncl; j <= nch; j++)
        {
            m[i][j] = (char*)srt_calloc(sz, sizeof(char));
            /* who did this ???  sprintf(m[i][j], "HELLO"); */
            if (!m[i][j])
            {
                free_smatrix_size(m, nrl, nrh, ncl, nch, sz);
                nrerror("allocation failure 2 in smatrix_size");
            }
        }
    }

    return m;
}

int** imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long  i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    int** m;

    /* allocate pointers to rows */
    m = (int**)srt_calloc(nrow + NR_END, sizeof(int*));
    if (!m)
        nrerror("allocation failure 1 in matrix()");

    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl] = (int*)srt_calloc(nrow * ncol + NR_END, sizeof(int));
    if (!m[nrl])
    {
        free((FREE_ARG)(m + nrl - NR_END));
        m = NULL;
        nrerror("allocation failure 2 in matrix()");
    }
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

long** lngmatrix(long nrl, long nrh, long ncl, long nch)
/* KNL: ansi-friendly version added for grf_f_mid_conv_2020.c 	    */
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long   i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    long** m;

    /* allocate pointers to rows */
    m = (long**)srt_calloc(nrow + NR_END, sizeof(long*));
    if (!m)
        nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl] = (long*)srt_calloc(nrow * ncol + NR_END, sizeof(long));
    if (!m[nrl])
    {
        free((FREE_ARG)(m + nrl - NR_END));
        m = NULL;
        nrerror("allocation failure 2 in matrix()");
    }
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

double** submatrix(
    double** a, long oldrl, long oldrh, long oldcl, long oldch, long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
    long     i, j, nrow = oldrh - oldrl + 1, ncol = oldcl - newcl;
    double** m;

    if (newrl != oldrh - oldrl + 1)
        nrerror("wrong no. of rows to submatrix");
    if (newcl != oldch - oldcl + 1)
        nrerror("wrong no. of cols to submatrix");

    /* allocate array of pointers to rows */
    m = (double**)srt_calloc(nrow + NR_END, sizeof(double*));
    if (!m)
        nrerror("allocation failure in submatrix()");
    m += NR_END;
    m -= newrl;

    /* set pointers to rows */
    for (i = oldrl, j = newrl; i <= oldrh; i++, j++)
        m[j] = a[i] + ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

double** convert_matrix(double* a, long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
    long     i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    double** m;

    /* allocate pointers to rows */
    m = (double**)srt_calloc(nrow + NR_END, sizeof(double*));
    if (!m)
        nrerror("allocation failure in convert_matrix()");
    m += NR_END;
    m -= nrl;

    /* set pointers to rows */
    m[nrl] = a - ncl;
    for (i = 1, j = nrl + 1; i < nrow; i++, j++)
        m[j] = m[j - 1] + ncol;
    /* return pointer to array of pointers to rows */
    return m;
}

double*** f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
    long      i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
    double*** t;

    /* allocate pointers to pointers to rows */
    t = (double***)srt_calloc(nrow + NR_END, sizeof(double**));
    if (!t)
        nrerror("allocation failure 1 in f3tensor()");
    t += NR_END;
    t -= nrl;

    /* allocate pointers to rows and set pointers to them */
    t[nrl] = (double**)srt_calloc(nrow * ncol + NR_END, sizeof(double*));
    if (!t[nrl])
    {
        free((FREE_ARG)(t + nrl - NR_END));
        nrerror("allocation failure 2 in f3tensor()");
    }
    t[nrl] += NR_END;
    t[nrl] -= ncl;

    /* allocate rows and set pointers to them */
    t[nrl][ncl] = (double*)srt_calloc(nrow * ncol * ndep + NR_END, sizeof(double));
    if (!t[nrl][ncl])
    {
        free((FREE_ARG)(t[nrl] + ncl - NR_END));
        free((FREE_ARG)(t + nrl - NR_END));
        nrerror("allocation failure 3 in f3tensor()");
    }
    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ndl;

    for (j = ncl + 1; j <= nch; j++)
        t[nrl][j] = t[nrl][j - 1] + ndep;
    for (i = nrl + 1; i <= nrh; i++)
    {
        t[i]      = t[i - 1] + ncol;
        t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
        for (j = ncl + 1; j <= nch; j++)
            t[i][j] = t[i][j - 1] + ndep;
    }

    /* return pointer to array of pointers to rows */
    return t;
}

long*** l3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
    long    i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
    long*** t;

    /* allocate pointers to pointers to rows */
    t = (long***)srt_calloc(nrow + NR_END, sizeof(long**));
    if (!t)
        nrerror("allocation failure 1 in l3tensor()");
    t += NR_END;
    t -= nrl;

    /* allocate pointers to rows and set pointers to them */
    t[nrl] = (long**)srt_calloc(nrow * ncol + NR_END, sizeof(long*));
    if (!t[nrl])
    {
        free((FREE_ARG)(t + nrl - NR_END));
        nrerror("allocation failure 2 in l3tensor()");
    }
    t[nrl] += NR_END;
    t[nrl] -= ncl;

    /* allocate rows and set pointers to them */
    t[nrl][ncl] = (long*)srt_calloc(nrow * ncol * ndep + NR_END, sizeof(long));
    if (!t[nrl][ncl])
    {
        free((FREE_ARG)(t[nrl] + ncl - NR_END));
        free((FREE_ARG)(t + nrl - NR_END));
        nrerror("allocation failure 3 in l3tensor()");
    }
    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ndl;

    for (j = ncl + 1; j <= nch; j++)
        t[nrl][j] = t[nrl][j - 1] + ndep;
    for (i = nrl + 1; i <= nrh; i++)
    {
        t[i]      = t[i - 1] + ncol;
        t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
        for (j = ncl + 1; j <= nch; j++)
            t[i][j] = t[i][j - 1] + ndep;
    }

    /* return pointer to array of pointers to rows */
    return t;
}

/* Allocates f4tensor */
double**** f4tensor(
    long m_min, long m_max, long n_min, long n_max, long o_min, long o_max, long p_min, long p_max)
/* allocate a double 4tensor with range t[m_min..m_max][n_min..n_max][o_min..o_max][p_min..p_max] */
{
    double**** t;
    int        i, j, k;
    long       m = m_max - m_min + 1;
    long       n = n_max - n_min + 1;
    long       o = o_max - o_min + 1;
    long       p = p_max - p_min + 1;

    t = (double****)srt_calloc(m + NR_END, sizeof(double***));
    if (!t)
    {
        return NULL;
    }
    t += NR_END;
    t -= m_min;

    t[m_min] = (double***)srt_calloc(m * n + NR_END, sizeof(double**));
    if (!t[m_min])
    {
        t += m_min;
        t -= NR_END;
        srt_free(t);
        return NULL;
    }
    t[m_min] += NR_END;
    t[m_min] -= n_min;

    t[m_min][n_min] = (double**)srt_calloc(m * n * o + NR_END, sizeof(double*));
    if (!t[m_min][n_min])
    {
        t[m_min] += n_min;
        t[m_min] -= NR_END;
        srt_free(t[m_min]);
        t += m_min;
        t -= NR_END;
        srt_free(t);
        return NULL;
    }
    t[m_min][n_min] += NR_END;
    t[m_min][n_min] -= o_min;

    t[m_min][n_min][o_min] = (double*)srt_calloc(m * n * o * p + NR_END, sizeof(double));
    if (!t[m_min][n_min][o_min])
    {
        t[m_min][n_min] += o_min;
        t[m_min][n_min] -= NR_END;
        srt_free(t[m_min][n_min]);
        t[m_min] += n_min;
        t[m_min] -= NR_END;
        srt_free(t[m_min]);
        t += m_min;
        t -= NR_END;
        srt_free(t);
        return NULL;
    }
    t[m_min][n_min][o_min] += NR_END;
    t[m_min][n_min][o_min] -= p_min;

    for (k = o_min + 1; k <= o_max; k++)
    {
        t[m_min][n_min][k] = t[m_min][n_min][k - 1] + p;
    }

    for (j = n_min + 1; j <= n_max; j++)
    {
        t[m_min][j] = t[m_min][j - 1] + o;
        for (k = o_min; k <= o_max; k++)
        {
            t[m_min][j][k] = t[m_min][j - 1][k] + o * p;
        }
    }

    for (i = m_min + 1; i <= m_max; i++)
    {
        t[i]        = t[i - 1] + n;
        t[i][n_min] = t[i - 1][n_min] + n * o;
        for (k = o_min; k <= o_max + 1; k++)
        {
            t[i][n_min][k] = t[i - 1][n_min][k] + n * o * p;
        }
        for (j = n_min + 1; j <= n_max; j++)
        {
            t[i][j]        = t[i][j - 1] + o;
            t[i][j][o_min] = t[i][j - 1][o_min] + o * p;
            for (k = o_min + 1; k <= o_max; k++)
            {
                t[i][j][k] = t[i][j][k - 1] + p;
            }
        }
    }
    return t;
}

void free_vector(double* v, long nl, long nh)
/* free a double vector allocated with vector() */
{
    free((FREE_ARG)(v + nl - NR_END));
    v = NULL;
}

void free_ivector(int* v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
    free((FREE_ARG)(v + nl - NR_END));
    v = NULL;
}

void free_cvector(unsigned char* v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
    free((FREE_ARG)(v + nl - NR_END));
    v = NULL;
}

void free_lvector(unsigned long* v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
    free((FREE_ARG)(v + nl - NR_END));
    v = NULL;
}

void free_lngvector(long* v, long nl, long nh)
/* KNL: ansi-friendly version for longs, 5/10/94 */
/* free a long vector allocated with lngvector() */
{
    free((FREE_ARG)(v + nl - NR_END));
    v = NULL;
}

void free_dvector(double* v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
    free((FREE_ARG)(v + nl - NR_END));
    v = NULL;
}

void free_svector(String* v, long nl, long nh)
/* free a String vector allocated with svector() */
{
    free((FREE_ARG)(v + nl - NR_END));
    v = NULL;
}

void free_svector_size(String* v, long nl, long nh, long sz)
/* free a String vector allocated with svector() */
{
    int i;

    for (i = nl; i <= nh; i++)
    {
        free((FREE_ARG)(v[i]));
        v[i] = NULL;
    }
    free((FREE_ARG)(v + nl - NR_END));
    v = NULL;
}

void free_matrix(double** m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by matrix() */
{
    free((FREE_ARG)(m[nrl] + ncl - NR_END));
    free((FREE_ARG)(m + nrl - NR_END));
    m = NULL;
}

void free_dmatrix(double** m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
    free((FREE_ARG)(m[nrl] + ncl - NR_END));
    free((FREE_ARG)(m + nrl - NR_END));
    m = NULL;
}

void free_smatrix(char*** m, long nrl, long nrh, long ncl, long nch)
/* free a string matrix allocated by smatrix() */
{
    free((FREE_ARG)(m[nrl] + ncl - NR_END));
    free((FREE_ARG)(m + nrl - NR_END));
    m = NULL;
}

void free_smatrix_size(char*** m, long nrl, long nrh, long ncl, long nch, long sz)
/* free a string matrix allocated by smatrix_size() */
{
    int i, j;
    for (i = nrl; i <= nrh; i++)
    {
        for (j = ncl; j <= nch; j++)
        {
            free((FREE_ARG)(m[i][j]));
            m[i][j] = NULL;
        }
    }
    free_smatrix(m, nrl, nrh, ncl, nch);
    m = NULL;
}

void free_imatrix(int** m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
    free((FREE_ARG)(m[nrl] + ncl - NR_END));
    free((FREE_ARG)(m + nrl - NR_END));
    m = NULL;
}

void free_lngmatrix(long** m, long nrl, long nrh, long ncl, long nch)
/* KNL: ansi-friendly version for longs        */
/* free a long matrix allocated by lngmatrix() */
{
    free((FREE_ARG)(m[nrl] + ncl - NR_END));
    free((FREE_ARG)(m + nrl - NR_END));
    m = NULL;
}

void free_submatrix(double** b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
    free((FREE_ARG)(b + nrl - NR_END));
    b = NULL;
}

void free_convert_matrix(double** b, long nrl, long nrh, long ncl, long nch)
/* free a matrix allocated by convert_matrix() */
{
    free((FREE_ARG)(b + nrl - NR_END));
    b = NULL;
}

void free_f3tensor(double*** t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* free a double f3tensor allocated by f3tensor() */
{
    free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
    free((FREE_ARG)(t[nrl] + ncl - NR_END));
    free((FREE_ARG)(t + nrl - NR_END));
    t = NULL;
}

void free_l3tensor(long*** t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* free a double f3tensor allocated by f3tensor() */
{
    free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
    free((FREE_ARG)(t[nrl] + ncl - NR_END));
    free((FREE_ARG)(t + nrl - NR_END));
    t = NULL;
}

/* Frees f4tensor */
void free_f4tensor(
    double**** t,
    long       m_min,
    long       m_max,
    long       n_min,
    long       n_max,
    long       o_min,
    long       o_max,
    long       p_min,
    long       p_max)
{
    free((FREE_ARG)(t[m_min][n_min][o_min] + p_min - NR_END));
    free((FREE_ARG)(t[m_min][n_min] + o_min - NR_END));
    free((FREE_ARG)(t[m_min] + n_min - NR_END));
    free((FREE_ARG)(t + m_min - NR_END));

    t = NULL;
}

#undef NR_END
#undef FREE_ARG
#undef nerror
