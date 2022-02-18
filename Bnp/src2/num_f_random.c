/* ==========================================================================
   FILE_NAME:	num_f_random.c

   PURPOSE:     random numbers generation
   ========================================================================== */

#include "math.h"
#include "num_h_allhdr.h"
#include "num_h_random.h"

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
#define NTAB 32
#define NDIV (1 + IMM1 / NTAB)
#define MC_GENRAND_EPS 1.2e-7
#define RNMX (1.0 - MC_GENRAND_EPS)

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0 / MBIG)

double uniform(long* idum)
/* ran2 from numerical recipes in C */
{
    long        j;
    long        k;
    static long idum2 = 123456789;
    static long iy    = 0;
    static long iv[NTAB];
    double      temp;

    if (*idum <= 0)
    {
        if (-(*idum) < 1)
            *idum = 1;
        else
            *idum = -(*idum);
        idum2 = (*idum);
        for (j = NTAB + 7; j >= 0; j--)
        {
            k     = (*idum) / IQ1;
            *idum = IA1 * (*idum - k * IQ1) - k * IR1;
            if (*idum < 0)
                *idum += IM1;
            if (j < NTAB)
                iv[j] = *idum;
        }
        iy = iv[0];
    }
    k     = (*idum) / IQ1;
    *idum = IA1 * (*idum - k * IQ1) - k * IR1;
    if (*idum < 0)
        *idum += IM1;
    k     = idum2 / IQ2;
    idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
    if (idum2 < 0)
        idum2 += IM2;
    j     = iy / NDIV;
    iy    = iv[j] - idum2;
    iv[j] = *idum;
    if (iy < 1)
        iy += IMM1;

    if ((temp = AM * iy) > RNMX)
        return RNMX;
    else
        return temp;
}

/* ----------------------------------------------------------------------------- */
/* Not so precise but MUCH MUCH faster... */
double uniform_fast(long* idum)
{
    static long inext, inextp;
    static long ma[56];
    static long iff = 0;
    long        mj, mk;
    long        i, ii, k;

    if (*idum < 0 || iff == 0)
    {
        iff = 1;
        mj  = labs(MSEED - labs(*idum));
        mj %= MBIG;
        ma[55] = mj;
        mk     = 1;
        for (i = 1; i <= 54; i++)
        {
            ii     = (21 * i) % 55;
            ma[ii] = mk;
            mk     = mj - mk;
            if (mk < MZ)
                mk += MBIG;
            mj = ma[ii];
        }
        for (k = 1; k <= 4; k++)
            for (i = 1; i <= 55; i++)
            {
                ma[i] -= ma[1 + (i + 30) % 55];
                if (ma[i] < MZ)
                    ma[i] += MBIG;
            }
        inext  = 0;
        inextp = 31;
        *idum  = 1;
    }
    if (++inext == 56)
        inext = 1;
    if (++inextp == 56)
        inextp = 1;
    mj = ma[inext] - ma[inextp];
    if (mj < MZ)
        mj += MBIG;
    ma[inext] = mj;
    return mj * FAC;
}

/* --------------------------------------------------------------------------
   Returns a random integer from [0,n]
   Used to generate random permutations in BalSam
   -------------------------------------------------------------------------- */

int random_int(int n, long* seed)
{
    int    i;
    double x, z;

    z = (double)(n);
    x = (1 + z) * uniform_fast(seed);
    i = (int)x;
    if ((double)i > x)
        --i;
    return i;
}

/* -------------------------------------------------------------------------
   The following functions are used to generate independent random Gaussian
   variables, from a single value to vectors to matrixes
   ------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------
        see numerical recipes in C
        returns a sample from N(0,1)
   ------------------------------------------------------------------------- */

static int iset = 0;

Err gauss_box_muller_init(void)
{
    iset = 0;

    return NULL;
}

double gauss_box_muller(long* seed)
{
    static double gset;
    double        fac, rsq, v1, v2;

    if (iset == 0)
    {
        do
        {
            v1 = 2.0 * uniform_fast(seed) - 1.0;
            v2 = 2.0 * uniform_fast(seed) - 1.0;

            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        fac  = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1 * fac;
        iset = 1;
        return (v2 * fac);
    }
    else
    {
        iset = 0;
        return gset;
    }
}

/* -------------------------------------------------------------------------
  Returns a randon Gaussian variable
        (faster than Gauss Box Muller and more efficient)
   ------------------------------------------------------------------------- */

double gauss_sample(long* seed)
{
    return inv_cumnorm_fast(uniform(seed));
}

/* -------------------------------------------------------------------------
  Returns a vector of randon Gaussian variables
                                        v[vl..vh]
   ------------------------------------------------------------------------- */

Err gauss_vector(double* v, long vl, long vh, long* seed)
{
    Err  err = NULL;
    long i;

    err = gauss_box_muller_init();

    for (i = vl; i <= vh; i++)
        v[i] = gauss_box_muller(seed);
    /* 		v[i] = gauss_sample(seed); */

    return err;
}

/* -------------------------------------------------------------------------
   Returns a matrix of random Gaussian variables :
                                        v[pl..ph][sl..sh]
  (Allocation is done outside)
  -------------------------------------------------------------------------- */
Err gauss_matrix(double** v, long pl, long ph, long sl, long sh, long* seed)
{
    Err  err = NULL;
    long i, j;

    err = gauss_box_muller_init();

    for (i = pl; i <= ph; i++)
    {
        for (j = sl; j <= sh; j++)
            v[i][j] = gauss_box_muller(seed);
        /*			v[i][j] = gauss_sample(seed); */
    }
    return err;
}

/* -------------------------------------------------------------------------
   Returns a cube of random Gaussian variables :
                                        v[pl..ph][bl..bh][sl..sh]
  (Allocation is done outside)
  -------------------------------------------------------------------------- */
Err gauss_cube(double*** v, long pl, long ph, long bl, long bh, long sl, long sh, long* seed)
{
    Err  err = NULL;
    long i, j, k;

    err = gauss_box_muller_init();

    for (i = pl; i <= ph; i++)
    {
        for (j = bl; j <= bh; j++)
        {
            for (k = sl; k <= sh; k++)
                v[i][j][k] = gauss_box_muller(seed);
            /*				v[i][j][k] = gauss_sample(seed);*/
        }
    }

    return err;

} /* END Err gauss_cube(...) */

/* -------------------------------------------------------------------------------------- */

/* Generation of a set of antithetic Gausssian samples  */

Err gauss_anti_cube(double*** v, long pl, long ph, long bl, long bh, long sl, long sh, long* seed)
{
    Err  err = NULL;
    long i, j, k;
    long num_path;
    long middle;
    long extra;

    /* See if there is an even or odd number of paths */
    num_path = pl - ph + 1;
    if (num_path % 2 != 0)
        extra = 1;
    else
        extra = 0;

    /* Offset from which the symmetry will be performed */
    middle = (num_path + extra) / 2;

    /* Generates the first half ( + 1 point if odd) of rand[i] */
    err = gauss_cube(v, pl, pl + middle, bl, bh, sl, sh, seed);
    if (err)
        return err;

    /* Makes the symetry (path wise) */
    for (i = pl; i < pl + middle - extra; i++)
    {
        for (j = bl; j <= bh; j++)
        {
            for (k = sl; k <= sh; k++)
            {
                v[i + middle][j][k] = -v[i][j][k];

            } /* END of loop on steps */

        } /* END of loop on Brownian */

    } /* END of loop on paths */

    /* Return a success message */
    return NULL;

} /* END Err gauss_anti_cube(...) */

/* -------------------------------------------------------------------------
   Returns a matrix of balanced Gaussian samples.
   The matrix is then described as follows:
                                        m[path][step]
   where:
                - pathindexstart <= path <= pathindexend
                - stepindexstart <= step <= stepindexend
  -------------------------------------------------------------------------- */

Err BalSampMatrix(
    double** m,
    int      pathindexstart,
    int      pathindexend,
    int      stepindexstart,
    int      stepindexend,
    long*    seed)
{
    Err    err      = NULL;
    int    path_num = pathindexend - pathindexstart + 1;
    int    step_num = stepindexend - stepindexstart + 1;
    int    i, j, k;
    int    numb_seg;
    double rand_point, rand_val;
    double dx, *sample = NULL, *temp_samp = NULL;

    sample = dvector(0, path_num);
    if (!sample)
        return serror("Allocation error in BalSampMatrix");
    temp_samp = dvector(0, path_num);
    if (!temp_samp)
        return serror("Allocation error in BalSampMatrix");
    numb_seg = path_num - pathindexstart + 1;

    /* Only splits the [0.5;1] segment and then takes the symetric **/
    if ((numb_seg % 2) == 1)
    {
        numb_seg -= 1;
        sample[numb_seg] = 0;
    }
    numb_seg = (int)floor(numb_seg / 2);
    dx       = 0.5 / (double)numb_seg;
    rand_val = 0;

    for (i = 0; i < numb_seg; i++)
    {
        rand_point = 0.5 + (double)(i + 0.5) * dx;
        /*Take middle point of segment */
        rand_val = inv_cumnorm_fast(rand_point);

        sample[(int)(2 * i)]     = rand_val;
        sample[(int)(2 * i + 1)] = -rand_val;
    }

    for (k = stepindexstart; k <= stepindexend; k++)
    {
        for (i = 0; i < path_num; i++)
            temp_samp[i] = sample[i];

        for (i = 0; i < path_num; i++)
        {
            /* uniform randomint for remaining numbers not choosen */
            j = random_int(path_num - 1 - i, seed);

            m[i + pathindexstart][k] = temp_samp[j];
            temp_samp[j]             = temp_samp[path_num - 1 - i];
        }
    }

    if (temp_samp)
        free_dvector(temp_samp, 0, path_num);
    temp_samp = NULL;
    if (sample)
        free_dvector(sample, 0, path_num);
    sample = NULL;

    return err;

} /* END Err BalSampMatrix(...) */

/* -------------------------------------------------------------------------
   The generation of a Cube of Balanced Gaussian samples
   The cube is described as follows:
                          rand[path][brownian][step]
   where:
                - pathindexstart     <= path     <= pathindexstart
                - stepindexstart     <= step     <= stepindexend
                - browindexstart     <= brow     <= browindexend
   Memory allocation for the cube is done OUTSIDE
  -------------------------------------------------------------------------- */

Err BalSampCube(
    double*** rand,
    long      pathindexstart,
    long      pathindexend,
    long      browindexstart,
    long      browindexend,
    long      stepindexstart,
    long      stepindexend,
    long*     seed)
{
    Err    err      = NULL;
    int    path_num = pathindexend - pathindexstart + 1;
    int    step_num = stepindexend - stepindexstart + 1;
    int    i, j, k;
    int    numb_seg;
    double rand_point, rand_val;
    double dx, *sample = NULL, *temp_samp = NULL;
    long   ibr;

    /* Memory allocation for local variables */
    sample    = dvector(0, path_num);
    temp_samp = dvector(0, path_num);
    numb_seg  = path_num;

    /* If ODD number of paths, make the last sample as a draw of 0 (central symmetry) */
    if ((numb_seg % 2) == 1)
    {
        numb_seg -= 1;
        sample[numb_seg] = 0;
    }

    /* Makes a set of uniformaly distributed Gaussian sample (on [.5;1] and symmetric on [0; .5] )
     */
    numb_seg = (int)floor(numb_seg / 2);
    dx       = 0.5 / (double)numb_seg;

    for (i = 0; i < numb_seg; i++)
    {
        /* Take middle point of segment for probability calculation */
        rand_point = 0.5 + (double)(i + 0.5) * dx;
        rand_val   = inv_cumnorm_fast(rand_point);

        /* Fills in the uniform Gaussian sample in a symmetrical manner */
        sample[(int)(2 * i)]     = rand_val;
        sample[(int)(2 * i + 1)] = -rand_val;
    }

    /* Loop on the number of Brownians */
    for (ibr = browindexstart; ibr <= browindexend; ibr++)
    {
        /* Loop on all the time steps */
        for (k = stepindexstart; k <= stepindexend; k++)
        {
            /* Initialises temp_samp as sample (uniform Gaussian mesh) */
            memcpy(temp_samp, sample, path_num * sizeof(double));

            /* Loop on all paths */
            for (i = 0; i < path_num; i++)
            {
                /* Chooses randomly the sample that will be affected (uniform choice on remaining )
                 */
                j = random_int(path_num - 1 - i, seed);

                /* Affects the randomly selected sample to the current path (storage in the cube) */
                rand[i + pathindexstart][ibr][k] = temp_samp[j];

                /* Copies the last sample value to the chosen sample (for the proper uniform choice)
                 */
                temp_samp[j] = temp_samp[path_num - 1 - i];

            } /* END of loop on paths */

        } /* END of loop on time steps */

    } /* END of loop on Brownians */

    /* Free memory for the local variables */
    if (temp_samp)
        free_dvector(temp_samp, 0, path_num);
    temp_samp = NULL;
    if (sample)
        free_dvector(sample, 0, path_num);
    sample = NULL;

    return err;
}

/* ========================================================================

        The following two functions are used to generate correlated
        random numbers from non correlated ones, inverting the
        correlation matrix using the method given in Numerical Recipes
        in C

   ========================================================================= */

/** compute from a correlation matrix the linear coefficients necessary
        to build a set of correlated random variables
    we write that
        random[i] = sum(j<=i) alpha[i][j]*gauss[j]
    the formulae used to compute the coefficients is thus:
        (1) alpha[i][j<i] = (correl[i][j]-sum(k<j) alpha[i][k]*alpha[j][k])
                        / alpha[j][j]
        (2) alpha[i][i] = sqrt( 1 - sum(j<i) alpha[j]^2 )
    in a second step, in order to avoid allocating a new ***double,
    we write
        random[i] = sum(j<i) coeff[i][j]*random[j] + alpha[i][i]*gauss[i]
    and compute the coeff[i][j] using the invert of a triangle matrix

**/

Err compute_coeff_from_correl(double** correl, long nbr, double** coeff)
{
    long     i, j, k, dim;
    double   norm, sum;
    double **alpha, **invert;

    alpha  = dmatrix(0, nbr - 1, 0, nbr - 1);
    invert = dmatrix(0, nbr - 1, 0, nbr - 1);

    sum = 0;
    for (i = 0; i < nbr; i++)
    {
        for (j = 0; j < i; j++)
            sum += fabs(correl[i][j]);
    }
    for (i = 0; i < nbr; i++)
    {
        for (j = 0; j < nbr; j++)
            coeff[i][j] = 0;
        coeff[i][i] = 1;
    }
    if (sum == 0)
        return NULL;
    /* First compute the alpha(dim;dim) matrix */
    for (i = 0; i < nbr; i++)
    {
        norm = 1;
        for (j = 0; j < i; j++)
        {
            sum = 0;
            for (k = 0; k < j; k++)
                sum += alpha[i][k] * alpha[j][k];
            if (alpha[j][j] == 0.0)
            {
                alpha[i][j] = 0.0;
                if (correl[i][j] != sum)
                    return (serror("Singularity in correlation matrix"));
            }
            else
            {
                alpha[i][j] = (correl[i][j] - sum) / alpha[j][j];
            }
            norm -= alpha[i][j] * alpha[i][j];
        }
        if (norm < 0)
            return serror("Correlation matrix is improperly defined...");
        alpha[i][i] = sqrt(norm);
    }

    for (dim = 0; dim < nbr; dim++)
    {
        /* Now invert the alpha(dim-1;dim-1) matrixes */
        for (i = 0; i < dim; i++)
        {
            for (j = 0; j < i; j++)
            {
                sum = 0;
                if (alpha[i][i] == 0.0)
                {
                    invert[i][j] = 0.0;
                }
                else
                {
                    for (k = j; k < i; k++)
                        sum += alpha[i][k] * invert[k][j];
                    invert[i][j] = -sum / alpha[i][i];
                }
            }
            if (alpha[i][i] == 0.0)
            {
                invert[i][i] = 0.0;
            }
            else
            {
                invert[i][i] = 1 / alpha[i][i];
            }
        }
        /* And compute the product with the alpha coefficients*/
        for (i = 0; i < dim; i++)
        {
            sum = 0;

            /*	bug-fix: used to be j<=dim which is not set */

            for (j = i; j < dim; j++)
                sum += alpha[dim][j] * invert[j][i];
            coeff[dim][i] = sum;
        }
        coeff[dim][dim] = alpha[dim][dim];
    }
    free_dmatrix(alpha, 0, nbr - 1, 0, nbr - 1);
    free_dmatrix(invert, 0, nbr - 1, 0, nbr - 1);

    return NULL;
}

/* ======================================================================= */

/* -----------------------------------------------------------------------
   From a set of independent gausssian samples **r, makes a set of
   correlated gaussian samples, using the linear coefficients **coeff
   derived from the correlation matrix
   The formula used is the following one:
    gauss_cor[i] = sum(j<i) coeff[i][j]*gauss_cor[j] + alpha[i][i]*gauss[j]
   Please note that
                (1) alpha[i][i] = coeff[i][i]
                (2) alpha[0][0] = 1  and is not in the loop...
   The r matrix is described as follows:
                        r[brown][step]
   and the function will modify the following:
                                        r[0..nbr-1][sl..sh]
   ----------------------------------------------------------------------- */

Err correl_random(double** r, long sl, long sh, long nbr, double** coeff)
{
    Err    err = NULL;
    long   i, j, k;
    double sum;

    for (i = 1; i < nbr; i++)
    {
        for (j = sl; j <= sh; j++)
        {
            sum = 0;
            for (k = 0; k <= i; k++)
                sum += coeff[i][k] * r[k][j];
            r[i][j] = sum;
        }
    }
    return err;
}

/* ======================================================================== */

/** compute from a correlation matrix the linear coefficients necessary
        to build a set of correlated random variables
    we use the Numerical Recipes in C subroutine modified by O. Van Eyseren,
    diagonalise_symmetric_matrix (in utdiagonalise.c), which is called from
        the subroutine below. **/

Err compute_eigen_from_correl(double** correl, long nbr, double** coeff)
{
    Err     err = NULL;
    long    i, j;
    double  sum;
    double* eigenval;

    eigenval = dvector(0, nbr - 1);
    sum      = 0;
    for (i = 0; i < nbr; i++)
    {
        for (j = 0; j < i; j++)
            sum += fabs(correl[i][j]);
    }
    for (i = 0; i < nbr; i++)
    {
        for (j = 0; j < nbr; j++)
            coeff[i][j] = 0;
        coeff[i][i] = 1;
    }
    if (sum == 0)
        return NULL;

    if (err = diagonalise_symmetric_matrix(correl, nbr, eigenval, coeff))
        return err;

    /* Checks the eigen values and make sure they are positive ! */
    for (j = 0; j < nbr; j++)
    {
        if (eigenval[j] < 0)
        {
            free_dvector(eigenval, 0, nbr - 1);
            return serror("Correlation matrix is not DEFINITE POSITIVE");
        }
    }

    /* And compute the product of the square root of the eigenval
        with the  coefficients, coeff*/

    for (j = 0; j < nbr; j++)
    {
        sum = 0;
        /* Square root calculation fails if the values are too small */
        if (eigenval[j] > 1e-16)
        {
            sum = sqrt(eigenval[j]);
        }
        for (i = 0; i < nbr; i++)
            coeff[i][j] *= sum;
    }

    free_dvector(eigenval, 0, nbr - 1);

    return NULL;
}

/* ------------------------------------------------------------------------ */

/* -----------------------------------------------------------------------
   From a set of independent gaussian samples **r, makes a set of
   correlated gaussian samples, using the eigenvalues and eigenvectors of
   correlation matrix.. We know that X=MR where M'M = C, the correlation
   matrix and R is the random matrix. Rewriting C=PDP' and M=sqrt(D)P.
   So:
    gauss_cor[i]= sum(k) m[i][k]*gauss_noncor[k]
   Now we should put m[i][k]=sqrt(d[i])*p[i][k] and this should be the value in
   coeff.
   A modification to the original programme below is needed as we cannot put
   the new  r matrix back into the old, so we create rr.
   The r matrix is described as follows:
                        r[brown][step]
   and the function will modify the following:
                                        r[0..nbr-1][sl..sh]
   ----------------------------------------------------------------------- */

Err correl_random_eigen(double** r, long sl, long sh, long nbr, double** coeff)
{
    Err      err = NULL;
    long     i, j, k;
    double   sum;
    double** rr;

    rr = dmatrix(0, nbr - 1, sl, sh);
    for (i = 0; i < nbr; i++)
    {
        for (j = sl; j <= sh; j++)
        {
            sum = 0;
            for (k = 0; k < nbr; k++)
                sum += coeff[i][k] * r[k][j];
            rr[i][j] = sum;
        }
    }
    for (i = 0; i < nbr; i++)
    {
        for (j = sl; j <= sh; j++)
        {
            r[i][j] = rr[i][j];
        }
    }

    free_dmatrix(rr, 0, nbr - 1, sl, sh);
    return err;
}

/* ======================================================================== */

/* ========================================================================= */

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef MC_GENRAND_EPS
#undef RNMX
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
