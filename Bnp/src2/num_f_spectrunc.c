/* ======================================================================
   FILENAME		: num_f_spectrunc.c

   PURPOSE		: generate randon numbers that sample properly the main
                  components of a Brownian motion (in terms of covariance)
                                  This function requires that memory has already been
                                  allocated and set to 0 for rand[path][brow][step]
                                  For several Brownian's, the same covar matrix is used,
                                  and the generation is made in the following order for
                                  the Sobol dimensions :
                                        - Eigen Vector 1:
                                                - Brownian 1
                                                - Brownian 2
                                                - ...
                                                - Brownian nbr

                                        - Eigen Vector 2:
                                                - Brownian 1
                                                - Brownian 2
                                                - ...
                                                - Brownian nbr
                  The i th step of the j th path for the k th Brownian motion
                                  is affected by all the factors as follows:

                                rand[j][k][i] = sum[l/factors] ( rand[j][l][i]* fac_coord[l] )

  ====================================================================== */

#include "math.h"
#include "num_h_allhdr.h"

static Err build_Brownian_covar_matrix(
    double* time_at_step, long sl, long sh, double** covar_matrix, long size);

/* -------------------------------------------------------------------------- */
/* --------------------------------------------------------------------------
   Given a covariance matrix, returns Brownian paths that have been generated
   buy a spectral decomposition (i.e. a PCA). The number of factors used to
   generate the Brownian paths are determined according to the precision
   parameter, that corresponds to the residual allowed variance (0.05% seems a
   good choice)
   The random realisations of the factors are generated using Sobol.
   The cube is filled as follows: rand[path][brow][step]
   -------------------------------------------------------------------------- */
Err SpecTruncCube(
    double*** rand,
    double*   time_at_steps, /* From [first_step] to [last_step] at least */
    double    precision,
    long      first_path,
    long      last_path,
    long      first_brow,
    long      last_brow,
    long      first_step,
    long      last_step,
    long*     seed)
{
    double*** temp_rand;
    double**  covar_matrix;
    double*   eigen_val;
    double**  eigen_vec;
    double    stdev;
    double    current_cum_var;
    double    total_cum_var;
    double    random;
    long      nfact;
    long      npath;
    long      nbrow;
    long      nstep;
    long      i;
    long      j;
    long      k;
    long      l;
    long      size;
    int       nrot;
    Err       err = NULL;
    long      sobol_fact;

    /* Set up a few usefull values */
    npath = last_path - first_path + 1;
    nbrow = last_brow - first_brow + 1;
    nstep = last_step - first_step + 1;

    /* Memory allocation for the covariance matrix */
    size         = nstep;
    covar_matrix = dmatrix(0, size - 1, 0, size - 1);

    /* Make the covariance matrix for the required time steps */
    err = build_Brownian_covar_matrix(
        time_at_steps, /* From [0] to [size -1] (sorted) */
        first_step,
        last_step,
        covar_matrix,
        size);
    if (err)
        return err;

    /* Memory allocation for the eigen values and eigen vectors */
    eigen_val = dvector(0, size - 1);
    eigen_vec = dmatrix(0, size - 1, 0, size - 1);

    /* Diagonalise the covariance matrix to extract eigen vectors and values */
    if (size < 10)
        err = jacobi_diagonalisation(covar_matrix, size, eigen_val, eigen_vec, &nrot);
    else
        err = diagonalise_symmetric_matrix(covar_matrix, size, eigen_val, eigen_vec);

    if (err)
    {
        free_dvector(eigen_val, 0, size - 1);
        free_dmatrix(eigen_vec, 0, size - 1, 0, size - 1);
        free_dmatrix(covar_matrix, 0, size - 1, 0, size - 1);
        return err;
    }

    /* Free the covariance matrix (not needed anymore) */
    free_dmatrix(covar_matrix, 0, size - 1, 0, size - 1);

    /* Determines how many factors are required to lock precision in variance */
    total_cum_var = 0.0;
    for (i = 0; i < size; i++)
        total_cum_var += eigen_val[i];
    current_cum_var = eigen_val[0];
    i               = 1;
    while ((1.0 - (current_cum_var / total_cum_var) > precision))
    {
        current_cum_var += eigen_val[i];
        i++;
    }
    nfact = i;

    /* Multiplies each of the used eigen vectors by the sqrt of its eigen value */
    for (i = 0; i < nfact; i++)
    {
        stdev = sqrt(eigen_val[i]);
        for (j = 0; j < size; j++)
            eigen_vec[j][i] *= stdev;
    }

    /* Create a temporary cube [path][EIGENVEC][BROW] for random numbers*/
    temp_rand = dcube(0, npath - 1, 0, nbrow - 1, 0, nfact - 1);

    /* Selects the number of dimensions for Sobol's cube (So far: all) */
    sobol_fact = nfact;

    /* Initialises Sobol to be able to fill in the first dimensions (factors) of the cube with
     * random numbers */
    err = sobol_init(0, npath - 1, 0, nbrow - 1, 0, sobol_fact - 1);
    if (err)
        return err;

    /* Fills the first dimension of the cube temp[0..np-1][0..nbr-1][0..sob-1] with Sobol */
    err = sobol_cube(temp_rand, 0, npath - 1, 0, nbrow - 1, 0, sobol_fact - 1);
    if (err)
        return err;

    /* Fills in the rest of the cube with an antithetic Gaussian noise (not worse than Sobol...) */
    err = gauss_anti_cube(temp_rand, 0, npath - 1, 0, nbrow - 1, sobol_fact, nfact - 1, seed);
    if (err)
        return err;

    /* Finishes with Sobol */
    err = sobol_free();
    if (err)
        return err;

    /* Generates the Brownian motions according to these eigen vectors */
    for (i = first_path; i <= last_path; i++)
    {
        /* Start loop on Brownian */
        for (j = first_brow; j <= last_brow; j++)
        {
            /* Start loop on eigen vectors */
            for (k = 0; k < nfact; k++)
            {
                /* The random number that will affect this rescaled factor */
                random = temp_rand[i - first_path][j - first_brow][k];

                /* Start loop on time steps (coordinate of eigen vectors) affected by this draw */
                for (l = first_step; l <= last_step; l++)
                {
                    rand[i][j][l] += random * eigen_vec[l - first_step][k];
                }

            } /* END loop on eigen vectors or factors */

            /* Transforms these absolute Brownian values into Gaussian increments */
            for (l = last_step; l > first_step; l--)
            {
                rand[i][j][l] -= rand[i][j][l - 1];
                rand[i][j][l] /= sqrt(time_at_steps[l] - time_at_steps[l - 1]);
            }
            /* Do not forget the first time step : divide by time */
            rand[i][j][first_step] /= sqrt(time_at_steps[first_step]);

        } /* END loop on Brownians */

    } /* END of loop on paths */

    /* Free all remaining allocated memory */
    free_dvector(eigen_val, 0, size - 1);
    free_dmatrix(eigen_vec, 0, size - 1, 0, size - 1);
    free_dcube(temp_rand, 0, npath - 1, 0, nbrow - 1, 0, nfact - 1);

    /* Return a success message */
    return NULL;

} /* END SpecTrunc(...) */

/* ------------------------------------------------------------------------------ */

/* Make the Brownian motion covariance matrix as inf(s,t) (allocation done outside) */

static Err build_Brownian_covar_matrix(
    double*  time_at_step, /* From [sl] to [sh] at least (sorted) */
    long     sl,
    long     sh,
    double** covar_matrix,
    long     size)
{
    long i;
    long j;

    /* First loop on all time steps from sl to sh*/
    for (i = 0; i < size; i++)
    {
        /* Second loop on all time steps with a lower time step */
        for (j = 0; j <= i; j++)
        {
            covar_matrix[i][j] = time_at_step[j + sl];
            covar_matrix[j][i] = covar_matrix[i][j];
        }
    }

    /* Return a success message */
    return NULL;

} /* END srt_f_build_Brownian_covar_matrix(...) */

/* -------------------------------------------------------------------------------- */
