/* ------------------------------------------------------------------------
        FILENAME: 	num_f_annealing.c

        AUTHORS:	O. Van Eyseren
                                A. Savine

        FUNCTION:   simulated_annealing

        PURPOSE:	provide a full function that minimises an error
                                criteria using the Simulated Annealing algorithm.
                                This function is intended as a generic wrapper to
                                the amebsa function as described in Numerical Recipes in C.
        INPUTS:
                        p:            initial simplex
                        y:            values of funk at the simplex points
                        ndim:         dimension of the space (== number of parameters )
                        ftol:         fractionnal convergence tolerance
                        niter:        maximum number of iterations
                        funcs():      the f(x) function to minimise
                        init_temp:    inital temperature in the system
                        decr_fact:    parameter to specify the speed of "cooling"
                        niter_per_t:  number of moves for each temperature
                        method:       method used for the annealation scheme
                        restart:      restart the search after a big drop (0 or 1)

        OUTPUTS:
                        p[1][...]:    best point ever
                        y[1]:         best func value ever reached
   ------------------------------------------------------------------------ */

#include "math.h"
#include "num_h_allhdr.h"
#include "num_h_annealing.h"

long sa_seed;

Err interp_annealing_method(String sMethod, AnnealMethod* pMethod)
{
    Err err = NULL;
    if ((!strcmp(sMethod, "ALPHA")) || (!strcmp(sMethod, "POWER")) ||
        (!strcmp(sMethod, "POWER_ALPHA")))
    {
        *pMethod = POWER_ALPHA;
    }
    else if (
        (!strcmp(sMethod, "EPSILON")) || (!strcmp(sMethod, "PROPORTIONAL_EPSILON")) ||
        (!strcmp(sMethod, "PROPORTIONAL")))
    {
        *pMethod = PROPORTIONAL_EPSILON;
    }
    else
    {
        return serror("Unknown Annealing method: %s ", sMethod);
    }

    return NULL;
}

/* ------------------------------------------------------------------------- */

Err simulated_annealing(
    double** p, /* From [1] to [ndim+1]  * [1] to [ndim] */
    double*  y, /* From [1] to [ndim+1] */
    long     ndim,
    double   ftol,
    long     niter,
    double (*funcs)(double*),
    double max_temp,
    double min_temp,
    double dec_fact)

{
    Err    err = NULL;
    double temperature;
    long   i;
    int    calls;

    double* pb; /* From [1] to [ndim] */
    double  yb = 1.0e+20;

    /* Initialise the pointer needed by  amebsa */
    pb = dvector(1, ndim);

    /* For communication with amebsa/amotsa */
    sa_seed = RANDINIT;

    /* Warn the user: Simulated Annealing algorithm */
    smessage(" initialising...");
    smessage(" please wait...");
    smessage("");

    /* Warn the user of the annalation scheme used */
    smessage(" starting annealing with :");
    smessage(" .initial temperature: %f", max_temp);
    smessage(" .niter: %d", niter);
    smessage(" .dec factor: %.8f", dec_fact);
    smessage("");

    /* Start iterations */
    for (i = 1; i <= niter; i++)
    {
        temperature =
            max_temp * pow(min_temp / max_temp,
                           ((double)floor(i * dec_fact * log((double)niter) / niter) - 1.00) /
                               (dec_fact * log((double)niter) - 1.00));

        smessage(" iteration: %d", (int)i);
        smessage(" temperature: %f", temperature);

        calls = (int)i;

        /* Call the main Simulated Annealing routine from Numerical Recipes */
        amebsa(p, y, ndim, pb, &yb, ftol, funcs, &calls, temperature);

        smessage(" current criteria: %.2f", 100 * sqrt(y[1]));
        smessage(" best criteria: %.2f", 100 * sqrt(yb));
        smessage("");

    } /* END for (i=1; i<=niter; i++) */

    /* Copy pb[1..ndim] to p[1][..] and yb to y[1] */
    for (i = 1; i <= ndim; i++)
    {
        p[1][i] = pb[i];
    }
    y[1] = yb;

    /* Free whatever has to be freed */
    free_dvector(pb, 1, ndim);

    /* Return a success message */
    return NULL;
}
