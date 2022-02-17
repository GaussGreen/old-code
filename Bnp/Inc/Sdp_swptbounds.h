#ifndef SWPBOUNDS
#define SWPBOUNDS



//#include "stdio.h"
//#include "math.h"
#include "utallhdr.h"

#ifdef __cplusplus
extern "C" {
#endif

Err smooth_caplets(double ***const_mat, double *const_val, int dim, int num_const, int print_level);

Err build_instrument_with_cvgs(double *vecs, double **matres, int dim, double *coverages);

/*  The main SDP function....... */
Err compute_bounds(double **calib_inst, double *market_vars, int ninst, double *vec_swpt, int dim, int minmax, double *result, int *return_error, double *error, int printlevel, double toler, int niter);
/*  ------------------------------ */

#ifdef __cplusplus
}
#endif


#endif

