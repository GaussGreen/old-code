// Author: AdA
// Date:   22/06/1999

#ifndef _SDP_H
#define _SDP_H

// Associated code: sdp.c
// Uses: some nr routines and some BLAS routines.

// BLAS from the Intel Math Kernel library can be linked directly,
// or the MiniBLAS provided for use in griffin can be used instead.
// Use sdplib.h for lib choices.

// This defines the main routine for Semi Deinite Programming
// It takes a set of matrix and doubles as constraints and
// a matrix for the optimisation direction.

// ------------------------------------------------------------------
// If in debug, include the Monitor header for convergence tracing.
// ------------------------------------------------------------------

#include "uterror.h"

Err srt_sdp(
    double*** const_mat,
    double*   const_val,
    double**  cmat,
    int       dim,
    int       num_const,
    double**  xmat,
    double*   obj_val,
    double*   error,
    int*      return_error,
    double    toler,
    int       printlevel,
    int       max_iter);

Err srt_bgm_sdp(
    double*** const_mat,
    double*   const_val,
    double**  cmat,
    int       dim,
    int       num_const,
    double**  xmat,
    double*   obj_val,
    double*   error,
    int*      return_error,
    double    toler,
    int       printlevel,
    int       max_iterations);

#endif
