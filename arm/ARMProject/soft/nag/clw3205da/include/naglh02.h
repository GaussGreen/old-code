#ifndef NAGLH02
#define NAGLH02  

/* <naglh02.h>
 *
 * Copyright 1997 Numerical Algorithms Group.
 *
 * Private header file for the h02 chapter of the NAG C Library.
 *
 * Mark 5, 1997
 *
 */

/* Define a 'magic number' denoting 'unset' value of the     */
/* int_obj_bnd option member.                                */
/* Need this because it is difficult to absolutely define a  */
/* valid range for the option (can be positive or negative). */
/* However a large negative value would be unreasonable so */
/* that is what we use.                                      */ 
#define INT_OBJ_BND_NOT_SET  -1.1e+25

/* Defines of field_code values for option structure fields. */
#define FC_NOCODE 0
#define FC_BEGIN 1
#define FC_BND_NAME 2
#define FC_BRANCH_DIR 3
#define FC_COL_LO_DEFAULT 4
#define FC_COL_UP_DEFAULT 5
#define FC_CRNAMES 6
#define FC_END 7
#define FC_FEAS_TOL 8
#define FC_FIRST_SOLN 9
#define FC_HROWS 10
#define FC_INF_BOUND 11
#define FC_INT_OBJ_BOUND 12
#define FC_INT_TOL 13
#define FC_LAMBDA 14
#define FC_LIST 15
#define FC_LOWER 16
#define FC_MAX_DEPTH 17
#define FC_MAX_DF 18
#define FC_MAX_ITER 19
#define FC_MAX_NODES 20
#define FC_NCOL_APPROX 21
#define FC_NODSEL 22
#define FC_NROW_APPROX 23
#define FC_OBJ_NAME 24
#define FC_OUTFILE 25
#define FC_OUTPUT_LEVEL 26
#define FC_PRINT_FUN 27
#define FC_PRINT_LEVEL 28
#define FC_PRIORITY 29
#define FC_PROB 30
#define FC_PROB_NAME 31
#define FC_RANGE_NAME 32
#define FC_RANK_TOL 33
#define FC_RHS_NAME 34
#define FC_SOLN_TOL 35
#define FC_STATE 36
#define FC_UPPER 37
#define FC_VARSEL 38

#endif /* not NAGLH02 */






