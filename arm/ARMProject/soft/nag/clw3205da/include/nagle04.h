#ifndef NAGLE04
#define NAGLE04  

/* <nagle04.h>
 *
 * Copyright 1991 Numerical Algorithms Group.
 *
 * Private header file for the e04 chapter of the NAG C Library.
 *
 * Mark 2, 1991
 *
 * Mark 4 revised. IER-1757 (Aug 1995).
 * Mark 5 revised. IER-2164 (Feb 1998).
 * Various new field codes for new e04 functions.
 */

/* Denotes a Boolean in options is not set */
#define BOOL_NOT_SET -111

/* Defines of field_code values for option structure fields. */
#define FC_NOCODE 0
#define FC_END 1
#define FC_LIST 2
#define FC_F_EST 3
#define FC_MAX_NF 4
#define FC_F_PREC 5
#define FC_MAX_ITER 6
#define FC_DEFAULTS 7  /* Not used */
#define FC_OPTIM_TOL 8
#define FC_PRINT_LEVEL 9
#define FC_VERIFY_GRAD 10
#define FC_MAX_LINE_STEP 11
#define FC_LINESEARCH_TOL 12
#define FC_OBJ_CHECK_STOP 13
#define FC_CON_CHECK_STOP 14
#define FC_OBJ_CHECK_START 15
#define FC_CON_CHECK_START 16
#define FC_HESSIAN 17
#define FC_INVALID 18
#define FC_DEBUG 19
#define FC_DEBUG_ITER 20
#define FC_OUTFILE 21
#define FC_PRINT_GCHECK 22
#define FC_STEP_MAX 23
#define FC_MINLIN 24
#define FC_DERIV_CHECK 25
#define FC_INIT_STATE 26
#define FC_LOCAL_SEARCH 27
#define FC_S 28
#define FC_V 29
#define FC_TDV 30
#define FC_HESD 31
#define FC_HESL 32
#define FC_DELTA 33
#define FC_STATE 34
#define FC_PRINT_FUN 35
#define FC_BEGIN 36
#define FC_AX 37
#define FC_PROB 38
#define FC_LAMBDA 39
#define FC_START 40
#define FC_FTOL 41
#define FC_FCHECK 42
#define FC_INF_STEP 43
#define FC_INF_BOUND 44
#define FC_RESET_FTOL 45
#define FC_CRASH_TOL 46
#define FC_FMAX_ITER 47
#define FC_HROWS 48
#define FC_RANK_TOL 49
#define FC_MAX_DF 50
#define FC_MINOR_PRINT_LEVEL 51
#define FC_OBJ_DERIV 52
#define FC_CON_DERIV 53
#define FC_PRINT_DERIV 54
#define FC_F_DIFF_INT 55
#define FC_C_DIFF_INT 56
#define FC_MINOR_MAX_ITER 57
#define FC_LIN_FEAS_TOL 58
#define FC_NONLIN_FEAS_TOL 59
#define FC_STEP_LIMIT 60
#define FC_CONF 61
#define FC_CONJAC 62
#define FC_R 63
#define FC_H 64
#define FC_NPROB 65
#define FC_N 66
#define FC_NCLIN 67
#define FC_NCNLIN 68
#define FC_MIN_INFEAS 69
#define FC_H_RESET_FREQ 70
#define FC_H_UNIT_INIT 71
#define FC_BND_NAME 72
#define FC_COL_LO_DEFAULT 73
#define FC_COL_UP_DEFAULT 74
#define FC_CRASH 75
#define FC_CRNAMES 76
#define FC_DERIV_WANT 77
#define FC_EST_DENSITY 78
#define FC_FACTOR_FREQ 79
#define FC_F_PREC_USED 80
#define FC_LU_FACTOR_TOL 81
#define FC_LU_SING_TOL 82
#define FC_LU_UPDATE_TOL 83
#define FC_MAX_SB 84
#define FC_MINIMIZE 85
#define FC_NCOL_APPROX 86
#define FC_NROW_APPROX 87
#define FC_NSB 88
#define FC_OBJ_NAME 89
#define FC_OUTPUT_LEVEL 90
#define FC_PARTIAL_PRICE 91
#define FC_PIVOT_TOL 92
#define FC_PROB_NAME 93
#define FC_RANGE_NAME 94
#define FC_RHS_NAME 95
#define FC_SCALE 96
#define FC_SCALE_TOL 97
#define FC_USE_HFWD_INIT 98

#endif /* not NAGLE04 */
