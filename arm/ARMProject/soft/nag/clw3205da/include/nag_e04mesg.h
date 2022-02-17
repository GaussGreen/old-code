/* <nag_e04mesg.h>
 *
 * Copyright 1990 Numerical Algorithms Group
 *
 * NAG C Library
 *
 * Mark 1, 1990
 *
 * Mark 5 revised. IER-2135 (Feb 1998).
 */
 
#ifndef NAG_E04MESG
#define NAG_E04MESG  

/* Output messages for Chapter e04 */

#define NM_NO_MESG 0
#define NM_END_LINE 1
#define NM_OPT_RES 2
#define NM_ITER_RES 3
#define NM_ITER_HEADINGS 4
#define NM_ITER_HEADINGS_LSQ_1 5
#define NM_ITER_HEADINGS_LSQ_2 6
#define NM_ITER_HEADINGS_QN_1 7
#define NM_ITER_HEADINGS_QN_2 8
#define NM_INIT_ITER_VALS 9
#define NM_INIT_ITER_VALS_LSQ_1 10
#define NM_INIT_ITER_VALS_QN_1 11
#define NM_ITER_VALS 12
#define NM_ITER_VALS_LSQ_1 13
#define NM_ITER_VALS_LSQ_2 14
#define NM_ITER_VALS_QN_1 15
#define NM_ITER_VALS_QN_2 16
#define NM_FREE 17
#define NM_UPPER 18
#define NM_LOWER 19
#define NM_CONST 20
#define NM_EXIT_ITER_NUM 21
#define NM_SOLN_RES 22
#define NM_SOLN_HEADINGS 23
#define NM_SOLN_HEADINGS_LSQ 24
#define NM_SOLN_VALS 25
#define NM_SSQ_SOLN_VALS 26
#define NM_1_D_VAL 27
#define NM_OBJ_VALUE 28
#define NM_SSQ_VALUE 29
#define NM_SOLN_FOUND 30
#define NM_WARNING 31
#define NM_LOCAL_SEARCH_DONE 32
#define NM_ITER_HEADINGS_SIMPLEX 33
#define NM_ITER_HEADINGS_SIMPLEX_1 34
#define NM_ITER_HEADINGS_SIMPLEX_2 35
#define NM_ITER_VALS_SIMPLEX 36
#define NM_ITER_VALS_SIMPLEX_1 37
#define NM_ITER_VALS_SIMPLEX_2 38
#define NM_ITER_VALS_SIMPLEX_3 39
#define NM_ITER_VALS_SIMPLEX_4 40
#define NM_SOLN_HEADINGS_SIMPLEX 41
#define NM_FUN_VAL_SIMPLEX 42
#define NM_RANK_OBFUN_MATRIX 43
#define NM_HESSIAN_INDEF 44
#define NM_MULT_TITLE 45
#define NM_CONSTRAINT_1 46
#define NM_CONSTRAINT_2 47
#define NM_CONSTRAINT_3 48
#define NM_CONSTRAINT_4 49
#define NM_CONSTRAINT_5 50
#define NM_CONSTRAINT_6 51
#define NM_CONSTRAINT_7 52
#define NM_CONSTRAINT_8 53
#define NM_CONSTRAINT_9 54
#define NM_CONSTRAINT_10 55
#define NM_CONSTRAINT_LINE_1 56
#define NM_CONSTRAINT_LINE_2 57
#define NM_CONSTRAINT_LINE_3 58
#define NM_CONSTRAINT_LINE_4 59
#define NM_CONSTRAINT_LINE_5 60
#define NM_CONSTRAINT_LINE_6 61
#define NM_CONSTRAINT_LINE_7 62
#define NM_CONSTRAINT_LINE_8 63
#define NM_CONSTRAINT_LINE_9 64
#define NM_CONSTRAINT_LINE_10 65
#define NM_MULTIPLIER_1 66
#define NM_MULTIPLIER_2 67
#define NM_MULTIPLIER_3 68
#define NM_MULTIPLIER_4 69
#define NM_MULTIPLIER_5 70
#define NM_MULTIPLIER_6 71
#define NM_MULTIPLIER_7 72
#define NM_MULTIPLIER_8 73
#define NM_MULTIPLIER_9 74
#define NM_MULTIPLIER_10 75
#define NM_MULT_VALUES_1 76
#define NM_MULT_VALUES_2 77
#define NM_ITER_HEAD_SHORT 78
#define NM_ITER_HEAD_FP_SHORT 79
#define NM_ITER_HEAD_QP_LONG 80
#define NM_ITER_HEAD_LP_LONG 81
#define NM_ITER_HEAD_FP_LONG 82
#define NM_ITER_VALUES_SHORT 83
#define NM_ITER_VALUES_QP_LONG 84
#define NM_ITER_VALUES_LP_LONG 85
#define NM_X_TITLE_1 86
#define NM_X_TITLE_2 87
#define NM_X_TITLE_3 88
#define NM_X_TITLE_4 89
#define NM_X_TITLE_5 90
#define NM_X_TITLE_6 91
#define NM_X_VALUES_1 92
#define NM_X_VALUES_2 93
#define NM_X_VALUES_3 94
#define NM_X_VALUES_4 95
#define NM_X_VALUES_6 96
#define NM_ITERATION_PROB 97
#define NM_STATUS 98
#define NM_NOT_SATISFY 99
#define NM_ITER_FEASIBLE 100
#define NM_MOVED 101
#define NM_REFACTOR 102
#define NM_QP_VAR_HEAD 103
#define NM_QP_LCON_HEAD 104
#define NM_QP_NLCON_HEAD 105
#define NM_FINAL_STATUS 106
#define NM_QP_V 107
#define NM_QP_L 108
#define NM_QP_N 109
#define NM_NO_BOUND 110
#define NM_UP_BOUND 111
#define NM_LOW_BOUND 112
#define NM_BOTH_BOUND 113
#define NM_EXIT_ITER 114
#define NM_FEASIBLE 115
#define NM_OPT_SOLN 116
#define NM_FINAL_OBJ 117
#define NM_MIN_SUMINF 118
#define NM_FINAL_SUMINF 119
#define NM_OBJ_GRAD 120
#define NM_ALL_G_OBJ 121
#define NM_SOME_G_OBJ 122
#define NM_SOME_CONST_G_OBJ 123
#define NM_ALL_CONST_G_OBJ 124
#define NM_SIMPLE_CHECK 125
#define NM_NO_SIMPLE_CHECK 126
#define NM_WRONG_GRAD 127
#define NM_OK_GRAD 128
#define NM_JAC_OK 129
#define NM_JAC_NOK 130
#define NM_SIMPLE_G_RES 131
#define NM_SUBFUN_REL_ERR 132
#define NM_COMP_CHECK 133
#define NM_SUBFUN_NO 134
#define NM_JAC_ELEM_HEAD 135
#define NM_OBJ_COMP_HEADINGS 136
#define NM_CON_COMP_HEADINGS 137
#define NM_COMP_VALS 138
#define NM_ASSUME_ZERO 139
#define NM_ZERO_ELEMENT 140
#define NM_ZERO_ELEMENTS 141
#define NM_ZERO_ELEMENT_IS 142
#define NM_ZERO_ELEMENTS_ARE 143
#define NM_ALL_OBJ_GOOD 144
#define NM_SOME_OBJ_WRONG 145
#define NM_ALL_JAC_GOOD 146
#define NM_SOME_JAC_WRONG 147
#define NM_SMALL_DERIV 148
#define NM_SIZE_ERROR 149
#define NM_SIZE_JAC_ERROR 150
#define NM_CHANGE_DER_LEVEL 151
#define NM_CON_GRAD 152
#define NM_ALL_G_CON 153
#define NM_SOME_G_CON 154
#define NM_CON_REL_ERR 155
#define NM_CON_NO 156
#define NM_SOME_CONST_G_CON 157
#define NM_ALL_CONST_G_CON 158
#define NM_FIN_DIFF_HEAD 159
#define NM_FIN_DIFF_INT 160
#define NM_PARAM_TITLE 161
#define NM_NUM_VAR 162
#define NM_NUM_VAR_NONEWL 163
#define NM_NUM_RESIDUALS 164
#define NM_NUM_MAT_ROWS 165
#define NM_MAX_LINE_STEP 166
#define NM_MINLIN 167
#define NM_MACH_PREC 168
#define NM_MACH_PREC_L 169
#define NM_OPTIM_TOL 170
#define NM_LINESEARCH_TOL 171
#define NM_F_EST 172
#define NM_F_PREC 173
#define NM_F_PREC_L 174
#define NM_VERIFY_GRAD 175
#define NM_MAX_ITER 176
#define NM_MAX_ITER_1 177
#define NM_PRINT_LEVEL 178
#define NM_PRINT_GCHECK 179
#define NM_OUTFILE 180
#define NM_DEBUG 181
#define NM_STEP_MAX 182
#define NM_PRINT_ITER 183
#define NM_DERIV_CHECK 184
#define NM_DERIV_CHECK_2 185
#define NM_INIT_STATE 186
#define NM_LOCAL_SEARCH 187
#define NM_MEM_ALLOC 188
#define NM_S 189
#define NM_V 190
#define NM_DELTA 191
#define NM_STATE 192
#define NM_HESL 193
#define NM_HESD 194
#define NM_NUM_LIN_CON 195
#define NM_PROB 196
#define NM_START 197
#define NM_NUM_NONLIN_CON 198
#define NM_INF_BOUND 199
#define NM_INF_STEP 200
#define NM_RANK_TOL 201
#define NM_RANK_TOL_R 202
#define NM_CRASH_TOL 203
#define NM_CRASH_TOL_L 204
#define NM_FTOL 205
#define NM_RESET_FTOL 206
#define NM_FCHECK 207
#define NM_MAX_DF 208
#define NM_HROWS 209
#define NM_MAX_FREE 210
#define NM_MAX_ACTIVE 211
#define NM_FMAX_ITER 212
#define NM_AX 213
#define NM_LAMBDA 214
#define NM_STEP_LIMIT 215
#define NM_LIN_FEAS_TOL 216
#define NM_NONLIN_FEAS_TOL 217
#define NM_MINOR_MAX_ITER 218
#define NM_HESSIAN 219
#define NM_HESSIAN_R 220
#define NM_F_DIFF_INT 221
#define NM_F_DIFF_INT_AUTO 222
#define NM_C_DIFF_INT 223
#define NM_C_DIFF_INT_AUTO 224
#define NM_OBJ_DERIV 225
#define NM_CON_DERIV 226
#define NM_PRINT_DERIV 227
#define NM_MINOR_PRINT_LEVEL 228
#define NM_MIN_INFEAS 229
#define NM_H_RESET_FREQ 230
#define NM_H_UNIT_INIT 231
#define NM_LAMBDA_L 232
#define NM_PROB_NAME 233
#define NM_OBJ_NAME 234
#define NM_RHS_NAME 235
#define NM_BND_NAME 236
#define NM_RANGE_NAME 237
#define NM_NAMES 238
#define NM_NEW_LINE 239
#define NM_BIGNUM_LIN_CON 240
#define NM_SPARSE_PROB_TYPE 241
#define NM_HESSIAN_COLUMNS 242
#define NM_MIN_OR_MAX 243
#define NM_SCALE 244
#define NM_CRASH 245
#define NM_PARTIAL_PRICE 246
#define NM_SUPERBASICS_LIMIT 247
#define NM_LU_FACTOR_TOL 248
#define NM_LU_SING_TOL 249
#define NM_STATE_L 250
#define NM_BIGNUM_VAR 251
#define NM_LU_UPDATE_TOL 252
#define NM_FACTOR_FREQ 253
#define NM_PIVOT_TOL 254
#define NM_SCALE_TOL 255
#define NM_DERIV_WANT 256
#define NM_USE_HFWD_INIT 257
#define NM_PRINT_DERIV_L 258
#define NM_COL_LO_DEFAULT 259
#define NM_COL_UP_DEFAULT 260
#define NM_INFINITY 261
#define NM_NCOL_APPROX 262
#define NM_NROW_APPROX 263
#define NM_MPS_SPARSE 264
#define NM_MPS_EST_DENSITY 265
#define NM_MPS_OUTLEVEL 266
#define NM_SET_TEXT 267
#define NM_SET_DOUBLE 268
#define NM_SET_INT 269
#define NM_RESET_DEF 270
#define NM_RESET 271
#define NM_TEXT_VALID_RANGE 272
#define NM_SET_DOUBLE_MIN 273
#define NM_INTRO 274
#define NM_QPNL_EX_MAJ_MIN 275
#define NM_QPNL_EX_MAJOR 276
#define NM_QPNL_HEAD_V 277
#define NM_QPNL_HEAD_L 278
#define NM_QPNL_VAR_CON_RES 279
#define NM_QPNL_HEAD_N 280
#define NM_QPNL_ncj_iter 281
#define NM_QPNL_ncj_head1 282
#define NM_QPNL_ncj_head2 283
#define NM_QPNL_ncj_iter_line_1 284
#define NM_QPNL_ncj_iter_line_2 285
#define NM_QPNL_ncj_head3 286
#define NM_QPNL_ncj_iter_line_3 287
#define NM_QPNL_ncj_head4 288
#define NM_QPNL_ncj_state_1 289
#define NM_QPNL_ncj_head5 290
#define NM_QPNL_ncj_state_2 291
#define NM_QPNL_ncj_work_set_fac1 292
#define NM_QPNL_ncj_work 293
#define NM_QPNL_ncj_diag_triang_r 294
#define NM_QPNL_ncj_val_diag_elem_r 295
#define NM_QPNL_ncj_dash_1 296
#define NM_QPNL_ncj_head6 297
#define NM_QPNL_ncj_gq_val 298
#define NM_QPNL_ncj_head7 299
#define NM_QPNL_ncj_kx_free_k 300
#define NM_QPNL_ncj_rlamda_activ 301
#define NM_QPNL_ncj_head8 302
#define NM_QPNL_ncj_kactiv 303
#define NM_QPNL_ncj_rlamda_k 304
#define NM_QPNL_ucc_lin_con_feas 305
#define NM_QPNL_ucc_bad_init_hess_r_refact 306
#define NM_QPNL_ucc_r_line_end 307
#define NM_QPNL_ucc_r_elem 308
#define NM_QPNL_ucc_w_line_end 309
#define NM_QPNL_ucc_w_lwrk1_j 310
#define NM_QPNL_ucc_opt_found 311
#define NM_QPNL_ucc_final_obj_val 312
#define NM_QPNL_ucc_min_sum_infeas 313
#define NM_QPNL_ucc_final_sum_infeas 314
#define NM_QPNL_ucc_maj_iter 315
#define NM_QPNL_uct_maj_iter_head 316
#define NM_QPNL_uct_prthdr_nlncon 317
#define NM_QPNL_uct_prthdr_Nlncon 318
#define NM_QPNL_uct_nlncon_line1 319
#define NM_QPNL_uct_Nnlncon_line1 320
#define NM_QPNL_uct_prthdr_nlncon_def 321
#define NM_QPNL_uct_prthdr_Nnlncon_def 322
#define NM_QPNL_uct_nlncon_def 323
#define NM_QPNL_uct_Nnlncon_def 324
#define NM_QPNL_uct_ncnln_zero_obj 325
#define NM_QPNL_uct_ncnln_Nzero_obj 326
#define NM_QPNL_uct_ncnln_Nzero_cons 327
#define NM_QPNL_uct_cons_head 328
#define NM_QPNL_uct_cons_istate 329
#define NM_QPNL_uct_gen_cons_head 330
#define NM_QPNL_uct_gen_cons_istate 331
#define NM_QPNL_uct_nlin_head 332
#define NM_QPNL_uct_nlin_state 333
#define NM_QPNL_uct_diag_t_head 334
#define NM_QPNL_uct_diag_t_val 335
#define NM_QPNL_uct_diag_r_head 336
#define NM_QPNL_uct_diag_r_val 337
#define NM_QPNL_uct_dashed 338
#define NM_QPNL_ucu_maj_itn_head 339
#define NM_QPNL_ucz_minor_itn_head 340
#define NM_REALLOC_BASIS_WKSP 341
#define NM_PARP_REDUCED 342
#define NM_PART_PRICE_INFO 343
#define NM_COLD_START 344
#define NM_WARM_START 345
#define NM_ITERATIONS 346
#define NM_NAME_ITER_OBJ 347
#define NM_INFEAS_SUMM 348
#define NM_FEAS_SPQP_SUMM 349
#define NM_SBASIC_SUMM 350
#define NM_DEGEN_SUMM 351
#define NM_X_PI_NORM_SCALED 352
#define NM_X_PI_NORM_UNSCALED 353
#define NM_MAX_INFEAS_SCALED 354
#define NM_MAX_INFEAS_UNSCALED 355
#define NM_CRASH_OPTION 356
#define NM_CRASH_LIN_EQ 357
#define NM_CRASH_LIN_INEQ 358
#define NM_CRASH_NONLIN 359
#define NM_CRASH_DETAILS 360
#define NM_SCALE_HEAD 361
#define NM_SCALE_PASS 362
#define NM_SCALE_PHASE_SUMM 363
#define NM_SCALE_NORMS 364
#define NM_SCALE_REDUCED 365
#define NM_ROW_SCALES_HEAD 366
#define NM_COL_SCALES_HEAD 367
#define NM_SCALES 368
#define NM_FACT_ITN 369
#define NM_BS_FACT_NSWAP 370
#define NM_BFACT_STATS 371
#define NM_ROW_CHECK 372
#define NM_BFACT_DETAILS 373
#define NM_BSFACT_DETAILS 374
#define NM_BFACT_LUFAC_RED 375
#define NM_BFACT_LUUPD_RED 376
#define NM_COL_REPL_SLACK 377
#define NM_COL_REPL_SLACK_ETC 378
#define NM_MAT_STATS_HEAD 379
#define NM_MAT_STATS 380
#define NM_MAT_DETAILS 381
#define NM_OBCOEF_MINMAX 382
#define NM_NLN_FAILED 383
#define NM_NO_MAX_PIV 384
#define NM_BIGDJ 385
#define NM_NORMS_RG_PI 386
#define NM_LU_FILE_GROW 387
#define NM_HESS_INDEF_ITN 388
#define NM_REDHESS_INDEF_ITN 389
#define NM_HESS_DIAGS 390
#define NM_ITN_INFEAS 391
#define NM_ITN_INFEAS_LONG 392
#define NM_WRONG_NO_BASICS 393
#define NM_ITN_FEAS_EQ 394
#define NM_SMALL_DIR_DERIV 395
#define NM_KSAVE_NOT_FOUND 396
#define NM_MAXPIV_TOO_SMALL 397
#define NM_BASICS_RECOMP 398
#define NM_HESS_SING 399
#define NM_VAR_INDEF 400
#define NM_SPQP_ITN_SHD_QP 401
#define NM_SPQP_ITN_SHD_FPLP 402
#define NM_ITN_SPQP_QP_SHORT 403
#define NM_ITN_SPQP_FPLP_SHORT 404
#define NM_SPQP_ITN_LHD_QP 405
#define NM_ITN_SPQP_QP_LONG 406
#define NM_SPQP_ITN_LHD_FPLP 407
#define NM_ITN_SPQP_FPLP_LONG 408
#define NM_ITN_FEAS_STRING 409
#define NM_PROB_DESCR 410
#define NM_SPQP_SOL_HEAD 411
#define NM_SPQP_SOL_LINE 412
#define NM_SPQP_COL_SOL_HEAD_LONG 413
#define NM_SPQP_ROW_SOL_HEAD_LONG 414
#define NM_SPQP_SOL_LINE_LONG 415
#define NM_ITER_HEADINGS_MODN_1 416
#define NM_INIT_ITER_VALS_MODN_1 417
#define NM_ITER_VALS_MODN_1 418
#define NM_F_PREC_TOO_LARGE 419
#define NM_F_PREC_TOO_SMALL 420
#define NM_HESS_DERIV_HEAD 421
#define NM_GRAD_HESS_DERIV_HEAD 422
#define NM_HESS_DERIV_RES 423
#define NM_GRAD_HESS_DERIV_RES 424

/* END OF DEFINES */


#ifdef NAG_MESG
char *nag_e04mesg[] =
{
": Dummy message for Chapter e04\n",

/* Result output */
"NM_END_LINE: \n",
"NM_OPT_RES: \nResults from %s:\n-------------------\n",
"NM_ITER_RES: \nIteration results:\n",
"NM_ITER_HEADINGS: \n  Itn  Nfun   Objective   Norm g   Norm x   \
Norm (x(k-1)-x(k))   Step\n",
"NM_ITER_HEADINGS_LSQ_1: \n  Itn  Nfun   Objective   Norm g   Norm x   \
Norm (x(k-1)-x(k))   Step    Grade\n",
"NM_ITER_HEADINGS_LSQ_2: \n       x                  g            Singular values\n",
"NM_ITER_HEADINGS_QN_1: \n  Itn  Nfun  Objective    Norm g    Norm x  \
Norm(x(k-1)-x(k))  Step    Cond H\n",
"NM_ITER_HEADINGS_QN_2: \n Variable         x                g          Status\n",
"NM_INIT_ITER_VALS: %4ld %5ld  %12.4e %8.1e %8.1e\n",
"NM_INIT_ITER_VALS_LSQ_1: %4ld %5ld  %12.4e %8.1e %8.1e                               %3ld\n",
"NM_INIT_ITER_VALS_QN_1: %4ld %5ld %12.4e %9.1e %9.1e                          %9.1e\n",
"NM_ITER_VALS: %4ld %5ld  %12.4e %8.1e %8.1e     %9.1e      %9.1e\n",
"NM_ITER_VALS_LSQ_1: %4ld %5ld  %12.4e %8.1e %8.1e     %9.1e      %9.1e  %3ld\n",
"NM_ITER_VALS_LSQ_2: %13.5e       %12.4e       %12.4e\n",
"NM_ITER_VALS_QN_1: %4ld %5ld %12.4e %9.1e %9.1e   %9.1e    %9.1e %9.1e\n",
"NM_ITER_VALS_QN_2: %5ld       %12.4e     %12.4e",
"NM_FREE:     Free\n",
"NM_UPPER:     Upper Bound\n",
"NM_LOWER:     Lower Bound\n",
"NM_CONST:     Constant\n",
"NM_EXIT_ITER_NUM: \nExit after %1ld iterations and %1ld \
function evaluations.\n",
"NM_SOLN_RES: \nFinal solution:\n",
"NM_SOLN_HEADINGS: \n Variable         x                g\n",
"NM_SOLN_HEADINGS_LSQ: \n       x               g            Residuals\n",
"NM_SOLN_VALS: %5ld       %12.4e     %12.4e\n",
"NM_SSQ_SOLN_VALS: %13.5e    %12.4e     %12.4e\n",
"NM_1_D_VAL:                                   %12.4e\n",
"NM_OBJ_VALUE: \nFinal objective function value = %12.7e.\n",
"NM_SSQ_VALUE: \nThe sum of squares is %11.4e.\n",
"NM_SOLN_FOUND: \nOptimal solution found.\n",
"NM_WARNING: \nN.B. Warning state occured during solution.\n",
"NM_LOCAL_SEARCH_DONE: \nLocal search performed.\n",

/* NAG_SIMPLEX messages */
"NM_ITER_HEADINGS_SIMPLEX: \n   Itn   Nfun    Fmin         Fmax\n",
"NM_ITER_HEADINGS_SIMPLEX_1: \n   Current estimate of solution x\n",
"NM_ITER_HEADINGS_SIMPLEX_2: \n   Vertex            Position vectors for the current simplex\n",
"NM_ITER_VALS_SIMPLEX: %5ld %5ld %12.4e %12.4e\n",
"NM_ITER_VALS_SIMPLEX_1:      %12.4e\n",
"NM_ITER_VALS_SIMPLEX_2:    %5ld       ",
"NM_ITER_VALS_SIMPLEX_3:   %12.4e%s",
"NM_ITER_VALS_SIMPLEX_4: \n",
"NM_SOLN_HEADINGS_SIMPLEX: \nVector x\n",
"NM_FUN_VAL_SIMPLEX: \nFinal Function value is %12.4e\n",

/* Convex QP/LS (e04ncc) messages */
"NM_RANK_OBFUN_MATRIX: \nRank of the objective function data matrix = %1ld\n",
"NM_HESSIAN_INDEF: XXX Hessian appears to be indefinite.\nXXX \
Maximum diagonal and off-diagonal ignored in the Cholesky factorization:\
\n %22.14e %22.12e",

/* NAG_QP_START messages */
"NM_MULT_TITLE: \n Lagrange multipliers for the constraints\n",
"NM_CONSTRAINT_1: \n Bound Constraints    Linear Constraints       \
Artificial Constraints\n",
"NM_CONSTRAINT_2: \n Bound Constraints    Linear Constraints\n",
"NM_CONSTRAINT_3: \n Bound Constraints                             \
Artificial Constraints\n",
"NM_CONSTRAINT_4: \n                      Linear Constraints       \
Artificial Constraints\n",
"NM_CONSTRAINT_5: \n Bound Constraints\n",
"NM_CONSTRAINT_6: \n                                               \
Artificial Constraints\n",
"NM_CONSTRAINT_7: \n                      Linear Constraints\n",
"NM_CONSTRAINT_8: \n Bound Constraints                      \
Artificial Constraints\n",
"NM_CONSTRAINT_9: \n Bound Constraints\n",
"NM_CONSTRAINT_10: \n                                        \
Artificial Constraints\n",
"NM_CONSTRAINT_LINE_1:  -----------------    ------------------       \
----------------------\n",
"NM_CONSTRAINT_LINE_2:  -----------------    ------------------\n",
"NM_CONSTRAINT_LINE_3:  -----------------                             \
----------------------\n",
"NM_CONSTRAINT_LINE_4:                       ------------------       \
----------------------\n",
"NM_CONSTRAINT_LINE_5:  -----------------\n",
"NM_CONSTRAINT_LINE_6:                                                \
----------------------\n",
"NM_CONSTRAINT_LINE_7:                       ------------------\n",
"NM_CONSTRAINT_LINE_8:  -----------------                      \
----------------------\n",
"NM_CONSTRAINT_LINE_9:  -----------------\n",
"NM_CONSTRAINT_LINE_10:                                         \
----------------------\n",
"NM_MULTIPLIER_1:  Index  Multiplier    Index  Multiplier             \
Multiplier\n",
"NM_MULTIPLIER_2:  Index  Multiplier    Index  Multiplier\n",
"NM_MULTIPLIER_3:  Index  Multiplier                                  \
Multiplier\n",
"NM_MULTIPLIER_4:                       Index  Multiplier             \
Multiplier\n",
"NM_MULTIPLIER_5:  Index  Multiplier\n",
"NM_MULTIPLIER_6:                                                     \
Multiplier\n",
"NM_MULTIPLIER_7:                       Index  Multiplier\n",
"NM_MULTIPLIER_8:  Index  Multiplier                           Multiplier\n",
"NM_MULTIPLIER_9:  Index  Multiplier\n",
"NM_MULTIPLIER_10:                                              Multiplier\n",

"NM_MULT_VALUES_1: %s   %s            %s\n",
"NM_MULT_VALUES_2: %s                %s\n",
"NM_ITER_HEAD_SHORT: \n  Itn Jdel  Jadd   Step    Ninf  Sinf/Obj    Bnd\
  Lin  Nart  Nrz  Norm Gz\n\n",
"NM_ITER_HEAD_FP_SHORT: \n  Itn Jdel  Jadd   Step    Ninf   Suminf     Bnd\
  Lin  Nart  Nrz  Norm Gz\n\n",
"NM_ITER_HEAD_QP_LONG: \n  Itn Jdel  Jadd   Step    Ninf  Sinf/Obj    Bnd\
  Lin  Nart  Nrz  Norm Gz   NOpt   Min Lm    Cond T  Cond Rz     Rzz\n\n",
"NM_ITER_HEAD_LP_LONG: \n  Itn Jdel  Jadd   Step    Ninf  Sinf/Obj    Bnd\
  Lin  Nart  Nrz  Norm Gz   NOpt   Min Lm    Cond T\n\n",
"NM_ITER_HEAD_FP_LONG: \n  Itn Jdel  Jadd   Step    Ninf   Suminf     Bnd\
  Lin  Nart  Nrz  Norm Gz   NOpt   Min Lm    Cond T\n\n",
"NM_ITER_VALUES_SHORT: %4ld %3ld %.1s %3ld %.1s %8.1e %3ld %12.4e \
%4ld %4ld %4ld %5ld %10.2e\n",
"NM_ITER_VALUES_QP_LONG: %4ld %3ld %.1s %3ld %.1s %8.1e %3ld %12.4e \
%4ld %4ld %4ld %5ld %10.2e %15s %8.1e %.8s %.9s\n",
"NM_ITER_VALUES_LP_LONG: %4ld %3ld %.1s %3ld %.1s %8.1e %3ld %12.4e \
%4ld %4ld %4ld %5ld %10.2e %15s %8.1e\n",
"NM_X_TITLE_1: \n Index  Value of x   State     Value of Ax  State\
   Diagonal T     Diagonal Rz\n\n",
"NM_X_TITLE_2: \n Index  Value of x   State     Value of Ax  State\
   Diagonal T\n\n",
"NM_X_TITLE_3: \n Index  Value of x   State     Value of Ax  State\n\n",
"NM_X_TITLE_4: \n Index  Value of x   State\n\n",
"NM_X_TITLE_5: \n Index  Value of x   State     Value of Ax  State\
                Diagonal Rz\n\n",
"NM_X_TITLE_6: \n Index  Value of x   State    Diagonal Rz\n\n",
"NM_X_VALUES_1:  %3ld %s    %s  %s %s\n",
"NM_X_VALUES_2:  %3ld %s    %s  %s\n",
"NM_X_VALUES_3:  %3ld %s    %s\n",
"NM_X_VALUES_4:  %3ld %s\n",
"NM_X_VALUES_6:  %3ld %s   %s\n",
"NM_ITERATION_PROB: \n\n%s iteration %1ld\n--------------\n",
/* this string is a set of letters for designation of constraint status */
"NM_STATUS:   L U E F A ",
"NM_NOT_SATISFY: \n XXX Warning. Cannot satisfy the constraints to the \
accuracy requested.\n",
"NM_ITER_FEASIBLE: \n Itn %1ld  -- Feasible point found.\n",
"NM_MOVED: \n Itn %1ld -- %1ld variables moved to their bounds\n",
"NM_REFACTOR: XXX  Iterative refinement.\n\
     The maximum violation is %14.2e in constraint %1ld\n",
"NM_QP_VAR_HEAD: \n  Varbl State    Value      Lower Bound  Upper Bound    \
Lagr Mult    Residual\n\n",
"NM_QP_LCON_HEAD: \n   LCon  State    Value     Lower Bound  Upper Bound    \
Lagr Mult    Residual\n\n",
"NM_QP_NLCON_HEAD: \n   NLCon State    Value     Lower Bound  Upper Bound    \
Lagr Mult    Residual\n\n",
"NM_FINAL_STATUS: --++FRLLULEQTF",
"NM_QP_V:  V ", /* these three letters, V, L, N, couple with 1 of 4 lines below */
"NM_QP_L:  L ",
"NM_QP_N:  N ",
"NM_NO_BOUND: %3ld    %.2s %13.5e     None          None     %11.3e %11.3e\n",
"NM_UP_BOUND: %3ld    %.2s %13.5e     None     %12.4e  %11.3e %11.3e\n",
"NM_LOW_BOUND: %3ld    %.2s %13.5e %12.4e   None     %11.3e %11.3e\n",
"NM_BOTH_BOUND: %3ld    %.2s %13.5e %12.4e %12.4e  %11.3e %11.3e\n",
"NM_EXIT_ITER: \nExit after %1ld iterations.\n",
"NM_FEASIBLE: \nFeasible point found.\n",
"NM_OPT_SOLN: \nOptimal %s solution found.\n",
"NM_FINAL_OBJ: \nFinal %s objective value = %15.7e\n\n",
"NM_MIN_SUMINF: \nMinimum sum of infeasibilities = %15.7e\n",
"NM_FINAL_SUMINF: \nFinal sum of infeasibilities = %15.7e\n",

/* Gradient check output */
"NM_OBJ_GRAD: \n\nVerification of the objective gradients.\n\
----------------------------------------\n",
"NM_ALL_G_OBJ: \nAll objective gradient elements have been set.\n",
"NM_SOME_G_OBJ: \nThe user sets %1ld out of %1ld objective gradient elements.\n\
Each iteration %1ld gradient elements will be estimated numerically.\n",
"NM_SOME_CONST_G_OBJ: \n%1ld missing objective gradient elements assigned as \
constant elements.\n",
"NM_ALL_CONST_G_OBJ: \nAll missing objective gradient elements found and assigned \
to be constants.\n",
"NM_SIMPLE_CHECK: \nSimple Check:\n",
"NM_NO_SIMPLE_CHECK: \nSimple check not possible: every column contains a constant\n\
or missing element.\n",
"NM_WRONG_GRAD: \nXXX The objective gradient seems to be incorrect.\n",
"NM_OK_GRAD: \nThe objective gradient seems to be ok.\n",
"NM_JAC_OK: \nThe Jacobian seems to be ok.\n",
"NM_JAC_NOK: \nXXX The Jacobian seems to be incorrect.\n",
"NM_SIMPLE_G_RES: \nDirectional derivative of the objective %18.8e\n\
Difference approximation                %18.8e\n",
"NM_SUBFUN_REL_ERR: \nThe largest relative error was %7.2e in subfunction %-3ld\n",
"NM_COMP_CHECK: \n\nComponent-wise check:\n",
"NM_SUBFUN_NO: \nSubfunction %1ld\n",
"NM_JAC_ELEM_HEAD:   Elem.   x[j]     dx[j]     Jacobian value  Difference approxn.  Itns.\n",
"NM_OBJ_COMP_HEADINGS: \n   i     x[i]       dx[i]         g[i]       \
Difference approxn.   Itns.\n",
"NM_CON_COMP_HEADINGS: DUMMY\n",
"NM_COMP_VALS: %4ld %10.2e %10.2e %16.8e %16.8e %-4s  %2ld %s\n",
"NM_ASSUME_ZERO: %4ld %10.2e              g element not assigned, \
assumed zero.\n",
"NM_ZERO_ELEMENT: \nOne zero assigned gradient element has been checked and\n\
confirmed as constant. Element is",
"NM_ZERO_ELEMENTS: \nThe %1ld zero assigned gradient elements have been checked and\n\
confirmed as constants. Elements are:\n",
"NM_ZERO_ELEMENT_IS: \nThe zero element is",
"NM_ZERO_ELEMENTS_ARE: \nThe zero elements are:\n",
"NM_ALL_OBJ_GOOD: \n%1ld objective gradient elements out of the %1ld assigned,\n\
set in columns %1ld through %1ld, seem to be ok.\n",
"NM_SOME_OBJ_WRONG: \nXXX There seem to be %1ld incorrect objective \
gradient elements out of\nthe %1ld assigned set in columns %1ld through %1ld.\n",
"NM_ALL_JAC_GOOD: \n%1ld Jacobian elements out of the %1ld\n\
set in cols %1ld through %1ld seem to be ok.\n",
"NM_SOME_JAC_WRONG: \nXXX There seem to be %1ld incorrect Jacobian elements out\n    of the %1ld set in cols %1ld through %1ld.\n",
"NM_SMALL_DERIV: \n%1ld objective \
gradient elements out of the %1ld assigned may be incorrect.\n\
However, if the derivatives are small the finite difference approximations\n\
could be in error. The check should be repeated using a different\n\
starting point.\n",
"NM_SIZE_ERROR: \nThe largest relative error was %7.2e in element %-3ld\n",
"NM_SIZE_JAC_ERROR: \nThe largest relative error was %7.2e in row %1ld column %1ld\n", 
"NM_CHANGE_DER_LEVEL: Derivative level changed to %s.\n",
"NM_CON_GRAD: \n\nVerification of the constraint gradients.\n\
-----------------------------------------\n",
"NM_ALL_G_CON: \nAll constraint gradient elements have been set.\n",
"NM_SOME_G_CON: \nThe user sets %1ld out of %1ld constraint gradient elements.\n\
Each iteration, %1ld gradient elements will be estimated numerically.\n",
"NM_CON_REL_ERR: \nThe largest relative error was %7.2e in constraint %-3ld\n", 
"NM_CON_NO: \nConstraint %1ld\n", 
"NM_SOME_CONST_G_CON: \n%1ld constant constraint gradient elements assigned.\n",
"NM_ALL_CONST_G_CON: \nAll missing Jacobian elements are constants,\n",
"NM_FIN_DIFF_HEAD: \n\nFinite difference intervals.\n\
----------------------------\n    j     x[j]       Forward dx[j]    Central dx[j]    Error Est.\n",
"NM_FIN_DIFF_INT: %5ld %10.2e %16.6e %16.6e %16.6e\n",

/* Option field list output */
"NM_PARAM_TITLE: \nParameters to %s\n--------------------\n\n",
"NM_NUM_VAR: Number of variables...........%3ld\n\n",
"NM_NUM_VAR_NONEWL: Number of variables...........%3ld\n",
"NM_NUM_RESIDUALS: Number of residuals...........%3ld    ",
"NM_NUM_MAT_ROWS: Objective matrix rows.........%3ld    ",
"NM_MAX_LINE_STEP: max_line_step...........%9.2e",
"NM_MINLIN: minlin............%s",
"NM_MACH_PREC:     machine precision.......%9.2e\n",
"NM_MACH_PREC_L: machine precision.......%9.2e\n",
"NM_OPTIM_TOL: optim_tol...............%9.2e",
"NM_LINESEARCH_TOL:     linesearch_tol..........%9.2e\n",
"NM_F_EST: f_est...................%9.2e",
"NM_F_PREC:     f_prec..................%9.2e\n",
"NM_F_PREC_L: f_prec..................%9.2e",
"NM_VERIFY_GRAD: verify_grad......%s",
"NM_MAX_ITER:     max_iter................%9ld\n",
"NM_MAX_ITER_1: max_iter................%9ld",
"NM_PRINT_LEVEL: print_level...%s",
"NM_PRINT_GCHECK:     print_gcheck............    %s\n",
"NM_OUTFILE: outfile.................%9s\n",
"NM_DEBUG: debug...................     %s    \
debug_iter..............%9ld\n",
"NM_STEP_MAX: step_max................%9.2e",
"NM_PRINT_ITER: print_iter..............%9ld",
"NM_DERIV_CHECK:     deriv_check.............    %s\n",
"NM_DERIV_CHECK_2: deriv_check.............    %s\n",
"NM_INIT_STATE: init_state....%s",
"NM_LOCAL_SEARCH:     local_search............    %s\n",
"NM_MEM_ALLOC: \nMemory allocation:\n",
"NM_S: s.......................     %s\n",
"NM_V: v.......................     %s    tdv.....................%9ld\n",
"NM_DELTA: delta...................     %s\n",
"NM_STATE: state...................     %s\n",
"NM_HESL: hesl....................     %s",
"NM_HESD:     hesd...................      %s\n",
 
"NM_NUM_LIN_CON: Linear constraints............%3ld    ",
"NM_PROB: prob....................%s    ",
"NM_START: start...................%s\n",
"NM_NUM_NONLIN_CON: Nonlinear constraints.........%3ld\n",
"NM_INF_BOUND: inf_bound...............%9.2e    ",
"NM_INF_STEP: inf_step................%9.2e\n",
"NM_RANK_TOL: rank_tol................%9.2e    ",
"NM_RANK_TOL_R: rank_tol................%9.2e\n",
"NM_CRASH_TOL: crash_tol...............%9.2e\n",
"NM_CRASH_TOL_L: crash_tol...............%9.2e",
"NM_FTOL: ftol....................%9.2e    ",
"NM_RESET_FTOL: reset_ftol..............%9ld\n",
"NM_FCHECK: fcheck..................%9ld    ",
"NM_MAX_DF: max_df..................%9ld\n",
"NM_HROWS: hrows...................%9ld",
"NM_MAX_FREE: Max. free variables.......%3ld    ", /* not documented */
"NM_MAX_ACTIVE:     Max. active constraints...%3ld\n", /* not documented */
"NM_FMAX_ITER: fmax_iter...............%9ld",
"NM_AX: ax......................     %s",
"NM_LAMBDA:     lambda.................      %s\n",
"NM_STEP_LIMIT: step_limit..............%9.2e",
"NM_LIN_FEAS_TOL: lin_feas_tol............%9.2e",
"NM_NONLIN_FEAS_TOL:     nonlin_feas_tol.........%9.2e\n",
"NM_MINOR_MAX_ITER:     minor_max_iter..........%9ld\n",
"NM_HESSIAN: hessian.................    %s",
"NM_HESSIAN_R:     hessian.................    %s\n",
"NM_F_DIFF_INT: f_diff_int..............     %9.2e",
"NM_F_DIFF_INT_AUTO: f_diff_int..............Automatic",
"NM_C_DIFF_INT:     c_diff_int..............     %9.2e\n",
"NM_C_DIFF_INT_AUTO:     c_diff_int..............Automatic\n",
"NM_OBJ_DERIV: obj_deriv...............    %s",
"NM_CON_DERIV:     con_deriv...............    %s\n",
"NM_PRINT_DERIV:     print_deriv.........%s\n",
"NM_MINOR_PRINT_LEVEL:     minor_print_level...%s\n",
"NM_MIN_INFEAS: min_infeas..............    %s\n",
"NM_H_RESET_FREQ:     h_reset_freq............%9ld\n",
"NM_H_UNIT_INIT: h_unit_init.............    %s\n",
"NM_LAMBDA_L: lambda..................     %s",
"NM_PROB_NAME: prob_name............... %8.8s\n",
"NM_OBJ_NAME: obj_name................ %8.8s    ",
"NM_RHS_NAME: rhs_name................ %8.8s\n",
"NM_BND_NAME: bnd_name................ %8.8s\n",
"NM_RANGE_NAME: range_name.............. %8.8s    ",
"NM_NAMES: crnames..............%s\n",
"NM_NEW_LINE: \n",
"NM_BIGNUM_LIN_CON: Linear constraints......%9ld    ",
"NM_SPARSE_PROB_TYPE: Problem type............sparse %s    ",
"NM_HESSIAN_COLUMNS: Hessian columns.........%9ld    ",
"NM_MIN_OR_MAX: minimize................    %s    ",                        
"NM_SCALE: scale.............%s    ",
"NM_CRASH: crash.............%s    ",
"NM_PARTIAL_PRICE: partial_price...........%9ld    ",
"NM_SUPERBASICS_LIMIT: max_sb..................%9ld    ",
"NM_LU_FACTOR_TOL: lu_factor_tol...........%9.2e    ",
"NM_LU_SING_TOL: lu_sing_tol.............%9.2e    ",
"NM_STATE_L: state...................     %s    ",
"NM_BIGNUM_VAR: Number of variables.....%9ld\n",
"NM_LU_UPDATE_TOL: lu_update_tol...........%9.2e\n",
"NM_FACTOR_FREQ: factor_freq.............%9ld\n",
"NM_PIVOT_TOL: pivot_tol...............%9.2e\n",
"NM_SCALE_TOL: scale_tol...............%9.2e\n",
"NM_DERIV_WANT: deriv_want......%s    ",
"NM_USE_HFWD_INIT: use_hfwd_init...........    %s    ",
"NM_PRINT_DERIV_L: print_deriv.........%s    ",
/* e04mzc - sparse MPS reader */
"NM_COL_LO_DEFAULT: col_lo_default..........%9.2e",
"NM_COL_UP_DEFAULT:     col_up_default..........%9.2e\n",
"NM_INFINITY: infinity................%9.2e\n",
"NM_NCOL_APPROX: ncol_approx.............%9ld",
"NM_NROW_APPROX:     nrow_approx.............%9ld\n",
"NM_MPS_SPARSE: sparse..................    %s    ",
"NM_MPS_EST_DENSITY: est_density.............%9.2e\n",
"NM_MPS_OUTLEVEL: output_level......%s",

/* Option range checking output */
"NM_SET_TEXT: %s set to %s\n",
"NM_SET_DOUBLE: %s set to %7.2e\n",
"NM_SET_INT: %s set to %1ld\n",
"NM_RESET_DEF: This change will cause the default value of %s\n\
to be reset to %7.2e.\n",
"NM_RESET: The value of %s has been reset from %7.2e to %7.2e\n\
because the valid range for %s is: %s.\n",
"NM_TEXT_VALID_RANGE: The valid range for %s is: %s.\n",
"NM_SET_DOUBLE_MIN: %s set to its minimum value %7.2e\n",

/* e04xyc introduction */
"NM_INTRO: \nOptional parameter setting for %s.\n\
--------------------------------------\n\n\
Option file: %s\n\n",

/* Set selected option values by typing structure member name \
* and required value,\ne.g. to select the option \
* for the printing out of the result of each iteration\n\
* and the final solution type:\n\n\
* print_level = Nag_Soln_Iter\n\n\
* To terminate the interactive session type end.\n\n"
*/

/* e04ucc */
"NM_QPNL_EX_MAJ_MIN: Exit from NP problem after %1ld major iterations, %1ld minor iterations.\n",
"NM_QPNL_EX_MAJOR: Exit from %s problem after %1ld iterations.\n",
"NM_QPNL_HEAD_V: \n\nVarbl State     Value      Lower Bound   Upper Bound    Lagr Mult    Residual\n",
"NM_QPNL_HEAD_L: \n\nL Con State     Value      Lower Bound   Upper Bound    Lagr Mult    Residual\n",
"NM_QPNL_VAR_CON_RES: %1s%4ld %1s %2s %13.5e %13s %13s %12s %12s\n",
"NM_QPNL_HEAD_N: \n\nN Con State     Value      Lower Bound   Upper Bound    Lagr Mult    Residual\n",
"NM_QPNL_ncj_iter: \n\n %s iteration %5ld\n",
"NM_QPNL_ncj_head1: \n\n  Itn Jdel    Jadd      Step  Ninf  Sinf/Objective  Bnd  Lin  Art   Zr   Norm Gz   Norm Gf    Cond T\n",
"NM_QPNL_ncj_head2: \n\n  Itn Jdel    Jadd      Step  Ninf  Sinf/Objective  Bnd  Lin  Art   Zr   Norm Gz   Norm Gf    Cond T   Cond Rz\n",
"NM_QPNL_ncj_iter_line_1: %4ld %4ld %s %4ld %s %8.1e %3ld %16.8e %4ld %4ld %4ld %4ld %9.1e %9.1e %9.1e\n",
"NM_QPNL_ncj_iter_line_2: %4ld %4ld %s %4ld %s %8.1e %3ld %16.8e %4ld %4ld %4ld %4ld %9.1e %9.1e %9.1e %9.1e\n",
"NM_QPNL_ncj_head3: \n\n  Itn     Step  Ninf  Sinf/Objective  Norm Gz\n",
"NM_QPNL_ncj_iter_line_3: %4ld %9.1e %4ld %15.6e %9.1e\n",
"NM_QPNL_ncj_head4: \n Values and status of the %s constraints\n ---------------------------------------\n\n Variables...\n",
"NM_QPNL_ncj_state_1:  %15.6e %4ld\n",
"NM_QPNL_ncj_head5: \n General linear constraints...\n",
"NM_QPNL_ncj_state_2:  %15.6e %4ld\n",
"NM_QPNL_ncj_work_set_fac1: \n Diagonals of %s working set factor T\n",
"NM_QPNL_ncj_work:  %15.6e",
"NM_QPNL_ncj_diag_triang_r:  Diagonals of %s triangle R\n",
"NM_QPNL_ncj_val_diag_elem_r:  %15.6e",
"NM_QPNL_ncj_dash_1:  --------------------------------------------------\n",
"NM_QPNL_ncj_head6: Multipliers for the artificial constraints\n\n",
"NM_QPNL_ncj_gq_val:      %11.2e",
"NM_QPNL_ncj_head7: Multipliers for the %s bound constraints\n",
"NM_QPNL_ncj_kx_free_k:  %1ld",
"NM_QPNL_ncj_rlamda_activ:  %11.2e ",
"NM_QPNL_ncj_head8: Multipliers for the %s linear constraints\n",
"NM_QPNL_ncj_kactiv:  %1ld",
"NM_QPNL_ncj_rlamda_k:  %11.2e ",
"NM_QPNL_ucc_lin_con_feas: \nThe linear constraints are feasible.\n",
"NM_QPNL_ucc_bad_init_hess_r_refact:  XXX  Bad initial Hessian,   R  refactorized.\n",
"NM_QPNL_ucc_r_line_end: R = \n",
"NM_QPNL_ucc_r_elem: %10.5f ",
"NM_QPNL_ucc_w_line_end: W = \n",
"NM_QPNL_ucc_w_lwrk1_j: %10.5f ",
"NM_QPNL_ucc_opt_found: \nOptimal solution found.\n",
"NM_QPNL_ucc_final_obj_val: \nFinal objective value = %16.7e\n",
"NM_QPNL_ucc_min_sum_infeas: \nMinimum sum of infeasibilities = %16.7e\n",
"NM_QPNL_ucc_final_sum_infeas: \nFinal sum of infeasibilities = %16.7e\n",
"NM_QPNL_ucc_maj_iter: \n\nMajor iteration %1ld\n",
"NM_QPNL_uct_maj_iter_head: \n\nMajor iteration %1ld\n===================\n\n",
"NM_QPNL_uct_prthdr_nlncon: \nMaj  Mnr  Step    Nfun  Merit function  Violtn  Norm Gz   Nz  Bnd  Lin  Nln  Penalty  Norm Gf  Cond H   Cond Hz  Cond T  Conv\n",
"NM_QPNL_uct_prthdr_Nlncon: \nMaj  Mnr  Step   Nfun       Objective  Norm Gz   Nz  Bnd  Lin  Norm Gf  Cond H  Cond Hz  Cond T  Conv\n",
"NM_QPNL_uct_nlncon_line1:  %3ld %3ld %8.1e %3ld %15.8e %8.1e %8.1e %4ld %4ld %4ld %4ld %8.1e %8.1e %8.1e %8.1e %8.1e %s %s %s %s\n",
"NM_QPNL_uct_Nnlncon_line1:  %3ld %3ld %8.1e %3ld %15.8e %8.1e %4ld %4ld %4ld %8.1e %8.1e %8.1e %8.1e %s %s %s %s\n",
"NM_QPNL_uct_prthdr_nlncon_def:\n\n   Maj  Mnr    Step   Merit function  Violtn  Norm Gz  Cond Hz\n",
"NM_QPNL_uct_prthdr_Nnlncon_def:\n  Maj  Mnr    Step      Objective Norm Gz Cond Hz\n",
"NM_QPNL_uct_nlncon_def:  %4ld %4ld %8.1e %15.6e %8.1e %8.1e %8.1e  %s\n",
"NM_QPNL_uct_Nnlncon_def:  %4ld %4ld %8.1e %15.6e %8.1e %8.1e  %s\n",
"NM_QPNL_uct_ncnln_zero_obj:  Nonlinear objective value = %15.6e\n",
"NM_QPNL_uct_ncnln_Nzero_obj: \n Nonlinear objective value = %15.6e\n",
"NM_QPNL_uct_ncnln_Nzero_cons:  Norm of the nonlinear constraint violations = %15.6e\n",
"NM_QPNL_uct_cons_head: \nValues of the constraints and their predicted status\n----------------------------------------------------\nVariables\n",
"NM_QPNL_uct_cons_istate: %15.6e %4ld",
"NM_QPNL_uct_gen_cons_head: General linear constraints\n",
"NM_QPNL_uct_gen_cons_istate: %15.6e %4ld",
"NM_QPNL_uct_nlin_head: Nonlinear constraints\n",
"NM_QPNL_uct_nlin_state: %15.6e %4ld",
"NM_QPNL_uct_diag_t_head: Diagonals of  T  = \n",
"NM_QPNL_uct_diag_t_val: %15.6e ",
"NM_QPNL_uct_diag_r_head: Diagonals of  R  = \n",
"NM_QPNL_uct_diag_r_val: %15.6e ",
"NM_QPNL_uct_dashed: ===============================================\n\n",
"NM_QPNL_ucu_maj_itn_head: \n-------------------------------\nStart of major itn. %1ld\n",
"NM_QPNL_ucz_minor_itn_head: Minor itn %1ld -- Re-solve QP subproblem.\n",

/* e04nkc - sparse convex QP */
"NM_REALLOC_BASIS_WKSP: XXX Basis factors too large for workspace.\n\
    Reallocating and restarting....\n",
"NM_PARP_REDUCED: XXX Partial price reduced from %1ld to %1ld\n",
"NM_PART_PRICE_INFO: \nScale option %-20s Partial price%7ld\n\
Partial price section size (A)               %9ld\n\
Partial price section size (I)               %9ld\n",
"NM_COLD_START: \nCold Start\n----------\n",
"NM_WARM_START: \nWarm Start\n----------\n",
"NM_ITERATIONS: \n\nIterations\n----------\n",
"NM_NAME_ITER_OBJ: \nProblem name                  %8s\n\
No. of iterations      %15ld   Objective value  %27.10e\n",
"NM_INFEAS_SUMM: No. of infeasibilities %15ld   Sum of infeas    %27.10e\n",
"NM_FEAS_SPQP_SUMM: No. of Hessian-vector products%8ld   Linear objective value %21.10e\n\
                                         Quadratic objective value %18.10e\n",
"NM_SBASIC_SUMM: No. of superbasics  %18ld   No. of basic nonlinears %20ld\n",
"NM_DEGEN_SUMM: No of degenerate steps  %14ld   Percentage      %28.2f\n",
"NM_X_PI_NORM_SCALED: Norm of x     (scaled) %15.1e   Norm of pi  (scaled) %23.1e\n",
"NM_X_PI_NORM_UNSCALED: Norm of x  %27.1e   Norm of pi       %27.1e\n",
"NM_MAX_INFEAS_SCALED: Max Primal inf(scaled) %7ld%8.1e   Max Dual inf(scaled)       %9ld%8.1e\n",
"NM_MAX_INFEAS_UNSCALED: Max Primal infeas      %7ld%8.1e   Max Dual infeas            %9ld%8.1e\n",
"NM_CRASH_OPTION: \n\nCrash option  %s\n",
"NM_CRASH_LIN_EQ: \nCrash on linear E  rows (equalities):\n",
"NM_CRASH_LIN_INEQ: \nCrash on linear LG rows (inequalities):\n",
"NM_CRASH_NONLIN: \nCrash on nonlinear rows:\n",
"NM_CRASH_DETAILS: Slacks   %6ld   Free cols%6ld   Preferred%6ld\n\
Unit     %6ld   Double   %6ld   Triangle %6ld   Pad      %6ld\n",
"NM_SCALE_HEAD: \n\nScaling\n-------\n\
                   Min elem     Max elem       Max col ratio\n",
"NM_SCALE_PASS: After pass %3ld %12.2e %12.2e %19.2f\n",
"NM_SCALE_PHASE_SUMM: \n            Min Scale                       Max scale\n\
%3.3s%7ld %10.1e           %3.3s%7ld %10.1e\n\
No. of %3.3s scales between 0.5 and 2.0 = %1ld (or %1.1f %%)\n",
"NM_SCALE_NORMS: Norm of fixed columns and slacks    %20.1e\n\
 (before and after row scaling)     %20.1e\n",
"NM_SCALE_REDUCED: Scales are large --- reduced by %1.1e.\n",
"NM_ROW_SCALES_HEAD: \n\nRow scales  r(i)         a(i,j)  =   r(i)  *  scaled a(i,j)  /  c(j)\n\
----------------\n",
"NM_COL_SCALES_HEAD: \n\nCol scales  c(i)      x(j)    =   c(j)  *  scaled x(j)\n\
----------------\n",

"NM_SCALES: %6ld%16.5e",
"NM_FACT_ITN: \nFactorize%6ld   Demand   %6ld   Iteration%6ld\n",
"NM_BS_FACT_NSWAP: XXX BS Factorize. No. of (B S) changes = %1ld\n",
"NM_BFACT_STATS: Nonlinear%6ld   Linear   %6ld   Slacks   %6ld   Elems    %6ld   Density%8.2f\n",
"NM_ROW_CHECK: Row check at itn %1ld.  Max residual = %.2e on row %1ld.\n\
Norm x = %.2e",
"NM_BFACT_DETAILS: Compressns%5ld   Merit%10.2f   lenL%11ld   lenU%11ld   Increase%7.2f\
   m %6ld  Ut%6ld  d1%6ld\n\
Lmax%11.1e   Bmax%11.1e   Umax%11.1e   Umin%11.1e   Growth%9.1e   Lt%6ld\
  bp%6ld  d2%6ld\n",
"NM_BSFACT_DETAILS: Compressns%5ld   Merit%10.2f   lenL%11ld   lenU%11ld   Increase%7.2f\
   m %6ld  Ut%6ld  d1%6ld\n\
                  BSmax%10.1e                                                         \
Lt%6ld  bp%6ld  d2%6ld\n",
"NM_BFACT_LUFAC_RED: Itn %1ld -- LU Factor Tolerance reduced from %.2e to %.2e.\n",
"NM_BFACT_LUUPD_RED: Itn %1ld -- LU Update Tolerance reduced from %.2e to %.2e.\n",
"NM_COL_REPL_SLACK: Column %1ld replaced by slack %1ld\n",
"NM_COL_REPL_SLACK_ETC: ...and so on.  Total slacks inserted = %1ld.\n",
"NM_MAT_STATS_HEAD: \nMatrix statistics\n-----------------\n      \
        Total  Normal (lower or upper bound)  Free (no bounds)  \
Fixed (equal bounds)  Bounded (distinct bounds)\n",
"NM_MAT_STATS: %s%12ld                   %12ld      %12ld\
          %12ld               %12ld\n",
"NM_MAT_DETAILS: \nNo. of matrix elements %20ld\n\
Density %35.1f  (excluding fixed columns, free rows, and RHS)\n\
Biggest %35.4e\nSmallest %34.4e\n\n\
No. of objective coefficients %13ld\n",
"NM_OBCOEF_MINMAX: Biggest %35.4e  (excluding fixed columns)\n\
Smallest %34.4e\n",
"NM_NLN_FAILED: XXX e04mzh: e04nln failed.  kbsq = %1ld\n",
"NM_NO_MAX_PIV: XXX e04mzh: No maximum pivot. Problem abandoned.\n",
"NM_BIGDJ: Biggest dj = %11.3e (variable %1ld)\n",
"NM_NORMS_RG_PI: Norm rg    = %11.3e (reduced gradient norm)\n\
Norm pi    = %11.3e (norm of dual variables)\n",
"NM_LU_FILE_GROW:  LU file has increased by a factor of %6.1f\n",
"NM_HESS_INDEF_ITN: Itn %7ld -- Hessian numerically indefinite.\n",
"NM_REDHESS_INDEF_ITN: Itn %7ld -- Indefinite reduced Hessian.\n",
"NM_HESS_DIAGS:                Square of diag, min diag = %9.1e, %9.1e.\n",
"NM_ITN_INFEAS: Itn %1ld -- Infeasible\n",
"NM_ITN_INFEAS_LONG: Itn %1ld -- Infeasible.  No. of infeasibilities = %1ld. Sum of infeasibilities = %1.1e.\n",
"NM_WRONG_NO_BASICS: XXX WARNING - %1ld basics specified.\n\
    Preferably should have been %1ld.\n",
"NM_ITN_FEAS_EQ: Itn %1ld -- Feasible point found (for %1ld equality constraints).\n",
"NM_SMALL_DIR_DERIV: XXX Small directional derivative = %1.1e.\n",
"NM_KSAVE_NOT_FOUND: XXX e04nlm: ksave not found. jqsave = %1ld.\n",
"NM_MAXPIV_TOO_SMALL: XXX e04nln: Max pivot (= %1.1) is too small.\n",
"NM_BASICS_RECOMP: Itn %1ld -- %1ld nonbasics set on bound, basics recomputed.\n",
"NM_HESS_SING: Itn %7ld -- Factorize detected a singular Hessian.\n\
               Rank = %1ld, Diag, min diag = %1.1e, %1.1e.\n",
"NM_VAR_INDEF: Itn %1ld -- Variable %1ld can %s indefinitely.\n",
"NM_SPQP_ITN_SHD_QP: \n   Itn      Step   Ninf   Sinf/Objective   Norm rg\n",
"NM_SPQP_ITN_SHD_FPLP: \n   Itn      Step   Ninf   Sinf/Objective\n",
"NM_ITN_SPQP_QP_SHORT: %6ld %9.1e %6ld %16.6e %9.1e\n",
"NM_ITN_SPQP_FPLP_SHORT: %6ld %9.1e %6ld %16.6e\n",
"NM_SPQP_ITN_LHD_QP: \n  Itn  pp        dj    +S    -S    -B    -B      Step     \
Pivot Ninf   Sinf/Objective     L     U Ncp   Norm rg   Ns   Cond Hz\n",
"NM_ITN_SPQP_QP_LONG: %5ld %3ld %9.1e %5ld %5ld %5ld %5ld %9.1e %9.1e %4ld %16.8e \
%5ld %5ld %3ld %9.1e %4ld %9.1e\n",
"NM_SPQP_ITN_LHD_FPLP: \n  Itn  pp        dj    +S    -S    -B    -B      Step     \
Pivot Ninf   Sinf/Objective     L     U Ncp\n",
"NM_ITN_SPQP_FPLP_LONG: %5ld %3ld %9.1e %5ld %5ld %5ld %5ld %9.1e %9.1e %4ld %16.8e \
%5ld %5ld %3ld\n",
"NM_ITN_FEAS_STRING: Itn %1ld -- Feasible %s\n",
"NM_PROB_DESCR: \nName            %8s\n\nStatus      %12s\n\
\nObjective       %8s  (%s)\n\
RHS             %8s\nBounds          %8s\nRanges          %8s\n",
"NM_SPQP_SOL_HEAD: \n%8s State         Value   Lower Bound  Upper Bound   \
Lagr Mult    Residual\n",
"NM_SPQP_SOL_LINE: %8.8s %1.1s %3.3s %13.5e %13s%13s %11.3e %11s\n",
"NM_SPQP_COL_SOL_HEAD_LONG: \nSection 1 - Columns\n\n\
  Number   Column  State      Activity    Obj Gradient    \
Lower Bound   Upper Bound   Reduced Gradnt     m+j\n",
"NM_SPQP_ROW_SOL_HEAD_LONG: \nSection 2 - Rows\n\n\
  Number      Row  State      Activity  Slack Activity    \
Lower Bound   Upper Bound    Dual Activity       i\n",
"NM_SPQP_SOL_LINE_LONG:  %7ld %8.8s  %1.1s %3.3s %13.5e   %13.5e  %13s %13s\
    %13.5e %7ld\n",
/* e04lbc */
"NM_ITER_HEADINGS_MODN_1:   Itn  Nfun  Objective    Norm g    Norm x  \
 Norm step    Step     CondH PosDef\n",
"NM_INIT_ITER_VALS_MODN_1: %4ld %5ld %12.4e %9.1e %9.1e                     %9.1e  %s\n",
"NM_ITER_VALS_MODN_1: %4ld %5ld %12.4e %9.1e %9.1e %9.1e %9.1e %9.1e  %s\n",

/* e04xac */
"NM_F_PREC_TOO_LARGE: \nXXX WARNING: options.f_prec is too large. \
Using default value instead.\n",
"NM_F_PREC_TOO_SMALL: \nXXX WARNING: options.f_prec is too small. \
Using default value instead.\n",
"NM_HESS_DERIV_HEAD: \n  j    X(j)    Fwd diff int  Cent diff int   \
Error est  Hess diag est  Nfun Info\n",
"NM_GRAD_HESS_DERIV_HEAD: \n  j    X(j)    Fwd diff int  Cent diff int\
   Error est     Grad est    Hess diag est  Nfun Info\n",
"NM_HESS_DERIV_RES: %3ld %9.2e %13.6e %13.6e %13.6e %13.6e  %5ld %s\n",
"NM_GRAD_HESS_DERIV_RES: %3ld %9.2e %13.6e %13.6e %13.6e %13.6e %13.6e  %5ld %s\n",

/* END OF MESSAGE STRINGS */

""
};
#else
extern char *nag_e04mesg[];
#endif

#endif  /* not NAG_E04MESG */

