typedef enum {
	IMSL_LICENSE_NUMBER_ENTRY        =     50,
	IMSL_WRONG_NUM_ARGS              =    101,
	IMSL_MAJOR_VIOLATION             =    102,
	IMSL_UNKNOWN_OPTION              =    103,
	IMSL_OPTIONAL_ARG_NULL_1         =    104,
	IMSL_OPTIONAL_ARG_NULL_2         =    105,
	IMSL_INVALID_OPTION_VALUE        =    109,
	IMSL_INVALID_OPT_VALUE           =    110,
	IMSL_REAL_OUT_OF_RANGE           =    120,
	IMSL_ILLEGAL_INTEGER_2           =    122,
	IMSL_INTEGER_OUT_OF_RANGE        =    132,
	IMSL_UNEXPECTED_NULL_POINTER     =    150,
	IMSL_UNEXPECTED_NULL_POINT       =    151,
	IMSL_BAD_ERROR_TYPE              =    160,
	IMSL_OUT_OF_MEMORY               =    200,
	IMSL_OUT_OF_MEMORY_1             =    201,
	IMSL_OUT_OF_MEMORY_2             =    202,
	IMSL_OUT_OF_MEMORY_3             =    203,
	IMSL_OUT_OF_MEMORY_4             =    204,
	IMSL_OUT_OF_MEMORY_5             =    205,
	IMSL_OUT_OF_MEMORY_6             =    206,
	IMSL_OUT_OF_MEMORY_7             =    207,
	IMSL_OUT_OF_MEMORY_8             =    208,
	IMSL_OUT_OF_MEMORY_9             =    209,
	IMSL_WORKSPACE_ERROR             =    220,
	IMSL_ERROR_ACTIVE_ALLOC_SIZE     =    221,
	IMSL_SIGNAL_INT                  =    502,
	IMSL_SIGNAL_ILL                  =    504,
	IMSL_SIGNAL_ABRT                 =    506,
	IMSL_SIGNAL_FPE                  =    508,
	IMSL_SIGNAL_BUS                  =    510,
	IMSL_SIGNAL_SGEV                 =    511,
	IMSL_SIGNAL_TERM                 =    515,
	IMSL_INEQUALITY_VIOLATION        =    601,
	IMSL_INEQUALITY_VIOLATION_2      =    602,
	IMSL_INEQUALITY_VIOLATION_3      =    603,
	IMSL_LARGER_N_VALUE_NEEDED       =    604,
	IMSL_LARGER_LDA_VALUE_NEEDED     =    605,
	IMSL_LARGER_LDB_VALUE_NEEDED     =    606,
	IMSL_COMPLEX_DIVIDE_BY_ZERO      =    607,
	IMSL_LEADING_DIM_OF_A_TOO_SMALL  =    608,
	IMSL_NEED_ZERO_LE_NLCA_LT_N      =    609,
	IMSL_NEED_ZERO_LE_NUCA_LT_N      =    610,
	IMSL_LEADING_DIM_OF_B_TOO_SMALL  =    611,
	IMSL_NEED_NLCA_LE_NLCB_LT_N      =    612,
	IMSL_NEED_NUCA_LE_NUCB_LT_N      =    613,
	IMSL_LENGTH_OF_VECTORS           =    614,
	IMSL_LENGTH_OF_VECTORS_2         =    615,
	IMSL_LENGTH_OF_VECTORS_3         =    616,
	IMSL_TOO_MANY_CODIAGONALS        =    617,
	IMSL_LOWER_CODIAGONALS           =    618,
	IMSL_UPPER_CODIAGONALS           =    619,
	IMSL_NOT_ENOUGH_ROWS             =    620,
	IMSL_NOT_ENOUGH_COLUMNS          =    621,
	IMSL_NUM_OF_EQUATIONS            =   1001,
	IMSL_ILL_CONDITIONED             =   1003,
	IMSL_SINGULAR_MATRIX             =   1004,
	IMSL_SINGULAR_TRI_MATRIX         =   1009,
	IMSL_NEGATIVE_ORDER              =   1010,
	IMSL_LDA_LESS_ORDER              =   1011,
	IMSL_LDFAC_LESS_ORDER            =   1012,
	IMSL_LDAINV_LESS_ORDER           =   1013,
	IMSL_COL_DIM_LESS_ORDER          =   1014,
	IMSL_NEED_NUM_ROWS_LT_LEAD       =   1015,
	IMSL_NEED_NUM_ROWS_GT_ZERO       =   1016,
	IMSL_NEED_NUM_COLS_GT_ZERO       =   1017,
	IMSL_NONPOSITIVE_MATRIX          =   1018,
	IMSL_TRANS_MUST_EQUAL_N_T_OR_C   =   1019,
	IMSL_NEED_LDA_GE_M               =   1020,
	IMSL_NEED_COL_A_LT_DIM           =   1021,
	IMSL_NEED_QT_GT_NUM_COL          =   1022,
	IMSL_NEED_Q_GT_NUM_COL           =   1023,
	IMSL_COL_DIM_LESS_COL            =   1024,
	IMSL_SPECIFY_FACTOR_ONLY         =   1025,
	IMSL_SPECIFY_SOLVE_ONLY          =   1027,
	IMSL_BASIS_SOLVE_ONLY            =   1028,
	IMSL_OUT_OF_MEMORY_FOR_Q         =   1029,
	IMSL_NEED_POS_COL_ROWS           =   1030,
	IMSL_NRQR_GT_LDQR                =   1031,
	IMSL_NEED_ROWS_LE_LEAD_DIM       =   1032,
	IMSL_IPATH_RANGE                 =   1033,
	IMSL_IPATH_RANGE_4               =   1034,
	IMSL_NO_MEM_FOR_FAC              =   1035,
	IMSL_NO_MEM_FOR_SYS              =   1036,
	IMSL_HERMITIAN_DIAG_REAL_2       =   1037,
	IMSL_CONDITION_ONLY_SPECIFIER    =   1038,
	IMSL_INVERSE_ONLY_SPECIFIER      =   1039,
	IMSL_FACTOR_ONLY_SPECIFIER       =   1040,
	IMSL_NO_SPACE_FOR_INVERSE        =   1041,
	IMSL_NEED_LDU_GT_ZERO            =   1042,
	IMSL_NEED_LDV_GT_ZERO            =   1043,
	IMSL_CONTROL_FLAG_ERROR          =   1044,
	IMSL_SLOWCONVERGENT_MATRIX       =   1045,
	IMSL_BAD_MAX_ITERATIONS          =   1046,
	IMSL_OUT_OF_SPACE                =   1047,
	IMSL_TOL_OUT_OF_RANGE            =   1048,
	IMSL_NOT_NONNEG_DEFINITE         =   1049,
	IMSL_INCONSISTENT_EQUATIONS      =   1050,
	IMSL_INCONSISTENT_EQUATIONS_2    =   1051,
	IMSL_ROW_ELMNTS_MUST_BE_ZERO     =   1052,
	IMSL_NOT_POSITIVE_DEFINITE_M     =   1053,
	IMSL_SINGULAR_M_MATRIX           =   1054,
	IMSL_NOT_POSITIVE_DEFINITE_A     =   1055,
	IMSL_SINGULAR_A_MATRIX           =   1056,
	IMSL_NO_CONVERGENCE              =   1057,
	IMSL_BAD_JACOBI_PRECONDITION     =   1058,
	IMSL_NEED_NCODA_GE_ZERO          =   1059,
	IMSL_NEED_LDA_GT_NCODA           =   1060,
	IMSL_U_COLUMN_DIM_ERROR          =   1061,
	IMSL_COLUMN_DIM_ERROR            =   1062,
	IMSL_WRONG_IDO_VALUE             =   1063,
	IMSL_WRONG_NEVAL_VALUE           =   1064,
	IMSL_BAD_CODIAGONAL_VALUE        =   1065,
	IMSL_LDA_LARGER_THAN_NCODA       =   1066,
	IMSL_IDEF_ARGUMENT               =   1067,
	IMSL_BEST_ESTIMATE_RETURNED      =   1068,
	IMSL_NRA_EXCEEDS_LDA_VALUE       =   1069,
	IMSL_POSITIVE_NRB_NCB_VALUES     =   1070,
	IMSL_NRB_EXCEEDS_LDB_VALUE       =   1071,
	IMSL_BAD_DIMENSIONS_FOR_TRANS    =   1072,
	IMSL_BAD_NLCA_VALUE_SPECIFIED    =   1074,
	IMSL_BAD_NUCA_VALUE_SPECIFIED    =   1075,
	IMSL_TOO_MANY_ROWS_IN_MATRIX_A   =   1076,
	IMSL_TOO_MANY_ROWS_IN_FAC        =   1077,
	IMSL_INVALID_MULT_STRING         =   1078,
	IMSL_MATMUL_DIM_MISMATCH_2       =   1079,
	IMSL_MATMUL_DIM_MISMATCH_3       =   1080,
	IMSL_MAT_MUL_UNDEFINED           =   1090,
	IMSL_BAD_SOLVE_FACTOR_INVERSE    =   1091,
	IMSL_BAD_SOLVE_ONLY              =   1092,
	IMSL_BAD_FACTOR_ONLY             =   1093,
	IMSL_BAD_INVERSE_ONLY            =   1094,
	IMSL_BAD_RESULT_COL_DIM          =   1095,
	IMSL_BAD_SOLVE_FACTOR            =   1096,
	IMSL_NUM_OF_LS_EQUATIONS         =   1097,
	IMSL_NUM_OF_LS_VARIABLES         =   1098,
	IMSL_NUM_OF_LS_CONSTRAINTS       =   1099,
	IMSL_BAD_CONSTRAINT_BOUNDS       =   1100,
	IMSL_BAD_CONSTRAINT_TYPE         =   1101,
	IMSL_BAD_SCALING                 =   1102,
	IMSL_BAD_OFFSET                  =   1103,
	IMSL_BAD_SCALE_FACTOR            =   1104,
	IMSL_BAD_LS_OPTION               =   1105,
	IMSL_BAD_LS_ROW_ACCUMULATION     =   1106,
	IMSL_BAD_COLUMN_ORDER            =   1107,
	IMSL_BAD_POLARITY_FLAGS          =   1108,
	IMSL_BAD_OFFSET_FOR_RANK         =   1109,
	IMSL_BAD_OFFSET_FOR_BLOWUP       =   1110,
	IMSL_BAD_BLOWUP_RECIPROCAL       =   1111,
	IMSL_BAD_OFFSET_FOR_PRETRI       =   1112,
	IMSL_NEG_PRETRI_FACTOR           =   1113,
	IMSL_BAD_ROW_SEPARATOR           =   1114,
	IMSL_MATRIX_ORDER_TOO_SMALL      =   2001,
	IMSL_NUMBER_MAX_EIGENVALUES      =   2002,
	IMSL_LDA_VALUE_TOO_SMALL         =   2003,
	IMSL_NEED_ELOW_LESS_THAN_EHIGH   =   2004,
	IMSL_LDEVEC_VALUE_TOO_SMALL      =   2005,
	IMSL_NEVAL_MXEVAL_MISMATCH       =   2006,
	IMSL_LOST_ORTHOGONALITY          =   2007,
	IMSL_SLOW_CONVERGENCE_SYM        =   2008,
	IMSL_SLOW_CONVERGENCE_GEN        =   2009,
	IMSL_HERMITIAN_DIAG_REAL         =   2010,
	IMSL_HERMITIAN_DIAG_REAL_1       =   2011,
	IMSL_NEED_LDA_GE_N               =   2012,
	IMSL_NEED_LDB_GE_N               =   2013,
	IMSL_NEED_N_GE_ZERO              =   2014,
	IMSL_INCX_EQUALS_ZERO            =   2015,
	IMSL_INCY_EQUALS_ZERO            =   2016,
	IMSL_INVALID_UPLO_VALUE          =   2017,
	IMSL_NRA_MUST_BE_AT_LEAST_1      =   2018,
	IMSL_NCA_MUST_BE_AT_LEAST_1      =   2019,
	IMSL_NEED_NRB_EQUAL_TO_NRA       =   2020,
	IMSL_NEED_NCB_EQUAL_TO_NCA       =   2021,
	IMSL_NEED_LDA_GE_NRA             =   2022,
	IMSL_NEED_LDB_GE_NRB             =   2023,
	IMSL_INVALID_LDA_VALUE_GIVEN     =   2024,
	IMSL_NEED_A_AND_N_GT_ZERO        =   2025,
	IMSL_NEED_N_LE_LDA               =   2026,
	IMSL_NEED_LDB_GE_MATRIX_ORDER    =   2027,
	IMSL_MATRIX_B_NOT_POS_DEFINITE   =   2028,
	IMSL_MATRIX_ORDER_NOT_POSITIVE   =   2029,
	IMSL_NEED_MATRIX_ORDER_LE_LDA    =   2030,
	IMSL_NEED_N_LE_LDFAC             =   2031,
	IMSL_SUBMATRIX_NOT_POS_DEFINITE  =   2032,
	IMSL_NEED_M_GE_ZERO              =   2033,
	IMSL_NEED_LDB_GE_M_AND_GT_ZERO   =   2034,
	IMSL_SIDE_EQ_L_LDA_TOO_SMALL     =   2035,
	IMSL_NEED_NEW_LDA_SIDE_EQUALS_R  =   2036,
	IMSL_SIDE_MUST_EQUAL_L_OR_R      =   2037,
	IMSL_UPLO_MUST_EQUAL_U_OR_L      =   2038,
	IMSL_TRANSA_MUST_EQUAL_N_T_OR_C  =   2039,
	IMSL_DIAG_MUST_EQUAL_U_OR_N      =   2040,
	IMSL_NO_EIGENVALUES_RETURNED     =   2041,
	IMSL_SLOW_CONVERGENCE_2          =   2042,
	IMSL_LOST_ORTHOGONALITY_2        =   2043,
	IMSL_USE_IMSL_VECTORS_OPTION     =   2044,
	IMSL_SPLINE_ORDER_X              =   3001,
	IMSL_SPLINE_ORDER_Y              =   3002,
	IMSL_SPLINE_DERIV_X              =   3003,
	IMSL_SPLINE_DERIV_Y              =   3004,
	IMSL_SPLINE_COEFF_X              =   3005,
	IMSL_SPLINE_COEFF_Y              =   3006,
	IMSL_X_NOT_WITHIN_KNOTS          =   3007,
	IMSL_Y_NOT_WITHIN_KNOTS          =   3008,
	IMSL_SPLINE_ORDER_ARB            =   3009,
	IMSL_SPLINE_COEFF_XY             =   3010,
	IMSL_SPLINE_LEFT_ENDPT           =   3011,
	IMSL_SPLINE_RIGHT_ENDPT          =   3012,
	IMSL_SPLINE_RIGHT_ENDPT_1        =   3013,
	IMSL_SPLINE_LEFT_ENDPT_1         =   3014,
	IMSL_SPLINE_LEFT_ENDPT_2         =   3015,
	IMSL_SPLINE_RIGHT_ENDPT_2        =   3016,
	IMSL_SPLINE_RIGHT_ENDPT_3        =   3017,
	IMSL_SPLINE_LEFT_ENDPT_3         =   3018,
	IMSL_SPLINE_LIMITS_X             =   3019,
	IMSL_SPLINE_LIMITS_Y             =   3020,
	IMSL_SPLINE_DATA_X               =   3021,
	IMSL_SPLINE_DATA_Y               =   3022,
	IMSL_SPLINE_LD_FDATA             =   3023,
	IMSL_XDATA_NOT_INCREASING        =   3024,
	IMSL_YDATA_NOT_INCREASING        =   3025,
	IMSL_SPLINE_ORDER_X_OR_Y         =   3026,
	IMSL_SPLINE_BAD_XYDATA           =   3027,
	IMSL_KNOT_MULTIPLICITY           =   3028,
	IMSL_KNOT_NOT_INCREASING         =   3029,
	IMSL_ILL_COND_INTERP_PROB        =   3030,
	IMSL_SPLINE_BAD_ORDER            =   3031,
	IMSL_SPLINE_BAD_COEFFS           =   3032,
	IMSL_SPLINE_ORDER_DERIV          =   3033,
	IMSL_DUPLICATE_XDATA_VALUES      =   3034,
	IMSL_SPLINE_NEED_DATA_PTS        =   3035,
	IMSL_SPLINE_EQUAL_LIMITS         =   3036,
	IMSL_LIMITS_LOWER_TOO_SMALL      =   3037,
	IMSL_LIMITS_UPPER_TOO_SMALL      =   3038,
	IMSL_LIMITS_UPPER_TOO_BIG        =   3039,
	IMSL_LIMITS_LOWER_TOO_BIG        =   3040,
	IMSL_SPLINE_LD_FDATA_1           =   3041,
	IMSL_SPLINE_POS_DATA_PTS         =   3042,
	IMSL_SPLINE_ORDER_POSI           =   3043,
	IMSL_SPLINE_MORE_COEF_REQ        =   3044,
	IMSL_SPLINE_LOW_ACCURACY         =   3045,
	IMSL_NEGATIVE_WEIGHTS            =   3046,
	IMSL_SPLINE_NONDEC_DATA          =   3047,
	IMSL_SPLINE_SMLST_ELEMNT         =   3048,
	IMSL_SPLINE_LRGST_ELEMNT         =   3049,
	IMSL_SPLINE_NO_POS_ELMNT         =   3050,
	IMSL_SPLINE_BAD_COEFFS_2         =   3051,
	IMSL_XDATA_TOO_LARGE             =   3052,
	IMSL_XDATA_TOO_SMALL             =   3053,
	IMSL_SPLINE_BAD_ORDER_2          =   3054,
	IMSL_ILL_COND_LIN_SYS            =   3055,
	IMSL_NO_CONV_NEWTON              =   3056,
	IMSL_XGUESS_MULTIPLICITY         =   3057,
	IMSL_XGUESS_NOT_INCREASING       =   3058,
	IMSL_OPT_KNOTS_STACKED_1         =   3059,
	IMSL_OPT_KNOTS_STACKED_2         =   3060,
	IMSL_INTERVAL_TOO_SMALL          =   3061,
	IMSL_INIT_GUESS_OUT_BOUNDS       =   3062,
	IMSL_ROUNDING_ERRORS_IN_X        =   3063,
	IMSL_FINAL_VALUE_AT_BOUND        =   3064,
	IMSL_MAX_FCN_EVAL_EXCEEDED       =   3065,
	IMSL_NEED_AT_LEAST_2_PTS         =   3066,
	IMSL_WRONG_ILEFT_VALUE           =   3067,
	IMSL_WRONG_IRIGHT_VALUE          =   3068,
	IMSL_NEED_AT_LEAST_4_PTS         =   3069,
	IMSL_NOT_PERIODIC                =   3070,
	IMSL_NEED_AT_LEAST_3_PTS         =   3071,
	IMSL_NEED_LARGER_ITMAX           =   3072,
	IMSL_MAX_ITERATIONS_REACHED      =   3073,
	IMSL_DIAG_ELMNT_NEAR_ZERO        =   3074,
	IMSL_NEGATIVE_SMPAR_VALUE        =   3075,
	IMSL_INTCEP_SHOULD_BE_ZERO       =   3076,
	IMSL_NBASIS_FCNS_TOO_SMALL       =   3077,
	IMSL_NEED_AT_LEAST_1_PT          =   3078,
	IMSL_BAD_IWT_OPTION              =   3079,
	IMSL_LINEAR_DEPENDENCE           =   3080,
	IMSL_LINEAR_DEPENDENCE_CONST     =   3081,
	IMSL_NEED_NXOUT_GT_ZERO          =   3082,
	IMSL_NEED_NYOUT_GT_ZERO          =   3083,
	IMSL_LDSUR_TOO_SMALL             =   3084,
	IMSL_XOUT_NOT_STRICTLY_INCRSING  =   3085,
	IMSL_YOUT_NOT_STRICTLY_INCRSING  =   3086,
	IMSL_ALL_POINTS_COLLINEAR        =   3087,
	IMSL_DUPLICATE_XYDATA_VALUES     =   3088,
	IMSL_NEED_NRA_AND_NCA_GT_ZERO    =   3089,
	IMSL_NEED_SMALLER_NRA_VALUE      =   3090,
	IMSL_NRA_MUST_BE_POSITIVE        =   3091,
	IMSL_NRQR_GREATER_THAN_LDQR      =   3092,
	IMSL_KBASIS_IS_NEGATIVE          =   3093,
	IMSL_IPATH_RANGE_2               =   3094,
	IMSL_MATRIX_R_EXACTLY_SINGULAR   =   3095,
	IMSL_X_AND_XPERMU_LENGTH_LE_0    =   3096,
	IMSL_IPATH_RANGE_3               =   3097,
	IMSL_IPERMU_RANGE                =   3098,
	IMSL_SPLINE_ORDER_ARB_2          =   3099,
	IMSL_SPLINE_DATA_PTS_XY          =   3100,
	IMSL_KNOT_DATA_INTERLACING       =   3101,
	IMSL_DATA_TOO_LARGE              =   3102,
	IMSL_DATA_TOO_SMALL              =   3103,
	IMSL_SPLINE_CD_FDATA             =   3104,
	IMSL_COL_DIM_SUR                 =   3105,
	IMSL_NEGATIVE_WEIGHTS_2          =   3106,
	IMSL_SPLINE_BAD_ORDER_1          =   3107,
	IMSL_IDERIV_NOT_POSITIVE         =   3108,
	IMSL_NINTV_NOT_POSITIVE          =   3109,
	IMSL_KORDER_NOT_POSITIVE         =   3110,
	IMSL_X_NOT_POSITIVE              =   3111,
	IMSL_Y_NOT_POSITIVE              =   3112,
	IMSL_INVALID_RADIAL_STRUCT       =   3113,
	IMSL_DIM_LESS_THAN_ONE           =   3114,
	IMSL_NUM_CENTERS_LESS_THAN_ONE   =   3115,
	IMSL_INVALID_CENTERS_POINTER     =   3116,
	IMSL_INVALID_COEFF_POINTER       =   3117,
	IMSL_INVALID_RAD_FUNC_POINTER    =   3118,
	IMSL_NUM_POINTS_LESS_THAN_ONE    =   3119,
	IMSL_CENTERS_GT_POINTS           =   3120,
	IMSL_SEED_NEGATIVE               =   3121,
	IMSL_XVEC_NOT_INCREASING         =   3122,
	IMSL_XVEC_LENGTH                 =   3123,
	IMSL_NXVAL_POSITIVE              =   3124,
	IMSL_NHARD_GT_TOTAL              =   3125,
	IMSL_XVAL_WITHIN_KNOTS           =   3126,
	IMSL_DER_TOO_SMALL               =   3127,
	IMSL_BAD_CNSTR_TYPE              =   3128,
	IMSL_BAD_RANGE                   =   3129,
	IMSL_BAD_PERIODIC_DER            =   3130,
	IMSL_HALF_CONSTRAINT             =   3131,
	IMSL_BAD_TYPE_PAIR               =   3132,
	IMSL_BAD_INTGR_CNSTR             =   3133,
	IMSL_BAD_INTGR_RANGE             =   3134,
	IMSL_BAD_INTGR_DER               =   3135,
	IMSL_ONLY_HALF_CONSTR            =   3136,
	IMSL_DER_GE_ZERO                 =   3137,
	IMSL_CONSTR_REMOVED              =   3138,
	IMSL_NO_FIT_OBTAINED             =   3139,
	IMSL_WEIGHT_LE_ZERO              =   3140,
	IMSL_MAX_SUBINTER_SMALL          =   4001,
	IMSL_RULE_SMALL                  =   4002,
	IMSL_ERR_ABS_SMALL               =   4003,
	IMSL_ERR_REL_SMALL               =   4004,
	IMSL_ERR_TOL_ZERO                =   4005,
	IMSL_ERR_REL_BIG                 =   4006,
	IMSL_MAX_SUBINTERVALS            =   4007,
	IMSL_ROUNDOFF_CONTAMINATION      =   4008,
	IMSL_PRECISION_DEGRADATION       =   4009,
	IMSL_EXTRAPOLATION_ROUNDOFF      =   4010,
	IMSL_DIVERGENT                   =   4011,
	IMSL_INTERVAL_BOUNDS             =   4012,
	IMSL_NUM_BREAK_POINTS            =   4013,
	IMSL_MAX_MOMENTS_REACHED         =   4014,
	IMSL_WEIGHT_CHOICE               =   4015,
	IMSL_MAX_CYCLES_REACHED          =   4016,
	IMSL_BAD_INTEGRAND_BEHAVIOR      =   4017,
	IMSL_MAX_CYCLES_ACHIEVED         =   4018,
	IMSL_EXTRAPOLATION_PROBLEMS      =   4019,
	IMSL_BAD_WEIGHT_CHOICE           =   4020,
	IMSL_ALPHA_ARGUMENT              =   4021,
	IMSL_BETA_ARGUMENT               =   4022,
	IMSL_INTEGRATION_LIMITS          =   4023,
	IMSL_C_AND_A_DIFFERENT           =   4024,
	IMSL_C_AND_B_DIFFERENT           =   4025,
	IMSL_MAX_STEPS_ALLOWED           =   4026,
	IMSL_NUM_QUADRATURE_POINTS       =   4027,
	IMSL_DIM_OF_HYPER_REC            =   4028,
	IMSL_MAX_EIGEN_GT_ZERO           =   4029,
	IMSL_MAX_EVALS_TOO_LARGE         =   4030,
	IMSL_NOT_CONVERGENT              =   4031,
	IMSL_RECURRENCE_COEFF            =   4032,
	IMSL_IWEIGHT_OUT_OF_RANGE        =   4033,
	IMSL_WEIGHT_FCN_PARAMETER        =   4034,
	IMSL_WRONG_NFIX_VALUE            =   4035,
	IMSL_NEED_MORE_QUAD_POINTS       =   4036,
	IMSL_RECURRENCE_COEFF_C          =   4037,
	IMSL_FIXED_POINT_LOCATION        =   4038,
	IMSL_NO_CONVERGE_100_ITERATIONS  =   4039,
	IMSL_NFIX_CANNOT_EQUAL_TWO       =   4040,
	IMSL_BAD_FIXED_POINT_VALUE       =   4041,
	IMSL_NO_FIXED_POINTS_ALLOWED     =   4042,
	IMSL_A_MUST_BE_OUTSIDE_INTERVAL  =   4043,
	IMSL_B_MUST_BE_OUTSIDE_INTERVAL  =   4044,
	IMSL_A_AND_B_ON_OPPOSITE_SIDES   =   4045,
	IMSL_ODE_T_CHANGED               =   5001,
	IMSL_ODE_TEND_UNCHANGED          =   5002,
	IMSL_ODE_TOO_MANY_EVALS          =   5003,
	IMSL_ODE_TOO_MANY_STEPS          =   5004,
	IMSL_HMIN_GT_HMAX                =   5005,
	IMSL_ODE_FAIL                    =   5006,
	IMSL_ODE_NEG_NEQ                 =   5007,
	IMSL_ODE_NEG_TOL                 =   5010,
	IMSL_MIN_STEPSIZE_TOO_SMALL      =   5011,
	IMSL_MAX_STEPSIZE_TOO_SMALL      =   5012,
	IMSL_INVALID_IMSL_FLOOR_USAGE    =   5013,
	IMSL_INORM_OUTSIDE_OF_RANGE      =   5014,
	IMSL_INPUT_MATRIX_A_IS_SINGULAR  =   5015,
	IMSL_INCORRECT_LDA_VALUE_GIVEN   =   5016,
	IMSL_INCORRECT_LDA_GIVEN         =   5017,
	IMSL_ARGUMENT_X_CHANGED_VALUE    =   5018,
	IMSL_ARGUMENT_XEND_IS_UNCHANGED  =   5019,
	IMSL_REPEATED_ERR_TEST_FAILURE   =   5020,
	IMSL_INTEGRATION_HALTED_1        =   5021,
	IMSL_INTEGRATION_HALTED_2        =   5022,
	IMSL_TOL_TOO_SMALL_OR_STIFF      =   5023,
	IMSL_PARAM_CANNOT_BE_NEGATIVE    =   5024,
	IMSL_NEED_INORM_LESS_THAN_FOUR   =   5025,
	IMSL_WRONG_METHOD_INDICATOR      =   5026,
	IMSL_WORKSPACE_REQUIREMENT       =   5027,
	IMSL_BAD_NCODA_VALUE_GIVEN       =   5028,
	IMSL_NEED_ZERO_LE_NLCA_LT_M      =   5029,
	IMSL_WANT_ZERO_LE_NUCA_LT_N      =   5030,
	IMSL_INCREASE_LDA_VALUE          =   5031,
	IMSL_IDO_OUT_OF_RANGE            =   5032,
	IMSL_INVALID_IDO_VALUE           =   5033,
	IMSL_INVALID_IDO_VALUE_2         =   5034,
	IMSL_SEQUENCE_LENGTH             =   6001,
	IMSL_REQ_ARGUMENT_IS_NULL        =   6002,
	IMSL_NEED_NCA_GT_ZERO            =   6003,
	IMSL_NEED_NRA_GT_ZERO            =   6004,
	IMSL_NEED_NRCOEF_GT_ZERO         =   6005,
	IMSL_NEED_NCCOEF_GT_ZERO         =   6006,
	IMSL_LDA_IS_LT_NRCOEF            =   6007,
	IMSL_LDCOEF_IS_LT_NRCOEF         =   6008,
	IMSL_NO_CONVERGE_MAX_ITER        =   7001,
	IMSL_POLYNOMIAL_DEGREE           =   7002,
	IMSL_ZERO_COEFF_1                =   7003,
	IMSL_ZERO_COEFF                  =   7004,
	IMSL_FEWER_ZEROS_FOUND           =   7005,
	IMSL_POLYNOMIAL_DEGREE_50        =   7006,
	IMSL_NUMBER_EQUN_UNKN_LT_1       =   7008,
	IMSL_ERRREL_LESS_THAN_ZERO       =   7009,
	IMSL_ITMAX_LESS_THAN_ZERO        =   7010,
	IMSL_TOO_MANY_FCN_EVALS          =   7011,
	IMSL_NO_BETTER_POINT             =   7012,
	IMSL_NO_PROGRESS                 =   7013,
	IMSL_MIN_AT_BOUND                =   8001,
	IMSL_NEED_A_LE_XGUESS_LE_B       =   8003,
	IMSL_MIN_AT_LOWERBOUND           =   8004,
	IMSL_MIN_AT_UPPERBOUND           =   8005,
	IMSL_STEP_TOLERANCE              =   8006,
	IMSL_LINEAR_CONSTRAINT_VALUE     =   8007,
	IMSL_NEQ_CANNOT_BE_GT_NCON       =   8008,
	IMSL_NO_MORE_PROGRESS            =   8009,
	IMSL_SYSTEM_INCONSISTENT         =   8010,
	IMSL_NEG_CONSTRAINT_VALUE        =   8011,
	IMSL_M_IS_LARGER_THAN_LDA        =   8013,
	IMSL_BOUNDS_INCONSISTENT         =   8014,
	IMSL_WRONG_CONSTRAINT_TYPE       =   8015,
	IMSL_PROB_UNBOUNDED              =   8016,
	IMSL_TOO_MANY_ITN                =   8017,
	IMSL_PROB_INFEASIBLE             =   8018,
	IMSL_NUMERIC_DIFFICULTY          =   8019,
	IMSL_N_MUST_BE_POSITIVE          =   8020,
	IMSL_NONNEGATIVE_CONSTRAINTS     =   8021,
	IMSL_WRONG_IBTYPE_VALUE          =   8022,
	IMSL_XSCALE_DIAGONAL_LT_ZERO     =   8023,
	IMSL_UPHILL_DIRECTION            =   8024,
	IMSL_TOO_MANY_LINESEARCH         =   8025,
	IMSL_NO_PROGRESS_MADE            =   8026,
	IMSL_QP_INCONSISTENT             =   8027,
	IMSL_TOO_MANY_FCN_EVAL           =   8029,
	IMSL_TOO_MANY_GRAD_EVAL          =   8030,
	IMSL_TOO_MANY_HESSIAN_EVAL       =   8031,
	IMSL_UNBOUNDED                   =   8032,
	IMSL_INEFFICIENT_PROB_SIZE       =   8034,
	IMSL_SCALING_WITH_IDENT_MATRIX   =   8035,
	IMSL_MAXGRAD_VALUE_TOO_SMALL     =   8039,
	IMSL_MAXHES_VALUE_TOO_SMALL      =   8040,
	IMSL_FSCALE_VALUE_TOO_SMALL      =   8041,
	IMSL_NO_FURTHER_PROGRESS         =   8049,
	IMSL_REL_FCN_TOLERANCE           =   8050,
	IMSL_FALSE_CONVERGENCE           =   8051,
	IMSL_WRONG_EPSFCN_VALUE          =   8052,
	IMSL_POS_XSCALE_ELMNTS_NEEDED    =   8053,
	IMSL_NEED_POSITIVE_NUM_FCNS      =   8054,
	IMSL_NEED_LDFJAC_GE_M            =   8056,
	IMSL_TOO_MANY_VARIABLES          =   8057,
	IMSL_NEGATIVE_STEP_TOL           =   8058,
	IMSL_NEGATIVE_REL_FCN_TOL        =   8059,
	IMSL_NEGATIVE_FALSE_CONV_TOL     =   8060,
	IMSL_NEGATIVE_ABS_FCN_TOL        =   8061,
	IMSL_NEED_NONNEGATIVE_STEPMX     =   8062,
	IMSL_NEED_NONNEGATIVE_DELTA      =   8063,
	IMSL_LITTLE_FCN_CHANGE           =   8064,
	IMSL_TOO_MANY_JACOBIAN_EVAL      =   8065,
	IMSL_NEED_POSITIVE_XSCALE_ELEM   =   8066,
	IMSL_NEED_POSITIVE_FSCALE_ELEM   =   8067,
	IMSL_NEED_POSITIVE_NDIGIT        =   8068,
	IMSL_NEED_POSITIVE_MXITER        =   8069,
	IMSL_NEED_POSITIVE_MAXFCN        =   8070,
	IMSL_NEED_POSITIVE_MAXJAC        =   8071,
	IMSL_NEED_POSITIVE_GRADTL        =   8072,
	IMSL_IP_TOO_MANY_ITN             =   8073,
	IMSL_IP_OUT_OF_MEMORY            =   8074,
	IMSL_IP_NO_PARENT                =   8075,
	IMSL_IP_NO_INT_SOL_FOUND         =   8076,
	IMSL_WRONG_STRATEGY_TYPE         =   8077,
	IMSL_NEGATIVE_GAP                =   8078,
	IMSL_IP_NEED_POS_MAX_ITN         =   8079,
	IMSL_CHEBY_NEG_TERMS             =   9001,
	IMSL_CHEBY_TOO_MANY_TERMS        =   9002,
	IMSL_CHEBY_RANGE                 =   9003,
	IMSL_CHEBY_EVAL_TERMS            =   9004,
	IMSL_CHEBY_EVAL_TOL              =   9005,
	IMSL_NEGATIVE_INTEGER            =   9006,
	IMSL_F_INVERSE_OVERFLOW          =   9008,
	IMSL_LARGE_ARG_OVERFLOW          =   9009,
	IMSL_X_AND_A_ARE_TOO_LARGE       =   9010,
	IMSL_ZERO_ARG_OVERFLOW           =   9011,
	IMSL_SMALL_ARG_OVERFLOW          =   9012,
	IMSL_LARGE_ABS_ARG_OVERFLOW      =   9013,
	IMSL_SMALL_ARG_UNDERFLOW         =   9014,
	IMSL_LARGE_ARG_UNDERFLOW         =   9015,
	IMSL_LARGE_ABS_ARG_UNDERFLOW     =   9016,
	IMSL_SMALL_ABS_ARG_UNDERFLOW     =   9017,
	IMSL_BETA_UNDERFLOW              =   9018,
	IMSL_NORMAL_UNDERFLOW            =   9019,
	IMSL_NEAR_NEG_INT_WARN           =   9020,
	IMSL_NEAR_NEG_INT_FATAL          =   9021,
	IMSL_CANNOT_FIND_XMIN            =   9022,
	IMSL_CANNOT_FIND_XMAX            =   9023,
	IMSL_ARG_ZERO                    =   9024,
	IMSL_LARGE_ABS_ARG_WARN          =   9025,
	IMSL_LARGE_ABS_ARG_FATAL         =   9026,
	IMSL_NON_POSITIVE_ARGUMENT       =   9027,
	IMSL_LARGE_ARG_TERMINAL          =   9028,
	IMSL_LARGE_ARG_WARN              =   9029,
	IMSL_ERF_ALGORITHM               =   9030,
	IMSL_BETA_NEG_ARG                =   9031,
	IMSL_FIRST_ARG_LT_ZERO           =   9032,
	IMSL_SECOND_ARG_LT_ZERO          =   9033,
	IMSL_NO_CONV_200_TS_TERMS        =   9034,
	IMSL_NO_CONV_200_CF_TERMS        =   9035,
	IMSL_NEED_ZERO_LT_X_LE_A         =   9036,
	IMSL_LT_HALF_ACCURATE            =   9037,
	IMSL_BOTH_ARGS_ARE_LT_ZERO       =   9038,
	IMSL_X_IS_LESS_THAN_MINUS_1      =   9039,
	IMSL_X_IS_TOO_CLOSE_TO_NEG_1     =   9040,
	IMSL_NEED_ZERO_LT_P_LT_ONE       =   9041,
	IMSL_DF_MUST_BE_AT_LEAST_ONE     =   9042,
	IMSL_P_OUTSIDE_OPEN_INTERVAL     =   9043,
	IMSL_PIN_MUST_BE_POSITIVE        =   9044,
	IMSL_QIN_MUST_BE_POSITIVE        =   9045,
	IMSL_P_OUTSIDE_EXCLUSIVE_INT     =   9046,
	IMSL_BEST_BETIN_APPROXIMATION    =   9047,
	IMSL_DFN_OR_DFD_IS_NEGATIVE      =   9048,
	IMSL_NEED_DFN_AND_DFD_GT_ZERO    =   9049,
	IMSL_DIST_FCN_SET_TO_ZERO        =   9050,
	IMSL_DF_MUST_BE_GE_POINT_5       =   9051,
	IMSL_UNABLE_TO_BRACKET_VALUE     =   9052,
	IMSL_CHI_2_INV_CDF_CONVERGENCE   =   9053,
	IMSL_LESS_THAN_ZERO              =   9054,
	IMSL_ARG_LESS_THAN_ZERO          =   9055,
	IMSL_SHAPE_PARAMETER_NEGATIVE    =   9056,
	IMSL_NEED_N_GREATER_THAN_ZERO    =   9058,
	IMSL_BAD_PROBABILITY_VALUE       =   9059,
	IMSL_GREATER_THAN_N              =   9060,
	IMSL_NEED_ARGUMENT_GT_ZERO       =   9061,
	IMSL_LOT_SIZE_TOO_SMALL          =   9062,
	IMSL_K_GREATER_THAN_N            =   9063,
	IMSL_THETA_MUST_BE_POSITIVE      =   9064,
	IMSL_NO_CONV_300_CF_TERMS        =   9065,
	IMSL_FORMAT_INVALID              =  10001,
	IMSL_FORMAT_INVALID_CONV         =  10002,
	IMSL_FORMAT_NO_CONV              =  10003,
	IMSL_FORMAT_STAR                 =  10004,
	IMSL_FORMAT_WIDE_FIELD           =  10005,
	IMSL_FORMAT_TOO_NARROW           =  10006,
	IMSL_FORMAT_BAD_CONV             =  10007,
	IMSL_FORMAT_TOO_WIDE             =  10008,
	IMSL_NARROW_PAGE                 =  10009,
	IMSL_WIDE_PAGE                   =  10010,
	IMSL_ZERO_LINE                   =  10011,
	IMSL_LONG_TITLE                  =  10012,
	IMSL_MUT_EXC_PRINT_OPTIONS       =  10013,
	IMSL_MUT_EXC_ROW_LABEL_OPT       =  10014,
	IMSL_MUT_EXC_COL_LABEL_OPT       =  10015,
	IMSL_MUT_EXC_COV_OPTIONS         =  10016,
	IMSL_MUT_EXC_TIE_OPTION          =  10017,
	IMSL_MUT_EXC_SCORE_OPTION        =  10018,
	IMSL_ILLEGAL_WRITE_OPTION        =  10019,
	IMSL_INVALID_PAGE_OPTION         =  10020,
	IMSL_BAD_PAGE_WIDTH              =  10021,
	IMSL_BAD_PAGE_LENGTH             =  10022,
	IMSL_ILLEGAL_OPT_ARG             =  11001,
	IMSL_NOBS_LESS_THAN_ONE          =  11002,
	IMSL_POLY_DEGREE_LT_ZERO         =  11003,
	IMSL_LARGER_NOBS_REQUIRED        =  11004,
	IMSL_NEED_LDX_GE_NOBS            =  11005,
	IMSL_NCOL_MUST_BE_GE_ONE         =  11006,
	IMSL_NEED_MAXDEG_GE_ZERO         =  11007,
	IMSL_WRONG_ICRIT_VALUE           =  11008,
	IMSL_BAD_ICRIT_OR_LOF_VALU       =  11009,
	IMSL_BAD_CRIT_OR_ICRIT_VALU      =  11010,
	IMSL_BAD_CRIT_IF_ICRIT_EQ_1      =  11011,
	IMSL_CONSTANT_XVALUES            =  11012,
	IMSL_PERFECT_FIT_POLY            =  11013,
	IMSL_CONSTANT_YVALUES            =  11014,
	IMSL_FEW_DISTINCT_XVALUES        =  11015,
	IMSL_PERFECT_FIT                 =  11016,
	IMSL_WRONG_LOF_VALUE             =  11017,
	IMSL_MISSING_VALUES_IN_X         =  11018,
	IMSL_BAD_NDEG_VALUE              =  11019,
	IMSL_NEGATIVE_D_ELMNTS           =  11020,
	IMSL_WRONG_DFE_VALUE             =  11021,
	IMSL_WRONG_SSE_VALUE             =  11022,
	IMSL_BAD_DFPE_VALUE              =  11023,
	IMSL_BAD_SSPE_VALUE              =  11024,
	IMSL_WRONG_IPRINT_VALUE          =  11025,
	IMSL_NEED_LDSQSS_GE_NDEG         =  11026,
	IMSL_NEED_LDTLOF_GE_NDEG         =  11027,
	IMSL_NEED_LARGER_LDCOEF          =  11028,
	IMSL_NEGATIVE_DFE_VALUE          =  11029,
	IMSL_NEGATIVE_STDB_VALUE         =  11030,
	IMSL_STRICTLY_POS_TABLE_ELMNTS   =  11031,
	IMSL_NEED_POSITIVE_THETA         =  11032,
	IMSL_NEED_POSITIVE_PIN_VALUE     =  11033,
	IMSL_NEED_POSITIVE_QIN_VALUE     =  11034,
	IMSL_ELMNTS_SET_TO_NAN           =  11035,
	IMSL_SHAPE_PARAMETER_A           =  11036,
	IMSL_A_IS_TOO_LARGE              =  11037,
	IMSL_COVARIANCE_SPECIFIERS       =  11038,
	IMSL_COVARIANCE_MEMORY_REQ       =  11039,
	IMSL_NO_MEMORY_VARIANCE_STATS    =  11040,
	IMSL_VECTOR_OF_MEANS_X_MEMORY    =  11041,
	IMSL_RESIDUAL_VECTOR_MEMORY      =  11042,
	IMSL_WRONG_INTCEP_VALUE          =  11043,
	IMSL_MODEL_REGRESSOR_REQUIRED    =  11044,
	IMSL_NEED_LARGER_LDR_VALUE       =  11045,
	IMSL_NEED_NOBS_GE_ZERO           =  11046,
	IMSL_RANK_DEFICIENT              =  11047,
	IMSL_DOWNDATING_REQUESTED        =  11048,
	IMSL_NO_STAT_INFERENCE           =  11049,
	IMSL_WRONG_VALUE_OF_TOL          =  11050,
	IMSL_ISUB_SHOULD_BE_ZERO         =  11051,
	IMSL_NEED_ONE_LE_INDDEP_LE_NCOL  =  11052,
	IMSL_NEED_ONE_LE_INDIND_LE_NCOL  =  11053,
	IMSL_NONNEG_FREQ_REQUEST_1       =  11054,
	IMSL_NONNEG_FREQ_REQUEST_2       =  11055,
	IMSL_NONNEG_WEIGHT_REQUEST_1     =  11056,
	IMSL_NONNEG_WEIGHT_REQUEST_2     =  11057,
	IMSL_TOLERALNCE_INCONSISTENT     =  11058,
	IMSL_TOLERALNCE_INCONSISTENT_2   =  11059,
	IMSL_REMAINING_ELMNTS_NOT_ZERO   =  11060,
	IMSL_NEED_RARG_GE_ZERO           =  11061,
	IMSL_ROW_OF_X_CONTAINED_NAN      =  11062,
	IMSL_NONPOSITIVE_NROW_VALUE      =  11063,
	IMSL_CONPRM_VALUE_TOO_BIG        =  11064,
	IMSL_CONPRV_VALUE_TOO_BIG        =  11065,
	IMSL_NUM_NONMISS_OBS_LT_ZERO     =  11066,
	IMSL_LESS_THAN_TWO_VALID_OBS     =  11067,
	IMSL_ZERO_SUM_OF_WEIGHTS         =  11068,
	IMSL_SUM_OF_WEIGHTS_ZERO         =  11069,
	IMSL_MAX_LESS_THAN_MIN           =  11070,
	IMSL_MIN_GREATER_THAN_MAX        =  11071,
	IMSL_NOT_ENOUGH_OBSERVATIONS     =  11072,
	IMSL_VARIANCE_UNDERFLOW          =  11073,
	IMSL_HIGH_ORDER_UNDERFLOW        =  11074,
	IMSL_FOURTH_ORDER_UNDERFLOW      =  11075,
	IMSL_CONSTANT_OBSERVATIONS       =  11076,
	IMSL_VAR_AND_STD_ARE_NEGATIVE    =  11077,
	IMSL_COEFF_OF_VARIATION_NAN      =  11078,
	IMSL_NEGATIVE_STD_VALUE          =  11079,
	IMSL_ERROR_IN_T_STATISTIC        =  11080,
	IMSL_NEGATIVE_VARIANCE           =  11081,
	IMSL_CHI_SQUARED_STAT_ERROR      =  11082,
	IMSL_NOBS_MUST_BE_GE_TWO         =  11083,
	IMSL_NEED_LARGER_FUZZ_VALUE      =  11084,
	IMSL_WRONG_ITIE_OPTION_USED      =  11085,
	IMSL_WRONG_ISCORE_OPTION_USED    =  11086,
	IMSL_NEED_LARGER_SAMPLE_SIZE     =  11087,
	IMSL_RANK_MUST_BE_AT_LEAST_ONE   =  11088,
	IMSL_NEED_RANK_LE_SAMPLE_SIZE    =  11089,
	IMSL_NEED_N_CATEGORIES_GE_TWO    =  11090,
	IMSL_NEED_N_PARAMETERS_GE_ZERO   =  11091,
	IMSL_INCORRECT_CDF_1             =  11092,
	IMSL_INCORRECT_CDF_2             =  11093,
	IMSL_X_VALUE_OUT_OF_RANGE        =  11094,
	IMSL_MISSING_DATA_ELEMENT        =  11095,
	IMSL_INCORRECT_CDF_3             =  11096,
	IMSL_TOO_MANY_CELL_DELETIONS     =  11097,
	IMSL_INCORRECT_CDF_5             =  11098,
	IMSL_ALL_OBSERVATIONS_MISSING    =  11099,
	IMSL_INCORRECT_CDF_4             =  11100,
	IMSL_EXPECTED_VAL_LESS_THAN_5    =  11101,
	IMSL_EXPECTED_VAL_LESS_THAN_1    =  11102,
	IMSL_NO_BOUND_AFTER_100_TRYS     =  11103,
	IMSL_NO_UNIQUE_INVERSE_EXISTS    =  11104,
	IMSL_CONVERGENCE_ASSUMED         =  11105,
	IMSL_TOO_MANY_OBS_DELETED        =  11106,
	IMSL_MORE_OBS_DEL_THAN_ENTERED   =  11107,
	IMSL_DIFFERENT_OBS_DELETED       =  11108,
	IMSL_ZERO_SUM_OF_WEIGHTS_2       =  11109,
	IMSL_ZERO_SUM_OF_WEIGHTS_3       =  11110,
	IMSL_CONSTANT_VARIABLE           =  11111,
	IMSL_INSUFFICIENT_DATA           =  11112,
	IMSL_TOO_FEW_VALID_OBS_CORREL    =  11113,
	IMSL_LOF_COL_DIM_4               =  11114,
	IMSL_POLY_COL_DIM_4              =  11115,
	IMSL_CHOOSE_S1_GREATER_S2        =  11116,
	IMSL_MUT_EXCLUSIVE               =  11117,
	IMSL_MUT_EXC_TALLY_OPT           =  11118,
	IMSL_IOPT_NEED_2_INTERVALS       =  11119,
	IMSL_IOPT_NEED_3_INTERVALS       =  11120,
	IMSL_XHI_LT_XLO                  =  11121,
	IMSL_DIV_NOT_MONOTONIC           =  11122,
	IMSL_N_INTERVALS_LT_1            =  11123,
	IMSL_NAN_NOT_ALLOWED             =  11124,
	IMSL_NEED_ABS_RHO_LE_1           =  11130,
	IMSL_ABS_RHO_EQ_1                =  11131,
	IMSL_YEAR_0                      =  12001,
	IMSL_YEAR_1582                   =  12002,
	IMSL_YEAR_45                     =  12003,
	IMSL_BAD_DAY                     =  12004,
	IMSL_NEGATIVE_DAYS               =  12005,
	IMSL_DAYS_1582                   =  12006,
	IMSL_LEAP_YEAR_DAYS              =  12007,
	IMSL_DAYS_PER_YEAR               =  12008,
	IMSL_ARG_OUT_OF_RANGE            =  12009,
	IMSL_INDEX_VARIABLE_VALUE        =  12010,
	IMSL_BAD_CONST_NAME              =  12111,
	IMSL_INCOMPATIBLE_UNITS          =  12112,
	IMSL_MASS_TO_FORCE               =  12113,
	IMSL_ILLEGAL_UNIT                =  12114,
	IMSL_EXCEEDED_MAX_LOOP_COUNTER   =  12115,
	IMSL_EXCLUSIVE_SYMFAC_REQUEST    =  13000,
	IMSL_EXCLUSIVE_NUMFAC_REQUEST    =  13001,
	IMSL_BAD_SQUARE_ROOT             =  13002,
	IMSL_NZ_LESS_THAN_N              =  13003,
	IMSL_BAD_NONZERO                 =  13004,
	IMSL_NONZERO_IN_UP_TRIANGLE      =  13005,
	IMSL_BAD_DIAG_SIZE_REQUEST       =  13006,
	IMSL_BAD_RETURN_NUMNZ_REQUEST    =  13007,
	IMSL_BAD_ROW_INDEX               =  13008,
	IMSL_BAD_COLUMN_INDEX            =  13009,
	IMSL_GROWTH_FACTOR_EXCEEDED      =  13010,
	IMSL_BAD_HB_COL_PTR              =  13011,
	IMSL_BAD_HB_ROW_IND              =  13012,
	IMSL_FACTORED_MATRIX_SINGULAR    =  13013,
	IMSL_GMRES_STAGNATION            =  13014,
	IMSL_BAD_KDMAX                   =  13015,
	IMSL_RETURN_STRING_ONLY          =  50000,
	IMSL_NESTED_TOO_DEEP             =  60000,
	IMSL_ERROR_STACK_MISMATCH        =  60001,
	IMSL_SMALL_ARGUMENT              =  60002,
	IMSL_BAD_ARGUMENT_TO_N1RNOF      =  60003,
	IMSL_BAD_CHECKSUM                =  63000,
	IMSL_BAD_CHECKSUM_1              =  63001,
	IMSL_BAD_CHECKSUM_2              =  63002,
	IMSL_BAD_CHECKSUM_3              =  63003,
	IMSL_BAD_CHECKSUM_4              =  63004,
	IMSL_BAD_IOP_VALUE               =  63005,
	IMSL_BAD_IOP_VALUE_2             =  63006,
	IMSL_LAST_ERROR_MSG                  
} Imsl_code;
