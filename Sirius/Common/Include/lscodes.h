/* termination codes: Negatives are non-algorithmic terminations */
/* positive are algorithmic terminations                        */
/* ** 2/02 jcp ** add code for time limit exceeded       */
/*-------------------------------------------------------------------*/

  enum lsgrg_termination_codes {
    _LSGRG_PROBLEM_STRUCTURE      = -17,  /* run structure different from setup structure */
    _LSGRG_ANAJAC_BAD_COL         = -16,
    _LSGRG_BAD_USER_NNZ           = -15,
    _LSGRG_ALLOCTBL_OVERFLOW      = -14,
    _LSGRG_MISSING_GCOMPX         = -13,
    _LSGRG_MISSING_PARSH          = -12,
    _LSGRG_MISSING_GCOMP          = -11,
    _LSGRG_BAD_OPTION_VALUE       = -10,
    _LSGRG_BOUNDS_ERROR           = -9,
    _LSGRG_LINEAR_VARS_ERROR      = -8,
    _LSGRG_BAD_NOBJ               = -7,
    _LSGRG_DIMENSION_ERROR        = -6,
    _LSGRG_BAD_COMMAND            = -5,
    _LSGRG_INTERNAL_ERROR         = -4,
    _LSGRG_INVERT_FAILURE         = -3,
    _LSGRG_INSFMEMORY             = -2, _LSGRG_BADINPUT              = -1,
    _LSGRG_STATUS_NOT_SET         =  0, _LSGRG_KTC                   =  1,
    _LSGRG_FRACTCHG               =  2, _LSGRG_ALLREMEDIES           =  3,
    _LSGRG_ITERATIONS             =  4, _LSGRG_UNBOUNDED             =  5,
    _LSGRG_INFEASIBLE_KTC         =  6, _LSGRG_INFEASIBLE_FRACTCHG   =  7,
    _LSGRG_INFEASIBLE_ALLREMEDIES =  8, _LSGRG_INFEASIBLE_ITERATIONS =  9,
    _LSGRG_INFEASIBLE             = 10,
    _LSGRG_SETUP_SUCCESS          = 11, _LSGRG_SHUTDOWN_SUCCESS      = 12,

    _LSGRG_REDOBJ_CONSTVIOL       = 13,
    _LSGRG_REDGRA_NB_LE_0         = 14,
    _LSGRG_XDOT_COLLEN            = 15,
    _LSGRG_XSAXPY_COLLEN          = 16,
    _LSGRG_GETBAS_INSFMEM         = 17,
    _LSGRG_XPIVOT_COLLEN          = 18,
    _LSGRG_CHUZQ_BADPIVOT         = 19,
    _LSGRG_XPIVOT_BASIS_ILLCOND   = 20,
    _LSGRG_XPIVOT_BASIS_SING      = 21,
    _LSGRG_XPIVOT_INSFMEM         = 22,
    _LSGRG_XPIVOT_OTHER_ERR       = 23,
    _LSGRG_CONSBS_REINVERT_BC     = 24,
    _LSGRG_CONSBS_BASIS_STRUCTURE = 25,
    _LSGRG_CONSBS_NOINVERT_SEARCH = 26,
    _LSGRG_CONSBS_BASIC_SLACK     = 27,
    _LSGRG_CONSBS_JPIV_0          = 28,
    _LSGRG_CONSBS_NB_TOOBIG       = 29,
    _LSGRG_CONDNM_BAD_NBLOCKS     = 30,
    _LSGRG_CONDNM_BAD_COLLEN      = 31,
    _LSGRG_DIREC_UPDATE_ERR       = 32,
    _LSGRG_PH0FAC_INVERT_FAILURE  = 33,
    _LSGRG_PH0PIV_BAD_INDEX       = 34,
    _LSGRG_PH0PIV_XPIVOT_FAILURE  = 35,
    _LSGRG_PH0PIV_BAD_ICOLS1      = 36,
    _LSGRG_PH0PIV_BAD_ICOLS2      = 37,
    _LSGRG_OTHER_RUNTIME          = 38,
    _LSGRG_USER_TERMINATION       = 39,
    _LSGRG_JAC_OVERFLOW           = 40,
    _LSGRG_OPTQUEST_ERROR         = 41,
    _LSGRG_TIME_LIMIT_EXCEEDED    = 42

};
