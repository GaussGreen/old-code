
#ifndef Fx3FMultiMCH
#define Fx3FMultiMCH

#include "McEBOptimisation.h"

/*	Main function */
/*	------------- */

Err mc_main_multi_3dfx(
    /*	Time data */
    long    npaths,
    long    num_col,
    double* time,
    long*   date,
    long    nb_dates,
    /*  Do PECS adjustment */
    int do_pecs,
    /* for Optimisation of exercise boundary */
    int        do_optimisation,
    int*       optimise,
    MCEBPARAMS params,
    LINK_UND   link_und,
    void**     func_parm_tab,
    /*	Payoff function */
    Err (*payoff_func)(
        /* Event */
        double evt_date,
        double evt_time,
        void*  func_parm,
        /* Market data */
        LINK_UND link,
        double*  sv,
        /* Results */
        int     num_col,
        double* res),
    double** res);

#endif
