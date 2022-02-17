#ifndef CHEYBETAPRICINGH
#define CHEYBETAPRICINGH

Err cheybeta_pricing_pde(
    /*	Time data		*/
    int nstept, double *time, double *date,

    /*	Discretisation	*/
    int nstepx, int nstepphi,
    /*	Model data		*/
    SrtUndPtr und,

    /*	Product data */
    void **func_parm_tab, int *eval_evt,

    /*	Market data */
    double *ifr, char *yc,

    /*	Payoff function */
    Err (*payoff_func)(double evt_date, double evt_time, void *func_parm,

                       /* Market data	*/
                       void *yc, double ifr,
                       /* Model data	*/
                       TermStruct *ts,
                       /* Gride data	*/
                       int lx, int ux, int lphi, int uphi, double *x,
                       double *phi,

                       /* Vector of results to be updated */
                       int nprod, double ***prod_val),
    /*	Result */
    int nprod, double *res);

Err cheybeta_pricing_mc(
    /*	Time data		*/
    int nstept, double *time, double *date,

    /*	Discretisation	*/
    int numpaths, int method, /* 0 balantisam      , 1 balantisam adjusted */
    SrtMCSamType gen_method,

    /*	Model data		*/
    SrtUndPtr und,

    /*	Product data */
    void **func_parm_tab, int *eval_evt,

    /*	Market data */
    double *ifr, char *yc, int *stop_vol, int *stop_lamda,

    /*	Payoff function */
    Err (*payoff_func)(double evt_date, double evt_time, void *func_parm,

                       /* Market data	*/
                       void *yc,
                       /* var of model	*/
                       double r, double x, double phi,
                       /* Vector of results to be updated */
                       int ncols, double *cols_val),
    /*	Result */
    int nprod, double **res);

#endif
