#ifndef __3F_QUANTO_GRFN_H
#define __3F_QUANTO_GRFN_H

Err grfn_payoff_4_3fquanto_tree(
    /* Event */
    double evt_date,
    double evt_time,
    void*  func_parm,
    /* Market data */
    void*  lgm2F_yc,
    double r1_lambda,
    double r2_lambda,
    double r1_phi,
    double r2_phi,
    double r1_r2_phi,
    double rho,
    void*  lgm1F_yc,
    double r3_lambda,
    double r3_phi,
    int    lgm2F_is_dom_for,
    /* Nodes data */
    long n1,
    long n2,
    long n3,
    /* i: d1, j: d2, k: d3, l = {0: r1, 1: r2, 2: r3} */
    double**** sv,
    /* Vector of results to be updated */
    long       nprod,
    double**** prod_val);

#endif
