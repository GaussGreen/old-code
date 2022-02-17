#ifndef OPDELIV_H
#define OPDELIV_H

/*	Delivery options */

double opdeliv(
    int num_bonds, /*	Number of bonds in the basket  , index 0 = cheapest */
    double *fwds,  /*	Forwards */
    double *conv_facts, /*	Conversion factors */
    double **cov,       /*	Covariance matrix */
    double mat,         /*	Maturity */
    int log_or_norm,    /*	0: norm  , 1: log */
    long npth);         /*	Number of Sobol points */

#endif