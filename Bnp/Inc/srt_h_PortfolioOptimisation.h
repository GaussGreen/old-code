
#ifndef SRT_H_PORTFOLIOOPTIMISATION_H
#define SRT_H_PORTFOLIOOPTIMISATION_H

Err srt_f_PortfolioOpt(long size, double *mu, double *sigma, double **Correl,
                       long ndays, double DrawDown, int level,
                       double *PortfComp);

#endif