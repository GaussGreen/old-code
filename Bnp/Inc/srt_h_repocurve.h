#ifndef SRT_H_REPOCURVE_H
#define SRT_H_REPOCURVE_H

Err srt_f_init_EQRepoCurve(Date today, char *ccy, double **repo_curve,
                           long ncols, long nrows, String repo_crv_name);

#endif
