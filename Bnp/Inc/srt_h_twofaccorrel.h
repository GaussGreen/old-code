/* =======================================================================

        FILENAME:   srt_h_twofaccorrel.h

        PURPOSE:    relevant routines for a two factor correlation calibration

        AUTHOR:     O. Van Eyseren (also inspired from A. Savine)

   ======================================================================= */

#ifndef SRT_H_TWOFACCORREL_H
#define SRT_H_TWOFACCORREL_H

/* -----------------------------------------------------------------------
        Given model parameters and historical correlation matrix      ,
        returns the standard deviation of the error between model
        and historical correlation.
   ----------------------------------------------------------------------- */
Err stdev_correl(double **hist_corr, double *tenor_mat, long num_tenor,
                 double alpha, double beta, double rho, double *stdev);

/* ------------------------------------------------------------------------
        Returns in chisq(sum of errors squared )      , the difference
        between model and historical correlation matrix      , for given
        model parameters (alpha      , beta      , rho) and ifr tenors
   ------------------------------------------------------------------------ */
Err chisq_correl(double **hist_corr, double *tenor_mat, long num_tenor,
                 double alpha, double beta, double rho, double *chisq);

/* -----------------------------------------------------------------------
        Given model parameters (alpha      , beta      , rho)
        returns the full correlation matrix implied by the model      , for a
   given set of ifr tenor matrurities (in years)
   ----------------------------------------------------------------------- */
Err model_corr_matrix(double **model_corr, double *tenor_mat, long num_tenor,
                      double alpha, double beta, double rho);

/* -----------------------------------------------------------------------
        Given model parameters (alpha      , beta      , rho)
        returns the correlation implied by the model
   ----------------------------------------------------------------------- */
double model_correl(double mat1, double mat2, double alpha, double beta,
                    double rho);

/* -----------------------------------------------------------------------
        Given model parameters (alpha      , beta      , rho) and a set of dates
   (in Y) return the confidence level explain by the two eigen vectors
   ----------------------------------------------------------------------- */

Err TwoFactorEigenConf(double *plSetOfDates, int iNumberOfDates, double dAlpha,
                       double dGamma, double dRho,
                       /* return a vector of two elements */
                       double *pdConfidenceLevel);

/* -----------------------------------------------------------------------
        Given a tenor string (1m      , 1y      , 1y6m)
        returns the corresponding maturity for the model
   ----------------------------------------------------------------------- */
Err interp_corr_tenor(String fra_tenor, double *fra_mat);

/* -----------------------------------------------------------------------
        Euivalent to model_correl fro an interface use
   ----------------------------------------------------------------------- */
Err srt_f_twofac_ifr_correl(String tenor1, String tenor2, String und_name,
                            double *correl);

/* -----------------------------------------------------------------------
  Given any model      , returns the normal correl between swap rate 1 and swap
  rate 2 seen at some observation date through a grfn tableau.
   ----------------------------------------------------------------------- */

Err swap_correl(Date obs_date, Date start1, Date end1, String cmp1,
                String basis1, Date start2, Date end2, String cmp2,
                String basis2, String und, int mdl_rows, char **paramStrings,
                char **valueStrings, String n_ln, double *correl);

/*

  Given any model      , returns the normal correl between instr 1 and instr 2
  senn at some observation date through a grfn tableau.

*/

Err generic_correl(Date obs_date, String inst1, String inst2, String und,
                   int mdl_rows, char **paramStrings, char **valueStrings,
                   double *correl);

#endif
