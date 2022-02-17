#ifndef OPSABRCALIBBVM_H
#define OPSABRCALIBBVM_H

Err srt_f_optbvmsabrvol(double Forward, double Strike, double Maturity,
                        double Sigma, double Alpha, double Delay, double Rho,
                        SrtDiffusionType input, SrtDiffusionType output,
                        double *vol);

Err srt_f_optbvmsabr_mr_vol(double Forward, double Strike, double Maturity,
                            double Sigma, double Alpha, double Delay,
                            double Rho, double Lambda, double *numericalparams,
                            int nparams, SrtDiffusionType input,
                            SrtDiffusionType output, double *vol);

Err opBVMsabrcalib(double Fwd, double maturity, int nbr_strikes,
                   double *strikes, double *marketvols, double *ATMVol,
                   double *alpha,
                   int freeze_alpha, /* if 0      , alpha is not calibrated */
                   double *lambda,
                   int freeze_lambda, /* If 0      , beta is not calibrated */
                   double *rho,
                   int freeze_rho, /* if 0      , rho is not calibrated */
                   double *fitting_error);

Err alphaHardLimit(double forward,
                   /*	Option specs	*/
                   double strike, double maturity, double putSpreadWidth,
                   /*	Model specs	*/
                   double volinput, double alpha, double beta, double rho,

                   SrtDiffusionType input,

                   double *alphaHardLimit);

Err alphaHardLimitDicho(double forward,
                        /*	Option specs	*/
                        double strike, double maturity, double putSpreadWidth,
                        /*	Model specs	*/
                        double volinput, double alpha, double beta, double rho,

                        SrtDiffusionType input,

                        double *alphaHardLimit);

Err betaHardLimitDicho(double forward,
                       /*	Option specs	*/
                       double strike, double maturity, double putSpreadWidth,
                       /*	Model specs	*/
                       double volinput, double alpha, double beta, double rho,

                       SrtDiffusionType input,

                       double *betaHardLimit);

Err rhoHardLimitDicho(double forward,
                      /*	Option specs	*/
                      double strike, double maturity, double putSpreadWidth,
                      /*	Model specs	*/
                      double volinput, double alpha, double beta, double rho,

                      SrtDiffusionType input,

                      double *rhoHardLimit);

Err SABRGenHardLimitDichoGen(
    double forward,
    /*	Option specs	*/
    double strike, double maturity, double putSpreadWidth,
    /*	Model specs	*/
    double volinput, double Alpha, double a, double b, double c, double Rho,
    double (*vol_local)(double x, double a, double b, double c, int type),
    SrtDiffusionType input,
    int HardLimType, // 0:Alpha      , 1:a      ,2:b      , 3:c      , 4:Rho
    double *HardLimit);

Err alphaHardLimitDichoGen(
    double forward,
    /*	Option specs	*/
    double strike, double maturity, double putSpreadWidth,
    /*	Model specs	*/
    double volinput, double alpha, double delay, double b, double c, double rho,
    double (*vol_local)(double x, double a, double b, double c, int type),
    SrtDiffusionType input, double *alphaHardLimit);

Err delayHardLimitDichoGen(
    double forward,
    /*	Option specs	*/
    double strike, double maturity, double putSpreadWidth,
    /*	Model specs	*/
    double volinput, double alpha, double delay, double b, double c, double rho,
    double (*vol_local)(double x, double a, double b, double c, int type),

    SrtDiffusionType input,

    double *delayHardLimit);

Err rhoHardLimitDichoGen(double forward,
                         /*	Option specs	*/
                         double strike, double maturity, double putSpreadWidth,
                         /*	Model specs	*/
                         double volinput, double alpha, double delay, double b,
                         double c, double rho,
                         double (*vol_local)(double x, double a, double b,
                                             double c, int type),

                         SrtDiffusionType input,

                         double *rhoHardLimit);

/**********************************************************************
The following functions compute SABR vol and some derivatives as an
explicit function of the ATM volatility. They come from a cpp file
***********************************************************************/

/* BVMSABR lognormal Vol */
double BVMSabrLognormalVol(double vol, double alpha, double lambda, double rho,
                           double forward, double strike, double maturity,
                           SABR_VOL_TYPE cst_vol);

/* Derivative with respect to the vovol alpha */
double BVMSabrDsigmaDalpha(double vol, double alpha, double lambda, double rho,
                           double forward, double strike, double maturity,
                           SABR_VOL_TYPE cst_vol);

/* Derivative with respect to the delay parameter lambda */
double BVMSabrDsigmaDbeta(double vol, double alpha, double lambda, double rho,
                          double forward, double strike, double maturity,
                          SABR_VOL_TYPE cst_vol);

/* Derivative with respect to the correlation parameter rho */
double BVMSabrDsigmaDrho(double vol, double alpha, double lambda, double rho,
                         double forward, double strike, double maturity,
                         SABR_VOL_TYPE cst_vol);

/* Derivative with respect to the ATM lognormal volatility */
double BVMSabrDsigmaDvolatm(double vol, double alpha, double lambda, double rho,
                            double forward, double strike, double maturity);

/**************************************************************************
*************************************************************************/

#endif
