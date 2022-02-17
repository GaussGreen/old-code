#ifndef opsabrgenericinterpH
#define opsabrgenericinterpH

#include "opsabrgeneric.h"
#include "srt_h_all.h"

Err BMM_Option_Price(double Forward, double Strike, double Maturity,
                     double Sigma, double Alpha, double Beta, double Rho,
                     double Pi, SrtCallPutType optiontype, double *Price);

Err BMM_Option_Price_quick(double Forward, double Strike, double Maturity,
                           double Sigma, double Alpha, double Beta, double Rho,
                           double Pi, SrtCallPutType optiontype, double *Price);

Err BMM_Option_Price2(double forward, double strike, double maturity,
                      double sigma, double alpha, double beta, double rho,
                      double zeta, SrtCallPutType optiontype, double *price);

Err BMM_Option_Price3(double forward, double strike, double maturity,
                      double sigma, double alpha, double beta, double rho,
                      double veta, double pi0, SrtCallPutType optiontype,
                      double *price);

Err BMMGen_Option_Price(double forward, double strike, double maturity,
                        double sigma, double alpha, double a, double b,
                        double c, double rho, double Pi,
                        SrtDiffusionType optiontype, double *price,
                        SrtDiffusionType TypeVolLoc);

Err BMM_Option_Price_From_States(double Strike, double Maturity, double Beta,
                                 double Fwd1, double Fwd2, double Sig1,
                                 double Sig2, double Pi,
                                 SrtCallPutType optiontype, double *Price);

Err BMM_Option_Price3_From_States(double forward, double strike,
                                  double maturity, double beta, double fwd,
                                  double fwd1, double fwd2, double sig,
                                  double sig1, double sig2, double pi0,
                                  double pi, SrtCallPutType optiontype,
                                  double *price);

Err BMMCalibOnSabrStates(double forward, double maturity, double sigma,
                         double beta, double alpha, double rho, double pi,
                         double NStd, double *Fwd1, double *Fwd2, double *Sig1,
                         double *Sig2, double *Pi, SrtDiffusionType TypeInput,
                         double *Calibres);

Err BMM2CalibOnSabrStates(double forward, double maturity, double sigma,
                          double beta, double alpha, double rho, double zeta,
                          double NStd, double *Fwd1, double *Fwd2, double *Sig1,
                          double *Sig2, double *Pi, SrtDiffusionType TypeInput,
                          double *Calibres);

Err BMMCalibOnSabr(double forward, double maturity, double sigma, double alpha,
                   double beta, double rho, double pi, double NStd,
                   SrtDiffusionType TypeInput, double *CalibSigma,
                   double *CalibAlpha, double *CalibBeta, double *CalibRho,
                   double *CalibPi, double *Calibres);

Err BMM2CalibOnSabr(double forward, double maturity, double atmvol,
                    double alpha, double beta, double rho, double zeta,
                    double Nstd, SrtDiffusionType TypeInput, double *Calibsigma,
                    double *Calibalpha, double *Calibbeta, double *Calibrho,
                    double *Calibpi, double *Calibres);

Err BMM3CalibOnSabr(double forward, double maturity, double atmvol,
                    double alpha, double beta, double rho, double veta,
                    double pi0, double Nstd, SrtDiffusionType TypeInput,
                    double *Calibsigma, double *Calibalpha, double *Calibbeta,
                    double *Calibrho, double *Calibveta, double *Calibpi0,
                    double *Calibres);

Err BMMGetStates(double forward, double maturity, double sigma, double alpha,
                 double beta, double rho, double pi, double *pSig1,
                 double *pSig2, double *pFwd1, double *pFwd2);

Err BMM2GetStates(double forward, double maturity, double sigma, double alpha,
                  double beta, double rho, double zeta, double *pSig1,
                  double *pSig2, double *pFwd1, double *pFwd2, double *pPi);

Err BMM3GetStates(double forward, double maturity, double sigma, double alpha,
                  double beta, double rho, double veta, double pi0,
                  double *pSig, double *pSig1, double *pSig2, double *pFwd,
                  double *pFwd1, double *pFwd2, double *pPi0, double *pPi);

Err BMMGenGetStates(double forward, double maturity, double sigma, double alpha,
                    double a, double b, double c, double rho, double pi,
                    SrtDiffusionType TypeVolLoc, double *pSig1, double *pSig2,
                    double *pFwd1, double *pFwd2);

Err BetaPrice2(double Forward, double Strike, double Maturity, double Sigma,
               double Beta, double Rho, double *Price);

Err BetaPrice3(double Forward, double Strike, double Sigma, double Maturity,
               double Beta, double *Price);

Err BMMCalibOnPrices(double forward, double maturity, double atmvol,
                     double alpha, double beta, double rho, double pi,
                     int nCalibStrikes, double *CalibStrikes, double *CalibVols,
                     double *CalibWeights, double *CalibSigmabeta,
                     double *CalibBeta, double *CalibAlpha, double *CalibRho,
                     double *CalibPi, double *Calibres);

Err BMM2CalibOnPrices(double forward, double maturity, double atmvol,
                      double alpha, double beta, double rho, double zeta,
                      int nCalibStrikes, double *CalibStrikes,
                      double *CalibVols, double *CalibWeights,
                      double *CalibSigmabeta, double *CalibBeta,
                      double *CalibAlpha, double *CalibRho, double *CalibPi,
                      double *Calibres);

Err BMMGenCalibOnSabrWithStrikesAndFirstGuess(
    double forward, double maturity, double atmvol, double alpha, double beta,
    double a, double b, double c, double rho, double pi,
    SrtDiffusionType TypeVolLoc, int nCalibStrikes, double *CalibStrikes,
    double *CalibVols, double *CalibSigmabeta, double *CalibAlpha,
    double *CalibA, double *CalibB, double *CalibC, double *CalibRho,
    double *CalibPi, double *Calibres);

Err srt_f_optbmmprice(double Forward, double Strike, double Maturity,
                      double Sigma, double Alpha, double Beta, double Rho,
                      double Pi, SrtDiffusionType input,
                      SrtCallPutType call_put, double *price);

Err srt_f_optbmmvol(double Forward, double Strike, double Maturity,
                    double Sigma, double Alpha, double Beta, double Rho,
                    double Pi, SrtDiffusionType input, SrtDiffusionType output,
                    double *vol);

Err srt_f_optbmmcalvol(double Forward, double Strike, double Maturity,
                       double Sigma, double Alpha, double Beta, double Rho,
                       double Pi, double NStd, SrtDiffusionType input,
                       SrtDiffusionType output, double *vol);

Err srt_f_optbmm2vol(double Forward, double Strike, double Maturity,
                     double Sigma, double Alpha, double Beta, double Rho,
                     double Zeta, SrtDiffusionType input,
                     SrtDiffusionType output, double *vol);

Err srt_f_optbmm2calvol(double Forward, double Strike, double Maturity,
                        double Sigma, double Alpha, double Beta, double Rho,
                        double Zeta, double NStd, SrtDiffusionType input,
                        SrtDiffusionType output, double *vol);

Err srt_f_optbmm3vol(double Forward, double Strike, double Maturity,
                     double Sigma, double Alpha, double Beta, double Rho,
                     double veta, double pi0, SrtDiffusionType input,
                     SrtDiffusionType output, double *vol);

Err srt_f_optbmm3calvol(double Forward, double Strike, double Maturity,
                        double Sigma, double Alpha, double Beta, double Rho,
                        double veta, double pi0, double NStd,
                        SrtDiffusionType input, SrtDiffusionType output,
                        double *vol);

Err srt_f_optbmmvolfromstates(double Strike, double Maturity, double Beta,
                              double Fwd1, double Fwd2, double Sig1,
                              double Sig2, double Pi, SrtDiffusionType output,
                              double *vol);

Err srt_f_optbmm3volfromstates(double Forward, double Strike, double Maturity,
                               double Beta, double Fwd, double Fwd1,
                               double Fwd2, double Sig, double Sig1,
                               double Sig2, double Pi0, double Pi,
                               SrtDiffusionType output, double *vol);

double op_bmm_calib(double F, double K, double T,
                    double sigma, // ATMLOG vol
                    double alpha, double beta, double rho, double pi);

double op_BMM2_calib(double F, double K, double T, double sigma, double alpha,
                     double beta, double rho, double zeta);

double op_bmm3_calib(double F, double K, double T, double sigma, double alpha,
                     double beta, double rho, double veta, double pi0);

Err BMMGenPrice(double Forward, double Strike, double SigmaBeta,
                double Maturity, double a, double b, double c, double *Price,
                SrtDiffusionType TypeVolLoc);

double op_bmmgen_calib(double F, double K, double T,
                       double sigma, // ATMLOG vol
                       double alpha, double a, double b, double c, double rho,
                       double pi, SrtDiffusionType TypeVolLoc);

Err srt_f_optbmmgenvol(double Forward, double Strike, double Maturity,
                       double Sigma, double Alpha, double a, double b, double c,
                       double Rho, double Pi, SrtDiffusionType input,
                       SrtDiffusionType output, double *vol,
                       SrtDiffusionType TypeVolLoc);
#endif