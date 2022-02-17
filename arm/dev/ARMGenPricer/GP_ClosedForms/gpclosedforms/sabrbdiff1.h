/*----------------------------------------------------------------------------*

    FILE: sabrbdiff1.h   

    Description: Formulaes of SABR Model in case Beta<>1

 *----------------------------------------------------------------------------*/


#ifndef _ARM_CF_SABRDIFF1_H
#define _ARM_CF_SABRDIFF1_H

/*---- Some definitions ----*/


/*---- For SABR ----*/

#ifndef K_LD
#define K_LD          0
#endif

#ifndef K_SABR_GEO
#define K_SABR_GEO    1
#endif 

#ifndef K_SABR_ARITH
#define K_SABR_ARITH  2
#endif

/* case Beta < 1 */

#ifndef K_SABR_IMPLNVOL
#define K_SABR_IMPLNVOL 3
#endif

#ifndef K_SABR_WEIGHT
#define K_SABR_WEIGHT   4
#endif

#ifndef K_SABR_IMPLNVOL2
#define K_SABR_IMPLNVOL2 5
#endif

#ifndef K_SABR_ANALYTIC
#define K_SABR_ANALYTIC  6
#endif

#ifndef K_SABR_NINTEGRATION
#define K_SABR_NINTEGRATION 7
#endif

#ifndef K_SABR_AUTOMATIC
#define K_SABR_AUTOMATIC 8
#endif

#ifndef K_SABR_DIRECTEXACT
#define K_SABR_DIRECTEXACT 9
#endif

#ifndef K_SABR_DIRECTEXACTSTRIKE
#define K_SABR_DIRECTEXACTSTRIKE 10
#endif

#ifndef K_SABR_DIRECTGEOMETRIC
#define K_SABR_DIRECTGEOMETRIC 11
#endif

#ifndef K_SABR_DIRECTARITHMETIC
#define K_SABR_DIRECTARITHMETIC 12
#endif

#ifndef K_SABR_NORMALEXACT
#define K_SABR_NORMALEXACT 13
#endif

#ifndef K_SABR_NORMALGEOMETRIC
#define K_SABR_NORMALGEOMETRIC 14
#endif

#ifndef K_SABR_NORMALARITHMETIC
#define K_SABR_NORMALARITHMETIC 15
#endif

#ifndef K_SABR_ANALYTICZP0
#define K_SABR_ANALYTICZP0 16
#endif

#ifndef K_SABR_ANALYTICZP2
#define K_SABR_ANALYTICZP2 17
#endif

#ifndef K_SABR_DIRECT_SC1
#define K_SABR_DIRECT_SC1 18
#endif

#ifndef K_SABR_DIRECT_SC2
#define K_SABR_DIRECT_SC2 19
#endif


extern double CptSABR_BetEqOne_ImplVol(double f,double K,double tex, double alpha,double rho, double nu);

extern double CptSABR_implicit_vol_direct_Series(double f,double K, double T,
                                                 double alpha, double beta,
                                                 double rho, double nu,
                                                 int flag);


extern double CptSABR_implicit_vol_direct(double f, double K, double matu,
                                          double alpha, double beta,
                                          double rho, double nu,
                                          int flag);

extern double CptSABR_implicit_vol_normal_Series(double f, double K, double matu,
                                                 double alpha, double beta,
                                                 double rho, double nu,
                                                 int flag);

extern double CptSABR_implicit_vol_normal(double f, double K, double matu,
                                          double alpha, double beta, double rho,
                                          double nu,
                                          int flag);

extern double CptSABR_implicit_vol_direct_DerAlpha(double f, double K, double matu,
                                                   double alpha, double beta,
                                                   double rho, double nu,
                                                   int flag);

extern double CptSABR_implicit_vol_normal_DerAlpha(double f, double K, double matu,
                                            double alpha, double beta, double rho,
                                            double nu,
                                            int flag);

extern double PhiAlpha(double f, double K,
                       double maturity,
                       double sigmaATM,
                       double alpha, double rho, double nu, double beta,
                       double weight, int type);

extern double DerivPhiAlpha(double f, double K,
                            double mat,
                            double alpha, double rho, double nu, double beta,
                            double weight, int type,
                            double sigATM = 0.0);


extern double ComputeAlpha(double f, double K,
                           double mat,
                           double sigmaATM,
                           double rho, double nu, double beta,
                           double weight, int type);

extern double ARM_CptSABRVolBetaDiffOneNormalExactSABR_GEO(double maturity,
                                                    double volOrAlpha,
                                                    double f,
                                                    double K,
                                                    double rho,
                                                    double nu,
                                                    double beta,
                                                    int sigmaOrAlphaFLAG);

extern double ARM_CptSABRVolBetaDiffOneDirectExactSABR_IMPLNVOL(double maturity,
                                                         double volOrAlpha,
                                                         double f,
                                                         double K,
                                                         double rho,
                                                         double nu,
                                                         double beta,
                                                         int sigmaOrAlphaFLAG);

extern double ComputeAlphaBetaEqOne(double f, double K,
                                    double maturity, 
                                    double sigmaATM,
                                    double Rho_t, double Nu_t,
                                    int type);


/*---------------------------------------------------------------*/
/*         Principal Function                                    */
/*---------------------------------------------------------------*/

extern double ARM_CalcSABRVolBetaDiffOne(double maturity,
                                         double volOrAlpha,
                                         double f,
                                         double K,
                                         double rho,
                                         double nu,
                                         double beta,
                                         int sigmaOrAlphaFLAG,
                                         int SABR_TYPE,
                                         double SABRWeight = 0.0);


#endif

/*----------------------------------------------------------------------------*/
/*---- End Of File ----*/

