#ifndef OPSABRGENERICCALIB_H
#define OPSABRGENERICCALIB_H

typedef struct {
  int nb_strikes;
  double forward;
  double maturity;
  double *Vec_strikes;
  double *Vec_marketvols;
  double *Vec_weights;
  double Beta;
  double a;
  double b;
  double c;
  double Veta;
  double Pi;
  double Pi0;
  double VegaFact;
  int ExitFlag; // flag for the optimizer of BMM
  SrtDiffusionType TypeVolLoc;
  FuncVolLocType vol_loc;
} Opt_params;

// Returns BS(@K) see .c for I/O
Err srt_f_optsabrgenvol(double Forward, double Strike, double Maturity,
                        double Sigma, double Alpha, double a, double b,
                        double c, double Rho, SrtDiffusionType input,
                        SrtDiffusionType output, SrtDiffusionType voLloc,
                        double *vol);

/* Calibrates the generic SABR parameters to a market smile */
Err opsabrgencalib(double Fwd, double maturity, int nbr_strikes,
                   double *strikes, double *marketvols, double *ATMVol,
                   double *alpha, int freeze_alpha, double *a, int freeze_a,
                   double *b, int freeze_b, double *c, int freeze_c,
                   double *rho, int freeze_rho, SrtDiffusionType VolLoc,
                   double *fitting_error);

/* Calibrates the generic SABR parameters to a market smile */
Err opsabrgencalibpro(double forward, double maturity, int nbr_strikes,
                      double *strikes, double *market_vols, double *InitGuess,
                      double *lbounds, double *ubounds, SrtDiffusionType VolLoc,
                      double *fitting_error);
#endif
