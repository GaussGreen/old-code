#ifndef __BVM_PERTURBATION_H
#define __BVM_PERTURBATION_H

double bvm_local_vol(double x, double a, double b, double c, int log_or_norm);

double BVM_First_Order_Perturbation_Price(double today,
						double Fwd,
						double Vol,
						double Strike,
						double Maturity,
						double gamma,
						double alpha,
						double rho,
						double VolInfinity,
						double Lambda,
						double (*local_vol)(double x, double a, double b, double c, int log_or_norm),
						int	NQuad_Her,
						double *w_her,
						double *x_her,
						int	NQuad_Leg,
						double *w_leg,
						double *x_leg,
						int LOG_NORM);


double BVM_Second_Order_Perturbation_Price(double today,
						double Fwd,
						double Vol,
						double Strike,
						double Maturity,
						double gamma,
						double alpha,
						double rho,
						double VolInfinity,
						double Lambda,
						double (*local_vol)(double x, double a, double b, double c, int log_or_norm),
						int	NQuad_Her,
						double *w_her,
						double *x_her,
						int	NQuad_Leg,
						double *w_leg,
						double *x_leg);


/*
double BVM_Perturbation_Price2(double today,
						double Fwd,
						double Vol,
						double Strike,
						double Maturity,
						double alpha,
						double rho,
						int	NQuad_Her,
						double *w_her,
						double *x_her,
						int	NQuad_Leg,
						double *w_leg,
						double *x_leg);
*/

Err sabrBVMmontecarlo(
					double forw,
					double vovol,
					double beta,
					double rho,
					double num_paths,
					double num_steps,
					double sigma,
					double maturity,
					double *strike,
					int    num_strikes,
					int modeltype, 
					int sampletype,/*0: randsam, 1: abs, 2:sobol, 3: SPECTRUNC*/
					double *impvol
					);

Err sabrBVMmontecarlo2(	
				double							forward,
				double							*strike,
				int								nb_strike,
				double							maturity,
				double							sigma_beta,
				double							alpha,
				double							beta,
				double							rho,
				double							lambda,
				int								npaths,
				int								nsteps,
				int								do_balsam,
				double							**res);

Err sabrBVMmontecarlo3(
					double forw,
					double vovol,
					double beta,
					double rho,
					double lambda,
					double num_paths,
					double num_steps,
					double sigma,
					double maturity,
					double *strike,
					int    num_strikes,
					int modeltype, 
					int sampletype,/*0: randsam, 1: abs, 2:sobol, 3: SPECTRUNC*/
					double *impvol
					);


double SABRNormalVol(double forward, double strike, double maturity,
				 double betavol, double vovol, double rho, double beta);

double BVMLOGVol(double forward, double strike, double maturity,
				 double betavol, double vovol, double rho, double eta);

Err BVMApproxPrice(double Xinit, 
				 double XStar, double GStar, 
				 double Sigma, double Lambda, 
				 double Maturity, 
				 long NHermite,
				 long NStrikes,
				 double *Strikes,
				 double	*ImpliedVol,
				 double	*Fwd);

Err BVMApproxPriceFromFwd(double Fwd,
				   double XStar, double GStar, 
				   double Sigma, double Lambda, 
				   double Maturity, 
				   long NHermite,
				   long NStrikes,
				   double *Strikes,
				   double	*ImpliedVol,
				   double *FwdOutPut);

Err BVMApproxAccurate(double Xinit,
					  double Strike,
					  double Maturity,
					  double XStar, double GStar,
					  double sigma,
					  double lambda,
					  long Order,
					  double *Price);

Err BVMApproxAccurate2(double Xinit,
					  double Strike,
					  double Maturity,
					  double XStar, double GStar,
					  double sigma,
					  double lambda,
					  long NQuadrature,
					  double *Price);

double g_function(double x, double xStar, double gStar, double lambda, int nderiv);

#endif



