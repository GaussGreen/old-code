#ifndef __VOLBONDCALIB_H
#define __VOLBONDCALIB_H


Err VolBondCalibration(
	char *yc_name,						/*	Name of the yield curve */
	char *vol_curve_name,				/*	Name of the market vol curve */
	char *ref_rate_name,					/*	Name of the reference rate */
	long startDate,
	long longEndDate,
	long shortEndDate,
	char *swaption_freq,				/*	Frequency and basis of swaption */
	char *swaption_basis,
	double alpha,				/*	Alpha, Gamma, Rho */
	double gamma,
	double rho,
	double *lambda,						/*	Result : Calibrated Lambda and Sigma */
	double *sig);

Err VolBondCalibrationOld(
	char *yc_name,						/*	Name of the yield curve */
	char *vol_curve_name,				/*	Name of the market vol curve */
	char *ref_rate_name,					/*	Name of the reference rate */
	long startDate,
	long longEndDate,
	long shortEndDate,
	char *swaption_freq,				/*	Frequency and basis of swaption */
	char *swaption_basis,
	double alpha,				/*	Alpha, Gamma, Rho */
	double gamma,
	double rho,
	double *lambda,						/*	Result : Calibrated Lambda and Sigma */
	double *sig);

#endif //__VOLBONDCALIB_H