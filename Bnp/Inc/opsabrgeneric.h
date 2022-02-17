#ifndef opsabrgenericH
#define	opsabrgenericH

#include "srt_h_all.h"

typedef double(*FuncVolLocType)(double , double , double , double , int );
void GetLocVolFromDiffusionType(SrtDiffusionType TypeVolLoc, FuncVolLocType* vol_loc);
double	vol_shifted_log	(double x, double a, double b, double c, int type);
double	vol_log_quadra	(double x, double a, double b, double c, int type);
double	vol_quadra		(double x, double a, double b, double c, int type);
double	vol_sabr		(double x, double a, double b, double c, int type);
double	vol_bvm			(double x, double a, double b, double c, int type);
double	vol_bvm2		(double x, double a, double b, double c, int type);
double	vol_bvmh		(double x, double a, double b, double c, int type);
double	vol_bvmh2		(double x, double a, double b, double c, int type);
double	vol_bvmc		(double x, double a, double b, double c, int type);
double	vol_parameters_check_from_SrtDiffusionType(double a, double b, double c, int type, SrtDiffusionType VolLocType);

double op_sabrgen(
					double							F,
					double							K,				
					double							T,
					double							sigma,
					double							alpha,
					double							a,
					double							b,
					double							c,
					double							rho,
					double							(*vol_local)(double x, double a, double b, double c, int type));


double op_sabrgen_calib(	
					double							F,
					double							K,
					double							T,
					double							sigma,
					double							alpha,
					double							a,
					double							b,
					double							c,
					double							rho,
					double							(*vol_local)(double x, double a, double b, double c, int type));

double op_sabr_distri(	
					double							F,
					double							K,
					double							T,
					double							sigma,
					double							alpha,
					double							a,
					double							b,
					double							c,
					double							rho,
					double							perc,
					int								type,	/* 0: density, 1: cumul */
					double							(*vol_local)(double x, double a, double b, double c, int type));

double op_sabr_find_limit(
					double	dist_tgt,
					double	F,
					double	T,
					double	sigma,
					double	alpha,
					double	a,
					double	b,
					double	c,
					double	rho,
					double	perc,
					double	(*vol_local)(double x, double a, double b, double c, int type));


#endif