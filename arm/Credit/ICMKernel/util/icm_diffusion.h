#ifndef _ICM_DIFFUSION_H_
#define _ICM_DIFFUSION_H_

#include <vector>
#include <map>
#include <math.h>
#include "ICMKernel/util./icm_macro.h" 
#include "ARMKernel/util/refvalue.h"
// 17783 using namespace std;

  
// dXt = -a.Xt.dt + sigmaTaux.dZt
// rt = Xt + b
// Donc drt = a*(b-rt)dt + sigmaTaux.dZt;
double TauxVasicek(double t,
				   double T,
				   double a,
				   double sigma,
				   double IRVol,
				   vector<double> ZCCurve);

double ZC(double time, vector<double> ZCCurve);

double Forward(double T,
			   double dt,
			   vector<double> ZCCurve);

double Diff(double Vinit,
			double Mean,
			double SpeedReversion,
			double Vol,
			double dt,
			double Norm01);

#endif // _ICM_DIFFUSION_H_