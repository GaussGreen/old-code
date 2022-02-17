
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/ostringstream.h"
#include "gpbase/stringmanip.h"
#include "gpbase/gpmatrixtriangular.h"
#include "gpbase/utilityport.h"
#include "gpclosedforms/smile_sabr.h"
#include "gpclosedforms/vanilla_normal.h"

#include "gpclosedforms/berm2datesproxy.h"

#include "gpclosedforms/gaussian_integrals.h"

#include "gpnumlib/ran2.h"
#include "gpnumlib/sobol.h"
#include "gpnumlib/gaussiananalytics.h"
#include "gpclosedforms/sabrvanilla.h"
#include "gpcalib/densityfunctors.h"
#include "gpclosedforms/sabrbdiff1.h"

CC_BEGIN_NAMESPACE(ARM)

void ARM_Berm2DatesProxy::price()
{
	int n, size = itsResetDates.size();
	double vol, nextopt, swap;
	double t = (itsResetDates[0] - AsOfDate) / 365.;
	double T = (itsResetDates[1] - AsOfDate) / 365.;
	double dt = T - t;

	// ARM_RandUniform_NRRan2 gen(-156);
	ARM_GP_T_Vector<size_t> nbsimul(1,itsNbSimul);
	ARM_RandomGenerator * gen = new ARM_Sobol(0);
	gen->reset(3,nbsimul,3);

	ARM_GP_Vector x(3);
	ARM_GP_Vector x1(itsNbSimul), x2(itsNbSimul), z(itsNbSimul);
	ARM_GP_Vector u1(itsNbSimul), u2(itsNbSimul);
	
	double rho1v = itsRho[0];
	double rho2v = itsRho[1];
	double rho12 = itsRateCorrel;

	double c1z = rho1v;
	double c11 = sqrt(1. - rho1v*rho1v);
	double c2z = rho2v;
	double c21 = rho12 > 0.999 ? sqrt(1.-rho2v*rho2v) : (rho12 - rho1v*rho2v) / c11;
	double c22 = rho12 > 0.999 ? 0. : sqrt(1.- rho2v*rho2v - c21*c21);

	double w1, w2;
	for(n = 0; n < itsNbSimul; n++)
	{
		gen->draw(x);
		w1 = ARM_GaussianAnalytics::cdfNormal_Inv(x[0]);
		w2 = ARM_GaussianAnalytics::cdfNormal_Inv(x[1]);
		z[n] = ARM_GaussianAnalytics::cdfNormal_Inv(x[2]);

		x1[n] = c1z * z[n] + c11 * w1;
		x2[n] = c2z * z[n] + c21 * w1 + c22 * w2;

		u1[n] = NormalCDF(x1[n]);
		u2[n] = NormalCDF(x2[n]);
	}
	
	delete gen;

	ARM_GP_Vector tx1(itsNbSimul), tx2(itsNbSimul);

	// diffusion du premier taux
	if(itsDiffType == 2 || itsDiffType == 3)
	{
		sabrFwdRate(tx1, u1, itsFwdRates[0], t, itsFwdRatesVol[0], itsRho[0], itsNu[0], 201);
		sabrFwdRate(tx2, u2, itsFwdRates[1], t, itsFwdRatesFwdVol[1], itsRho[1], itsNu[1], 201);
	}

	itsBerm = 0.;
	itsEuropean[0] = itsEuropean[1] = 0.;

	for(n = 0; n < itsNbSimul; n++)
	{
		if(itsDiffType == 0)
		{
			tx1[n] = itsFwdRates[0] + itsFwdRatesVol[0] * x1[n];
			tx2[n] = itsFwdRates[1] + itsFwdRatesFwdVol[1] * x2[n];
		}
		else if(itsDiffType == 1)
		{
			tx1[n] = itsFwdRates[0] * exp(-0.5*itsFwdRatesVol[0]*itsFwdRatesVol[0]*t+itsFwdRatesVol[0]*sqrt(t)*x1[n]);
			tx2[n] = itsFwdRates[1] * exp(-0.5*itsFwdRatesFwdVol[1]*itsFwdRatesFwdVol[1]*t+itsFwdRatesFwdVol[1]*sqrt(t)*x2[n]);
		}

		vol = itsVols[1] * exp(-0.5*itsNu[1]*itsNu[1]*t + itsNu[1]*sqrt(t)*z[n]);

		if(itsDiffType == 3)
		{
			vol = CptSABR_implicit_vol_direct(tx2[n],itsStrike,dt,vol,1.,itsRho[1],itsNu[1],ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL);
		}

		// la prochaine option
		if(itsDiffType == 0)
			nextopt = VanillaOption_N (tx2[n],vol,itsStrike,dt,1) * itsAnnuity[1];
		else
			nextopt = BS(tx2[n],itsStrike,dt,vol,1) * itsAnnuity[1];

		itsEuropean[1] += nextopt;
		itsEuropean[0] += tx1[n] - itsStrike > 0. ? (tx1[n] - itsStrike) * itsAnnuity[0] : 0.;

		swap = (tx1[n] - itsStrike) * itsAnnuity[0];

		itsBerm += swap > nextopt ? swap : nextopt;
	}

	itsBerm /= itsNbSimul;
	itsEuropean[0] /= itsNbSimul;
	itsEuropean[1] /= itsNbSimul;
}

void ARM_Berm2DatesProxy::sabrFwdRate(ARM_GP_Vector& fwdrate, const ARM_GP_Vector& u, double fwd, double t, 
									double alpha, double rho, double nu, int nbPoints)
{
	fwdrate.resize(itsNbSimul);

	if(nbPoints == itsNbSimul)
	{
		for(int k = 0; k < itsNbSimul; k++)
		{
			fwdrate[k] = SABR_smile::inverse_distribution(fwd,NormalCDF(u[k]),t,alpha,1.,rho,nu,ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL,40,fwd/4.,1.5,fwd/2.);
		}

		return;
	}

	ARM_GP_Vector x(nbPoints), y(nbPoints);

	double lbound = -6, ubound = 6;
	double h = (ubound - lbound) / (nbPoints + 1);
	double x0 = lbound, stddev = alpha*sqrt(t);

	ARM_GP_Vector put(nbPoints);
	for(int k = 0; k < nbPoints; k++, x0 += h)
	{		
		x[k] = fwd * exp(-0.5*stddev*stddev + stddev*x0);
		put[k] = x[k] - fwd + SABR_smile::call_option(fwd,x[k],t,alpha,1.,rho,nu,ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL,40,fwd/4.,1.5,fwd/2.);
	}

	double pd = 1. - SABR_smile::digital_call_option(fwd,x[0],t,alpha,1.,rho,nu,ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL,40,fwd/4.,1.5,fwd/2.);

	for(k = 0; k < nbPoints; k++)
	{
		y[k] = k == 0 ? pd : k == nbPoints - 1 ? (put[k] - put[k-1]) / (x[k] - x[k-1]) : (put[k+1] - put[k-1]) / (x[k+1] - x[k-1]);
	}

	for(k = 0; k < itsNbSimul; k++)
	{
		if(u[k] < y[0] + K_DOUBLE_TOL) 
			fwdrate[k] = x[0];
		else if(u[k] > y[nbPoints-1] - K_DOUBLE_TOL)
			fwdrate[k] = x[nbPoints-1];
		else
		{
			int i = 0;
			while(y[i] < u[k]) i++;

			fwdrate[k] = x[i-1] + (u[k] - y[i-1])*(x[i] - x[i-1])/(y[i] - y[i-1]);
		}
	}
}


CC_END_NAMESPACE()
