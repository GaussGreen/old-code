
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/ostringstream.h"
#include "gpbase/stringmanip.h"
#include "gpbase/gpmatrixtriangular.h"
#include "gpbase/utilityport.h"
#include "gpclosedforms/smile_sabr.h"


#include "gpclosedforms/SpreadVolBondProxy.h"

#include "gpclosedforms/gaussian_integrals.h"

#include "gpnumlib/ran2.h"
#include "gpnumlib/sobol.h"
#include "gpnumlib/gaussiananalytics.h"
#include "gpclosedforms/sabrvanilla.h"
#include "gpcalib/densityfunctors.h"
#include "gpclosedforms/sabrbdiff1.h"
#include "gpclosedforms/spreadoption_lognormal_formula.h"

CC_BEGIN_NAMESPACE(ARM)

void ARM_SpreadVBProxy::price()
{
	int n, k, size = itsResetDates.size();
	double t, T, dt, var, sfvar, efvar, vvar, vol1, vol2, K, corr;

	// ARM_RandUniform_NRRan2 gen(-156);
	ARM_GP_T_Vector<size_t> nbsimul(1,itsNbSimul);
	ARM_RandomGenerator * gen = new ARM_Sobol(0);
	gen->reset(4,nbsimul,4);

	std::vector<double> x(4);
	std::vector<double> x1(itsNbSimul), x2(itsNbSimul), z(itsNbSimul), z1(itsNbSimul), z2(itsNbSimul), zc(itsNbSimul);
	std::vector<double> u1(itsNbSimul), u2(itsNbSimul);
	std::vector<double> sfwd1(itsNbSimul), sfwd2(itsNbSimul), efwd1(itsNbSimul), efwd2(itsNbSimul);
	std::vector<double> svol1(itsNbSimul), svol2(itsNbSimul);

	double type2atm, adjust, cpn3;

	for(n = 0; n < itsNbSimul; n++)
	{
		gen->draw(x);
		z1[n] = ARM_GaussianAnalytics::cdfNormal_Inv(x[0]);
		z2[n] = ARM_GaussianAnalytics::cdfNormal_Inv(x[1]);
		z[n] = ARM_GaussianAnalytics::cdfNormal_Inv(x[2]);
		zc[n] = ARM_GaussianAnalytics::cdfNormal_Inv(x[3]);
	}
	
	delete gen;
	int nbLeg = 12;

	double rho12, rho1v, rho2v;
	double c1z, c11, c2z, c21, c22;
	
	std::vector<double> spreadopt(itsNbSimul);

	for(k = 0; k < size - 1; k++)
	{
		t = (itsResetDates[k] - AsOfDate)/365.;
		T = (itsResetDates[k+1] - AsOfDate)/365.;
		dt = T - t;
			
		rho12 = itsCorrel12[k+1];
		rho1v = itsRho1[k];
		rho2v = itsRho2[k];
		
		c1z = rho1v;
		c11 = sqrt(1. - rho1v*rho1v);
		c2z = rho2v;
		c21 = rho12 > 0.999 ? sqrt(1.-rho2v*rho2v) : (rho12 - rho1v*rho2v) / c11;
		c22 = rho12 > 0.999 ? 0. : sqrt(1.- rho2v*rho2v - c21*c21);

		for(n = 0; n < itsNbSimul; n++)
		{
			x1[n] = c1z * z[n] + c11 * z1[n];
			x2[n] = c2z * z[n] + c21 * z1[n] + c22 * z2[n];

			u1[n] = NormalCDF(x1[n]);
			u2[n] = NormalCDF(x2[n]);
		}

		if(itsSABRDiff)
		{
			sabrFwdRate(sfwd1, u1, itsFwdRates1[k], t, itsFwdRatesVol1[k], itsRho1[k], itsNu1[k],201);
			sabrFwdRate(efwd1, u1, itsFwdRates1[k+1], t, itsFwdRatesFwdVol1[k+1], itsRho1[k], itsNu1[k],201);

			sabrFwdRate(sfwd2, u2, itsFwdRates2[k], t, itsFwdRatesVol2[k], itsRho2[k], itsNu2[k],201);
			sabrFwdRate(efwd2, u2, itsFwdRates2[k+1], t, itsFwdRatesFwdVol2[k+1], itsRho2[k], itsNu2[k],201);
		}

		for(n = 0; n < itsNbSimul; n++)
		{
			if(itsSABRDiff == false)
			{
				sfvar		= itsFwdRatesVol1[k] * itsFwdRatesVol1[k] * t;
				efvar		= itsFwdRatesFwdVol1[k+1] * itsFwdRatesFwdVol1[k+1] * t;

				sfwd1[n]	= itsFwdRates1[k] * exp(-0.5*sfvar + sqrt(sfvar) * x1[n]);
				efwd1[n]	= itsFwdRates1[k+1] * exp(-0.5*efvar + sqrt(efvar) * x1[n]);

				sfvar		= itsFwdRatesVol2[k] * itsFwdRatesVol2[k] * t;
				efvar		= itsFwdRatesFwdVol2[k+1] * itsFwdRatesFwdVol2[k+1] * t;

				sfwd2[n]	= itsFwdRates2[k] * exp(-0.5*sfvar + sqrt(sfvar) * x2[n]);
				efwd2[n]	= itsFwdRates2[k+1] * exp(-0.5*efvar + sqrt(efvar) * x2[n]);
			}

			svol1[n] = itsVols1[k] * exp(-0.5*itsNu1[k]*itsNu1[k]*t + itsNu1[k]*sqrt(t) * z[n]);
			svol2[n] = itsVols2[k] * exp(-0.5*itsNu2[k]*itsNu2[k]*t + itsNu2[k]*sqrt(t) * z[n]);

			corr = zc[n] * itsCorrVol * itsCorrVol * (fabs(itsCorrMeanRev) < K_DOUBLE_TOL ? t : (1. - exp(-2.*itsCorrMeanRev*t)) / (2.*itsCorrMeanRev));

			corr = fabs(corr) < K_DOUBLE_TOL ? itsCorrel12[k+1] : tanh(atanh(itsCorrel12[k+1]) + corr);
			
			if(itsPriceType1)
			{
				K = itsfwdLeverage[k] * (sfwd1[n] - sfwd2[n]) + itsfwdStrikes[k];

				vol1 = itsSABRDiff == false ? svol1[n] : CptSABR_implicit_vol_direct(efwd1[n],efwd2[n] + K,dt,svol1[n],1.,itsRho1[k],itsNu1[k],ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL);
				vol2 = itsSABRDiff == false ? svol2[n] : CptSABR_implicit_vol_direct(efwd2[n],efwd1[n] - K,dt,svol2[n],1.,itsRho2[k],itsNu2[k],ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL);

				if(itsType1sens == 1 || itsType1sens == -1)
				{
					itsStdVB1[k] += Export_LogNormal_SpreadOption(efwd1[n],efwd2[n],vol1,vol2,corr,K,dt,itsType1sens,ARM_CF_SpreadDigitalOption_Formula::SPREADOPTION,nbLeg);
				}
				else
				{
					itsStdVB1[k] += Export_LogNormal_SpreadOption(efwd1[n],efwd2[n],vol1,vol2,corr,K,dt,1,ARM_CF_SpreadDigitalOption_Formula::SPREADOPTION,nbLeg)
								  + Export_LogNormal_SpreadOption(efwd1[n],efwd2[n],vol1,vol2,corr,K,dt,-1,ARM_CF_SpreadDigitalOption_Formula::SPREADOPTION,nbLeg);

					if(fabs(itsLevier[k]) > K_DOUBLE_TOL)
					{
						double put1 = Export_LogNormal_SpreadOption(efwd1[n],efwd2[n],vol1,vol2,corr,itsFix[k]/itsLevier[k]+K,dt,-1,ARM_CF_SpreadDigitalOption_Formula::SPREADOPTION,nbLeg);
						double put2 = Export_LogNormal_SpreadOption(efwd1[n],efwd2[n],vol1,vol2,corr,-itsFix[k]/itsLevier[k]+K,dt,-1,ARM_CF_SpreadDigitalOption_Formula::SPREADOPTION,nbLeg);
						double put3 = Export_LogNormal_SpreadOption(efwd1[n],efwd2[n],vol1,vol2,corr,K,dt,-1,ARM_CF_SpreadDigitalOption_Formula::SPREADOPTION,nbLeg);

						itsCoupon[k] += itsLevier[k] * (put1 + put2 - 2 * put3);
					}
				}
			}

			if(itsPriceType2)
			{
				K = itsfwdLeverage[k] * (efwd1[n] - efwd2[n]) + itsfwdStrikes[k];

				vol1 = itsSABRDiff == false ? svol1[n] : CptSABR_implicit_vol_direct(efwd1[n],efwd2[n] + K,dt,svol1[n],1.,itsRho1[k],itsNu1[k],ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL);
				vol2 = itsSABRDiff == false ? svol2[n] : CptSABR_implicit_vol_direct(efwd2[n],efwd1[n] - K,dt,svol2[n],1.,itsRho2[k],itsNu2[k],ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL);

				if(itsType2sens == 1 || itsType2sens == -1)
				{
					spreadopt[n] = Export_LogNormal_SpreadOption(efwd1[n],efwd2[n],vol1,vol2,corr,K,dt,itsType2sens,ARM_CF_SpreadDigitalOption_Formula::SPREADOPTION,nbLeg);
					itsStdVB2[k] += spreadopt[n];
				}
				else
					itsStdVB2[k] += Export_LogNormal_SpreadOption(efwd1[n],efwd2[n],vol1,vol2,corr,K,dt,1,ARM_CF_SpreadDigitalOption_Formula::SPREADOPTION,nbLeg)
								  + Export_LogNormal_SpreadOption(efwd1[n],efwd2[n],vol1,vol2,corr,K,dt,-1,ARM_CF_SpreadDigitalOption_Formula::SPREADOPTION,nbLeg);
			}
		}

		if(itsPriceType1) itsStdVB1[k] /= itsNbSimul;
		if(itsPriceType2) itsStdVB2[k] /= itsNbSimul;
		itsCoupon[k] /= itsNbSimul;
	}
}

void ARM_SpreadVBProxy::sabrFwdRate(std::vector<double>* fwdrate, const std::vector<double>* u, double fwd, double t, 
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

	std::vector<double> x(nbPoints), y(nbPoints);

	double lbound = -6, ubound = 6;
	double h = (ubound - lbound) / (nbPoints + 1);
	double x0 = lbound, stddev = alpha*sqrt(t);

	std::vector<double> put(nbPoints);
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


string ARM_SpreadVBProxy::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "ARM_SpreadVBProxy\n";

	return os.str();	
}

CC_END_NAMESPACE()
