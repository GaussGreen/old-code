
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/ostringstream.h"
#include "gpbase/stringmanip.h"
#include "gpbase/gpmatrixtriangular.h"
#include "gpbase/utilityport.h"
#include "gpclosedforms/smile_sabr.h"


#include "gpclosedforms/VolBondMinMaxProxy.h"

#include "gpclosedforms/gaussian_integrals.h"

#include "gpnumlib/ran2.h"
#include "gpnumlib/sobol.h"
#include "gpnumlib/gaussiananalytics.h"
#include "gpclosedforms/sabrvanilla.h"
#include "gpcalib/densityfunctors.h"
#include "gpclosedforms/sabrbdiff1.h"

CC_BEGIN_NAMESPACE(ARM)

void ARM_VBMinMaxProxy::price()
{
	int n, k, size = itsResetDates.size();
	double t, T, dt, var, sfvar, efvar, vvar, vol, K;

	// ARM_RandUniform_NRRan2 gen(-156);
	ARM_GP_T_Vector<size_t> nbsimul(1,itsNbSimul);
	ARM_RandomGenerator * gen = new ARM_Sobol(0);
	gen->reset(2,nbsimul,2);

	ARM_GP_Vector x(2);
	ARM_GP_Vector z1(itsNbSimul), z2(itsNbSimul);
	ARM_GP_Vector u(itsNbSimul);

	ARM_GP_Vector smaxfwd(itsNbSimul), sminfwd(itsNbSimul), svol(itsNbSimul), emax(itsNbSimul), emin(itsNbSimul);
	ARM_GP_Vector sfwd(itsNbSimul), efwd(itsNbSimul);
	
	double type2atm, adjust, cpn3;

	for(n = 0; n < itsNbSimul; n++)
	{
		gen->draw(x);
		u[n] = x[0];
		z1[n] = ARM_GaussianAnalytics::cdfNormal_Inv(x[0]);
		z2[n] = ARM_GaussianAnalytics::cdfNormal_Inv(x[1]);
	}
	
	delete gen;

	for(k = 0; k < size - 1; k++)
	{
		t = (itsResetDates[k] - AsOfDate)/365.;
		T = (itsResetDates[k+1] - AsOfDate)/365.;
		dt = T - t;

		sfvar = itsTotalVol[k] * itsTotalVol[k] * t;
		efvar = itsLeftVol[k+1] * itsLeftVol[k+1] * t;

		vvar = itsNu[k] * itsNu[k] * t;


		if(itsSABRDiff)
		{
			sabrFwdRate(sfwd, u, itsFwdRates[k], t, itsTotalVol[k], itsRho[k], itsNu[k],201);
			sabrFwdRate(efwd, u, itsFwdRates[k+1], t, itsLeftVol[k+1], itsRho[k], itsNu[k],201);
		}

		for(n = 0; n < itsNbSimul; n++)
		{
			if(itsSABRDiff == false)
			{
				sfwd[n]	= itsFwdRates[k] * exp(-0.5*sfvar + sqrt(sfvar) * z1[n]);
				efwd[n]	= itsFwdRates[k+1] * exp(-0.5*efvar + sqrt(efvar) * z1[n]);
			}

			smaxfwd[n] = itsMaxChoice == 0 ? sfwd[n] : itsMaxChoice == 1 ? efwd[n] : sfwd[n] > efwd[n] ? sfwd[n] : efwd[n];
			sminfwd[n] = itsMinChoice == 0 ? sfwd[n] : itsMinChoice == 1 ? efwd[n] : sfwd[n] < efwd[n] ? sfwd[n] : efwd[n];

			svol[n] = itsRightVol[k] * exp(-0.5*vvar + sqrt(vvar) * (itsRho[k] * z1[n] + sqrt(1. - itsRho[k]*itsRho[k])*z2[n]));

			if(itsPriceType1)
			{
				K = itsStartLev[k] * sfwd[n] + itsStartAdd[k];

				vol = itsSABRDiff == false ? svol[n] : CptSABR_implicit_vol_direct(efwd[n],K,dt,svol[n],1.,itsRho[k],itsNu[k],ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL);
					
				if(itsType1sens == 1 || itsType1sens == -1)
					itsStdVB1[k] += BS(efwd[n],K,dt,vol,itsType1sens);
				else
					itsStdVB1[k] += BS(efwd[n],K,dt,vol,1)
								  + BS(efwd[n],K,dt,vol,-1);

				if(itsVBLev[k] > K_DOUBLE_TOL)
				{
					double stk	= (itsCapStartLev[k] * sfwd[n] + itsCapStartAdd[k]);
					double put1 = K + stk/itsVBLev[k] < 0. ? 0. : BS(efwd[n],K + stk/itsVBLev[k], dt, vol, -1);
					double put2 = K - stk/itsVBLev[k] < 0. ? 0. : BS(efwd[n],K - stk/itsVBLev[k], dt, vol, -1);
					double put3 = K < 0. ? 0. : BS(efwd[n],K, dt, vol, -1);

					itsCoupon[k] += itsVBLev[k] * (put1 + put2 - 2.*put3);
				}
				else if(itsVBLev[k] < K_DOUBLE_TOL)
				{
					double stk	= - (itsCapStartLev[k] * sfwd[n] + itsCapStartAdd[k]);
					double call	= BS(efwd[n],K + stk/fabs(itsVBLev[k]), dt, vol, 1);
					double put	= BS(efwd[n],K - stk/fabs(itsVBLev[k]), dt, vol, -1);

					itsCoupon[k] += fabs(itsVBLev[k]) * (call + put);
				}
			}

			if(itsPriceType2)
			{
				K = itsStartLev[k] * efwd[n] + itsStartAdd[k];

				vol = itsSABRDiff == false ? svol[n] : CptSABR_implicit_vol_direct(efwd[n],K,dt,svol[n],1.,itsRho[k],itsNu[k],ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL);

				if(itsType2sens == 1 || itsType2sens == -1)
					itsStdVB2[k] += BS(efwd[n],K,dt,vol,itsType2sens);
				else
					itsStdVB2[k] += BS(efwd[n],K,dt,vol,1)
								  + BS(efwd[n],K,dt,vol,-1);
			}

			if(itsPriceType3)
			{

				vol = itsSABRDiff == false ? svol[n] : CptSABR_implicit_vol_direct(efwd[n],efwd[n],dt,svol[n],1.,itsRho[k],itsNu[k],ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL);
				var = vol * vol * dt;

				type2atm	= 2 * BS(efwd[n],efwd[n],dt,vol,1);
				
				emax[n]		= E_Max(smaxfwd[n], var);
				emin[n]		= E_Min(sminfwd[n], var);
				cpn3		= emax[n] - emin[n];

				// ajustement
				adjust			= (1. - (1. - type2atm/cpn3) * sqrt(itsMinMaxFreq));

				itsStdVB3[k]	+= cpn3 * adjust;
				itsMaxRate[k]	+= emax[n] * adjust;
				itsMinRate[k]	+= emin[n] * adjust;
			}
		}

		if(itsPriceType1) itsStdVB1[k] /= itsNbSimul;
		if(itsPriceType2) itsStdVB2[k] /= itsNbSimul;
		if(itsPriceType3)
		{
			itsMaxRate[k]	/= itsNbSimul;
			itsMinRate[k]	/= itsNbSimul;
			itsStdVB3[k]	/= itsNbSimul;
		}
		itsCoupon[k] /= itsNbSimul;
	}
}

void ARM_VBMinMaxProxy::sabrFwdRate(ARM_GP_Vector& fwdrate, const ARM_GP_Vector& u, double fwd, double t, 
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


string ARM_VBMinMaxProxy::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "ARM_VBMinMaxProxy\n";

	return os.str();	
}

/*
void ARM_VBMinMaxProxy::price()
{
	int k, size = itsResetDates.size();
	double t, T, dt, maxfwd, minfwd, var;

	for(k = 0; k < size - 1; k++)
	{
		t = itsResetDates[k];
		T = itsResetDates[k+1];
		dt = (T - t) / 365.;
		var = itsVols[k] * itsVols[k] * dt;

		if(itsFwdRates[k] > itsFwdRates[k+1])
		{
			maxfwd = itsFwdRates[k];
			minfwd = itsFwdRates[k+1];
		}
		else
		{
			maxfwd = itsFwdRates[k+1];
			minfwd = itsFwdRates[k];
		}
		
		itsMaxRate[k] = E_Max(maxfwd, var);
		itsMinRate[k] = E_Min(minfwd, var);

		itsCpn[k] = (itsMaxRate[k] - itsMinRate[k]) * itsIT[k] * itsDF[k];
		
		itsPrice += itsCpn[k];
	}
}

double ARM_VBMinMaxProxy::E_Max(double fwd, double var)
{
	CalcBound func(this, fwd, sqrt(var));
	double bbound;
	double bound = brentSolve(func, 1e-8,fwd,10*fwd,1e-8,100,0,&bbound);

	return Integrale(fwd, var, fwd, bound);
}

double ARM_VBMinMaxProxy::E_Min(double fwd, double var)
{
	CalcBound func(this, fwd, sqrt(var));
	double bbound;
	double bound = brentSolve(func, 1e-8,1e-8,fwd,1e-8,100,0,&bbound);

	return Integrale(fwd, var, bound, fwd);
}

double ARM_VBMinMaxProxy::Integrale(double fwd, double var, double lbound, double ubound)
{
	GaussLegendre_Coefficients c(itsNbLegPts);

	double scale = ubound - lbound, sum = 0., x;
	for(int i = 0; i < itsNbLegPts; i++)
	{
		x = lbound + scale * (0.5 + 0.5 * c.get_point(i));
		sum += density(fwd,var,x) * 0.5 * c.get_weight(i);
	}
	
	sum *= scale;
	
	return sum;
}

string ARM_VBMinMaxProxy::toString(const string& indent, const string& nextIndent) const
{
	return ARM_RootObject::toString(indent,nextIndent);
}
*/
CC_END_NAMESPACE()
