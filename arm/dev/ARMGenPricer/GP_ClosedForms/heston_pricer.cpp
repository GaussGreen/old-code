#include "gpclosedforms/heston_pricer.h"

#include <gpbase\cloneutilityfunc.h>

CC_BEGIN_NAMESPACE(ARM)

double ARM_HestonOptionPricer::price()
{
	if(m_shstrike < K_DOUBLE_TOL)
	{
		return m_sens == 1 ? m_forward : 0.;
	}

	double p = m_shfwd - m_shstrike * value() / ARM_NumericConstants::ARM_PI;

	if(m_shift > 0.)
		return m_sens == 1 ? p : p - m_shfwd + m_shstrike;
	else
		return m_sens == 1 ? p - m_shfwd + m_shstrike : p;
}

void ARM_HestonOptionPricer::ComputeCoeff(const complex<double>& u, complex<double>& A, complex<double>& B, 
										  double kappa, double theta, double rho, double nu)
{
	double deltat;
	
	if (m_reset+K_NEW_DOUBLE_TOL<m_times[m_idx])
	{
		if (m_idx>0)
		{
			deltat = m_reset - m_times[m_idx-1];
			m_idx--;
		}
		else
			deltat = m_reset;
	}
	else
		deltat = m_reset - m_times[m_idx];

	ComputeCoeffs(m_idx , deltat, u, A, B, kappa, theta, rho, nu);

	for(int k = m_idx; k >= 0; k--)
	{
		deltat = k == 0 ? m_times[k] : m_times[k] - m_times[k-1];

		ComputeCoeffs(k, deltat , u, A, B, kappa, theta, rho, nu);

		if (k != 0)
			B *= m_levels[k] * m_levels[k] / m_levels[k-1] / m_levels[k-1];
	}
}

void ARM_HestonOptionPricer::ComputeCoeffs(int k, double deltat, const complex<double>& u, complex<double>& A, complex<double>& B, 
										   double kappa, double theta, double rho, double vvol)
{
	if(deltat < K_DOUBLE_TOL) return;

	if(fabs(m_levels[k]) < K_DOUBLE_TOL) return;

	vvol *= m_levels[k];
	theta *= SQR(m_levels[k]);

	complex<double> I(0.,1.), un(1.,0.);
	// beta = kappa - u x i x rho x vvol
	complex<double> beta(kappa + rho * vvol * u.imag(), - rho * vvol * u.real());
	complex<double> ab = u * (u + I);
	ab *= vvol * vvol;
	// D = sqrt(beta^2 + u(u+i)x vvol^2)
	complex<double> d = std::sqrt(beta*beta + ab);
	if(fabs(d.real() + beta.real()) < K_DOUBLE_TOL && fabs(d.imag() + beta.imag()) < K_DOUBLE_TOL) d *= -1.;
	complex<double> ddt(-d.real()*deltat,-d.imag()*deltat);
	complex<double> dexp = std::exp(ddt);
	complex<double> reccB(B.real()*SQR(vvol),B.imag()*SQR(vvol));
	// G = (beta - D - reccB)/(beta + D - reccB)
	complex<double> G = (beta - d - reccB) / (beta + d - reccB);
	// B = (beta - D) / (vvol^2) x (1 - exp(-D x tau)) / (1 - G exp(-D x tau))
	complex<double> addB = (beta - d - reccB)*(un - dexp)/(un - G*dexp);
	addB /= SQR(vvol);
	B += addB;
	// A = (kappa x theta / vvol^2) x ( (beta - D) x tau - 2 log((G exp(-D x tau) - 1)/(G - 1))
	complex<double> betamoinsD = beta - d;
	betamoinsD *= deltat;
	complex<double> addA = std::log((G * dexp - un) / (G - un));
	addA *= -2.;
	addA += betamoinsD;
	addA *= kappa * theta / SQR(vvol);
	A += addA;
}

void ARM_HestonOptionPricerVVolt::ComputeCoeff(const complex<double>& u, complex<double>& A, complex<double>& B,
											   double kappa, double theta, double rho, double nu)
{
	double dt = m_idxNu == 0 ? m_reset : m_reset - m_nut.GetAbscisses()[m_idxNu-1]/365.;

	ComputeCoeffs(0,dt,u,A,B,kappa,theta,rho,m_nut[m_idxNu]);

	for(int k = m_idxNu-1; k >= 0; k--)
	{
		dt = k == 0 ? m_nut.GetAbscisses()[k]/365. : (m_nut.GetAbscisses()[k] - m_nut.GetAbscisses()[k-1])/365.;

		ComputeCoeffs(0,dt,u,A,B,kappa,theta,rho,m_nut[k]);
	}
}

CC_END_NAMESPACE()
