
#ifndef _GP_CF_HESTONPRICER_H
#define _GP_CF_HESTONPRICER_H


#include "firsttoinc.h"
#include "gpbase/port.h"

#include "gpnumlib/levmarq.h"
#include "gpnumlib/solver.h"
#include "gpclosedforms/inverse.h"

#include <cmath>
#include <complex>

#include "gpclosedforms/basic_distributions.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/oscillatory_integrals.h"
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/vanilla_normal.h"
#include "gpbase/numericconstant.h"
#include "gpbase/curve.h"

using std::sqrt;
using std::log;
using std::real;
using std::imag;
using std::exp;
using std::complex;

CC_BEGIN_NAMESPACE(ARM)


///////////////////////////////////////////////////////////////////////////////
///
///			Process :
///			dS= S(level*V^(1/2) dW1
///			dV=kappa*(theta-V)dt +vvol*V^(1/2) dW2 
///			dW1.dW2=rho*dt
///
///			Jump J : probability lambda, volatility sigmaJ, log-size muJ
///
/////////////////////////////////////////////////////////////////////////////:
class ARM_HestonOptionPricer : public OscillatoryIntegral
{
protected:

	double				m_reset;
	double				m_forward;
	double				m_strike;
	std::vector<double>		m_strikes;
	int					m_sens;

	double				m_v0;
	double				m_kappa;
	double				m_theta;
	double				m_rho;
	double				m_nu;
	double				m_shift;

	std::vector<double>		m_times;
	std::vector<double>		m_levels;

	bool				m_isDigital;

	double				m_shfwd;
	double				m_shstrike;
	complex<double>		m_logFK;

	std::vector<double>		m_shstrikes;
	vector<complex<double> > m_logFKs;

	int					m_idx;

public:

	ARM_HestonOptionPricer(double reset, 
						   double forward, 
						   double strike, 
						   int sens, 
						   double v0, 
						   double kappa, 
						   double theta, 
						   double rho, 
						   double nu, 
						   double shift,
						   const std::vector<double>& times, 
						   const std::vector<double>& levels,
						   bool isDigital = false)
	:
		m_reset(reset), m_forward(forward), m_strike(strike), m_sens(sens),
	m_v0(v0), m_kappa(kappa), m_theta(theta), m_rho(rho), m_nu(nu), m_shift(shift),
	m_times(times), m_levels(levels), m_isDigital(isDigital), 
	OscillatoryIntegral(0, 80)
	{
		m_shfwd	= m_forward / fabs(m_shift);
		m_shstrike = m_shift < 0. ? m_shfwd * (1. + fabs(m_shift)) - m_strike : m_strike + (1. - m_shift) * m_shfwd;

		m_logFK = std::log(complex<double>(m_shfwd,0.0)/complex<double>(m_shstrike,0.0));

		if(times.size() == 1 || times.size() == 0) 
			m_idx = 0;
		else
		{
			if(m_times.size() != m_levels.size())
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_HestonOptionPricer::times and levels have different size");
			}

			m_idx = 0;

			while(m_idx < times.size() && times[m_idx] < m_reset) m_idx ++;
			if (m_idx == times.size())
				m_idx--;
		}
		m_strikes.push_back(strike);
	}

	ARM_HestonOptionPricer(double reset, 
						   double forward, 
						   double strike, 
						   int sens, 
						   double v0, 
						   double kappa, 
						   double theta, 
						   double rho, 
						   double nu, 
						   double shift,
						   double level,
						   bool isDigital = false,
						   bool isManual = false,
						   double bound = 0.0) :
	m_reset(reset), m_forward(forward), m_strike(strike), m_sens(sens),
	m_v0(v0), m_kappa(kappa), m_theta(theta), m_rho(rho), m_nu(nu), m_shift(shift), m_isDigital(isDigital),
	OscillatoryIntegral(0, 80)
	{
		m_shfwd	= m_forward / fabs(m_shift);
		m_shstrike = m_shift < 0. ? m_shfwd * (1. + fabs(m_shift)) - m_strike : m_strike + (1. - m_shift) * m_shfwd;

		m_logFK = std::log(complex<double>(m_shfwd,0.0)/complex<double>(m_shstrike,0.0));

		m_idx = 0;

		m_levels.resize(1, level);
		m_strikes.push_back(strike);
	}

public:

	virtual int getNbIntegrals() {return m_strikes.size();};
	double	price();
	void prices(std::vector<double>& prices);
	double  proba(const std::vector<double>& strikes);

	double oscillatorycomponant(double k)	{return 0.;};
	double mapper(double k)					{return k;};
	double inversemapper(double S)			{return S;};
	double mapperjacobian(double k)			{return 1;};

	void setStrikes(const std::vector<double>& strikes);

	virtual double integrand(double x)
	{
		double value;
		if (m_isDigital)
		{
			complex<double> A(0.,0.), B(0.,0.), u(x,0);

			ComputeCoeff(u, A, B, m_kappa, m_theta, m_rho, m_nu);

			double real = -x * m_logFK.imag() + A.real() + B.real() * m_v0 * SQR(m_levels[0]);
			double im	= x * m_logFK.real() + A.imag() + B.imag() *  m_v0 *  SQR(m_levels[0]);

			value = exp(real) * sin(im) / (x);
		}
		else
		{
			complex<double> A(0.,0.), B(0.,0.), u(x, - 0.5);

			ComputeCoeff(u, A, B, m_kappa, m_theta, m_rho, m_nu);

			double real = 0.5 * m_logFK.real() - x * m_logFK.imag() + A.real() + B.real() * m_v0 * SQR(m_levels[0]);
			double im	= x * m_logFK.real() + 0.5 * m_logFK.imag()  + A.imag() + B.imag() *  m_v0 *  SQR(m_levels[0]);

			value = exp(real) * cos(im) / (x * x + 0.25);
		}

		return value;
	}

	virtual void integrands(double x,std::vector<double>& values)
	{
		if (m_isDigital)
		{
			complex<double> A(0.,0.), B(0.,0.), u(x,0);
			ComputeCoeff(u, A, B, m_kappa, m_theta, m_rho, m_nu);

			int i;
			double real,im;
			for (i = 0; i < m_strikes.size(); ++i)
			{
				real = -x * m_logFKs[i].imag() + A.real() + B.real() * m_v0 * SQR(m_levels[0]);
				im	= x * m_logFKs[i].real() + A.imag() + B.imag() *  m_v0 *  SQR(m_levels[0]);

				values[i] = exp(real) * sin(im) / (x);
			}
		}
		else
		{
			complex<double> A(0.,0.), B(0.,0.), u(x, - 0.5);

			ComputeCoeff(u, A, B, m_kappa, m_theta, m_rho, m_nu);

			int i;
			double real,im;
			for (i = 0; i < m_strikes.size(); ++i)
			{
				real = 0.5 * m_logFKs[i].real() - x * m_logFKs[i].imag() + A.real() + B.real() * m_v0 * SQR(m_levels[0]);
				im	= x * m_logFKs[i].real() + 0.5 * m_logFKs[i].imag() + A.imag() + B.imag() *  m_v0 *  SQR(m_levels[0]);

				values[i] = exp(real) * cos(im) / (x * x + 0.25);
			}
		}
	}

protected:

	virtual void ComputeCoeff(const complex<double>& u, complex<double>& A, complex<double>& B, double kappa, double theta, double rho, double nu);

	void ComputeCoeffs(int k, double deltat, const complex<double>& u, complex<double>& A, complex<double>& B, double kappa, double theta, double rho, double nu);
};


///////////////////////////////////////////////////////////////////////////////
///
///			Process :
///			dS= S(level*V^(1/2) dW1
///			dV=kappa*(theta-V)dt +vvol*V^(1/2) dW2 
///			dW1.dW2=rho*dt
///
///			
///
/////////////////////////////////////////////////////////////////////////////:
class ARM_HestonFFTPricer : public ARM_HestonOptionPricer
{
private:
	double itsParityCoef;

public:

	ARM_HestonFFTPricer(double reset, 
						   double forward, 
						   double strike, 
						   int sens, 
						   double v0, 
						   double kappa, 
						   double theta, 
						   double rho, 
						   double nu, 
						   double beta,
						   const std::vector<double>& times, 
						   const std::vector<double>& levels,
						   bool isDigital = false,
						   double parityCoef = 1.0):
	ARM_HestonOptionPricer( reset, forward, 
						    strike, 
						    sens, 
						    v0, 
						    kappa, 
						    theta, 
						    rho, 
						    nu, 
						    beta,
						    times, 
						    levels,
						    isDigital),
							itsParityCoef(parityCoef)
	{
		if(fabs(beta) < 0.01) 
			beta = 0.01;
		double tmp = (1.0- beta)*itsParityCoef/beta;

		m_shstrike = m_strike + tmp;
		m_shfwd = m_forward + tmp;

		m_shfwd	= m_forward / fabs(m_shift);
		m_shstrike = m_shift < 0. ? m_shfwd * (1. + fabs(m_shift)) - m_strike : m_strike + (1. - m_shift) * m_shfwd;

		m_logFK = std::log(complex<double>(m_shfwd,0.0)/complex<double>(m_shstrike,0.0));

		if(times.size() == 1 || times.size() == 0) 
			m_idx = 0;
		else
		{
			if(m_times.size() != m_levels.size())
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_HestonOptionPricer::times and levels have different size");
			}

			m_idx = 0;

			while(m_idx < times.size() && times[m_idx] < m_reset) m_idx ++;
			if (m_idx == times.size())
				m_idx--;
		}
		m_strikes.push_back(strike);
	}
	/*ARM_HestonFFTPricer(double reset, 
						   double forward, 
						   double strike, 
						   int sens, 
						   double v0, 
						   double kappa, 
						   double theta, 
						   double rho, 
						   double nu, 
						   double shift,
						   double level,
						   bool isDigital = false,
						   bool isManual = false,
						   double bound = 0.0) :
	m_reset(reset), m_forward(forward), m_strike(strike), m_sens(sens),
	m_v0(v0), m_kappa(kappa), m_theta(theta), m_rho(rho), m_nu(nu), m_shift(shift), m_isDigital(isDigital),
	OscillatoryIntegral(0, 80)
	{
		m_shfwd	= m_forward / fabs(m_shift);
		m_shstrike = m_shift < 0. ? m_shfwd * (1. + fabs(m_shift)) - m_strike : m_strike + (1. - m_shift) * m_shfwd;

		m_logFK = std::log(complex<double>(m_shfwd,0.0)/complex<double>(m_shstrike,0.0));

		m_idx = 0;

		m_levels.resize(1, level);
		m_strikes.push_back(strike);
	}*/

public:

	virtual int getNbIntegrals() {return m_strikes.size();};
	double	price();
	void prices(std::vector<double>* prices);
	double  proba(std::vector<double> strikes);

	double oscillatorycomponant(double k)	{return 0.;};
	double mapper(double k)					{return k;};
	double inversemapper(double S)			{return S;};
	double mapperjacobian(double k)			{return 1;};

	void setStrikes(const std::vector<double>* strikes);

	virtual double integrand(double x)
	{
		double value;
		if (m_isDigital)
		{
			complex<double> A(0.,0.), B(0.,0.), u(x,0);

			ComputeCoeff(u, A, B, m_kappa, m_theta, m_rho, m_nu);

			double real = -x * m_logFK.imag() + A.real() + B.real() * m_v0 * SQR(m_levels[0]);
			double im	= x * m_logFK.real() + A.imag() + B.imag() *  m_v0 *  SQR(m_levels[0]);

			value = exp(real) * sin(im) / (x);
		}
		else
		{
			complex<double> A(0.,0.), B(0.,0.), u(x, - 0.5);

			ComputeCoeff(u, A, B, m_kappa, m_theta, m_rho, m_nu);

			double real = 0.5 * m_logFK.real() - x * m_logFK.imag() + A.real() + B.real() * m_v0 * SQR(m_levels[0]);
			double im	= x * m_logFK.real() + 0.5 * m_logFK.imag()  + A.imag() + B.imag() *  m_v0 *  SQR(m_levels[0]);

			value = exp(real) * cos(im) / (x * x + 0.25);
		}

		return value;
	}

	virtual void integrands(double x,std::vector<double>& values)
	{
		if (m_isDigital)
		{
			complex<double> A(0.,0.), B(0.,0.), u(x,0);
			ComputeCoeff(u, A, B, m_kappa, m_theta, m_rho, m_nu);

			int i;
			double real,im;
			for (i = 0; i < m_strikes.size(); ++i)
			{
				real = -x * m_logFKs[i].imag() + A.real() + B.real() * m_v0 * SQR(m_levels[0]);
				im	= x * m_logFKs[i].real() + A.imag() + B.imag() *  m_v0 *  SQR(m_levels[0]);

				values[i] = exp(real) * sin(im) / (x);
			}
		}
		else
		{
			complex<double> A(0.,0.), B(0.,0.), u(x, - 0.5);

			ComputeCoeff(u, A, B, m_kappa, m_theta, m_rho, m_nu);

			int i;
			double real,im;
			for (i = 0; i < m_strikes.size(); ++i)
			{
				real = 0.5 * m_logFKs[i].real() - x * m_logFKs[i].imag() + A.real() + B.real() * m_v0 * SQR(m_levels[0]);
				im	= x * m_logFKs[i].real() + 0.5 * m_logFKs[i].imag() + A.imag() + B.imag() *  m_v0 *  SQR(m_levels[0]);

				values[i] = exp(real) * cos(im) / (x * x + 0.25);
			}
		}
	}

protected:

	virtual void ComputeCoeff(const complex<double>& u, complex<double>& A, complex<double>& B, double kappa, double theta, double rho, double nu);

	void ComputeCoeffs(int k, double deltat, const complex<double>& u, complex<double>& A, complex<double>& B, double kappa, double theta, double rho, double nu);
};

























class ARM_MixteHestonOptionPricer : public ARM_HestonOptionPricer
{
protected:
	double	m_sigma;

public:

	ARM_MixteHestonOptionPricer(double reset, 
						   double forward, 
						   double strike, 
						   int sens, 
						   double sigma,
						   double v0, 
						   double kappa, 
						   double theta, 
						   double rho, 
						   double nu, 
						   double shift,
						   const std::vector<double>& times, 
						   const std::vector<double>& levels,
						   bool isDigital = false) :
	ARM_HestonOptionPricer(reset, forward, strike, sens, v0, kappa, theta, rho, nu, shift, times, levels, isDigital),
		m_sigma(sigma)
	{}

	ARM_MixteHestonOptionPricer(double reset, 
						   double forward, 
						   double strike, 
						   int sens, 
						   double sigma,
						   double v0, 
						   double kappa, 
						   double theta, 
						   double rho, 
						   double nu, 
						   double shift,
						   double level,
						   bool isDigital = false) :
	ARM_HestonOptionPricer(reset, forward, strike, sens, v0, kappa, theta, rho, nu, shift, level, isDigital),
		m_sigma(sigma)
	{}

	virtual double integrand(double x)
	{
		complex<double> A(0.,0.), B(0.,0.), u(x, - 0.5);

		double value;

		if (m_isDigital)
		{
			complex<double> A(0.,0.), B(0.,0.), u(x,0);

			ComputeCoeff(u, A, B, m_kappa, m_theta, m_rho, m_nu);

			double real =  - x * m_logFK.imag() + A.real() + B.real() * m_v0 * SQR(m_levels[0])
						  -0.5* m_sigma * m_sigma * m_reset* x * x;
			double im	= x * m_logFK.real() + A.imag() + B.imag() *  m_v0 *  SQR(m_levels[0])
						  -0.5* m_sigma * m_sigma * m_reset * x;

			value = exp(real) * sin(im) / (x);
		}
		else
		{
			complex<double> A(0.,0.), B(0.,0.), u(x, - 0.5);

			ComputeCoeff(u, A, B, m_kappa, m_theta, m_rho, m_nu);

			double real = 0.5 * m_logFK.real() - x * m_logFK.imag() + A.real() + B.real() * m_v0 * SQR(m_levels[0])
						- 0.5 * (x * x + 0.25) * m_sigma * m_sigma * m_reset;
			double im	= x * m_logFK.real() + 0.5 * m_logFK.imag() +  A.imag() + B.imag() *  m_v0 *  SQR(m_levels[0]);

			value = exp(real) * cos(im) / (x * x + 0.25);
		}

		return value;
	}

	virtual void integrands(double x,std::vector<double>& values)
	{
		complex<double> A(0.,0.), B(0.,0.), u(x, - 0.5);

		if (m_isDigital)
		{
			complex<double> A(0.,0.), B(0.,0.), u(x,0);

			ComputeCoeff(u, A, B, m_kappa, m_theta, m_rho, m_nu);

			int i;
			double real,im;
			for (i = 0; i < m_strikes.size(); ++i)
			{
				real = - x * m_logFKs[i].imag() + A.real() + B.real() * m_v0 * SQR(m_levels[0])
						  -0.5* m_sigma * m_sigma * m_reset* x * x;
				im	= x * m_logFKs[i].real() + A.imag() + B.imag() *  m_v0 *  SQR(m_levels[0])
						  -0.5* m_sigma * m_sigma * m_reset * x;
				values[i] = exp(real) * sin(im) / (x);
			}
		}
		else
		{
			complex<double> A(0.,0.), B(0.,0.), u(x, - 0.5);

			ComputeCoeff(u, A, B, m_kappa, m_theta, m_rho, m_nu);

			int i;
			double real,im;
			for (i = 0; i < m_strikes.size(); ++i)
			{
				real = 0.5 * m_logFKs[i].real() - x * m_logFKs[i].imag() + A.real() + B.real() * m_v0 * SQR(m_levels[0])
						- 0.5 * (x * x + 0.25) * m_sigma * m_sigma * m_reset;
				im	= x * m_logFKs[i].real() + 0.5 * m_logFKs[i].imag() + A.imag() + B.imag() *  m_v0 *  SQR(m_levels[0]);

				values[i] = exp(real) * cos(im) / (x * x + 0.25);
			}
		}
	}

};

class ARM_Heston2BOptionPricer : public ARM_HestonOptionPricer
{
protected:
	double		m_v02;
	double		m_kappa2;
	double		m_theta2;
	double		m_rho2;
	double		m_nu2;

public:
	ARM_Heston2BOptionPricer(double reset, 
						   double forward, 
						   double strike, 
						   int sens, 
						   double v01, 
						   double kappa1, 
						   double theta1, 
						   double rho1, 
						   double nu1, 
						   double v02,
						   double kappa2,
						   double theta2,
						   double rho2,
						   double nu2,
						   double shift,
						   const std::vector<double>& times, 
						   const std::vector<double>& levels) :
	ARM_HestonOptionPricer(reset, forward, strike, sens, v01, kappa1, theta1, rho1, nu1, shift, times, levels),
		m_v02(v02), m_kappa2(kappa2), m_theta2(theta2), m_rho2(rho2), m_nu2(nu2)
	{}

	ARM_Heston2BOptionPricer(double reset, 
						   double forward, 
						   double strike, 
						   int sens, 
						   double v01, 
						   double kappa1, 
						   double theta1, 
						   double rho1, 
						   double nu1, 
						   double v02,
						   double kappa2,
						   double theta2,
						   double rho2,
						   double nu2,
						   double shift,
						   double level) :
	ARM_HestonOptionPricer(reset, forward, strike, sens, v01, kappa1, theta1, rho1, nu1, shift, level),
		m_v02(v02), m_kappa2(kappa2), m_theta2(theta2), m_rho2(rho2), m_nu2(nu2)
	{}

	virtual double integrand(double x)
	{
		complex<double> A1(0.,0.), B1(0.,0.), A2(0.,0.), B2(0.,0.), u(x, - 0.5);

		ComputeCoeff(u, A1, B1, m_kappa, m_theta, m_rho, m_nu);
		ComputeCoeff(u, A2, B2, m_kappa2, m_theta2, m_rho2, m_nu2);
		
		double real = 0.5 * m_logFK.real() - x * m_logFK.imag() + A1.real() + B1.real() * m_v0 * SQR(m_levels[0])
					+ A2.real() + B2.real() * m_v02 * SQR(m_levels[0]);
		double im	= x * m_logFK.real() + 0.5 * m_logFK.imag() + A1.imag() + B1.imag() *  m_v0 *  SQR(m_levels[0])
					+ A2.imag() + B2.imag() * m_v02 * SQR(m_levels[0]);

		return exp(real) * cos(im) / (x * x + 0.25);
	}
};

class ARM_HestonOptionPricerVVolt : public ARM_HestonOptionPricer
{
protected:
	ARM_Curve	m_nut;
	int			m_idxNu;

public:
	ARM_HestonOptionPricerVVolt(double reset, 
						   double forward, 
						   double strike, 
						   int sens, 
						   double v0, 
						   double kappa, 
						   double theta, 
						   double rho, 
						   const ARM_Curve& nut,
						   double shift,
						   double level) : 
	ARM_HestonOptionPricer(reset, forward, strike, sens, v0, kappa, theta, rho, 1., shift, level), m_nut(nut)
	{
		m_idxNu = 0;
		int size = m_nut.GetAbscisses().size();
		while(m_nut.GetAbscisses()[m_idxNu]/365. < reset) 
		{
			if(m_idxNu == size-1) break;
			m_idxNu++;
		}
	}

public:
	void	SetLevel(double level) {m_levels[0] = level;};

protected:

	virtual void ComputeCoeff(const complex<double>& u, complex<double>& A, complex<double>& B, double kappa, double theta, double rho, double nu);

};

CC_END_NAMESPACE()

#endif
