
#ifndef _GP_CF_SMILECALIB_H
#define _GP_CF_SMILECALIB_H

#include "firsttoinc.h"
#include "gpbase/port.h"
#include "gpbase/gpvector.h"
#include "gpbase/typedef.h"
#include "gpbase/removenagwarning.h"
#include "gpclosedforms/basic_distributions.h"
#include "gpnumlib/levmarq.h"
#include "gpnumlib/solver.h"
#include "gpclosedforms/vanilla_normal.h"
#include "gpclosedforms/inverse.h"
#include "gpclosedforms/sabrvanilla.h"
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/heston_pricer.h"
#include "gpclosedforms/bisabr_spreadoption.h"
#include "gpclosedforms/sabrbdiff1.h"
#include "gpclosedforms/normal_heston.h"
#include "gpclosedforms/sabrimpliedvol.h"
#include "gpclosedforms/vanilla_shifted_lognormal.h"
#include "gpclosedforms/merton.h"

CC_BEGIN_NAMESPACE(ARM)

class ARM_SmileCalibration_Params : public ARM_RootObject
{
public:

	enum {SABR, Heston, NormHeston, BiSABR};

protected:
	
	int	m_type;

public:

	ARM_SmileCalibration_Params()
	{
	}

	ARM_SmileCalibration_Params(int type)
	{
		m_type = type;
	}

	ARM_SmileCalibration_Params(const ARM_SmileCalibration_Params& rhs) :
	m_type(rhs.m_type)
	{
	}

	virtual ~ARM_SmileCalibration_Params()
	{
	}

public:

	virtual ARM_Object*		Clone() const { return new ARM_SmileCalibration_Params(*this); };
	virtual string			ExportShortName() const {return "LSCPX";};
	virtual	string			toString(const string& indent="",const string& nextIndent="") const;

	int		type() const	{return m_type;};
};

inline string ARM_SmileCalibration_Params::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "Smile Calibration Param Structure \n";

	return os.str();	
}

class ARM_SmileCalibration
{
protected:

	std::vector<double>		m_resets;
	std::vector<double>		m_forwards;

	ARM_VectorVector	m_mktvols;
	ARM_VectorVector	m_strikes;
	std::vector<double>		m_weights;
	
	bool				m_calibConstraint;
	std::vector<double>		m_constraintStrikes;
	std::vector<double>		m_constraintVols;
	double				m_constrainedParamLBound;
	double				m_constrainedParamUBound;

	int					m_idx;
	int					m_locidx;

	double				VolToBP;

	bool				m_calibsuccess;
	double				m_quaderror;

public:

	ARM_SmileCalibration() 
	{
		VolToBP = sqrt(250.)/10000.;
	}

	virtual ~ARM_SmileCalibration();

public:

	void	Init(double reset, double forward, 
				 const std::vector<double>& mktvols, const std::vector<double>& strikes,
				 bool calibConstraint, double constraintK, double constraintVol,
				 ARM_SmileCalibration_Params * params);

	void	Init(const std::vector<double>& resets, const std::vector<double>& forwards,
				 const ARM_VectorVector& mktvols, const ARM_VectorVector& strikes,
				 bool calibConstraint, const std::vector<double>& constraintStrikes, const std::vector<double>& constraintVols,
				 const std::vector<double>& weights,
				 ARM_SmileCalibration_Params * params);


	virtual bool	Calibrate();

	bool	IsCalibSuccess() const	{return m_calibsuccess;};
	double	QuadraticErr() const	{return m_quaderror;};

protected:

	virtual void	InitOneReset(ARM_SmileCalibration_Params * params) = 0;

	virtual void	InitMultipleResets(ARM_SmileCalibration_Params * params) = 0;

	virtual double	ImpliedVol(double strike, int sens) = 0;
	
	double			ImpliedVol(double strike) {return strike > m_forwards[m_idx] ? ImpliedVol(strike, 1) : ImpliedVol(strike, -1);};

	virtual void	SetParam(double param[]) = 0;

	virtual void	SetConstrainedParam(double x) = 0;

	void			CalibConstraint();

	virtual void	SetLocalParam(double param[]) {};

	virtual void	SetIdx(int k);

	virtual void	InitfGuess(std::vector<double>& fguess) = 0;

	class CalibConstraintFunc : public DoubleToDoubleFunc
	{
	private:
		ARM_SmileCalibration * _this;

	public:
		CalibConstraintFunc(ARM_SmileCalibration * calib)
		{
			_this = calib;
		}

	public:

		double operator()(double x) const
		{
			_this->SetConstrainedParam(x);

			return _this->ImpliedVol(_this->m_constraintStrikes[_this->m_idx]);
		}

	};
	
	class CalibFunc : public ARM_LEVMARQFunc 
	{
	private:
		ARM_SmileCalibration * _this;

	public:
		CalibFunc(ARM_SmileCalibration * calib)
		{
			_this = calib;
		}

	public:

		void operator()(double p[], double hx[], int m, int n, void * adata = NULL) const
		{
			_this->SetParam(p);

			for(int k = 0, ind = 0; k < _this->m_resets.size(); k++)
			{
				_this->SetIdx(k);

				if(fabs(_this->m_weights[k]) < K_DOUBLE_TOL) continue;

				_this->SetLocalParam(p);

				if(_this->m_calibConstraint) _this->CalibConstraint();

				for(int i = 0; i < (*_this->m_strikes[k]).size(); i++, ind++)
				{
					hx[ind] = (*_this->m_mktvols[k])[i] - _this->ImpliedVol((*_this->m_strikes[k])[i]);
				}
			}
		}
	};

	friend class CalibConstraintFunc;
	friend class CalibFunc;

protected:
	
	void AdjustVolToBP()
	{
		for(int k = 0; k < m_mktvols.size(); k++) *m_mktvols[k]  *= 1/VolToBP;
		
		std::transform(m_constraintVols.begin(), m_constraintVols.end(), m_constraintVols.begin(),
               std::bind1st(std::multiplies<double>(),1/VolToBP));
	}

protected:

	class Local_VolInverse : public DoubleToDoubleFunc
	{
	private:
		double	m_fwd;
		double	m_strike;
		int		m_sens;
		double	m_mat;
		
	public:
		Local_VolInverse(double fwd, double K, int sens, double T)
		{
			m_fwd		= fwd;
			m_strike	= K;
			m_sens		= sens;
			m_mat		= T;
		}

	public:
		double operator()(double x) const
		{
			return BS(m_fwd, m_strike, m_mat, x, m_sens);
		}
	};

	double	LocalBSImpliedVol(double fwd, double K, double T, int sens, double price);
};

inline double ARM_SmileCalibration::LocalBSImpliedVol(double fwd, double K, double T, int sens, double price)
{
	double best;
	
	Local_VolInverse func(fwd, K, sens, T);

	double vol = brentSolve(func, price, 1e-8, 2., 1e-8, 100, 0, &best);

	return vol;
}

class ARM_SmileCalibration_SABR;

class ARM_SmileCalibration_Params_SABR : public ARM_SmileCalibration_Params
{
protected:

	double	m_alpha;
	double	m_beta;
	double	m_rho;
	double	m_nu;
	int		m_sabrflag;
	double	m_shift;

	std::vector<double>	m_levels;
	std::vector<double>	m_shifts;
	std::vector<double>	m_rhos;

	bool	m_calibalpha;
	bool	m_calibbeta;
	bool	m_calibnu;
	bool	m_calibshift;
	bool	m_calibrho;
	bool	m_calibeachrho;

public:

	ARM_SmileCalibration_Params_SABR(double alpha, double beta, double rho, double nu, int sabrflag, double shift,
		bool calibalpha, bool calibbeta, bool calibrho, bool calibnu, bool calibshift, bool calibeachrho = false)
	{
		m_alpha			= alpha;
		m_beta			= beta;
		m_rho			= rho;
		m_nu			= nu;
		m_sabrflag		= sabrflag;
		m_shift			= m_beta < 1. ? 1. : shift;

		m_calibalpha	= calibalpha;
		m_calibbeta		= calibbeta;
		m_calibrho		= calibrho;
		m_calibnu		= calibnu;
		m_calibshift	= m_beta < 1. ? false : calibshift;
		m_calibeachrho	= calibrho && calibeachrho ? true : false;
	}

	ARM_SmileCalibration_Params_SABR(const ARM_SmileCalibration_Params_SABR& rhs) : 
		ARM_SmileCalibration_Params(rhs), m_alpha(rhs.m_alpha), m_beta(rhs.m_beta), m_rho(rhs.m_rho), m_nu(rhs.m_nu),
		m_sabrflag(rhs.m_sabrflag), m_shift(rhs.m_shift), m_calibalpha(rhs.m_calibalpha), m_calibbeta(rhs.m_calibbeta),
		m_calibrho(rhs.m_calibrho), m_calibnu(rhs.m_calibnu), m_calibshift(rhs.m_calibshift), m_calibeachrho(rhs.m_calibeachrho),
		m_levels(rhs.m_levels), m_shifts(rhs.m_shifts), m_rhos(rhs.m_rhos)
	{
	}

public:

	double					alpha() const	{return m_alpha;};
	double					beta() const	{return m_beta;};
	double					rho() const		{return m_rhos[0];};
	double					nu() const		{return m_nu;};
	int						flag() const	{return m_sabrflag;};

	const std::vector<double>&	shifts() const	{return m_shifts;};
	const std::vector<double>&	levels() const	{return m_levels;};
	const std::vector<double>&	rhos() const	{return m_rhos;};

	bool					calibalpha() const		{return m_calibalpha;};
	bool					calibbeta() const		{return m_calibbeta;};
	bool					calibnu() const			{return m_calibnu;};
	bool					calibshift() const		{return m_calibshift;};
	bool					calibrho() const		{return m_calibrho;};
	bool					calibeachrho() const	{return m_calibeachrho;};

	bool&					calibalpha()		{return m_calibalpha;};
	bool&					calibbeta()			{return m_calibbeta;};
	bool&					calibnu()			{return m_calibnu;};
	bool&					calibshift() 		{return m_calibshift;};
	bool&					calibrho()			{return m_calibrho;};
	bool&					calibeachrho()		{return m_calibeachrho;};

	friend class ARM_SmileCalibration_SABR;
};

class ARM_SmileCalibration_SABR : public ARM_SmileCalibration
{
protected:

	ARM_SmileCalibration_Params_SABR *	m_params;

protected:

	void	InitOneReset(ARM_SmileCalibration_Params * params);

	void	InitMultipleResets(ARM_SmileCalibration_Params * params);

	double	ImpliedVol(double strike, int sens);

	void	SetParam(double param[]);

	void	SetConstrainedParam(double x);

	void	SetLocalParam(double param[]);

	void	InitfGuess(std::vector<double>& fguess);
};

inline double ARM_SmileCalibration_SABR::ImpliedVol(double strike, int sens)
{
	double alpha = m_params->m_levels.size() == 0 ? m_params->m_alpha : m_params->m_alpha * m_params->m_levels[m_idx];
	double rho = m_params->m_calibeachrho ? m_params->m_rhos[m_idx] : m_params->m_rhos[0];

	if(fabs(m_params->m_beta - 1.) < K_DOUBLE_TOL && m_params->m_shifts[m_idx] < 1.)
	{
		double fwd	= m_forwards[m_idx] / fabs(m_params->m_shifts[m_idx]);
		double K	= m_params->m_shifts[m_idx] < 0. ? fwd * (1. + fabs(m_params->m_shifts[m_idx])) - strike : strike - (m_params->m_shifts[m_idx] - 1.) * fwd; 

		double vol	= CptSABR_implicit_vol_direct(fwd, K, m_resets[m_idx], alpha, 1., rho, m_params->m_nu, 3);

		double opt	= BS(fwd, K, m_resets[m_idx], vol, sens * (m_params->m_shifts[m_idx] < 0. ? -1 : 1));

		return LocalBSImpliedVol(m_forwards[m_idx], strike, m_resets[m_idx], sens, opt);
	}
	else
	{
		return CptSABR_implicit_vol_direct(m_forwards[m_idx], strike, m_resets[m_idx], alpha, m_params->m_beta, rho, m_params->m_nu, m_params->m_sabrflag);
	}
}


class ARM_SmileCalibration_Heston;
class ARM_SmileCalibration_NormalHeston;

class ARM_SmileCalibration_Params_Heston : public ARM_SmileCalibration_Params 
{
protected:

	double			m_v0;
	double			m_kappa;
	double			m_theta;
	double			m_nu;
	double			m_rho;
	double			m_shift;
	double			m_level;

	std::vector<double>	m_shifts;
	std::vector<double>	m_levels;
	std::vector<double>	m_rhos;
	std::vector<double>	m_sigmas;

	bool			m_calibkappa;
	bool			m_calibtheta;
	bool			m_calibrho;
	bool			m_calibnu;
	bool			m_calibshift;
	bool			m_calibeachrho;

	bool			m_bootstraplevel;
	
	bool			m_normalHeston;

public:
	ARM_SmileCalibration_Params_Heston()
	{}

	ARM_SmileCalibration_Params_Heston(double v0, double kappa, double theta, double rho, double nu, double shift, double level,
		bool calibkappa, bool calibtheta, bool calibrho, bool calibnu, bool calibshift, bool calibeachrho = false, bool bootstraplevel = false, bool normalHeston = false) 
	{
		build(v0, kappa, theta, rho, nu, shift, level, std::vector<double>(1,0.), calibkappa, calibtheta, calibrho, calibnu, calibshift, calibeachrho, bootstraplevel, normalHeston);
	}

	ARM_SmileCalibration_Params_Heston(double v0, double kappa, double theta, double rho, double nu, double shift, double level, const std::vector<double>& sigmas,
		bool calibkappa, bool calibtheta, bool calibrho, bool calibnu, bool calibshift, bool calibeachrho = false, bool bootstraplevel = false, bool normalHeston = false) 
	{
		build(v0, kappa, theta, rho, nu, shift, level, sigmas, calibkappa, calibtheta, calibrho, calibnu, calibshift, calibeachrho, bootstraplevel, normalHeston);
	}

	ARM_SmileCalibration_Params_Heston(const ARM_SmileCalibration_Params_Heston& rhs) : 
		ARM_SmileCalibration_Params(rhs), m_v0(rhs.m_v0), m_kappa(rhs.m_kappa), m_theta(rhs.m_theta), 
		m_nu(rhs.m_nu), m_shift(rhs.m_shift), m_level(rhs.m_level), 
		m_calibkappa(rhs.m_calibkappa), m_calibtheta(rhs.m_calibtheta), m_calibrho(rhs.m_calibrho), m_calibnu(rhs.m_calibnu),
		m_calibshift(rhs.m_calibshift), m_calibeachrho(rhs.m_calibeachrho), m_bootstraplevel(rhs.m_bootstraplevel), m_normalHeston(rhs.m_normalHeston),
		m_shifts(rhs.m_shifts), m_levels(rhs.m_levels), m_rhos(rhs.m_rhos), m_sigmas(rhs.m_sigmas)
	{
	}

public:

	void		build(double v0, double kappa, double theta, double rho, double nu, double shift, double level,
					bool calibkappa, bool calibtheta, bool calibrho, bool calibnu, bool calibshift, 
					bool calibeachrho = false, bool bootstraplevel = false, bool normalHeston = false)
	{
		build(v0, kappa, theta, rho, nu, shift, level, std::vector<double>(1,0.), calibkappa, calibtheta, calibrho, calibnu,
			calibshift, calibeachrho, bootstraplevel, normalHeston);
	}
					
	void		build(double v0, double kappa, double theta, double rho, double nu, double shift, double level,
					const std::vector<double>& sigmas,
					bool calibkappa, bool calibtheta, bool calibrho, bool calibnu, bool calibshift, 
					bool calibeachrho = false, bool bootstraplevel = false, bool normalHeston = false)
	{
		m_v0				= v0;
		m_kappa				= kappa;
		m_theta				= theta;
		m_rho				= rho;
		m_nu				= nu;
		m_sigmas			= sigmas;

		m_shift				= normalHeston ? 1. : shift;
		m_level				= level;

		m_calibkappa		= calibkappa;
		m_calibtheta		= calibtheta;
		m_calibrho			= calibrho;
		m_calibnu			= calibnu;
		m_calibshift		= normalHeston ? false : calibshift;
		m_calibeachrho		= calibrho && calibeachrho ? true : false;

		m_bootstraplevel	= normalHeston ? false : bootstraplevel;

		m_normalHeston		= normalHeston;
	}

	double					kappa() const	{return m_kappa;};
	double					v0() const		{return m_v0;};
	double					theta() const	{return m_theta;};
	double					rho() const		{return m_rhos[0];};
	double					nu() const		{return m_nu;};
	double					level() const	{return m_levels[0];};
	double					shift() const	{return m_shifts[0];};

	bool					normalHeston() const	{return m_normalHeston;};

	bool					calibkappa() const		{return m_calibkappa;};
	bool					calibtheta() const		{return m_calibtheta;};
	bool					calibrho() const		{return m_calibrho;};
	bool					calibnu() const			{return m_calibnu;};
	bool					calibshift() const		{return m_calibshift;};
	bool					calibeachrho() const	{return m_calibeachrho;};

	bool&					calibkappa()			{return m_calibkappa;};
	bool&					calibtheta() 			{return m_calibtheta;};
	bool&					calibrho()				{return m_calibrho;};
	bool&					calibeachrho()			{return m_calibeachrho;};
	bool&					calibnu()				{return m_calibnu;};
	bool&					calibshift()			{return m_calibshift;};

	const std::vector<double>&	levels() const	{return m_levels;};
	const std::vector<double>&	shifts() const	{return m_shifts;};
	const std::vector<double>&	rhos() const	{return m_rhos;};

	friend class ARM_SmileCalibration_Heston;
	friend class ARM_SmileCalibration_NormalHeston;

};

class ARM_SmileCalibration_Heston : public ARM_SmileCalibration
{
protected:

	ARM_SmileCalibration_Params_Heston *	m_params;
	
	double								m_thetaLBound;
	double								m_nuLBound;

protected:

	void	InitOneReset(ARM_SmileCalibration_Params * params);

	void	InitMultipleResets(ARM_SmileCalibration_Params * params);

	double	ImpliedVol(double strike, int sens);

	void	SetParam(double param[]);

	void	SetConstrainedParam(double x);

	void	SetLocalParam(double param[]);

	void	InitfGuess(std::vector<double>& fguess);
};

inline double ARM_SmileCalibration_Heston::ImpliedVol(double strike, int sens)
{
	double price = 0.;
	double rho = m_params->m_calibeachrho ? m_params->m_rhos[m_idx] : m_params->m_rhos[0];

	if(m_params->m_bootstraplevel)
	{
		ARM_MixteHestonOptionPricer pricer(m_resets[m_idx], m_forwards[m_idx], strike, sens, m_params->m_sigmas[m_idx],
			m_params->m_v0, m_params->m_kappa, m_params->m_theta, rho, m_params->m_nu, m_params->m_shifts[m_idx],
			m_resets, m_params->m_levels);

		price = pricer.price();
	}
	else
	{
		ARM_MixteHestonOptionPricer pricer(m_resets[m_idx], m_forwards[m_idx], strike, sens, m_params->m_sigmas[m_idx],
			m_params->m_v0, m_params->m_kappa, m_params->m_theta, rho, m_params->m_nu, m_params->m_shifts[m_idx], 
			m_params->m_levels[m_idx]);

		price = pricer.price();
	}

	return LocalBSImpliedVol(m_forwards[m_idx], strike, m_resets[m_idx], sens, price);
}

class ARM_SmileCalibration_NormalHeston : public ARM_SmileCalibration_Heston
{
protected:

	void	InitOneReset(ARM_SmileCalibration_Params * params);

	void	InitMultipleResets(ARM_SmileCalibration_Params * params);

	double	ImpliedVol(double strike, int sens);

	void	InitfGuess(std::vector<double>& fguess);
};

inline double ARM_SmileCalibration_NormalHeston::ImpliedVol(double strike, int sens)
{
	double v0 = m_params->m_v0 * m_params->m_levels[m_idx] * m_params->m_levels[m_idx];
	double theta = m_params->m_theta * m_params->m_levels[m_idx] * m_params->m_levels[m_idx];
	double nu = m_params->m_nu * m_params->m_levels[m_idx];
	double rho = m_params->m_calibeachrho ? m_params->m_rhos[m_idx] : m_params->m_rhos[0];

	double lambda = -1.; //fabs(strike - m_forwards[m_idx]) < K_DOUBLE_TOL ? -1. : 0.1;
	double opt = NormalHeston(rho, m_params->m_kappa, theta, nu, v0, m_forwards[m_idx], strike, m_resets[m_idx], sens,lambda, 20,20);

	bool success;
	return  opt < K_DOUBLE_TOL ? 0. : VanillaImpliedVol_N(m_forwards[m_idx], opt, strike, m_resets[m_idx], sens, NULL, &success) / VolToBP;	
}

class ARM_SmileCalibration_BiSABR;

class ARM_SmileCalibration_Params_BiSABR : public ARM_SmileCalibration_Params
{
protected:

	double	m_rhoS1S2;
	double	m_rhoS1V2;
	double	m_rhoS2V1;
	double	m_rhoV1V2;

	double	m_forward1;
	double	m_alpha1;
	double	m_beta1;
	double	m_rho1;
	double	m_nu1;

	double	m_forward2;
	double	m_alpha2;
	double	m_beta2;
	double	m_rho2;
	double	m_nu2;

	bool	m_calibrhoSS;
	bool	m_calibrhoSV;
	bool	m_calibrhoVV;

public:
	
	ARM_SmileCalibration_Params_BiSABR(double fwd1, double alpha1, double beta1, double rho1, double nu1,
						 double fwd2, double alpha2, double beta2, double rho2, double nu2,
						 double rhoS1S2, double rhoS1V2, double rhoS2V1, double rhoV1V2,
						 bool calibrhoS1S2, bool calibrhoSV, bool calibrhoVV)
	{
		m_forward1	= fwd1;
		m_alpha1	= alpha1;
		m_beta1		= beta1;
		m_rho1		= rho1;
		m_nu1		= nu1;

		m_forward2	= fwd2;
		m_alpha2	= alpha2;
		m_beta2		= beta2;
		m_rho2		= rho2;
		m_nu2		= nu2;

		m_rhoS1S2	= rhoS1S2;
		m_rhoS1V2	= rhoS1V2;
		m_rhoS2V1	= rhoS2V1;
		m_rhoV1V2	= rhoV1V2;

		m_calibrhoSS	= calibrhoS1S2;
		m_calibrhoSV	= calibrhoSV;
		m_calibrhoVV	= calibrhoVV;		
	}

public:

	double	rhoS1S2() const	{return m_rhoS1S2;};
	double	rhoS1V2() const	{return m_rhoS1V2;};
	double	rhoS2V1() const {return m_rhoS2V1;};
	double	rhoV1V2() const {return m_rhoV1V2;};

	friend class ARM_SmileCalibration_BiSABR;
};

class ARM_SmileCalibration_BiSABR : public ARM_SmileCalibration
{
protected:

	ARM_SmileCalibration_Params_BiSABR *	m_params;

protected:

	void	InitOneReset(ARM_SmileCalibration_Params * params);

	void	InitMultipleResets(ARM_SmileCalibration_Params * params);

	double	ImpliedVol(double strike, int sens);

	void	SetParam(double param[]);

	void	SetConstrainedParam(double x);

	void	InitfGuess(std::vector<double>& fguess);
};

inline double ARM_SmileCalibration_BiSABR::ImpliedVol(double strike, int sens)
{
	double opt = BiSABR_SpreadOption2(m_params->m_forward1, m_params->m_alpha1, m_params->m_beta1, m_params->m_rho1, m_params->m_nu1,
									 m_params->m_forward2, m_params->m_alpha2, m_params->m_beta2, m_params->m_rho2, m_params->m_nu2,
									 strike, m_resets[m_idx], m_params->m_rhoS1S2, m_params->m_rhoV1V2, m_params->m_rhoS1V2, m_params->m_rhoS2V1);

	bool success;

	return VanillaImpliedVol_N(m_params->m_forward1 - m_params->m_forward2, opt, strike, m_resets[m_idx], sens, NULL, &success) / VolToBP;
}

class ARM_SmileCalibration_Heston2b;

class ARM_SmileCalibration_Params_Heston2b : public ARM_SmileCalibration_Params 
{
protected:

	double			m_v01;
	double			m_kappa1;
	double			m_theta1;
	double			m_nu1;
	double			m_rho1;
	double			m_v02;
	double			m_kappa2;
	double			m_theta2;
	double			m_nu2;
	double			m_rho2;
	double			m_shift;
	double			m_level;

	std::vector<double>	m_shifts;
	std::vector<double>	m_levels;
	std::vector<double>	m_rhos1;
	std::vector<double>	m_rhos2;

	bool			m_calibkappa1;
	bool			m_calibtheta1;
	bool			m_calibrho1;
	bool			m_calibnu1;
	bool			m_calibeachrho1;
	bool			m_calibkappa2;
	bool			m_calibtheta2;
	bool			m_calibrho2;
	bool			m_calibnu2;
	bool			m_calibeachrho2;
	bool			m_calibshift;

	bool			m_bootstraplevel;

public:
	ARM_SmileCalibration_Params_Heston2b()
	{}

	ARM_SmileCalibration_Params_Heston2b(double v01, double kappa1, double theta1, double rho1, double nu1, 
		double v02, double kappa2, double theta2, double rho2, double nu2,
		double shift, double level,
		bool calibkappa1, bool calibtheta1, bool calibrho1, bool calibnu1, 
		bool calibkappa2, bool calibtheta2, bool calibrho2, bool calibnu2, 
		bool calibshift, bool calibeachrho1 = false, bool calibeachrho2 = false, bool bootstraplevel = false) 
	{
		build(v01, kappa1, theta1, rho1, nu1, v02, kappa2, theta2, rho2, nu2, 
			shift, level, calibkappa1, calibtheta1, calibrho1, calibnu1, 
			calibkappa2, calibtheta2, calibrho2, calibnu2, 
			calibshift, calibeachrho1, calibeachrho2, bootstraplevel);
	}


public:

	void		build(double v01, double kappa1, double theta1, double rho1, double nu1, 
					  double v02, double kappa2, double theta2, double rho2, double nu2,
					  double shift, double level,
					  bool calibkappa1, bool calibtheta1, bool calibrho1, bool calibnu1, 
					  bool calibkappa2, bool calibtheta2, bool calibrho2, bool calibnu2, 
					  bool calibshift, bool calibeachrho1 = false, bool calibeachrho2 = false, bool bootstraplevel = false) 
	{
		m_v01				= v01;
		m_kappa1			= kappa1;
		m_theta1			= theta1;
		m_rho1				= rho1;
		m_nu1				= nu1;
		m_v02				= v02;
		m_kappa2			= kappa2;
		m_theta2			= theta2;
		m_rho2				= rho2;
		m_nu2				= nu2;

		m_shift				= shift;
		m_level				= level;

		m_calibkappa1		= calibkappa1;
		m_calibtheta1		= calibtheta1;
		m_calibrho1			= calibrho1;
		m_calibnu1			= calibnu1;
		m_calibeachrho1		= calibrho1 && calibeachrho1 ? true : false;

		m_calibkappa2		= calibkappa2;
		m_calibtheta2		= calibtheta2;
		m_calibrho2			= calibrho2;
		m_calibnu2			= calibnu2;
		m_calibeachrho2		= calibrho2 && calibeachrho2 ? true : false;

		m_calibshift		= calibshift;
		m_bootstraplevel	= bootstraplevel;

	}

	double					kappa1() const	{return m_kappa1;};
	double					v01() const		{return m_v01;};
	double					theta1() const	{return m_theta1;};
	double					rho1() const	{return m_rhos1[0];};
	double					nu1() const		{return m_nu1;};
	double					kappa2() const	{return m_kappa2;};
	double					v02() const		{return m_v02;};
	double					theta2() const	{return m_theta2;};
	double					rho2() const	{return m_rhos2[0];};
	double					nu2() const		{return m_nu2;};

	double					level() const	{return m_levels[0];};
	double					shift() const	{return m_shifts[0];};

	bool					calibkappa1() const		{return m_calibkappa1;};
	bool					calibtheta1() const		{return m_calibtheta1;};
	bool					calibrho1() const		{return m_calibrho1;};
	bool					calibnu1() const		{return m_calibnu1;};
	bool					calibeachrho1() const	{return m_calibeachrho1;};

	bool					calibkappa2() const		{return m_calibkappa2;};
	bool					calibtheta2() const		{return m_calibtheta2;};
	bool					calibrho2() const		{return m_calibrho2;};
	bool					calibnu2() const		{return m_calibnu2;};
	bool					calibeachrho2() const	{return m_calibeachrho2;};

	bool					calibshift() const		{return m_calibshift;};

	bool&					calibkappa1()			{return m_calibkappa1;};
	bool&					calibtheta1() 			{return m_calibtheta1;};
	bool&					calibrho1()				{return m_calibrho1;};
	bool&					calibeachrho1()			{return m_calibeachrho1;};
	bool&					calibnu1()				{return m_calibnu1;};

	bool&					calibkappa2()			{return m_calibkappa2;};
	bool&					calibtheta2() 			{return m_calibtheta2;};
	bool&					calibrho2()				{return m_calibrho2;};
	bool&					calibeachrho2()			{return m_calibeachrho2;};
	bool&					calibnu2()				{return m_calibnu2;};

	bool&					calibshift()			{return m_calibshift;};

	const std::vector<double>&	levels() const	{return m_levels;};
	const std::vector<double>&	shifts() const	{return m_shifts;};
	const std::vector<double>&	rhos1() const	{return m_rhos1;};
	const std::vector<double>&	rhos2() const	{return m_rhos2;};

	friend class ARM_SmileCalibration_Heston2b;

};

class ARM_SmileCalibration_Heston2b : public ARM_SmileCalibration
{
protected:

	ARM_SmileCalibration_Params_Heston2b *	m_params;
	
	double								m_thetaLBound;
	double								m_nuLBound;

protected:

	void	InitOneReset(ARM_SmileCalibration_Params * params);

	void	InitMultipleResets(ARM_SmileCalibration_Params * params);

	double	ImpliedVol(double strike, int sens);

	void	SetParam(double param[]);

	void	SetConstrainedParam(double x);

	void	SetLocalParam(double param[]);

	void	InitfGuess(std::vector<double>& fguess);
};

inline double ARM_SmileCalibration_Heston2b::ImpliedVol(double strike, int sens)
{
	double price = 0.;
	double rho1 = m_params->m_calibeachrho1 ? m_params->m_rhos1[m_idx] : m_params->m_rhos1[0];
	double rho2 = m_params->m_calibeachrho2 ? m_params->m_rhos2[m_idx] : m_params->m_rhos2[0];

	if(m_params->m_bootstraplevel)
	{
		ARM_Heston2BOptionPricer pricer(m_resets[m_idx], m_forwards[m_idx], strike, sens, 
			m_params->m_v01, m_params->m_kappa1, m_params->m_theta1, rho1, m_params->m_nu1, 
			m_params->m_v02, m_params->m_kappa2, m_params->m_theta2, rho2, m_params->m_nu2, 
			m_params->m_shifts[m_idx],m_resets, m_params->m_levels);

		price = pricer.price();
	}
	else
	{
		ARM_Heston2BOptionPricer pricer(m_resets[m_idx], m_forwards[m_idx], strike, sens, 
			m_params->m_v01, m_params->m_kappa1, m_params->m_theta1, rho1, m_params->m_nu1, 
			m_params->m_v02, m_params->m_kappa2, m_params->m_theta2, rho2, m_params->m_nu2, 
			m_params->m_shifts[m_idx], m_params->m_levels[m_idx]);

		price = pricer.price();
	}

	return LocalBSImpliedVol(m_forwards[m_idx], strike, m_resets[m_idx], sens, price);
}

class ARM_SmileCalibration_SABR2beta;

class ARM_SmileCalibration_Params_SABR2beta : public ARM_SmileCalibration_Params
{
protected:

	double	m_alpha;
	double	m_beta1;
	double	m_beta2;
	double	m_rho;
	double	m_nu;
	double	m_zero;
	double  m_lambda;

	
	bool	m_calibalpha;
	bool	m_calibbeta1;
	bool	m_calibbeta2;
	bool	m_calibrho;
	bool	m_calibnu;
	bool	m_calibzero;
	bool	m_caliblambda;

public:

	ARM_SmileCalibration_Params_SABR2beta(double alpha, double beta1, double beta2, double rho, double nu, double zero, double lambda,
		bool calibalpha, bool calibbeta1, bool calibbeta2, bool calibrho, bool calibnu, bool calibzero, bool caliblambda)
	{
		m_alpha			= alpha;
		m_beta1			= beta1;
		m_beta2			= beta2;
		m_rho			= rho;
		m_nu			= nu;
		m_zero			= zero;
		m_lambda		= lambda;

		m_calibalpha	= calibalpha;
		m_calibbeta1	= calibbeta1;
		m_calibbeta2	= calibbeta2;
		m_calibrho		= calibrho;
		m_calibnu		= calibnu;
		m_calibzero		= calibzero;
		m_caliblambda	= caliblambda;
	}

	ARM_SmileCalibration_Params_SABR2beta(const ARM_SmileCalibration_Params_SABR2beta& rhs) : 
		ARM_SmileCalibration_Params(rhs), m_alpha(rhs.m_alpha), m_beta1(rhs.m_beta1), m_beta2(rhs.m_beta2), m_rho(rhs.m_rho), m_nu(rhs.m_nu),
		m_zero(rhs.m_zero), m_lambda(rhs.m_lambda), m_calibalpha(rhs.m_calibalpha), m_calibbeta1(rhs.m_calibbeta1), m_calibbeta2(rhs.m_calibbeta2),
		m_calibrho(rhs.m_calibrho), m_calibnu(rhs.m_calibnu), m_calibzero(rhs.m_calibzero), m_caliblambda(rhs.m_caliblambda)
	{
	}

public:

	double					alpha() const	{return m_alpha;};
	double					beta1() const	{return m_beta1;};
	double					beta2() const	{return m_beta2;};
	double					rho() const		{return m_rho;};
	double					nu() const		{return m_nu;};
	double					zero() const	{return m_zero;};
	double					lambda() const	{return m_lambda;};
	
	friend class ARM_SmileCalibration_SABR2beta;
};

class ARM_SmileCalibration_SABR2beta : public ARM_SmileCalibration
{
protected:

	ARM_SmileCalibration_Params_SABR2beta *	m_params;

protected:

	void	InitOneReset(ARM_SmileCalibration_Params * params);

	void	InitMultipleResets(ARM_SmileCalibration_Params * params);

	double	ImpliedVol(double strike, int sens);

	void	SetParam(double param[]);

	void	SetConstrainedParam(double x);

	void	SetLocalParam(double param[]);

	void	InitfGuess(std::vector<double>& fguess);

public:

	static double vol(double fwd, double strike, double mat, double alpha, double beta1, double beta2, double rho, double nu, double zero, double lambda)
	{
		double sens = fwd>strike?-1.:1.;

		if (fabs(zero)>K_NEW_DOUBLE_TOL)
		{
			double vol	= sabr2b_implicit_vol(fwd,strike,mat,alpha,beta1,beta2,rho,nu,zero,lambda);

			double opt	= BS(fwd-zero, strike-zero, mat, vol, sens);

			double best;
	
			Local_VolInverse func(fwd, strike, sens, mat);

			double out = brentSolve(func, opt, 1e-8, 2., 1e-8, 100, 0, &best);

			return out;
		}
		else
		{
			return sabr2b_implicit_vol(fwd,strike,mat,alpha,beta1,beta2,rho, nu,0.,lambda);
		}
	}
};

inline double ARM_SmileCalibration_SABR2beta::ImpliedVol(double strike, int sens)
{
	double fwd = m_forwards[m_idx];
	double mat = m_resets[m_idx];
	if (fabs(m_params->m_zero)>K_NEW_DOUBLE_TOL)
	{
		double vol	= sabr2b_implicit_vol(fwd,strike,mat,m_params->m_alpha,m_params->m_beta1,m_params->m_beta2,m_params->m_rho, m_params->m_nu,m_params->m_zero,m_params->m_lambda);

		double opt	= BS(fwd-m_params->m_zero, strike-m_params->m_zero, mat, vol, sens);

		return LocalBSImpliedVol(fwd, strike, mat, sens, opt);
	}
	else
	{
		return sabr2b_implicit_vol(fwd,strike,mat,m_params->m_alpha,m_params->m_beta1,m_params->m_beta2,m_params->m_rho, m_params->m_nu,0.,m_params->m_lambda);
	}
}

class ARM_SmileCalibration_Merton;

class ARM_SmileCalibration_Params_Merton : public ARM_SmileCalibration_Params 
{
private:
	double		m_sigma;
	double		m_lambda1;
	double		m_U1;
	double		m_lambda2;
	double		m_U2;

	bool		m_calibsigma;
	bool		m_caliblambda1;
	bool		m_calibU1;
	bool		m_caliblambda2;
	bool		m_calibU2;

public:
	ARM_SmileCalibration_Params_Merton()
	{}

	ARM_SmileCalibration_Params_Merton(double sigma, double lambda1, double U1, double lambda2, double U2, 
		bool calibsigma, bool caliblambda1, bool calibU1, bool caliblambda2, bool calibU2) : 
		m_sigma(sigma), m_lambda1(lambda1), m_U1(U1), m_lambda2(lambda2), m_U2(U2), 
		m_calibsigma(calibsigma), m_caliblambda1(caliblambda1), m_calibU1(calibU1),
		m_caliblambda2(caliblambda2), m_calibU2(calibU2)
	{}

	~ARM_SmileCalibration_Params_Merton()
	{}

public:

	double		sigma() const			{return m_sigma;};
	double		lambda1() const			{return m_lambda1;}
	double		U1() const				{return m_U1;};
	double		lambda2() const			{return m_lambda2;}
	double		U2() const				{return m_U2;};

	bool		calibsigma() const		{return m_calibsigma;};
	bool		caliblambda1() const	{return m_caliblambda1;};
	bool		calibU1() const			{return m_calibU1;};
	bool		caliblambda2() const	{return m_caliblambda2;};
	bool		calibU2() const			{return m_calibU2;};

	bool&		calibsigma()			{return m_calibsigma;};
	bool&		caliblambda1()			{return m_caliblambda1;};
	bool&		calibU1()				{return m_calibU1;};
	bool&		caliblambda2()			{return m_caliblambda2;};
	bool&		calibU2()				{return m_calibU2;};

	friend class ARM_SmileCalibration_Merton;
};

class ARM_SmileCalibration_Merton : public ARM_SmileCalibration
{
protected:

	ARM_SmileCalibration_Params_Merton *	m_params;

protected:

	void	InitOneReset(ARM_SmileCalibration_Params * params);

	void	InitMultipleResets(ARM_SmileCalibration_Params * params);

	double	ImpliedVol(double strike, int sens);

	void	SetParam(double param[]);

	void	SetConstrainedParam(double x);

	void	SetLocalParam(double param[]);

	void	InitfGuess(std::vector<double>& fguess);
};


class ARM_SmileCalibration_Spread2Heston;

class ARM_SmileCalibration_Params_Spread2Heston : public ARM_SmileCalibration_Params
{
private:
	ARM_SmileCalibration_Params_Heston	m_params1;
	ARM_SmileCalibration_Params_Heston	m_params2;
	double								m_shift1;
	double								m_shift2;
	double								m_correl;

public:
	ARM_SmileCalibration_Params_Spread2Heston(double v0, double kappa, double theta, bool calibtheta) : 
	  m_params1(v0, kappa, theta, 0., 0., 1., 1., false, calibtheta, true, true, false),
	  m_params2(v0, kappa, theta, 0., 0., 1., 1., false, true, true, true, false)
	{
	}

	double		v0() const	{return m_params1.v0();};
	double		kappa() const {return m_params1.kappa();};
	double		theta() const {return m_params1.theta();};
	double		nu() const {return m_params1.nu();};
	double		rho1() const {return m_params1.rho();};
	double		rho2() const {return m_params2.rho();};
	double		shift1() const {return m_shift1;};
	double		shift2() const {return m_shift2;};
	double		level1() const {return m_params1.level();};
	double		level2() const {return m_params2.level();};
	double		correl() const {return m_correl;};

	friend class ARM_SmileCalibration_Spread2Heston;
};

inline double Spread2HestonVanilla(double reset, double fwd1, double fwd2, double strike, int callPut, 
								   double v0, double kappa, double theta, 
								   double nu, double rho1, double rho2, double shift1, double shift2,
								   double level1, double level2, double correl,
								   double index1lev = 1., double index2lev = 1.)
{
	double shfwd1	= index1lev * fwd1 / shift1;
	double shfwd2	= index2lev * fwd2 / shift2;
	double lev1		= level1 * shfwd1;
	double lev2		= level2 * shfwd2;
	double lev		= sqrt(lev1 * lev1 + lev2 * lev2 - 2. * lev1 * lev2 * correl);
	double nhrho	= (rho1 * lev1 - rho2 * lev2) / lev;
	
	if(nhrho > 0.999) nhrho = 0.999;
	if(nhrho < -0.999) nhrho = -0.999;

	double nhv0		= v0 * lev * lev;
	double nhtheta	= theta * lev * lev;
	double nhnu		= nu * lev;
	double fwd		= shfwd1 - shfwd2;
	double K		= strike + (1. - shift1) * shfwd1 - (1. - shift2) * shfwd2;

	double lambda	= -1.; //calib ? -1. : 0.1; //fabs(strike - (fwd1 - fwd2)) < K_DOUBLE_TOL ? -1. : 0.1;

	double price	= NormalHeston(nhrho, kappa, nhtheta, nhnu, nhv0, fwd, K, reset, callPut,lambda);

	return price < 0. ? 0. : price;
}

class ARM_SmileCalibration_Spread2Heston : public ARM_SmileCalibration 
{
private:
	double			m_reset;
	double			m_fwd1;
	double			m_fwd2;

	std::vector<double>	m_mktvols1;
	std::vector<double>	m_strikes1;
	double			m_constrvol1;
	double			m_constrK1;
	std::vector<double>	m_mktvols2;
	std::vector<double>	m_strikes2;
	double			m_constrvol2;
	double			m_constrK2;

	std::vector<double>	m_mktvolspread;
	std::vector<double>	m_strikespread;
	double			m_constrVolSpread;
	double			m_constrKSpread;

	ARM_SmileCalibration_Params_Spread2Heston * m_params;

public:
	ARM_SmileCalibration_Spread2Heston()
	{}

	~ARM_SmileCalibration_Spread2Heston()
	{}

public:

	void	Init(double reset, 
				double fwd1, const std::vector<double>& mktvols1, const std::vector<double>& strikes1, double constrvol1, double constrK1, 
				double fwd2, const std::vector<double>& mktvols2, const std::vector<double>& strikes2, double constrvol2, double constrK2,
				const std::vector<double>& mktvolspread, const std::vector<double>& strikespread, double constrvolSpread, double constrKSpread,
				ARM_SmileCalibration_Params_Spread2Heston * params);

	bool	Calibrate();

	void	CalibrateCorrel(double reset, double fwd1, double fwd2, double mktvolspread, double strikespread,
				ARM_SmileCalibration_Params_Spread2Heston * params);

protected:

	virtual void	InitOneReset(ARM_SmileCalibration_Params * params) {};

	virtual void	InitMultipleResets(ARM_SmileCalibration_Params * params) {};

	virtual double	ImpliedVol(double strike, int sens) {return 0.;};
	
	virtual void	SetParam(double param[]) {};

	virtual void	SetConstrainedParam(double x) {};

	virtual void	InitfGuess(std::vector<double>& fguess) {};

private:
	
	void	SetParamsForFwdCalib(double p[]);
	void	SetConstParamForFwd(double x, int i);
	void	SetCorrel(double x) {m_params->m_correl = x;};
	void	SetParamsForSpreadCalib(double p[]);

	double	ImpliedVolFwd(double strike, int i);
	double	ImpliedVolSpread(double strike);

	void	CalibFwd();
	void	CalibFwdConstraint(int i);
	void	CalibSpreadConstraint();

	class CalibSpreadFunc : public ARM_LEVMARQFunc
	{
	private:
		ARM_SmileCalibration_Spread2Heston * _this;

	public:
		CalibSpreadFunc(ARM_SmileCalibration_Spread2Heston * calib)
		{
			_this = calib;
		}

	public:
		void operator()(double p[], double hx[], int m, int n, void * adata = NULL) const
		{
			_this->SetParamsForSpreadCalib(p);

			_this->CalibFwd();

			_this->CalibSpreadConstraint();

			for(int i = 0; i < _this->m_strikespread.size(); i++)
			{
				hx[i] = _this->m_mktvolspread[i] - _this->ImpliedVolSpread(_this->m_strikespread[i]);
			}
		}
	};

	class CalibSpreadConstrFunc : public DoubleToDoubleFunc
	{
	private:
		ARM_SmileCalibration_Spread2Heston * _this;

	public:
		CalibSpreadConstrFunc(ARM_SmileCalibration_Spread2Heston * calib)
		{
			_this = calib;
		}

	public:
		double operator()(double x) const
		{
			_this->SetCorrel(x);

			return _this->ImpliedVolSpread(_this->m_constrKSpread);
		}
	};

	friend class CalibSpreadFunc;
	friend class CalibSpreadConstrFunc;
};

inline double ARM_SmileCalibration_Spread2Heston::ImpliedVolSpread(double strike)
{
	int sens	= strike > m_fwd1 - m_fwd2 ? 1 : -1;

	double opt	= Spread2HestonVanilla(m_reset, m_fwd1, m_fwd2, strike, sens, m_params->v0(), m_params->kappa(), 
						m_params->theta(), m_params->nu(), m_params->rho1(), m_params->rho2(), 
						m_params->shift1(), m_params->shift2(), m_params->level1(), m_params->level2(),
						m_params->correl());
	
	bool success;

	return  opt < K_DOUBLE_TOL ? 0. : VanillaImpliedVol_N(m_fwd1-m_fwd2, opt, strike, m_reset, sens, NULL, &success) / VolToBP;	
}

class ARM_SmileCalibration_Spread2HestonCAP;

class ARM_Spread2Heston_CAPCalibration_Params
{
private:
	std::vector<double>							m_fwdReset;
	std::vector<double>							m_fwd1;
	std::vector<double>							m_fwd2;
	ARM_VectorVector						m_mktvols1;
	ARM_VectorVector						m_strikes1;
	std::vector<double>							m_atm1;
	ARM_VectorVector						m_mktvols2;
	ARM_VectorVector						m_strikes2;
	std::vector<double>							m_atm2;
	ARM_IntVector							m_isCalibrated;
	ARM_IntVector							m_currFwdCalib;

	ARM_SmileCalibration_Params_Heston *	m_params1;
	ARM_SmileCalibration_Params_Heston *	m_params2;
	std::vector<double>							m_shift1;
	std::vector<double>							m_shift2;

	double									m_correl;
	
	bool									m_constrainedCorrel;
	bool									m_calibAllShift;

public:
	ARM_Spread2Heston_CAPCalibration_Params()
	{}

	~ARM_Spread2Heston_CAPCalibration_Params();

	ARM_Spread2Heston_CAPCalibration_Params(const std::vector<double>& fwdReset, const std::vector<double>& fwd1, 
		const ARM_GP_Matrix& mktvols1, const ARM_GP_Matrix& strikes1, const std::vector<double>& fwd2,
		const ARM_GP_Matrix& mktvols2, const ARM_GP_Matrix& strikes2, double v0, double kappa, double theta, bool calibtheta, 
		bool caliballShift, bool constrainedCorrel);


	double		v0() const	{return m_params1[0].v0();};
	double		kappa() const {return m_params1[0].kappa();};
	double		theta(int i) const {return m_params1[i].theta();};
	double		nu(int i) const {return m_params1[i].nu();};
	double		rho1(int i) const {return m_params1[i].rho();};
	double		rho2(int i) const {return m_params2[i].rho();};
	double		shift1(int i) const {return m_shift1[i];};
	double		shift2(int i) const {return m_shift2[i];};
	double		level1(int i) const {return m_params1[i].level();};
	double		level2(int i) const {return m_params2[i].level();};
	double		correl() const {return m_correl;};

// FIXMEFRED: mig.vc8 (22/05/2007 18:45:14): doesnt exist
	//double		theta(double reset);
	//double		nu(double reset);
	//double		rho1(double reset);
	//double		rho2(double reset);
	//double		shift1(double reset);
	//double		shift2(double reset);
	//double		level1(double reset);
	//double		level2(double reset);

	void		getData(double reset, double& theta, double& nu, double& rho1, double& rho2, double& shift1,
					double& shift2, double& level1, double& level2);

	int			nbFwdToCalibrate(double firstCapReset, double lastCapReset);
	void		Calibrate(double firstCapReset, double lastCapReset);
	void		SetIsCalibrated(double firstCapReset, double lastCapReset);

	friend class ARM_SmileCalibration_Spread2HestonCAP;

private:

	void		getInterpolatePoints(double reset, int& kinf, int& ksup);
};

class ARM_SmileCalibration_Spread2HestonCAP
{
private:
	std::vector<double>	m_resets;
	std::vector<double>	m_dfs;
	std::vector<double>	m_fwd1;
	std::vector<double>	m_fwd2;

	std::vector<double>		m_capprices;
	std::vector<double>		m_capstrikes;
	double				m_constrCapPrice;
	double				m_constrCapStrike;

	ARM_Spread2Heston_CAPCalibration_Params * m_params;

public:

	ARM_SmileCalibration_Spread2HestonCAP()
	{}

	virtual ~ARM_SmileCalibration_Spread2HestonCAP()
	{
	}

public:

	void	Init(const std::vector<double>& resets, const std::vector<double>& dfs, const std::vector<double>& fwd1, const std::vector<double>& fwd2,
				 const std::vector<double>& capprices, const std::vector<double>& capstrikes, double constrCapPrice, double constrCapStrike,
				 ARM_Spread2Heston_CAPCalibration_Params * params);


	void	Calibrate();

	double	cap(double strike);

private:
	void	SetParamsForFwdCalib(double p[]);
	void	SetConstParamForFwd(double x, int i);
	void	SetCorrel(double x) {m_params->m_correl = x;};
	void	SetParamsForSpreadCalib(double p[]);

	double	ImpliedVolFwd(double strike, int i);
	double	ImpliedVolSpread(double strike);

	void	CalibFwd();
	void	CalibFwdConstraint(int i);
	void	CalibSpreadConstraint();
	
	class CalibSpreadFunc : public ARM_LEVMARQFunc
	{
	private:
		ARM_SmileCalibration_Spread2HestonCAP * _this;

	public:
		CalibSpreadFunc(ARM_SmileCalibration_Spread2HestonCAP * calib)
		{
			_this = calib;
		}

	public:
		void operator()(double p[], double hx[], int m, int n, void * adata = NULL) const
		{
			_this->SetParamsForSpreadCalib(p);

			_this->CalibFwd();

			_this->CalibSpreadConstraint();

			for(int i = 0; i < _this->m_capstrikes.size(); i++)
			{
				hx[i] = _this->cap(_this->m_capstrikes[i]) - _this->m_capprices[i];
			}
		}
	};

	class CalibSpreadConstrFunc : public DoubleToDoubleFunc
	{
	private:
		ARM_SmileCalibration_Spread2HestonCAP * _this;

	public:
		CalibSpreadConstrFunc(ARM_SmileCalibration_Spread2HestonCAP * calib)
		{
			_this = calib;
		}

	public:
		double operator()(double x) const
		{
			_this->SetCorrel(x);

			return _this->cap(_this->m_constrCapStrike);
		}
	};

	friend class CalibSpreadFunc;
	friend class CalibSpreadConstrFunc;
};

void ARM_Spread2Heston_TOTEMCalibration(const std::vector<double>& TOTEMMaturities, const std::vector<double>& TOTEMStrikes,
										const std::vector<double>& TOTEMPrices,
										const std::vector<double>& FullSchedReset,
										const std::vector<double>& FullSchedAnnuity,
										const std::vector<double>& FullSchedFwd1, const std::vector<double>& FullSchedFwd2,
										const std::vector<double>& FwdCalibReset,
										const std::vector<double>& LongFwds, const ARM_GP_Matrix& LongVols, const ARM_GP_Matrix& LongK,
										const std::vector<double>& ShortFwds, const ARM_GP_Matrix& ShortVols, const ARM_GP_Matrix& ShortK,
										double v0, double kappa, double theta, bool locShifts,
										std::vector<double>& Mat, std::vector<double>& Theta, std::vector<double>& Nu, 
										std::vector<double>& LongRho, std::vector<double>& LongShift, std::vector<double>& LongLevel,
										std::vector<double>& ShortRho, std::vector<double>& ShortShift, std::vector<double>& ShortLevel,
										std::vector<double>& Correl);

CC_END_NAMESPACE()

#endif