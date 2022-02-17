#ifndef _ICM_CUBIC_SPLINE_H_
#define _ICM_CUBIC_SPLINE_H_

#include "ARMKernel\glob\firsttoinc.h"
#include "ICMKernel\util\icm_macro.h"
#include "ARMKernel\glob\linalg.h"
#include "ICMKernel\glob\icm_enums.h"
#include "ICMKernel\util\icm_rootfinder1D.h"
#include "ICMKernel\util\icm_RootFinderND.h"
#include <vector>

/////////////////////////////
//
//		méthode d'interpolation par spline cubique
//
/////////////////////////////
#define DEPS 1.e-9

/**
static void __stdcall objfun_EvaluateDerivSecond(Integer n, double x[], double *objf,
                   double g[], Nag_Comm *comm);

static void __stdcall objfun_EvaluateC4cond(Integer n, double x[], double *objf,
                   double g[], Nag_Comm *comm);
**/ 

typedef enum ENUM_SPLINE_EMMODE {
		_emNone=0,
        _emC1=1,
		_emC2=2,
		_emC1_alpha=3,
		_emOptimized=4,
		_emOrder4=5,
		_emOrder4Convex=6,
		_emOptimizedC4=7,
		_emShumaker=8,
} eSPLINE_EMMODE; 

typedef enum ENUM_EMMODE {
        _emConstant,
		_emLinear,
		_emLinearTanh
} eEMMODE; 

class CCubicspline
{
public:
	vector<double>	m_alpha;
	double *		m_xvalues;
	double *		m_yvalues;
	int				m_size;

	eEMMODE			m_emmode;

	double			m_leftSlope;
	double			m_rightSlope;

	bool			m_needtorebuild;
	double *		m_mk;
	double *		m_ak;
	double *		m_bk;
	double *		m_ck;
	double *		m_dk;
	double *		m_ek;
	double *		m_fk;

	eSPLINE_EMMODE	m_Condition;
	int				m_cal_size;
	bool			m_calibrate;

public:
	CCubicspline(void)
	{
		Init();
	}

	void Init(void)
	{
		m_alpha.clear();
		m_xvalues		= 0;
		m_yvalues		= 0;
		m_size			= 0;

		m_emmode		= _emConstant;
		m_leftSlope		= 0.;
		m_rightSlope	= 0.;

		m_needtorebuild	= true;
		m_mk			= 0;
		m_ak			= 0;
		m_bk			= 0;
		m_ck			= 0;
		m_dk			= 0;
		m_ek			= 0;
		m_fk			= 0;

		m_Condition	= _emNone;
		m_cal_size		= 0;
		m_calibrate	=true;
	}

	~CCubicspline(void)
	{
		destroy();
	}

	void destroy()
	{
		if (m_xvalues) delete m_xvalues;
		if (m_yvalues) delete m_yvalues;
		if (m_mk) delete m_mk;
		if (m_ak) delete m_ak;
		if (m_bk) delete m_bk;
		if (m_ck) delete m_ck;
		if (m_dk) delete m_dk;
		if (m_ek) delete m_ek;
		if (m_fk) delete m_fk;

		m_xvalues		= 0;
		m_yvalues		= 0;
		m_size			= 0;

		m_emmode		= _emConstant;
		m_leftSlope		= 0.;
		m_rightSlope	= 0.;

		m_needtorebuild	= true;
		m_mk			= 0;
		m_ak			= 0;
		m_bk			= 0;
		m_ck			= 0;
		m_dk			= 0;
		m_ek			= 0;
		m_fk			= 0;

		m_Condition	= _emNone;
	}

	CCubicspline& operator = (const CCubicspline& spline)
	{
		if(&spline == this) return *this;

		destroy();
		m_alpha = spline.m_alpha;
		m_calibrate = spline.m_calibrate;

		m_size	= spline.m_size;

		if(m_size == 0) return *this;

		m_xvalues	= new double [m_size];
		m_yvalues	= new double [m_size];

		for(int k = 0; k < m_size; m_xvalues[k] = spline.m_xvalues[k], m_yvalues[k] = spline.m_yvalues[k++]);

		m_emmode		= spline.m_emmode;
		m_leftSlope		= spline.m_leftSlope;
		m_rightSlope	= spline.m_rightSlope;

		m_needtorebuild	= spline.m_needtorebuild;

		if(m_needtorebuild) return *this;

		m_ak			= new double [m_size];
		m_bk			= new double [m_size];

		for(k = 0; k < m_size; m_ak[k] = spline.m_ak[k], m_bk[k] = spline.m_bk[k++]);

		m_Condition	= spline.m_Condition;

		if(m_Condition)
		{
			m_ck	= new double [m_size];
			m_dk	= new double [m_size];
			m_ek	= new double [m_size];
			m_fk	= new double [m_size];
			m_mk	= new double [m_size];

			for(k = 0; k < m_size; m_mk[k] = spline.m_mk[k++]);

			for(k = 0; k < m_size; m_ck[k] = spline.m_ck[k], m_dk[k] = spline.m_dk[k],m_ek[k] = spline.m_ek[k],m_fk[k] = spline.m_fk[k++]);
		}
		else
		{
			m_mk	= new double [m_size];

			for(k = 0; k < m_size; m_mk[k] = spline.m_mk[k++]);
		}

		return *this;
	}

	CCubicspline(const double * xvalues,const double * yvalues,int size,eEMMODE emmode,eSPLINE_EMMODE condition = _emNone, double leftSlope = 0.,
			double rightSlope = 0.,vector<double>& alpha = vector<double>(),bool calibrate = true)
	{
		Init();
		build(xvalues,yvalues,size,emmode,condition,leftSlope,rightSlope,alpha,calibrate);
	}

public:
	// construction
	void build(const double * xvalues,const double * yvalues,int size,eEMMODE emmode,
		eSPLINE_EMMODE condition = _emNone, double leftSlope = 0.,double rightSlope = 0.,
		vector<double>& alpha = vector<double>(),bool calibrate = true);

	// interpolation / extrapolation
	double operator()(double x);

	// pente
	double derivative(double x);

private:
	// construction
	void buildSpline();

	void findMk(double * fk, double * hk);

	void buildSplineWithC1Condition();
	void buildSplineWithC2Condition();
	void buildSplineOptimized();
	void buildSplineOptimizedMulti();
	void buildSplineOrder4();
	void buildSplineOrder4CondDerivs();
	void buildSplineOptimizedMultiC4();
	void buildSplineShumaker();

	// extrapolation
	double leftExtrapol(double x);

	double rightExtrapol(double x);

	double EvaluateDerivSecond(const double& x);
};

inline double CCubicspline::operator ()(double x)
{
	if(m_size == 0) return CREDIT_DEFAULT_VALUE;

	if(m_needtorebuild) buildSpline();

	if(x < m_xvalues[0])		return leftExtrapol(x);
	if(x > m_xvalues[m_size-1])	return rightExtrapol(x);

	for(int k = 0; k < m_size - 1; k++)
	{
		if(x > m_xvalues[k] - DEPS && x < m_xvalues[k+1] + DEPS)
		{
			if ((m_Condition==_emOrder4)||(m_Condition==_emOrder4Convex))
			{
				return m_ek[k] * pow(x - m_xvalues[k],4.) + m_ak[k] * pow(x - m_xvalues[k],3.) + m_bk[k] * SQR(x - m_xvalues[k]) + m_ck[k] * (x - m_xvalues[k]) + m_dk[k];
			}
			else if(m_Condition == _emShumaker)
			{
				//Test : quelle partie de la decomposition de l'intervalle
				if (x <= m_mk[k+1])
					return m_ak[k+1] + m_bk[k+1] * (x - m_xvalues[k]) + m_ck[k+1] * (x - m_xvalues[k]) * (x - m_xvalues[k]);
				else
					return m_dk[k+1] + m_ek[k+1] * (x - m_mk[k+1]) + m_fk[k+1] * (x - m_mk[k+1]) * (x - m_mk[k+1]);
			}
			else if(m_Condition)
			{
				return m_ak[k] * pow(x - m_xvalues[k],3.) + m_bk[k] * SQR(x - m_xvalues[k]) + m_ck[k] * (x - m_xvalues[k]) + m_dk[k];
			}
			else
			{
				return m_mk[k+1] * pow(x - m_xvalues[k],3.) / (6. * (m_xvalues[k+1] - m_xvalues[k]))
					+ m_mk[k] * pow(m_xvalues[k+1] - x,3.) / (6. * (m_xvalues[k+1] - m_xvalues[k]))
					+ m_ak[k] * (x - m_xvalues[k]) + m_bk[k];
			}
		}
	}

	return CREDIT_DEFAULT_VALUE;
}

inline double CCubicspline::derivative(double x)
{
	if(m_size == 0) return 0.;

	if(m_needtorebuild) buildSpline();

	if(x < m_xvalues[0])
	{
		return m_emmode == _emLinear ? m_ak[0] : (leftExtrapol(x + 0.0001) - leftExtrapol(x - 0.0001)) / 0.0002;
	}

	if(x > m_xvalues[m_size-1])
	{
		return m_emmode == _emLinear ? m_ak[m_size-2] : (rightExtrapol(x + 0.0001) - rightExtrapol(x - 0.0001)) / 0.0002;
	}

	for(int k = 0; k < m_size - 1; k++)
	{
		if(x > m_xvalues[k] - DEPS && x < m_xvalues[k+1] + DEPS)
		{
			if(m_Condition)
			{
				return 3. * m_ak[k] * SQR(x - m_xvalues[k]) + 2. * m_bk[k] * (x - m_xvalues[k]) + m_ck[k];
			}
			else
			{
				return m_mk[k+1] * 3. * SQR(x - m_xvalues[k]) / (6. * (m_xvalues[k+1] - m_xvalues[k]))
					+ m_mk[k] * 3. * SQR(m_xvalues[k+1] - x) / (6. * (m_xvalues[k+1] - m_xvalues[k]))
					+ m_ak[k];
			}
		}
	}

	return 0.;
}

inline double CCubicspline::leftExtrapol(double x)
{
	switch(m_emmode)
	{
	case _emLinear:
		return m_ak[0] * (x - m_xvalues[0]) + m_bk[0];

	case _emLinearTanh:
		if(m_ak[0] < 0.)
			return m_ak[0] * (x - m_xvalues[0]) + m_bk[0];
		else
		{
			double slope = m_yvalues[0] == 0. ? 0. : m_ak[0] / m_yvalues[0];
			return m_yvalues[0] * (1. + tanh(slope * (x - m_xvalues[0])));
		}

	default: // constant
		return m_yvalues[0];
	}
}

inline double CCubicspline::rightExtrapol(double x)
{
	switch(m_emmode)
	{
	case _emLinear:
		return m_yvalues[m_size-1] + m_ak[m_size-2] * (x - m_xvalues[m_size-1]);

	case _emLinearTanh:
		if(m_ak[m_size-2] > 0)
			return m_yvalues[m_size-1] + m_ak[m_size-2] * (x - m_xvalues[m_size-1]);
		else
		{
			double slope = m_yvalues[m_size-1] == 0. ? 0. : m_ak[m_size-2] / m_yvalues[m_size-1];
			return m_yvalues[m_size-1] * (1. + tanh(slope * (x - m_xvalues[m_size-1])));
		}
		
	default: // constant
		return m_yvalues[m_size-1];
	}
}

class _spline_context
{
public:

	double* m_xvalues;
	double* m_yvalues;
	int		m_size;
	eEMMODE	m_emmode;
	eSPLINE_EMMODE	m_Condition;
	double	m_leftSlope;
	double	m_rightSlope;

	_spline_context(const _spline_context& in)
	{
		m_xvalues=in.m_xvalues;
		m_yvalues=in.m_yvalues;
		m_size=in.m_size;
		m_emmode=in.m_emmode;
		m_Condition=in.m_Condition;
		m_leftSlope=in.m_leftSlope;
		m_rightSlope=in.m_rightSlope;
	}

	void reset()
	{
	m_xvalues = NULL;
	m_yvalues = NULL;
	m_size=0;
	m_emmode=_emLinear;
	m_Condition=_emNone;;
	m_leftSlope=0.;
	m_rightSlope=0.;
	}

	_spline_context()
	{reset();}
};

#endif