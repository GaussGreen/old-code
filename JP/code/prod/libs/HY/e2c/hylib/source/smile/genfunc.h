// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 1999 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
//
// 1.0 10/19/99 Neil Yang
// 
//

#if ! defined _GENFUNC_
#define _GENFUNC_

#include "kplatform.h"

#include <exception>
/* #include "systemh.h "  */
#include "math.h"
#include "General/General.h"
#include "RootFindBrent.h"
#include "kexception.h"

/* Created by Neil Yang 10/18/99 */

class SoverFunc;

class BaseFunction :public CM::Object
{
	double m_amplitude;
	double m_stddev;
	double m_time;
public:
	BaseFunction(): m_amplitude(1), m_stddev(1){}
	BaseFunction(double amplitude, double stddev, double time):m_amplitude(amplitude),m_stddev(stddev) ,m_time(time) {}
	void set_amplitude(double amplitude)
	{
		m_amplitude=amplitude;
	}

	void set_stddev(double stddev)
	{
		m_stddev=stddev;
	}
	void set_time(double time)
	{
		m_time=time;
	}
	double get_amplitude() const {return m_amplitude;}
	double get_stddev() const {return m_stddev;}
	double get_time() const {return m_time;}

	//double calibrate(double x);

	//double integral(

	double deriv(double x) const;
	double inverse(double y) const;

	virtual double operator()(double x) const=0;

};


class SoverFunc
{
	double m_value;
	const BaseFunction* m_function;

public:
	SoverFunc(const BaseFunction* func,double value=0) : m_function(func), m_value(value) {}
	double operator() (double x) const 
	{
		
		double temp =(*m_function)(x);

		
		return (temp - m_value);
	}
};




class LogNormalFunction : public BaseFunction
{
	

public:
 
	double operator()(double x) const
	{
		get_amplitude();
		return get_amplitude()*exp(get_stddev()*x);
	}

};


class NormalLogNormalFunction : public BaseFunction
{
	

public:

	double operator()(double x) const
	{
		return log(1+get_amplitude()*exp(get_stddev()*x));
	}

};

class NormalLogNormalFunction2 : public BaseFunction
{
	

public:

	double operator()(double x) const
	{
		return get_amplitude()*log(1+exp(get_stddev()*x));
	}

};

class DualLogNormalFunction : public BaseFunction
{
	double m_turnPoint;
	double m_vol2;
	double m_a;

public:
	DualLogNormalFunction() : m_turnPoint(0), m_vol2(0), m_a(0){}
	DualLogNormalFunction(double turnPoint, double vol2, double a) : m_turnPoint(turnPoint), m_vol2(vol2), m_a(a){}
	double operator()(double x) const
	{
		double stddev2=m_vol2*sqrt(get_time());
		double e1=get_amplitude()*exp(get_stddev()*x);
		double e2 = get_amplitude()*exp(stddev2*x);
		double ec = exp(get_stddev()*x)/exp(get_stddev()*m_turnPoint);
		double weight = 1/(1+pow(ec,m_a));

		return e2*(1-weight)+e1*weight;
	}

};

class AnnuityIndex : public DualLogNormalFunction
{
	double m_maturity;
	double m_coupon;

public:
	AnnuityIndex(double maturity, double coupon, double turnPoint, double vol2, double a) : 
	  DualLogNormalFunction(turnPoint, vol2, a),m_maturity(maturity), m_coupon(coupon)
	  {
		
	  }

		double operator() (double x) const
		{
			double spread = DualLogNormalFunction::operator()(x);
			double r =0.0;

			return (1-1/pow((1+r+spread),m_maturity))/(r+spread);
		}

};

class AssetToStockMapping1 : public BaseFunction
{
	double   m_x1;
	double   m_level;
	double   m_debtPerShare;

	
public:
	AssetToStockMapping1(double x1, double level, double debtPerShare =1) 
	{
		if(level < x1/(1+x1))
		{
			throw KException("L<x/1+x.\n");
		}
		m_x1 = x1;
		m_level = level;
		m_debtPerShare = debtPerShare; 
	}

	double get_debtPerShare(){return m_debtPerShare;}
	double get_lim(){return m_level;}
	double get_x1(){return m_x1;}

	double operator() (double v) const
	{
		if (v> m_level)
		{
			return (v-1+ (1 - m_level)*pow(m_level/v,m_x1)) *m_debtPerShare;
		}
		else
		{
			return 0.0;
		}
	}

	double assetVolToStockVol(double stock, double assetVol)
	{
		double asset = inverse(stock);
		double stockVol = deriv(asset)*asset/stock*assetVol;
		return stockVol;
	}

	double stockVolToAssetVol(double stock, double stockVol)
	{
		double asset = inverse(stock);
		double assetVol = stock/asset/deriv(asset)*stockVol;

		return assetVol;
	}
};

class DummyMapping : public BaseFunction
{
	

	
public:
	DummyMapping() {}

	double operator() (double v) const
	{
		return v;
	}
};
#endif



