
#ifndef _KPROCESS_H_
#define	_KPROCESS_H_

#include "kplatform.h"

#include <ostream>
#include <fstream>
#include <map>
#include <vector>
#include "kratecurve.h"  
#include "kdate.h"
#include "kmatrix.h"
#include "kfunc.h"
#include "kplatdep.h"
#include "kconstant.h"



class	KNormalDistribution
{
public:
	KNormalDistribution(void){}
	KNormalDistribution(const KMatrix<double> &covar): _covar(covar){}
	~KNormalDistribution(void){}
	//assume m1_i and m2_j are independent
	void	operator +=(const KNormalDistribution &m){ _covar += m._covar;}
	void	operator -=(const KNormalDistribution &m){ _covar += m._covar;}
	KNormalDistribution	operator +(const KNormalDistribution &m)const{ KNormalDistribution	ans(*this);
																	ans += m;
																	return ans;}
	KNormalDistribution	operator -(const KNormalDistribution &m)const{ KNormalDistribution	ans(*this);
																	ans -= m;
																	return ans;}

	KNormalDistribution	operator *(const KValarray<double> &a)const{KNormalDistribution	ans(*this);
																	ans *= a;
																	return ans;}
	KNormalDistribution	operator *(double a)const{KNormalDistribution	ans(*this);
																	ans *= a;
																	return ans;}

	void	operator *=(const KValarray<double> &a);
	void	operator *=(double a);
	friend	KNormalDistribution	operator * (const KValarray<double> &a, const KNormalDistribution &m)
																	{return m * a;}
	friend	KNormalDistribution	operator * (double a, const KNormalDistribution &m){return m * a;}

	//variance of sum(E_i*X_i)
	double	variance(const KValarray<double> &coeff)const;
	//varaince of sum X_i
	double	variance()const;
   	KInput	print()const{return _covar.print();}
	friend std::ostream & operator<<(std::ostream &out, const KNormalDistribution &m)
	{out<<m.print(); return out;}
private:
	//covaraince matrix
	KMatrix<double>	_covar;
};

//KNormalProcess represents the following normal process
// X = sum{X_i|i=1, ..., N} and dX_i = -b_i*X_i*dt + s_i(t)*dZ_i(t)
// <dZ_i, dZ_j> = r_ij*dt, s_i(t) is a step function.
//In theory, spot volatilities can be any functions, but to achieve
//computation effeciancy, we restrict to only step functions. 
class	KNormalProcess
{
public:
	KNormalProcess(void):_beta(0., 1), _corr(1, 1, 1.){}
	KNormalProcess(const KValarray<double> &beta,
				const KMatrix<double> &corr):
				_beta(beta), _corr(corr){}
	KNormalProcess(	const KValarray<double> &beta,
					const KMatrix<double> &corr,
					const TableFunc<double, KValarray<double> >	&vol):
					_beta(beta), _corr(corr), _vol(vol){}
	~KNormalProcess(void){}
	int	get_dim()const{return _beta.size();}
	KValarray<double>	get_mean_reversion()const{return _beta;}
	TableFunc<double, KValarray<double> > get_vol()const{return _vol;}
	//get the average spot vol for the period of [t1, t2], so that 
	//it gives the correct variance for each factor over any period, which
	//has [t1, t2] as a subperiod.
	//the potential pitfall is that, if the spot vol is very volatile in the period
	//of [t1, t2], the interpolation may not give the correct average correlation
	//of the end distribution.
	//the handle is mostly used in approximating spot vol on tree, which 
	//should minimize potential errors as the period is in general very short.
	KValarray<double>  get_avg_vol(double t1, double t2);
	KMatrix<double>	get_correlation()const{return _corr;}
	void	set_vol(const TableFunc<double, KValarray<double> > &vol){_vol = vol;}
	//D(t2) = D(t2, t1) + exp(-beta *(t2-t1)) * D(t1)
	//D(0) = 0
	KNormalDistribution	get_distribution(double t2,	double	t1 = 0.); 
	friend	std::ostream &operator<<(std::ostream &, const KNormalProcess &);
protected:
	KValarray<double>	&beta(){return _beta;}
	KMatrix<double>		&corr(){return _corr;}
private:
	//mean reversion
	KValarray<double>	_beta;
	//correlation matrix
	KMatrix<double>	_corr;
	//spot vol info: _vol[i].second is the spot vol between
	//_vol[i-1].first and _vol[i].first
	TableFunc<double, KValarray<double> > _vol;
};

//A (a,t) = (1 - exp (at)) / (a*t)
double	A (double,double);

//IRVol class: assuming instantaneous fwd rates follow normal or lognormal process
//in HJM framwork.
class	 KIRProcess:public KNormalProcess    
{
public:
	KIRProcess(void){}
	KIRProcess(const KRateCurve &tcurve,
				const KValarray<double> &beta,
				const KMatrix<double> &corr,
				double	power = 1.):_tcurve(tcurve),
				KNormalProcess(beta, corr){_power = power;} 
	
	~KIRProcess(void){}
	//calibrate spot vol for given benchmark bs vol.
	void	calibrate(const KMap(KDate, KIndexRate) &index,
					const KMap(KDate, double) &bs_vol,
					const KValarray<double> &weight);
	//compute BS Vol for any index
	double	get_bs_vol(const KDate&, const KIndexRate &index);
	//convert dates to time in double
	KValarray<double>  get_avg_vol(const KDate &d1, const KDate &d2)
					{return KNormalProcess::get_avg_vol(GetYearDiff(_tcurve.valueDate(), d1),
														GetYearDiff(_tcurve.valueDate(), d2));}
	KNormalDistribution	get_distribution(const KDate &d2, const KDate &d1 = 0)
					{return KNormalProcess::get_distribution(GetYearDiff(_tcurve.valueDate(), d2),
														GetYearDiff(_tcurve.valueDate(), d1));}
 	
private:
	//_power = 1, lognormal; _power = 0. normal
	double	_power;
	KRateCurve	_tcurve;
	double	indexVolCoeff(const KDate	&fwdDate,
							double	beta,
							const KIndexRate &index)const;
};



class	KEqProcess:public KNormalProcess    
{
public:
	KEqProcess(void){}
	~KEqProcess(void){}
	
	double	&spot_price(void){return _spot;}
	double	spot_price(void)const{return _spot;}
	KString	&dividend_type(void){return _divType;}
	KString	dividend_type(void)const{return _divType;}
	
	double	get_dividend(const KDate &d)const
			{KMap(KDate, double)::const_iterator  iter;  
			if((iter = _div.find(d)) != _div.end()) return iter->second;
			else	return 0.;}
	void	set_dividend(const KDate &d, double a){_div[d] = a;}

	KTimeLine	get_timeline(void)const;
	//t is years from value Date
	double	get_fwd_price(const KDate &fwd, const KRateCurve &);
	//calibrate spot vol for given benchmark bs vol.
	void	calibrate(const KDate &valueDate, const KMap(KDate, double) &bs_vol);
	double	get_bs_vol(const KDate&);
	//convert dates to time in double
	KValarray<double>  get_avg_vol(const KDate &d1, const KDate &d2)
					{return KNormalProcess::get_avg_vol(GetYearDiff(_valueDate, d1),
														GetYearDiff(_valueDate, d2));}
	KNormalDistribution	get_distribution(const KDate &d2, const KDate &d1 = 0)
					{return KNormalProcess::get_distribution(GetYearDiff(_valueDate, d2),
														GetYearDiff(_valueDate, d1));}
 	
	friend	std::ostream &operator<<(std::ostream &, const KEqProcess &);
private:
	KDate	_valueDate;
	//spot price on _valueDate
	double	_spot;     
	//dividend info
	KString	_divType;
	KMap(KDate, double)	_div;

};

#endif

