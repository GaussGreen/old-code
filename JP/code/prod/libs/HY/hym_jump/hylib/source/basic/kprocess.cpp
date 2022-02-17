
#include <math.h>
#include "kprocess.h"
#include "kutility.h"
#include "kexception.h"
#include "kvalarray.h"

double	KNormalDistribution::variance(const KValarray<double> &coeff)const
{
	if(coeff.size() != _covar.row_size())
		KException("wrong size!");
	double	var = 0.;
	int	i, j, size = coeff.size();
	for(i = 0; i<size; i++)
	{
		var += coeff[i] * coeff[i] * _covar[i][i];
		for(j = 0; j<i; j++)
			var += 2* coeff[i] * coeff[j] * _covar[i][j];
	}
	return var;
}
double	KNormalDistribution::variance()const
{
	double	var = 0.;
	int	i, j, size = _covar.row_size();
	for(i = 0; i<size; i++)
	{
		var +=  _covar[i][i];
		for(j = 0; j<i; j++)
			var += 2 * _covar[i][j];
	}
	return var;
}


void	KNormalDistribution::operator *=(const KValarray<double> &a)
{
	int	i, j, size = _covar.row_size();
	for(i=0; i<size; i++)
		for(j=i; j<size; j++)
			_covar[i][j] *= a[i]*a[j];
	for(i=0; i<size; i++)
		for(j=0; j<i; j++)
			_covar[i][j] = _covar[j][i];
}
void	KNormalDistribution::operator *=( double a)
{
	int	i, j, size = _covar.row_size();
	for(i=0; i<size; i++)
		for(j=i; j<size; j++)
			_covar[i][j] *= a * a;
	for(i=0; i<size; i++)
		for(j=0; j<i; j++)
			_covar[i][j] = _covar[j][i];
}

KValarray<double>  KNormalProcess::get_avg_vol(double t1, double t2)
{
	if(t2 < t1)
		throw KException("Wrong order of time!");
	
	int	i, size = _beta.size();
	double	t, T;
	TableFunc<double, KValarray<double> >::const_iterator	iter, iter1, iter2;
	
	iter1 = _vol.lower_bound(t1);
	iter2 = _vol.lower_bound(t2);

	T = t1;
	KValarray<double>	avg_vol(0.,size);
	for(iter = iter1; iter != iter2; iter++)
	{
		t = T;
		T = iter->first;
		//get function value for the period, use average to avoid
		//potential end point problem for step function. 
		KValarray<double> &spot = _vol((t+T)/2);
		for(i=0; i<size; i++)
			avg_vol[i] += exp( -2 * _beta[i] * (t2 - T)) 
							* A( 2 * _beta[i], T - t) * ( T - t)
							* spot[i] * spot[i];
	}
	//last interval	(iter = iter1)
	t = T;
	T = t2;
	KValarray<double> &spot = _vol((t+T)/2);		
	for(i=0; i<size; i++)
		avg_vol[i] += A( 2 * _beta[i], T - t) * ( T - t)
						* spot[i] * spot[i];

	for(i=0; i<size; i++)
	{
		avg_vol[i] /= A( 2 * _beta[i], t2 - t1) * ( t2 - t1);
		avg_vol[i] = sqrt(avg_vol[i]);
	}
	return  avg_vol;	
}

KNormalDistribution	KNormalProcess::get_distribution(double	t2, double t1)
{
	if(t2 < t1)
		throw KException("wrong order of time!");
	//if t1 is negative, initialize at 0
	if(t1 < 0)
		t1 = 0;
	int	i, j, size = _beta.size();
	double	t, T;
	TableFunc<double, KValarray<double> >::const_iterator	iter, iter1, iter2;
	
	iter1 = _vol.lower_bound(t1);
	iter2 = _vol.lower_bound(t2);

	T = t1;
	KMatrix<double>	extra(size, size, 0.);
	for(iter = iter1; iter != iter2; iter++)
	{
		t = T;
		T = iter->first;
		//get function value for the period, use average to avoid
		//potential end point problem for step function. 
		KValarray<double> &spot = _vol((t+T)/2);
		for(i=0; i<size; i++)
			for(j=i; j<size; j++)
				extra[i][j] += exp( -(_beta[i] + _beta[j]) * (t2 - T)) 
							* A( (_beta[i] + _beta[j]), T - t) * ( T - t)
							* _corr[i][j] * spot[i] * spot[j];
	}
	//last interval	(iter = iter1)
	t = T;
	T = t2;
	KValarray<double> &spot = _vol((t+T)/2);		
	for(i=0; i<size; i++)
		for(j=i; j<size; j++)
			extra[i][j] += A( (_beta[i] + _beta[j]), T - t) * ( T - t)
							* _corr[i][j] * spot[i] * spot[j];

	for(i=0; i<size; i++)
		for(j=0; j<i; j++)
			extra[i][j] = extra[j][i];
	
	return  KNormalDistribution(extra);
}

std::ostream &operator<<(std::ostream &out, const KNormalProcess &process)
{
	out<<"Mean Reversion:\t"<<process._beta<<std::endl;
	out<<"Correlation Matrix"<<std::endl;
	out<<process._corr<<std::endl;
	out<<"Spot volatilities"<<std::endl;
	TableFunc<double, KValarray<double> >::const_iterator iter;
	for(iter = process._vol.begin(); iter != process._vol.end(); iter++)
		out<<iter->first<<"\t"<<iter->second<<std::endl;
	return out;
}

//*****  A  ******************************************************************
//       exponential decay function: A (a,t) = (1 - exp (at)) / (a*t)
//
double  A ( double	a, double    t)
{
	double	x;
	if (fabs (a * t) < 0.00000001)      
		return (1.);
	x = (1. - exp (- a * t)) / (a * t);
	return (x);
} 



double		KIRProcess::indexVolCoeff(const KDate	&fwdDate,
										double	beta,
										const KIndexRate &index) const
{
    double  expZero, zeroPrice, previousZero;           
	double	annuity, B, C, D,couponPeriod;		
	KDate	startDate, endDate;
	int		i, numCoupon = (int) (index.matInMonth() * index.freq() / 12);

    annuity = B = C = 0.;
	startDate = fwdDate;
	previousZero = _tcurve.get_pv(FLAT_FORWARD, fwdDate);
  	expZero = 	previousZero;
	for (i = 1; i <= numCoupon; i++)
	{
     	endDate = fwdDate + (int)(12/index.freq() * i) * ONE_MONTH;
		couponPeriod = GetYearDiff(startDate, endDate, index.dayCountConv());
		zeroPrice = _tcurve.get_pv(FLAT_FORWARD, endDate);
		annuity += couponPeriod * zeroPrice;
		//distribution switch
		if(_power == 1.)
			D = log(previousZero/zeroPrice) ;
		else if(_power == 0.)
			D = GetYearDiff(startDate, endDate);
		else
			throw KException("Not yet implemented!");
		C += exp (-beta * GetYearDiff(fwdDate, startDate)) 
				* D	* A (beta, GetYearDiff(startDate, endDate));
		B += C * couponPeriod * zeroPrice;
		//prepare for the next coupon
		previousZero = zeroPrice;
		startDate = endDate;
	}  
	B *=(expZero - zeroPrice) / annuity;   
	B += C * zeroPrice;
	B /= expZero - zeroPrice;
	return(B);
}

double	KIRProcess::get_bs_vol(const KDate &fwdDate, const KIndexRate &index)
{
	int	i, size = get_dim();
	double	t = GetYearDiff(_tcurve.valueDate(), fwdDate);
	KValarray<double>	B(size);
	KValarray<double>	&beta = get_mean_reversion(); 
	KNormalDistribution	td = get_distribution(t);
	for(i = 0; i<size; i++)
		B[i] = indexVolCoeff(fwdDate,beta[i], index);
	double	vol = td.variance(B)/t;
	return sqrt(vol);
}

void	KIRProcess::calibrate(const KMap(KDate, KIndexRate) &index,
								const KMap(KDate, double) &bs_vol,
								const KValarray<double> &weight)
{		
	int	j, size = get_dim();
	double	t1, t2, spotVol;
	KValarray<double>	B(size);
	KValarray<double>	mean_reversion(size);
	KValarray<double>	beta = get_mean_reversion();
	KMap(KDate, KIndexRate)::const_iterator	iter;
	TableFunc<double, KValarray<double> >::iterator	iter1;
	KMap(KDate, double)::const_iterator	iter2;
	KDate	valueDate = _tcurve.valueDate();
	
	//set spot vol
	KMap(double, KValarray<double> ) vol;
	for(iter = index.begin(); iter != index.end(); iter++)
		vol[GetYearDiff(valueDate, iter->first)] = weight;
	TableFunc<double, KValarray<double> >	spot_vol(vol.begin(), vol.end(),
											KInterp<double, KValarray<double> >::r_step());
	set_vol(spot_vol);

	KNormalDistribution	td(KMatrix<double>(size, size, 0.));
	KNormalDistribution	tmpDist;

	t2 = 0;
	for(iter = index.begin(), iter1 = spot_vol.begin(), iter2 = bs_vol.begin(); 
			iter != index.end(); iter++, iter1++, iter2++)
	{
		t1 = t2;
		t2 = GetYearDiff(valueDate, iter->first);
		for(j=0; j<size; j++)
			mean_reversion[j] = exp(-beta[j]* (t2 - t1));
		td *= mean_reversion;
		tmpDist = get_distribution(t2, t1);
	
		for(j = 0; j<size; j++)
			B[j] = indexVolCoeff( iter->first,beta[j], iter->second);
		//total variance
		spotVol = t2 * iter2->second * iter2->second;
		//incremental variance
		spotVol -= td.variance(B);
		spotVol /= tmpDist.variance();
		spotVol = sqrt(spotVol);
		iter1->second *= spotVol;
		//prepare for the next calculation
		td +=  spotVol * tmpDist;
	}
	//set spot vol
	set_vol(spot_vol);
}
	//calibrate spot vol for given benchmark bs vol.
void	KEqProcess::calibrate(const KDate &valueDate, const KMap(KDate, double) &bs_vol)
{
	_valueDate = valueDate;
	//set spot vol
	double	t1, t2, v1, v2, v;
	KMap(KDate, double)::const_iterator	iter;
	KMap(double, KValarray<double> ) vol;
	v2 = 0.;
	t2 = 0.;
	for(iter = bs_vol.begin(); iter != bs_vol.end(); iter++)
	{
		t1 = t2;
		v1 = v2;
		t2 = GetYearDiff(valueDate, iter->first);
		v2 = t2 * iter->second * iter->second;
		if( t2 == t1)
			v = 0;
		else
		{
			if( v2 < v1)
				throw KException("Error in total variance!"); 
			v = sqrt((v2 - v1)/(t2 - t1));
		}
		KValarray<double> spot(v,  1);
		vol[t2] = spot;
	}
	TableFunc<double, KValarray<double> >	spot_vol(vol.begin(), vol.end(),
											KInterp<double, KValarray<double> >::r_step());
	set_vol(spot_vol);
}

double	KEqProcess::get_bs_vol(const KDate &fwdDate)
{
	double	t = GetYearDiff(_valueDate, fwdDate);
	KNormalDistribution	td = get_distribution(t);
	double	vol = td.variance()/t;
	return sqrt(vol);
}

KTimeLine	KEqProcess::get_timeline()const
{
	KTimeLine	ans;
	KMap(KDate, double)::const_iterator	iter;
	for(iter = _div.begin(); iter!= _div.end(); iter++)
	{
		KTimePoint	t;
		t[__divDate] = 1;
		t[__eqDividend] = iter->second;
		ans.insert(iter->first, t);
	}
	return ans;
}
double	KEqProcess::get_fwd_price(const KDate &date, const KRateCurve &zero)
{
	double	fwd = _spot;
	if(_divType != __NonDiv)
	{
		KDate	valueDate = zero.valueDate();
		KMap(KDate, double)::const_iterator	iter, iter1, iter2;
		if(_divType == __ContYieldDiv)
		{	
			double	length;
			iter1 = _div.upper_bound(valueDate);
			if( iter1 == _div.end())
				throw KException("No dividend specification for the current period!");
			iter2 = _div.upper_bound(date);

			if(iter1 == iter2)
			{	 
				length = GetYearDiff(valueDate, date);
				fwd /= pow(1. + iter1->second, length);
			}
			else
			{
				length = GetYearDiff(valueDate, iter1->first);
				fwd /= pow(1. + iter1->second, length);
				iter = iter1;
				iter1++;
				while( iter1 != iter2)
				{
					length = GetYearDiff(iter->first, iter1->first);
					fwd /= pow(1. + iter1->second, length);
					iter = iter1;
					iter1++;
				}
				if(iter1 == iter2 && iter->first < date)
				{
					length = GetYearDiff(iter->first, date);
					fwd /= pow(1. + iter->second, length);
				}
			}
		}
		else if(_divType == __FixAmountDiv)
		{
			for(iter = _div.begin(); iter!= _div.end() && iter->first <= date; iter++)
				if(iter->first > valueDate)
					fwd -= iter->second * zero.get_pv(FLAT_FORWARD, iter->first);
		}
		else
		{
			for(iter = _div.begin(); iter!= _div.end() && iter->first <= date; iter++)
				if(iter->first > valueDate)		
					fwd *= 1. - iter->second;
		}
	}
	fwd /= zero.get_pv(FLAT_FORWARD,date);
	return fwd;
}

std::ostream &operator<<(std::ostream &out, const KEqProcess &s)
{
	out<<"Value Date\t"<<s._valueDate<<std::endl;
	out<<"SpotPrice\t"<<s._spot<<std::endl;
	out<<"DividendType\t"<<s._divType<<std::endl;
	KMap(KDate, double)::const_iterator	iter;
	for(iter = s._div.begin(); iter!= s._div.end(); iter++)
		out<<iter->first<<"\t"<<iter->second<<std::endl;
	out<<"Underlying process:"<<std::endl;
	out<< KNormalProcess(s);
	return out;
}
