#include <math.h>
#include "kratecurve.h"
#include "kexception.h"
#include "konedimsolver.h"

KRateCurve::KRateCurve(KInput &input)
{
	KInputMap &c_map = input.get_map();
	KDate	startDate(c_map.get("VALUE_DATE"));
	KInputMap &curve_map = c_map["RATE_CURVE"].get_map();

	_tcurve = GtoNewTCurve( startDate, curve_map.size(),  1.,  ACT_365);
	if( _tcurve == NULL)
		throw KException("KRateCurve constructor failed!");
	KInputMap::iterator iter;
	KMap(KDate, double)::iterator iter1;

	int	i;
	double	rate;
	//for sorting purpose
	KMap(KDate, double)	curve;
	for(iter = curve_map.begin(); iter!= curve_map.end(); iter++)
	{
		KDate	tmp(iter->first);
		curve_map.get(iter->first) >>rate;
		curve[tmp] = rate;
	}
	for(iter1 = curve.begin(), i=0; iter1!= curve.end(); iter1++, i++)
	{
		_tcurve->fArray[i].fDate = iter1->first;
		_tcurve->fArray[i].fRate = iter1->second;
	}
}

KRateCurve::KRateCurve(const KRateCurve& rhs)
{
	if (rhs._tcurve == NULL) _tcurve = NULL;	
	else _tcurve = GtoCopyCurve (rhs._tcurve);
}




KRateCurve::KRateCurve(TDate baseDate, TDate* dates, double* rates, int numPts,double basis,long dayCountConv)
{
	TCurve* temp = GtoMakeTCurve(baseDate,
								 dates,
								 rates,
								 numPts,
								 basis,
								 dayCountConv);
	if(temp== NULL)
	{
		throw KException(" TCurve build failed.");
	}
	else
	{
		_tcurve = temp;
	}

}



KRateCurve& KRateCurve::operator=(const KRateCurve& rhs)
{
	if(this != &rhs)
	{
		destroy();
		if (rhs._tcurve)
			_tcurve = GtoCopyCurve (rhs._tcurve);
		else 
			_tcurve = NULL;
	}
	return *this;
}
KInput KRateCurve::print()const
{
	KInputMap input;
	input.store("Value_Date", KDate(_tcurve->fBaseDate).str());

	KInputMap	curve;
	for(int i = 0; i < _tcurve->fNumItems;	i++)
	{
		KInput	tmp;
		_tcurve->fArray[i].fRate >> tmp;
		curve[KDate(_tcurve->fArray[i].fDate).str()] = tmp;
	}
	input["Rate_Curve"] = KInput(curve);
	return input;
}

void KRateCurve::shiftCurve(double x)
{
	for(int i = 0; i< size();  i++)
	{
	
		_tcurve->fArray[i].fRate += x;
	}
}

double	KRateCurve::get_rate(int interpType, const KDate &endDate, const KDate &startDate) const
{
	KDate	d = valueDate();	
	if(startDate >= endDate || endDate < d)
		throw KException("start date is after end date! get_rate failed");
	double	rate1, rate2;
	if( GtoInterpRate(endDate, _tcurve, interpType, &rate2) == FAILURE)
		throw KException("get_rate failed!");
	if(startDate <= d)
		return rate2;
	else
	{
		if( GtoInterpRate(startDate, _tcurve, interpType, &rate1) == FAILURE)
			throw KException("get_rate failed!");
		double	t1 = GetYearDiff(d, startDate);
		double	t2 = GetYearDiff(d, endDate);
		return pow( pow(1. + rate2, t2)/pow(1. + rate1, t1), 1./(t2 - t1)) - 1.;
	}
	return 0.;
}
double	KRateCurve::get_pv(int interpType, const KDate &d) const
{
	if(d < valueDate())
		throw KException("date is before value date! get_pv failed");
	double	pv;
	if( GtoInterpPV(d, _tcurve, interpType, &pv) == FAILURE)
		throw KException("get_rate failed!");
	return pv;
}

KIndexRate	KIndexRate::operator + (double s) const
{
	KIndexRate tmp(*this);
	tmp._spread += s;
	return tmp;
}

KIndexRate	KIndexRate::operator * (double w) const
{
	KIndexRate tmp(*this);
	tmp._weight *= w;
	return tmp;
}


	//present value of all cash flows from valuedate
double	KCashFlow::get_pv(const KRateCurve &zero)const
{
	KDate	valueDate = zero.valueDate();
	KCashFlow::const_iterator	iter = lower_bound(valueDate);
	double	pv = 0;
	for(;iter != end(); iter++)
		pv += iter->second * zero.get_pv(FLAT_FORWARD, iter->first);
	return pv;
}
//present value of all cash flows from valuedate
//discounted with flat yield 
double	KCashFlow::get_pv(const KDate &valueDate, double y)const
{
	KCashFlow::const_iterator	iter = lower_bound(valueDate);
	double	pv = 0;
	for(;iter != end(); iter++)
		pv += iter->second / pow(1. + y, GetYearDiff(valueDate,iter->first));
		
	return pv;
}
//flat yield for given pv
double	KCashFlow::get_yield(const KDate &valueDate, double p, double guess)
{
	//initialize start date
	_startDate = valueDate;
	KSecantRoot<KCashFlow>	yield(0.001);
	return yield.solve(*this,guess, p);
}
