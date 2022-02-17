
#ifndef	_K_RATECURVE_H_
#define	_K_RATECURVE_H_
#include <ostream>
#include <fstream>
#include "kinput.h"
#include "kdate.h"
#include "kcurve.h"
#include "General/General.h"
#include "fltrate.h"  /* GtoForwardRate */

extern "C" {
#include "tcurve.h"
#include "interp.h"

}

//interpolation methods
const	int	FLAT_FORWARD   = GTO_FLAT_FORWARDS;
const	int	PARABOLIC_FORWARD   = GTO_PARABOLIC_FORWARDS;
const   int LINEAR =GTO_LINEAR_INTERP;

//KRateCurve is not inheritated from KCurve.
//It only is a wrapper to GtoZeroCurve
class KRateCurve:public CM::Object
{
public:
	KRateCurve(){_tcurve = NULL;}
	KRateCurve(TCurve* tc) {_tcurve = tc;}
	KRateCurve(KInput &);
	KRateCurve(TDate baseDate, TDate* dates, double* rates, int numPts,double basis=1,
		       long daycountConv=GTO_ACT_365F);

	~KRateCurve(){destroy();}

	KRateCurve(const KRateCurve&);
	KRateCurve& operator=(const KRateCurve&);

	void shiftCurve(double x);
	int size() const{return _tcurve->fNumItems;}
	KDate valueDate() const{return KDate(_tcurve->fBaseDate);}
	KDate endDate() const {return KDate(_tcurve->fArray[_tcurve->fNumItems-1].fDate);}

	int   get_index(KDate date)
	{
		int i=0;
		while(date < KDate(_tcurve->fArray[i].fDate))
		{
			i++;
		}
		return i;
	}

	KDate iDate(int i){ return KDate(_tcurve->fArray[i].fDate);}
	double iRate(int i){return double(_tcurve->fArray[i].fRate);}
	//annual compounding rate [startDate, endDate] with daycountconvention ACT_365
	double  get_rate(int interpType, const KDate startDate, TFloatRate* fRateInfo) const
	{
		double rate;
		if(GtoForwardRate(_tcurve,
						  interpType,
						  fRateInfo,
						  startDate,
						  &rate) == FAILURE)
		{
			throw KException(" failed to get floating rate.");
		}

		return rate;
	}
	double	get_rate(int interpType, const KDate &endDate, const KDate &startDate = 0) const;
	double	get_rate(const KDate &endDate) const
	{
		return get_rate(FLAT_FORWARD,endDate);
	}
	double	get_pv(int interpType, const KDate &d) const;
	double	get_pv( const KDate &d) const
	{
		return get_pv(FLAT_FORWARD,d);
	}
	KInput	print()const;
	friend	std::ostream & operator<<(std::ostream & out, const KRateCurve& c)

	{// out<<c.print(); 
		return out;}
private:
	void	destroy(){	if (_tcurve)	GtoFreeTCurve (_tcurve);}
	TCurve* _tcurve;
};


class	KIndexRate
{	
	//Index rate _weight *(I + spread)
	int		_matInMonth;
    int		_freq;
	int		_dayCountConv;
	double	_spread;
	double	_weight;
public:
	KIndexRate(void){}
	KIndexRate(	int	mat,int	f,int dayCount){
						_matInMonth = mat; _freq = f; 
						_dayCountConv = dayCount;
						_spread = 0.; _weight = 1.;}
	~KIndexRate(void){}
	int	matInMonth()const{return _matInMonth;}
	int	freq()const{return _freq;}
	int		dayCountConv()const{return _dayCountConv;}
	double	spread()const{return _spread;}
	double	weight()const{return _weight;}

	KIndexRate	operator + (double s) const;	
	KIndexRate	operator * (double w) const;	
};

class KCashFlow:public KCurve<KDate, double> 
{
public:
	KCashFlow(){}
	KCashFlow(const KVector(KDate) &d, const KVector(double) &t):KCurve<KDate, double>(d, t){}

	~KCashFlow(){}
	KCashFlow	operator + (const KCashFlow &c)const
	{	KCashFlow	ans(*this); ans += c; return ans;}
	//present value of all cash flows from valuedate
	double	get_pv(const KRateCurve &zero)const;
	//present value of all cash flows from valuedate
	//discounted with flat yield 
	double	get_pv(const KDate &valueDate, double y)const;
	//flat yield for given pv
	double	get_yield(const KDate &valueDate, double p, double guess = 0.05);
//the only purpose of operator() is to be able to use solver
//use with caution!!!!
	double operator()(double y)const{return get_pv(_startDate, y);} 
private:
	//temp variable, in order to use solver
	KDate _startDate;
};
#endif

