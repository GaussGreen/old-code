/**********************************************************
*
*	zrbank.cpp
*
**********************************************************/
#include "tree.h"
#include <math.h>
#include "zrbank.h"

void	ZeroBank::update(BaseTree *tr)
{  
	KDate	&curDate = tr->get_curr_date();
	KMap(KDate, DS_Pair)::iterator	iter;
	//clean
	iter = _zrSlice.begin();
	while(iter != _zrSlice.end())
	{
		if(iter->second.first > curDate)
			iter = _zrSlice.erase(iter);
		else
			iter++;
	}
	//dev
	for(iter = _zrSlice.begin(); iter != _zrSlice.end(); iter++)
	{
		if(curDate < iter->first)
			tr->dev(iter->second.second);
	}
	//init
	iter = _zrSlice.find(curDate);
	if(iter != _zrSlice.end())
		iter->second.second.resize(tr, 1.);

}

void	ZeroBank::insert_date(BaseTree *tr)const
{
	KMap(KDate, DS_Pair)::const_iterator	iter;
	for(iter = _zrSlice.begin(); iter != _zrSlice.end(); iter++)
	{
		tr->insert_tp(iter->first);
		tr->insert_tp(iter->second.first);
	}
	 
}

void	ZeroBank::add_index(const KDate &setDate, const KIndexRate	&rate)
{
	KMap(KDate, DS_Pair)::iterator	iter;
	KDate	zrDate;
	int	numParCoupon = rate.matInMonth() * rate.freq() / 12;	
	for(int i = numParCoupon; i > 0; i--)
	{
		zrDate = setDate + int(i * 12 /rate.freq()) * ONE_MONTH;
		iter = _zrSlice.find(zrDate);
		if(iter != _zrSlice.end())
			iter->second.first = Min(setDate, iter->second.first); 
		else
		{
			_zrSlice.insert(KMap(KDate, DS_Pair)::value_type(
				zrDate, std::make_pair(setDate, DTimeSlice())));
		}
	}
}
 
DTimeSlice	ZeroBank::get_index_rate ( const KDate & curDate,  const KIndexRate	&rate)const 	
{
	double	couponPeriod;                       
	KDate	date1, date2;
	KMap(KDate, DS_Pair)::const_iterator	iter;
	
	DTimeSlice	annuity;
    int	nbReset = rate.matInMonth() * rate.freq() / 12;  
	date1 = curDate;
    for (int i = 1; i <= nbReset; i++)
    {
		date2 = curDate + (int) (12 * i / rate.freq()) * ONE_MONTH;
		iter = _zrSlice.find(date2);
		if(iter == _zrSlice.end())
			throw	KException("No zero slice with maturity\t")<<date2;
		couponPeriod = GetYearDiff(date1, date2,rate.dayCountConv());
		if(i == 1)
			annuity = couponPeriod * iter->second.second;
		else
			annuity += couponPeriod * iter->second.second;
		date1 = date2;
	} 
	return rate.weight() *((1. - iter->second.second)/annuity + rate.spread());
}   

