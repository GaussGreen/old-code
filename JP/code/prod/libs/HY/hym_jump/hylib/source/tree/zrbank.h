#ifndef	_ZR_BANK_H_
#define	_ZR_BANK_H_

#include "kplatform.h"
#include "kdate.h"
#include "kcurve.h"
#include "timeslice.h"
#include "kplatdep.h"
#include <map>


//ZeroBank is ordered according to ascending order of end date
class	ZeroBank
{
public:
	ZeroBank(void){}
	~ZeroBank(void){} 
	void	add_index(const KDate &setDate, const KIndexRate	&rate);
	void	insert_date(BaseTree *)const;
	//remove all slices whose beginDate > currentDate, dev all slices whose end date
	//>currentDate and initialize all slices whose end date is cuurent date
	void	update(BaseTree *tr);
	DTimeSlice	get_index_rate( const KDate & curDate,  const KIndexRate	&rate)const;
private:
	typedef	std::pair<KDate, DTimeSlice> DS_Pair;
	KMap(KDate, DS_Pair)	_zrSlice;
};


#endif