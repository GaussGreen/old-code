/****************************************************************************/
/*      Declaration file for structure of data.                             */
/****************************************************************************/

#ifndef	_BASE_TREE_H_
#define	_BASE_TREE_H_

#include "kplatform.h"

#include <iostream>
#include <vector>
#include "timeslice.h"
#include "kdate.h"
#include "kprocess.h"
#include "kcurve.h"
#include "kptr.h"
#include "kfunction.h"
#include "kutility.h"

//const	double	INTERVAL_RATIO = .1;
class BaseTree
{
	KRateCurve		_zero;
	KTimeLine		_tpLine;
	KTimeLine::iterator
					_tlPos;
protected:
	//tmp variables to speeding up dev operation
	DTimeSlice	_temp;
	FTimeSlice	_tempFunc;
protected:
	void	timeLineInit(int ppy);
	typedef	KTimeLine::iterator	iterator;
	typedef	KTimeLine::const_iterator	const_iterator;
	KTimeLine::iterator tl_valueDate(){ KDate valueDate = _zero.valueDate(); return _tpLine.find(valueDate); }
	KTimeLine::iterator tl_begin(){ return _tpLine.begin(); }
	KTimeLine::const_iterator tl_begin() const{ return _tpLine.begin(); }
	KTimeLine::iterator tl_end(){ return _tpLine.end(); }
	KTimeLine::const_iterator tl_end() const { return _tpLine.end(); }
	size_t		tl_size()const{return _tpLine.size();}
	KTimePoint	get_prev_tp()const{KTimeLine::iterator iter = _tlPos; iter--; return iter->second;}
	KTimePoint	get_next_tp()const{KTimeLine::iterator iter = _tlPos; iter++; return iter->second;}
public:
	BaseTree(void){}
	~BaseTree(void){}
	//time line handles
//	KDate   get_valueDate()const{return tl_begin()->first;}
	KDate   get_valueDate()const{return _zero.valueDate();}
	KDate	get_curr_date()const{return _tlPos->first;}
	KDate	get_next_date()const{KTimeLine::iterator iter = _tlPos; iter++; return iter->first;}
	KTimePoint	get_curr_tp()const{return _tlPos->second;}
	void	insert_tp(const KDate &d, const KTimePoint &x = KTimePoint())
				{_tpLine.insert(d, x);}
	void	insert_tp(const KDate &d, const KString &s, double x)
				{KTimePoint	t; t[s] = x; _tpLine.insert(d, t);}
	void	set_last_date(const KDate &d){insert_tp(d, __lastTreeDate, 1.);}
	void	printTimeLine(std::ostream &out){out<<_tpLine;}
		
	KRateCurve	get_index_zero()const{return _zero;}
	//set_index_zero needs to be set before calling init
	void	set_index_zero(const KRateCurve &zr){_zero = zr;}
	//tree initialization routine
	void	init(int ppy,bool rebuildTimeLine = true);
	//roll back tree return true if curTlPos is valid, false otherwise 
	bool	roll_back();
protected:
	//a handle which enables users to add more data before calling treeSpotVol
	virtual	void	addData(void){}
	virtual	void	treeSpotVol(void ) = 0;
	virtual	void	limits(void) = 0; 
	virtual	void	drift(void) = 0; 
	virtual	void	buildLattice() = 0;
	virtual	void	memInit(void) = 0;
	virtual	void	initDebtDrift(void) = 0;		//HY4.2
public:
	virtual	KValarray<int>	get_curr_limit(void) const=0;
	virtual	KValarray<int>	get_next_limit(void) const=0;
	virtual	KValarray<int>	get_max_limit(void) const=0;
//	void	dev(DTimeSlice &){throw KException("dev is not defined!");}
//	void	dev(FTimeSlice &){throw KException("dev is not defined!");}

	void dev(DTimeSlice &ts1, bool isVolShift = false)
				{throw KException("dev is not defined!");}	//HY3.4v

	void dev(FTimeSlice &ts1, bool isVolShift = false)
				{throw KException("dev is not defined!");}	//HY3.4v

	virtual	void	print(std::ostream &out)const
	{throw KException("print is not defined!");}
};




#endif




