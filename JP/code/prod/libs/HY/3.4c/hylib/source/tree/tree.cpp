
/**********************************************************
*
*	tree.cpp
*
***********************************************************/
#include "tree.h"
#include <math.h>

extern std::ofstream out;

void    BaseTree::timeLineInit(int ppy)
{
	KDate	valueDate,
			tmpDate,
			curDate,
			nextDate;		
    int		i, nodePerPeriod;
	double	period;

    valueDate = _zero.valueDate();
	//Add valueDate to time line
	_tpLine.insert(valueDate);

	KTimeLine::iterator	iter, nIter;
	// erase all dates before value date and dates after last date
	//add enough dates to have the right ppy
	iter = _tpLine.find(valueDate);
	_tpLine.erase(_tpLine.begin(), iter);
	//find last date
	for(iter =  _tpLine.begin(); iter!= _tpLine.end(); iter++)
		if(iter->second[__lastTreeDate] != 0)
			break;
	if(iter != _tpLine.end())
	{
		iter++;
		_tpLine.erase(iter, _tpLine.end());
	}
	
	iter = _tpLine.begin();

	while(iter != _tpLine.end())
	{
		nIter = iter;
		nIter++;
		if(nIter == _tpLine.end())
			break;

		period =  GetDaysDiff(iter->first, nIter->first);
		nodePerPeriod = (int)ceil(ppy * period / 365);
		
		for (i = 1; i < nodePerPeriod; i++)
			_tpLine.insert(iter->first.dFwd(i * period / nodePerPeriod));
		iter++;
	}  
	//add one more date
	nIter = _tpLine.end();
	nIter--;
	iter = nIter;
	iter--;
	_tpLine.insert(nIter->first.dFwd(GetDaysDiff(iter->first,nIter->first)));
	
	//fill zero rate and price
	int yy=1;
	for(iter = _tpLine.begin(); iter !=_tpLine.end(); iter++)
	{	
		iter->second[__timeFromValueDate] = GetYearDiff(valueDate, iter->first);
		iter->second[__zeroPrice] = _zero.get_pv(FLAT_FORWARD, iter->first);

	//	out<<"i="<<yy<<", "<<GtoFormatDate(iter->first)<<"zero="<<iter->second[__zeroPrice]<<"\n"<<std::endl;
	//	yy++;
	}
	/*fill fwd rate and time between tp for time points*/
	
	for(iter = _tpLine.begin(); iter !=_tpLine.end(); iter++)
	{
		nIter = iter;
		nIter++;
		if(nIter == _tpLine.end())
			break;
		iter->second[__timeTillNextTp] = GetYearDiff(iter->first, nIter->first);

//		out<<"i="<<yy<<", "<<GtoFormatDate(iter->first)<<" , "<<GtoFormatDate(nIter->first)<<"\n"<<std::endl;
		
		/* forward rate between tp and nTp*/
		iter->second[__fwdRateTillNextTp] = iter->second[__zeroPrice]/nIter->second[__zeroPrice] - 1.;
//		out<<"i="<<yy<<", zero1="<<iter->second[__zeroPrice]<<", zero2="<<nIter->second[__zeroPrice]<<", fwd="<<iter->second[__fwdRateTillNextTp]<<")\n";
//		yy++;
	}
	nIter = _tpLine.end();
	nIter--;
	iter = nIter;
	iter--;
	nIter->second[__timeTillNextTp] = iter->second[__timeTillNextTp];
	nIter->second[__fwdRateTillNextTp] = iter->second[__fwdRateTillNextTp];

}


void	BaseTree::init(int ppy, bool rebuildTimeLine )
{
	if(rebuildTimeLine == true)
	{
		timeLineInit(ppy);
	}
	//add additional data with default empty
	addData();
	//volatility calibration
	treeSpotVol();
	//compute limits for each time slices
	limits();
	//compute drift 
	drift();
	_tlPos = _tpLine.end();
}

/**********************************************************************
*
*	roll_back: Cleaning all all zero prices which are no longer
*	needed, compute all necessary discounting information, discount all
*	zero prices and add new zeros needed from now on.
*
***********************************************************************/
bool	BaseTree::roll_back()
{
	bool	status = false;
	if(_tlPos != _tpLine.begin())
	{
   		if(_tlPos == _tpLine.end())
		{
			_tlPos--;
			_temp.resize(this);
			_tempFunc.resize(this);
			memInit();
		
		}
		else
		{
			_tlPos--;
			/*compute discount information*/
			buildLattice();
		}
		status = true;
	}
	return status;
}


